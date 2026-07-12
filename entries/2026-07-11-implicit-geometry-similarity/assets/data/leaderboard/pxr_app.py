import io
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import boto3
import gradio as gr
import numpy as np
import pandas as pd
from config import (
    ACTIVITY_DATASET_SIZE,
    ACTIVITY_LEADERBOARD_S3,
    CURRENT_PHASE,
    HOURS_BETWEEN_SUBMISSIONS,
    PHASE_2_ACTIVITY_ENTRIES_S3,
    REQUIRED_ACTIVITY_COLUMNS,
    S3_BUCKET,
    STRUCTURE_DATASET_SIZE,
    STRUCTURE_LEADERBOARD_S3,
)
from gradio.themes.utils import sizes
from gradio_leaderboard import ColumnFilter, Leaderboard
from head_to_head import render_comparisons
from loguru import logger
from models import Submission, _safeify_username
from submission_store import _fetch_last_submission_date, upload_submission

# S3FS was causing 403 errors in testing on some machines, unclear why.
# Switching to boto3 client which works reliably
s3_client = boto3.client("s3", region_name="us-east-1")


def select_correct_leaderboard_to_display(track: str, phase: int) -> str:
    """Determine which leaderboard CSV to load based on track and phase."""
    track_paths = {
        "activity": ACTIVITY_LEADERBOARD_S3,
        "structure": STRUCTURE_LEADERBOARD_S3,
    }
    activity_phase_namespace = {1: "interim", 2: "all", 0: "final"}
    structure_phase_namespace = {1: "interim", 2: "interim", 0: "final"}
    try:
        base_path = track_paths[track]
    except KeyError as exc:
        raise ValueError(f"Invalid track: {track}") from exc
    try:
        suffix = (
            activity_phase_namespace
            if track == "activity"
            else structure_phase_namespace
        )[phase]
    except KeyError as exc:
        raise ValueError(f"Invalid phase: {phase}") from exc
    new_path = base_path.replace("/interim/", f"/{suffix}/")
    logger.info(
        f"Selected leaderboard path for track={track} phase={phase}: {new_path}"
    )
    return new_path


def sort_leaderboard_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Sort leaderboard columns in a consistent order."""
    start_cols = ["rank", "username", "Submitted", "model_report_link"]
    end_cols = ["Proprietary Data", "CLD"]
    start_cols = [c for c in start_cols if c in df.columns]
    end_cols = [c for c in end_cols if c in df.columns]
    middle_cols = [c for c in df.columns if c not in start_cols + end_cols]
    ordered_cols = start_cols + middle_cols + end_cols
    ordered_cols = [c for c in ordered_cols if "Unnamed" not in c]
    return df[ordered_cols]


def make_user_clickable(name: str) -> str:
    if not name.startswith("https://"):
        link = f"https://huggingface.co/{name}"
    else:
        link = name
    return f'<a target="_blank" href="{link}" style="color: var(--link-text-color); text-decoration: underline;text-decoration-style: dotted;">{name}</a>'


def make_tag_clickable(tag: str) -> str:
    if not tag or (tag == ""):
        return "Not submitted"
    if not tag.startswith(("http://", "https://")):
        return "Invalid link"
    return f'<a target="_blank" href="{tag}" style="color: var(--link-text-color); text-decoration: underline;text-decoration-style: dotted;">link</a>'


def hide_username_for_anonymous_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Replace usernames with aliases for anonymous entries."""
    df.loc[df["anonymous"], "username"] = df.loc[df["anonymous"], "user_alias"]
    return df.drop(columns=["user_alias", "anonymous"], errors="ignore")


_ACTIVITY_EMPTY = pd.DataFrame(
    columns=[
        "rank",
        "username",
        "Submitted",
        "model_report_link",
        "Proprietary Data",
        "MAE",
        "RAE",
        "R2",
        "Spearman ρ",
        "Kendall's τ",
    ]
)
_STRUCTURE_EMPTY = pd.DataFrame(
    columns=[
        "rank",
        "username",
        "Submitted",
        "model_report_link",
        "Proprietary Data",
        "LDDT-PLI",
        "BiSyRMSD",
        "LDDT-LP",
        "Coverage",
    ]
)


def _prepare_activity_df(df: pd.DataFrame, for_download: bool = False) -> pd.DataFrame:
    """Sort and rename activity leaderboard columns (no HTML).

    Args:
        df: Raw leaderboard DataFrame from S3.
        for_download: If True, keep mean and std as separate columns with full
            precision for the downloadable CSV. If False (default), drop std
            columns and round means to 4 dp for the live leaderboard.

    Returns:
        Prepared DataFrame.

    """
    df = df.sort_values("rank", ascending=True).reset_index(drop=True)
    rename_map = {
        "MAE_mean": "MAE",
        "RAE_mean": "RAE",
        "R2_mean": "R2",
        "Spearman_R_mean": "Spearman ρ",
        "Kendall's_Tau_mean": "Kendall's τ",
        "MAE_std": "MAE std",
        "RAE_std": "RAE std",
        "R2_std": "R2 std",
        "Spearman_R_std": "Spearman ρ std",
        "Kendall's_Tau_std": "Kendall's τ std",
        "submitted_at": "Submitted",
    }
    df = df.rename(columns=rename_map)
    if for_download:
        pass  # keep mean and std columns as-is with full precision
    else:
        std_cols = [c for c in df.columns if c.endswith(" std")]
        df = df.drop(columns=std_cols)
        for col in ["MAE", "RAE", "R2", "Spearman ρ", "Kendall's τ"]:
            if col in df.columns:
                df[col] = df[col].round(4)
    df["Submitted"] = pd.to_datetime(df["Submitted"], utc=True).dt.strftime(
        "%Y-%m-%d %H:%M UTC"
    )
    df["Proprietary Data"] = (
        df["used_proprietary_data"].fillna(False).map({True: "Yes", False: "No"})
        if "used_proprietary_data" in df.columns
        else "No"
    )
    df = df.drop(columns=["used_proprietary_data"], errors="ignore")
    return sort_leaderboard_columns(df)


def _prepare_structure_df(df: pd.DataFrame, for_download: bool = False) -> pd.DataFrame:
    """Sort and rename structure leaderboard columns (no HTML).

    Args:
        df: Raw leaderboard DataFrame from S3.
        for_download: If True, collapse mean/std pairs into 'XX±YY' strings for
            the downloadable CSV. If False (default), keep mean values as plain
            floats for numeric sorting in the live leaderboard.

    Returns:
        Prepared DataFrame.

    """
    df = df.sort_values("LDDT-PLI_mean", ascending=False).reset_index(drop=True)
    rename_map = {
        "LDDT-PLI_mean": "LDDT-PLI",
        "BiSyRMSD_mean": "BiSyRMSD",
        "LDDT-LP_mean": "LDDT-LP",
        "Ligand_RMSD_mean": "Ligand RMSD",
        "coverage_mean": "Coverage",
        "LDDT-PLI_std": "LDDT-PLI std",
        "BiSyRMSD_std": "BiSyRMSD std",
        "LDDT-LP_std": "LDDT-LP std",
        "Ligand_RMSD_std": "Ligand RMSD std",
        "coverage_std": "Coverage std",
    }
    df = df.rename(columns=rename_map)
    if for_download:
        pass  # keep mean and std columns as-is with full precision
    else:
        std_cols = [c for c in df.columns if c.endswith(" std")]
        df = df.drop(columns=std_cols)
        for col in ["LDDT-PLI", "BiSyRMSD", "LDDT-LP", "Ligand RMSD", "Coverage"]:
            if col in df.columns:
                df[col] = df[col].round(4)
    df = df.rename(columns={"submitted_at": "Submitted"})
    df["Submitted"] = pd.to_datetime(df["Submitted"], utc=True).dt.strftime(
        "%Y-%m-%d %H:%M UTC"
    )
    df["model_report_link"] = df["model_report_link"].fillna("")
    df["Proprietary Data"] = (
        df["used_proprietary_data"].fillna(False).map({True: "Yes", False: "No"})
        if "used_proprietary_data" in df.columns
        else "No"
    )
    df = df.drop(columns=["used_proprietary_data"], errors="ignore")
    return sort_leaderboard_columns(df)


def _load_csv_from_s3(key: str) -> pd.DataFrame:
    """Load a CSV file from S3 into a DataFrame."""
    obj = s3_client.get_object(Bucket=S3_BUCKET, Key=key)
    return pd.read_csv(io.BytesIO(obj["Body"].read()))


def load_activity_leaderboard() -> pd.DataFrame:
    """Load the activity leaderboard from S3."""
    logger.info("Refreshing activity leaderboard...")
    try:
        df = _load_csv_from_s3(
            select_correct_leaderboard_to_display("activity", CURRENT_PHASE)
        )
    except Exception as exc:
        logger.warning("Could not load activity leaderboard: {}", exc)
        return _ACTIVITY_EMPTY
    df = _prepare_activity_df(df)
    df["username"] = df["username"].map(make_user_clickable)
    df = hide_username_for_anonymous_entries(df)
    # Column was being loaded as a float, convert to str first
    df["model_report_link"] = df["model_report_link"].astype("string").fillna("")
    df["model_report_link"] = df["model_report_link"].map(make_tag_clickable)
    logger.info("Activity leaderboard loaded: {} entries.", len(df))
    return df


def load_activity_entries() -> pd.DataFrame:
    """Load the activity entries from S3."""
    logger.info("Refreshing activity entries...")
    try:
        df = _load_csv_from_s3(PHASE_2_ACTIVITY_ENTRIES_S3)
    except Exception as exc:
        logger.warning("Could not load activity entries: {}", exc)
        return pd.DataFrame(
            {
                "Username": ["OpenADMET_Baseline"],
                "Submitted": ["2024-04-01 00:00 UTC"],
                "Model Report Link": ["missing_link"],
            }
        )
    df = df[["display_username", "submitted_at", "model_report_link"]].copy()
    df = df.rename(
        columns={
            "display_username": "Username",
            "submitted_at": "Submitted",
            "model_report_link": "Model Report Link",
        }
    )
    df["Submitted"] = pd.to_datetime(df["Submitted"], utc=True).dt.strftime(
        "%Y-%m-%d %H:%M UTC"
    )
    df["Model Report Link"] = df["Model Report Link"].astype("string").fillna("")
    df["Model Report Link"] = df["Model Report Link"].map(make_tag_clickable)
    logger.info("Activity entries loaded: {} entries.", len(df))
    return df


def load_structure_leaderboard() -> pd.DataFrame:
    """Load the structure leaderboard from S3."""
    logger.info("Refreshing structure leaderboard...")
    try:
        df = _load_csv_from_s3(
            select_correct_leaderboard_to_display("structure", CURRENT_PHASE)
        )
    except Exception as exc:
        logger.warning("Could not load structure leaderboard: {}", exc)
        return _STRUCTURE_EMPTY
    df = _prepare_structure_df(df)
    df["username"] = df["username"].map(make_user_clickable)
    df = hide_username_for_anonymous_entries(df)
    df["model_report_link"] = df["model_report_link"].map(make_tag_clickable)
    logger.info("Structure leaderboard loaded: {} entries.", len(df))
    return df


def download_activity_leaderboard() -> str:
    """Write the activity leaderboard to a temp CSV and return the file path."""
    try:
        df = _load_csv_from_s3(
            select_correct_leaderboard_to_display("activity", CURRENT_PHASE)
        )
    except Exception as exc:
        logger.warning("Could not load activity leaderboard for download: {}", exc)
        df = _ACTIVITY_EMPTY
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", prefix="activity_leaderboard_", delete=False
        ) as f:
            df.to_csv(f, index=False)
            return f.name
    df = _prepare_activity_df(df, for_download=True)
    df = hide_username_for_anonymous_entries(df)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", prefix="activity_leaderboard_", delete=False
    ) as f:
        df.to_csv(f, index=False)
        return f.name


def download_structure_leaderboard() -> str:
    """Write the structure leaderboard to a temp CSV and return the file path."""
    try:
        df = _load_csv_from_s3(
            select_correct_leaderboard_to_display("structure", CURRENT_PHASE)
        )
    except Exception as exc:
        logger.warning("Could not load structure leaderboard for download: {}", exc)
        df = _STRUCTURE_EMPTY
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", prefix="structure_leaderboard_", delete=False
        ) as f:
            df.to_csv(f, index=False)
            return f.name
    df = _prepare_structure_df(df, for_download=True)
    df = hide_username_for_anonymous_entries(df)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", prefix="structure_leaderboard_", delete=False
    ) as f:
        df.to_csv(f, index=False)
        return f.name


def _format_submission_time_message(last_submission: datetime, track: str) -> str:
    """Format a message indicating when the user can next submit next."""
    track_name = "an activity" if track == "activity" else "a structure"
    next_submission_time = last_submission + pd.Timedelta(
        hours=HOURS_BETWEEN_SUBMISSIONS
    )
    time_remaining = next_submission_time - datetime.now(timezone.utc)
    seconds_left = max(0, int(time_remaining.total_seconds()))
    hours, rem = divmod(seconds_left, 3600)
    minutes, seconds = divmod(rem, 60)
    wait_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    return (
        f"Error: You submitted {track_name} prediction on "
        f"{last_submission.strftime('%Y-%m-%d %H:%M:%S (UTC)')}.\n"
        f"Please wait {wait_str} before submitting again."
    )


def submit_predictions(
    username,
    user_alias,
    anon_checkbox,
    participant_name,
    discord_username,
    email,
    affiliation,
    model_tag,
    paper_checkbox,
    proprietary_data_checkbox,
    track_select,
    file_input,
):
    """Handle a competition submission from the web UI or API.

    Validates required fields and the uploaded file, then records the submission.

    Args:
        username: HuggingFace username (required).
        user_alias: Optional alias for anonymous display on the leaderboard.
        anon_checkbox: If True, display alias instead of username.
        participant_name: Real name (private, not displayed publicly).
        discord_username: Discord handle (optional).
        email: Contact email (optional).
        affiliation: Institutional affiliation (optional).
        model_tag: Link to method report (optional).
        paper_checkbox: Opt-in for future publication inclusion.
        proprietary_data_checkbox: Whether proprietary data was used in training.
        track_select: "Activity Prediction" or "Structure Prediction".
        file_input: Path to the uploaded submission file.

    Returns:
        gr.update with a status message and visible=True.

    """
    # --- required field validation ---
    if not username or not username.strip():
        return gr.update(
            value="Error: Hugging Face username is required.", visible=True
        )
    if not track_select:
        return gr.update(value="Error: Please select a track.", visible=True)
    if file_input is None:
        return gr.update(value="Error: Please upload a submission file.", visible=True)

    file_path = Path(file_input)

    # --- file format validation ---
    if track_select == "Activity Prediction":
        suffix = file_path.suffix.lower()
        if suffix == ".parquet":
            try:
                df = pd.read_parquet(file_path)
            except Exception as exc:
                return gr.update(
                    value=f"Error: Could not read parquet file — {exc}", visible=True
                )
        elif suffix == ".csv":
            try:
                df = pd.read_csv(file_path)
            except Exception as exc:
                return gr.update(
                    value=f"Error: Could not read CSV file — {exc}", visible=True
                )
        else:
            return gr.update(
                value="Error: Activity submissions must be a .parquet or .csv file.",
                visible=True,
            )
        if len(df) != ACTIVITY_DATASET_SIZE:
            return gr.update(
                value=f"Error: Expected {ACTIVITY_DATASET_SIZE} rows, got {len(df)}.",
                visible=True,
            )
        missing = REQUIRED_ACTIVITY_COLUMNS - set(df.columns)
        if missing:
            return gr.update(
                value=f"Error: Missing required columns: {missing}", visible=True
            )
        if df["pEC50"].isnull().any():
            return gr.update(
                value="Error: pEC50 column contains NaN values.", visible=True
            )
        if not np.isfinite(df["pEC50"]).all():
            return gr.update(
                value="Error: pEC50 column contains infinite values.", visible=True
            )
        last_submission = _fetch_last_submission_date(
            "activity", _safeify_username(username.strip())
        )
        logger.info(
            f"Last submission date for user {username.strip()!r}: {last_submission}"
        )
        if (
            last_submission
            and (datetime.now(timezone.utc) - last_submission).total_seconds()
            < HOURS_BETWEEN_SUBMISSIONS * 3600
        ):
            return gr.update(
                value=_format_submission_time_message(
                    last_submission, track="activity"
                ),
                visible=True,
            )

    elif track_select == "Structure Prediction":
        if file_path.suffix.lower() != ".zip":
            return gr.update(
                value="Error: Structure submissions must be a .zip file.", visible=True
            )
        try:
            with zipfile.ZipFile(file_path) as zf:
                n_files = len(zf.namelist())
        except Exception as exc:
            return gr.update(
                value=f"Error: Could not read zip file — {exc}", visible=True
            )
        if n_files != STRUCTURE_DATASET_SIZE:
            return gr.update(
                value=f"Error: Expected {STRUCTURE_DATASET_SIZE} files in zip, got {n_files}.",
                visible=True,
            )
        last_submission = _fetch_last_submission_date(
            "structure", _safeify_username(username.strip())
        )
        logger.info(
            f"Last submission date for user {username.strip()!r}: {last_submission}"
        )
        if (
            last_submission
            and (datetime.now(timezone.utc) - last_submission).total_seconds()
            < HOURS_BETWEEN_SUBMISSIONS * 3600
        ):
            return gr.update(
                value=_format_submission_time_message(
                    last_submission, track="structure"
                ),
                visible=True,
            )

    # --- build submission model and persist to S3 ---
    submission = Submission(
        username=username.strip(),
        user_alias=user_alias.strip(),
        anonymous=anon_checkbox,
        participant_name=participant_name.strip(),
        discord_username=discord_username.strip(),
        email=email.strip(),
        affiliation=affiliation.strip(),
        model_report_link=model_tag.strip(),
        include_in_publication=paper_checkbox,
        used_proprietary_data=proprietary_data_checkbox,
        track=track_select,
        filename=file_path.name,
    )
    logger.info(
        f"Submission received: id={submission.submission_id!r} "
        f"user={submission.username!r} track={submission.track!r}"
    )
    upload_submission(submission, file_path)

    display_name = submission.display_name
    return gr.update(
        value=(
            f"Submission received from {display_name!r} for the {track_select} track. "
            "Your predictions are being processed and will appear on the leaderboard within 2 hours."
        ),
        visible=True,
    )


# Theme and Title setup
with gr.Blocks(
    title="OpenADMET PXR Blind Challenge",
    fill_height=False,
    theme=gr.themes.Ocean(text_size=sizes.text_lg),
) as demo:
    ### Header Section
    with gr.Row():
        with gr.Column(scale=7):
            gr.Markdown("""
                ## OpenADMET Blind Challenge: Predicting PXR Induction 🧬
                The next OpenADMET blind challenge focuses on predicting **Pregnane-X Receptor (PXR)** induction. PXR is a nuclear hormone receptor and master regulator of drug-metabolizing enzymes and transporters. Because compounds that induce PXR can derail drug discovery projects by causing adverse interactions, accurate prediction is critical.

                **Two Tracks are Available:**
                1. **Activity Prediction:** Predict pEC50 values for a test set of 513 compounds (split into two stages).
                2. **Structure Prediction:** Predict the bound structures of 184 ligands to PXR.

                Go to the **Leaderboard** to see current standings. To participate, head to the **Submit** tab to upload your predictions.

                **Challenge Timeline:** Submissions open April 1st. Phase 1 concludes May 25th; Phase 2 and the Structure Track conclude July 1st.
""")
        with gr.Column(scale=2):
            gr.Image(
                value="./_static/challenge_logo.png",  # Update this to your new logo path
                show_label=False,
                show_download_button=False,
                width="400px",
            )

    # --- Welcome Markdown Content ---
    welcome_md = """
# 💊 OpenADMET PXR Blind Challenge

## Background: Why PXR Matters

Evaluating PXR liabilities is a fundamental pillar of a late-stage ADMET cascade. PXR functions as a xenobiotic sensor, detecting foreign compounds and marshalling drug-metabolizing enzymes and transporters by activating their transcription. It primarily regulates **CYP3A4**, the enzyme responsible for metabolizing approximately **50% of all marketed drugs**.

Activation of PXR can lead to:
- **Drug-Drug Interactions (DDIs):** Accelerated metabolism can reduce co-administered drug concentrations to sub-optimal levels.
- **Hepatotoxicity:** Increased production of reactive, toxic metabolites.
- **Chemoresistance:** Enhanced clearance of chemotherapeutic agents in tumor cells.

Drug discovery teams face a unique challenge with PXR due to its large, flexible ligand-binding pocket, which accommodates a wide range of chemical structures. PXR is also relatively underrepresented in the literature, with only **~800 high-quality pEC50 values** from nearly 150 papers in ChEMBL.

## The Activity Dataset

At Octant, OpenADMET has generated a PXR induction dataset of more than **11,000 compounds** using a low-cost, high-fidelity in-house assay. Compounds were sourced primarily from two Enamine libraries (Discovery Diversity 10 set and FDA Approved Drugs set) along with subsequent orders of follow-on compounds, and profiled through a rigorous multi-step assay flow reminiscent of an on-target drug discovery program.

The dataset was built through the following stages:

- **Primary Screen:** 11,362 diverse compounds screened at a single concentration.
- **Dose-Response:** 4,325 compounds selected for an 8-concentration dose-response (with extensive counter-screening in a PXR-null cell line to evaluate specificity).
- **Refinement:** 211 compounds showed EC50 ≤ 1 µM (pEC50 ≥ 6).
- **Counter-Screen:** 63 compounds selected based on minimal activity in a PXR-null cell line to confirm on-target specificity.
- **Analog Expansion Set:** Similarity searches (ECFP4 Tanimoto > 0.4) of these 63 actives yielded the **513-compound test set**, ordered from the Enamine US on-demand catalog and fully assayed with dose-response curves.

This design mimics a **lead optimization scenario**, shifting from broad hit-finding to detailed exploration of Structure-Activity Relationships (SAR). The analog set contains detailed SAR and activity cliffs that should prove challenging for models. Cumulatively, this represents the **largest PXR activity dataset available in the literature**.

### Developing the Assay

The assay employs a cell-based reporter system to measure PXR induction using a well-established two-part chimeric design. This fusion protein is composed of the ligand-binding domain (LBD) of human PXR attached to a heterologous DNA-binding domain, which acts via a reporter construct containing the corresponding DNA response element upstream of a luciferase gene. This approach reduces crosstalk, provides superior specificity, and maximizes signal-to-noise ratio compared to using the native PXR promoter. The system is stably integrated into the selected cell line to ensure consistent reporter activity across screening runs.

When a test compound acts as an agonist, it binds the PXR-LBD and induces a conformational change that promotes recruitment of transcriptional activators, driving luciferase expression. Multiple cell lines and genetic constructs were evaluated and optimized against reference compounds to achieve robust, reproducible data consistent with literature pEC50 values.

A parallel **counter-screen** using an identical reporter system with nonsense mutations in the chimeric PXR gene eliminates false positives — filtering out general transcriptional activators and HDAC inhibitors from true PXR agonists.

### How the Dataset Was Constructed

To determine the optimal concentrations for the primary screen, a pilot screen was run at 10, 30, and 100 µM. The 10 and 30 µM concentrations were selected, as their respective hit rates (~17% and ~51%) yielded a manageable set for follow-up dose-response curves while biasing toward the most potent compounds and mitigating solubility issues. A few thousand compounds yielded enough activity to be promoted to a full 8-point concentration dose response curve (DRC). By fitting these DRCs, the EC50s for these few thousand compounds were estimated. These were combined with data from initial direct-to-DRC experiments conducted early in the program's development, yielding a training dataset of 4,140 EC50 values.

For hit-calling, a linear model was fitted on logged data with a fixed effect term for each compound contrasted against the negative control. This assigns standard errors, p-values, and confidence intervals based on replicated control conditions. The Benjamini-Hochberg method was applied to control the false discovery rate (FDR). A compound is classified as a hit when its log₂ fold-change exceeds 1 and FDR < 5%.

Hit expansion targeted compounds with EC50 ≤ 1 µM and at least 1.5 log-unit difference between the primary assay pEC50 and the counterassay pEC50. This yielded 63 selective hits, from which analogs were selected from the Enamine US on-demand catalog with ECFP4 Tanimoto similarity > 0.4, forming the 513-compound test set.

## The Structure Dataset

PXR's large, flexible binding pocket is highly dynamic and capable of recognising ligands of vastly different sizes and shapes — a structural plasticity that represents a significant challenge for structure-based design.

The structure dataset comprises **184 small molecules** — a mix of fragment-sized compounds and active compounds from the activity track — selected for the highest quality electron density. X-ray crystal structures have been determined at **UCSF (Fraser Lab)** but remain blinded until the challenge concludes. Fragments were soaked into apo crystals in the P2₁2₁2₁ crystal form at a nominal concentration of 10 mM. Data were collected at **NSLS-II** using the AMX and FMX beamlines. Data were reduced using Autoproc; electron density maps were analysed for fragment binding events using PanDDA. Ligands were modelled in COOT and refined with phenix.refine.

In addition, **68 structures from the PDB** have been re-refined and will be released as part of the training data package.

> **Note:** A small number of late-breaking crystal structures of Enamine deck compounds may be added to the structure track test set early in the challenge period. We will announce any additions in the `#pxr-challenge` channel on [Discord](https://discord.gg/MY5cEFHH3D). In the meantime, please submit as normal using the current 184-compound set.

## 🧪 The Challenge Tracks

Participants can compete in either or both tracks.

### 1. Activity Prediction Track

Participants predict **pEC50 values** for the 513-compound analog set. An extensive data package will be provided for the training set, including PXR pEC50 and Emax, null-line pEC50 and Emax, and supporting raw data.

The track proceeds in two phases:
- **Phase 1:** Predict activity for all 513 compounds. Analog Set 1 will serve as a **live leaderboard** during this period.
- **Phase 2:** EC50 values for Analog Set 1 are unblinded; participants then refine predictions for the remaining **Analog Set 2**. No live leaderboard — predictions are fully blinded until the deadline.

The **primary evaluation metric** is **RAE** (Relative Absolute Error) on pEC50. Extensive secondary metrics (MAE, R², Spearman ρ, Kendall's τ) and error estimation via bootstrapping are also reported.

### 2. Structure Prediction Track

Participants predict the **bound protein-ligand complex** for each of the 184 ligands, given their SMILES strings. Participants may use any computational approach, including protein structure prediction tools, docking algorithms, or existing PDB structures.

A **live leaderboard** is maintained using half of the 184 structures; the remaining half are held out and only scored at the final deadline.

#### Submission format
Submit a `.zip` archive containing exactly **184 PDB files**, one per ligand. The ligand must be present in each file with residue name **`LIG`** — this is how our scoring pipeline identifies the small molecule in the complex.

#### How scoring works
Each predicted complex is scored against the crystallographic reference using **[OpenStructure](https://openstructure.org/)** (OST). Two scorers are run per ligand pair:

- **SCRMSDScorer** — superimposes the binding site (Cα atoms within 8 Å of the reference ligand) and computes the symmetry-corrected ligand RMSD, accounting for all valid atomic symmetry mappings. Also yields LDDT-LP as a by-product.
- **LDDTPLIScorer** — computes LDDT-PLI without requiring superposition, evaluating whether predicted protein–ligand contacts are preserved relative to the reference structure.

If a structure contains multiple ligand assignments, the best-scoring pair (highest LDDT-PLI, then lowest BiSyRMSD) is used. All metrics are bootstrapped over 1000 resamples of the compound set to enable statistical comparison between submissions.

Structures where OST cannot form a valid ligand assignment (e.g. wrong atom connectivity) are penalised rather than ignored: LDDT-PLI and LDDT-LP are set to 0.0 (worst possible) and BiSyRMSD is set to 20.0 Å. This ensures a submission cannot inflate its score by omitting difficult cases.

| Metric | Direction | Description |
|--------|-----------|-------------|
| **LDDT-PLI** *(primary)* | ↑ | Fraction of reference protein–ligand contacts preserved; superposition-free and symmetry-aware |
| **BiSyRMSD** | ↓ | Binding-site symmetry-corrected RMSD of the ligand after binding-site superposition |
| **LDDT-LP** | ↑ | LDDT applied to the ligand pocket residues |
| **Coverage** | ↑ | Fraction of the 184 compounds for which a valid ligand assignment was found |

## ✅ How to Participate
1. **Register**: Create an account with Hugging Face.
2. **Get the Data**: The training and test sets are available on the [**HuggingFace dataset page**](https://huggingface.co/datasets/openadmet/pxr-challenge-train-test).
3. **Walk through the tutorial**: The [PXR Challenge Tutorial](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) contains track-specific notebooks:
   - **Activity track:** Step-by-step guide to importing pEC50 training data and training an LGBM baseline model, including how to use the single-concentration counter-screen data.
   - **Structure track:** Guide to importing the test set compounds and validating your submission, with pre-generated Boltz-2 structures as a reference starting point.

   The tutorial also includes a validation script to check your files before uploading.
4. **Join the Community**: Get support and coordinate in the `#pxr-challenge` channel on our [Discord](https://discord.gg/MY5cEFHH3D).
5. **Build and Refine**: Use the training data to build your models. In Phase 2, use the unblinded Analog Set 1 results to refine predictions for Analog Set 2.
6. **Submit**: Follow the instructions in the *Submit* tab.

## 📅 Timeline
| Date | Action |
|:--- |:--- |
| **March 17** | Challenge announced |
| **April 1** | Training/Test sets released; SMILES for structure prediction released; Submissions open |
| **May 25** | Phase 1 concludes; all Phase 1 submissions due; interim Activity leaderboard |
| **May 26** | Analog Set 1 unblinded |
| **July 1 (23:59:59 UTC)** | Phase 2 & Structure Track concludes; all remaining submissions due |

## Acknowledgements
We thank the experimentalists at **Octant** and **UCSF** (Fraser Lab) for all their hard work making this challenge possible. In particular we would like to thank Sam Sabaat, Scott Simpkins, Yuning Shen, Bryan Jiang, Henry Chan, Jeff Tang, Ayesha Ghazali, Theo Tarver, Steven Edgar, Dominic Ky and many others from Octant, and Galen Correy, Nikhil Gupta and the Fraser Lab at UCSF.
    """

    # --- Custom CSS for Tables ---
    gr.HTML("""
            <style>
            #welcome-md table {
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
            }
            #welcome-md th, #welcome-md td {
                padding: 10px;
                border: 1px solid rgba(0,0,0,0.1);
            }
            #welcome-md thead th {
                background: var(--panel-background-fill, #f5f5f7);
            }
            </style>
            """)

    with gr.Row():
        with gr.Column(scale=1):
            gr.HTML("""
                <div style="margin: 12px 0 20px 0; text-align: center;">
                    <a href="https://openadmet.ghost.io/announcing-the-next-openadmet-blind-challenge-predicting-pxr-induction"
                       target="_blank"
                       style="display: inline-block; padding: 12px 28px; background: #2563eb; color: #fff;
                              font-size: 1.05rem; font-weight: 600; border-radius: 8px; text-decoration: none;
                              box-shadow: 0 2px 8px rgba(37,99,235,0.25);">
                        📄 Challenge Announcement →
                    </a>
                </div>
            """)
        with gr.Column(scale=1):
            gr.HTML("""
                <div style="margin: 12px 0 20px 0; text-align: center;">
                    <a href="https://huggingface.co/datasets/openadmet/pxr-challenge-train-test"
                       target="_blank"
                       style="display: inline-block; padding: 12px 28px; background: #059669; color: #fff;
                              font-size: 1.05rem; font-weight: 600; border-radius: 8px; text-decoration: none;
                              box-shadow: 0 2px 8px rgba(5,150,105,0.25);">
                        🗄️ Training &amp; Test Dataset →
                    </a>
                </div>
            """)

    with gr.Tabs(elem_classes="tab-buttons"):
        with gr.TabItem("📖 About"):
            gr.Markdown(welcome_md, elem_id="welcome-md")

        with gr.TabItem("🚀 Leaderboard"):
            gr.Markdown("### Leaderboards")
            # Activity is live in phase 1; structure is live in phases 1 and 2.
            activity_leaderboard_timer = gr.Timer(value=30, active=(CURRENT_PHASE == 1))
            structure_leaderboard_timer = gr.Timer(
                value=30, active=(CURRENT_PHASE in (1, 2))
            )
            with gr.Tabs():
                with gr.TabItem("Activity Prediction (EC50)"):
                    leaderboard_headers = {
                        1: "#### Live Phase 1 Leaderboard",
                        2: "#### Interim Leaderboard",
                        0: "#### Final Leaderboard",
                    }
                    leaderboard_note = {
                        1: "The leaderboard is live during Phase 1, showing scores for the predictions on the Analog Set 1 compounds.",
                        2: "Phase 1 has concluded! The static interim leaderboard now shows scores for the predictions on all (Analog Sets 1 and 2) compounds.",
                        0: "The final leaderboard is now live, showing scores for the predictions on the Analog Set 2 compounds. Thank you to all participants!",
                        -1: "The final leaderboard will be available shortly, stay tuned!",
                    }
                    gr.Markdown(leaderboard_headers.get(CURRENT_PHASE, ""))
                    gr.Markdown(leaderboard_note.get(CURRENT_PHASE, ""))
                    if CURRENT_PHASE != -1:
                        gr.DownloadButton(
                            label="Download CSV",
                            value=download_activity_leaderboard,
                            size="sm",
                        )
                        activity_lb = Leaderboard(
                            value=load_activity_leaderboard(),
                            select_columns=[
                                "rank",
                                "username",
                                "Submitted",
                                "model_report_link",
                                "Proprietary Data",
                                "MAE",
                                "RAE",
                                "R2",
                                "Spearman ρ",
                                "Kendall's τ",
                            ],
                            search_columns=["username"],
                            filter_columns=[
                                ColumnFilter(
                                    "Proprietary Data",
                                    type="checkboxgroup",
                                    label="Proprietary Data",
                                    choices=[("Yes", "Yes"), ("No", "No")],
                                    default=[("Yes", "Yes"), ("No", "No")],
                                )
                            ],
                            datatype=[
                                "number",
                                "html",
                                "str",
                                "html",
                                "str",
                                "number",
                                "number",
                                "number",
                                "number",
                                "number",
                            ],
                        )
                with gr.TabItem("Structure Prediction (Pose)"):
                    leaderboard_headers = {
                        1: "#### Live Phase 1 Leaderboard",
                        2: "#### Live Phase 2 Leaderboard",
                        0: "#### Final Leaderboard",
                    }
                    leaderboard_note = {
                        1: "The leaderboard is live during Phase 1, showing scores for the predictions on the Analog Set 1 compounds.",
                        2: "The leaderboard is live during Phase 2, showing scores for the predictions on the Analog Set 1 compounds.",
                        0: "The final leaderboard is now live, showing scores for the predictions on the Analog Set 2 compounds. Thank you to all participants!",
                        -1: "The final leaderboard will be available shortly, stay tuned!",
                    }
                    gr.Markdown(leaderboard_headers.get(CURRENT_PHASE, ""))
                    gr.Markdown(leaderboard_note.get(CURRENT_PHASE, ""))
                    if CURRENT_PHASE != -1:
                        gr.DownloadButton(
                            label="Download CSV",
                            value=download_structure_leaderboard,
                            size="sm",
                        )
                        structure_lb = Leaderboard(
                            value=load_structure_leaderboard(),
                            select_columns=[
                                "rank",
                                "username",
                                "Submitted",
                                "model_report_link",
                                "Proprietary Data",
                                "LDDT-PLI",
                                "BiSyRMSD",
                                "LDDT-LP",
                                "Coverage",
                            ],
                            search_columns=["username"],
                            filter_columns=[
                                ColumnFilter(
                                    "Proprietary Data",
                                    type="checkboxgroup",
                                    label="Proprietary Data",
                                    choices=[("Yes", "Yes"), ("No", "No")],
                                    default=[("Yes", "Yes"), ("No", "No")],
                                )
                            ],
                            datatype=[
                                "number",
                                "html",
                                "str",
                                "html",
                                "str",
                                "number",
                                "number",
                                "number",
                                "number",
                            ],
                        )
            if CURRENT_PHASE != -1:
                activity_leaderboard_timer.tick(
                    fn=load_activity_leaderboard, outputs=[activity_lb]
                )
                structure_leaderboard_timer.tick(
                    fn=load_structure_leaderboard, outputs=[structure_lb]
                )

        with gr.TabItem("✉️ Submit"):
            if CURRENT_PHASE in (-1, 0):
                gr.Markdown("""
                    ## The challenge is now closed. Thank you to all participants! The final leaderboard will be available shortly in the **Leaderboard** tab.
                    
                    If you were unable to get your submission in before the deadline, please contact us through Discord and we *may* be able to help.
                
                    ---
                """)
            else:
                gr.Markdown("""
                    ### Submit your Predictions
                    New to the challenge? See the [PXR Challenge Tutorial](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) for a step-by-step walkthrough of data loading, model building, and preparing your submission. The tutorial repository also includes a validation script to check your files before uploading. Training and test sets are available on the [HuggingFace dataset page](https://huggingface.co/datasets/openadmet/pxr-challenge-train-test).

                    #### Activity Track
                    Submit a **`.parquet` or `.csv`** file with exactly **513 rows** and three columns:

                    | Column | Type | Description |
                    |--------|------|-------------|
                    | `SMILES` | string | SMILES string for the compound |
                    | `Molecule Name` | string | Compound identifier (e.g. `OADMET-00000`) |
                    | `pEC50` | float | Predicted pEC50 — no `NaN` or `inf` |

                    **Example (CSV):**
                    ```
                    SMILES,Molecule Name,pEC50
                    CCO,OADMET-00000,6.23
                    c1ccccc1,OADMET-00001,5.87
                    ...
                    ```

                    #### Structure Track
                    Submit a **`.zip`** archive containing exactly **184 `.pdb` files**, one per compound, named after the compound identifier (e.g. `x00011-1.pdb`). Each file must be a full protein–ligand complex with the ligand residue named **`LIG`**.
                """)

                with gr.Row():
                    # --- Column 1: Leaderboard identity ---
                    with gr.Column(scale=1):
                        gr.Markdown("### Leaderboard Identity")
                        username_input = gr.Textbox(
                            label="HuggingFace Username *",
                            placeholder="your-hf-username",
                            info="Required. Used to track submissions across phases.",
                        )
                        user_alias = gr.Textbox(
                            label="Alias",
                            placeholder="Optional display name for the leaderboard",
                        )
                        anon_checkbox = gr.Checkbox(
                            label="Submit anonymously (show alias instead of username on "
                            "leaderboard and Discord submission validation bot)",
                            value=False,
                        )

                    # --- Column 2: Contact & publication details ---
                    with gr.Column(scale=1):
                        gr.Markdown("### Contact & Publication *(optional, private)*")
                        participant_name = gr.Textbox(
                            label="Full Name",
                            placeholder="Jane Smith",
                            info="Not displayed publicly; used for publication tracking.",
                        )
                        discord_username = gr.Textbox(
                            label="Discord Username",
                            placeholder="janesmit#1234",
                        )
                        email = gr.Textbox(
                            label="Email",
                            placeholder="jane@example.com",
                        )
                        affiliation = gr.Textbox(
                            label="Affiliation",
                            placeholder="University / Company",
                        )
                        model_tag = gr.Textbox(
                            label="Method Report Link",
                            placeholder="https://...",
                            info="Required before the deadline to appear on the final leaderboard.",
                        )
                        paper_checkbox = gr.Checkbox(
                            label="Include me in a future Challenge publication",
                            value=False,
                        )
                        proprietary_data_checkbox = gr.Checkbox(
                            label="I used proprietary data (not publicly available) in training my model",
                            value=False,
                            info="Displayed publicly on the leaderboard.",
                        )

                    # --- Column 3: Track & file ---
                    with gr.Column(scale=1):
                        gr.Markdown("### Submission")
                        track_select = gr.Radio(
                            ["Activity Prediction", "Structure Prediction"],
                            label="Track *",
                        )
                        file_input = gr.File(label="Upload File *")

                    # --- Submit row ---
                    with gr.Row():
                        with gr.Column(scale=1):
                            pass
                        with gr.Column(scale=2):
                            submit_btn = gr.Button(
                                "Submit predictions", variant="primary", size="lg"
                            )
                            submit_msg = gr.Textbox(
                                label="Submission Status",
                                lines=2,
                                visible=False,
                                interactive=False,
                            )
                        with gr.Column(scale=1):
                            pass

                    submit_btn.click(
                        fn=submit_predictions,
                        inputs=[
                            username_input,
                            user_alias,
                            anon_checkbox,
                            participant_name,
                            discord_username,
                            email,
                            affiliation,
                            model_tag,
                            paper_checkbox,
                            proprietary_data_checkbox,
                            track_select,
                            file_input,
                        ],
                        outputs=[submit_msg],
                    )

        with gr.TabItem("🛠️ FAQ"):
            gr.Markdown("""
1. **What is the difference between Phase 1 and Phase 2?**
> Phase 1 is a blind test of all 513 compounds, with Analog Set 1 scored on a **live leaderboard**. In Phase 2, Analog Set 1 results are unblinded so you can incorporate them into your training pipeline and refine predictions for **Analog Set 2**, which is scored fully blind at the final deadline.

2. **Can I compete in one or both tracks?**
> Yes — you can compete in either or both tracks independently.

3. **I want to participate with colleagues. Should we submit together or separately?**
> Please submit under a single HuggingFace account or alias. Do not submit the same predictions from multiple accounts.

4. **Can I use AlphaFold or other protein prediction tools?**
> Yes! For the Structure Track you are encouraged to use any computational tools, including protein structure prediction and docking algorithms. Existing PDB structures may also be used. You may also use structural data to inform activity predictions, and vice versa.

5. **Is there a limit on submissions?**
> During the first week, submissions are capped at **once every 4 hours** to allow us to catch any unexpected scoring issues. After the first week this will likely move to once per day. Only your **latest** submission counts towards the leaderboard. Please validate your files using the [validation script in the tutorial repo](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) before uploading.

6. **When do final entries have to be submitted?**
> All final submissions must be received by **11:59 PM UTC on July 1**.

7. **Where do I ask technical questions?**
> For computational and infrastructure questions, use the `#pxr-challenge` channel on the [OpenADMET Discord](https://discord.gg/MY5cEFHH3D). For questions directed at the experimentalists (assay design, data interpretation, crystallography), post in Discord and the team will respond there.

8. **What file formats are required?**
> - **Activity Track:** `.parquet` or `.csv` with exactly 513 rows and columns `SMILES`, `Molecule Name`, `pEC50`.
> - **Structure Track:** `.zip` archive containing exactly 184 `.pdb` files with the ligand residue named `LIG`.
> See the full format details in the *Submit* tab and the [tutorial repo](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main).

9. **What should be in my submitted PDB file?**
> Each file must contain the **PXR monomer** (not the homodimer) plus your predicted ligand pose. The ligand must use residue name `LIG` exactly. Bond orders must match the SMILES provided in the challenge set — the scoring backend uses graph isomorphism, so any mismatch will result in a scoring failure. Use the [validation script](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) to check your files before uploading.

10. **What is the PXR sequence used in the challenge?**
> A FASTA file containing the PXR construct sequence is provided in the [tutorial repository](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main).

11. **How are submissions evaluated?**
> **Activity Track:** The primary metric is **RAE** (Relative Absolute Error) on pEC50. Secondary metrics — MAE, R², Spearman ρ, and Kendall's τ — are also reported. All metrics are bootstrapped over 1000 resamples.
>
> **Structure Track:** Scoring uses the **OpenStructure** pipeline (CASP15 reference implementation). Each predicted complex is scored pairwise against all reference ligands (NxM matrix), and the best-scoring pair is selected. The primary metric is **LDDT-PLI** (superposition-free). The secondary metric is **BiSyRMSD** (binding-site symmetry-corrected RMSD). Submit the **PXR monomer**, not the homodimer. Structures where no graph isomorphism match can be formed for the ligand receive a penalty score (LDDT-PLI = 0.0, BiSyRMSD = 20.0 Å). All metrics are bootstrapped over 1000 resamples.

12. **Where can I find my submission status?**
> Submission receipts and scoring updates are posted to the `#pxr-challenge-submissions` channel on the [OpenADMET Discord](https://discord.gg/MY5cEFHH3D). You will receive feedback for every entry, including error details if your submission failed.

13. **My entry hasn't appeared on the leaderboard. What should I do?**
> First, re-run the [validation script](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) to confirm your file meets all formatting requirements. Then check `#pxr-challenge-submissions` on [Discord](https://discord.gg/MY5cEFHH3D) for error messages. If nothing has appeared after 2 hours, post in `#pxr-challenge` for further assistance.

14. **My submission to the Structure Track was only partially scored. What should I do?**
> Partial scoring usually means one or more PDB files have bond orders that do not match the SMILES in the challenge set. Check each affected file and ensure the ligand connectivity exactly matches the corresponding SMILES. The [validation script](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main) can help identify mismatches.

15. **Can I use external data to train my models?**
> Yes, external data is allowed. You may also cross-use structural and activity data across tracks.

16. **What are the rules around proprietary data?**
> Proprietary data is allowed. When submitting, tick the "I used proprietary data" checkbox to disclose this. The flag is displayed publicly on the leaderboard and can be filtered by other users.

17. **Is a methodology report required?**
> Yes. A methodology report (code preferred; a written report is the minimum) must be submitted by the challenge close date (July 1) for your entry to appear on the final leaderboard. Submit the link via the Method Report Link field in the Submit tab.

18. **Where can I find the validation scripts?**
> Validation scripts for both tracks are in the [PXR Challenge Tutorial repository](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main).

19. **Will there be a summary paper?**
> Yes, we plan to publish a summary paper covering the challenge results. This may be combined with the Expansion challenge paper. Participants who opt in via the submission form will be considered for inclusion.

---

## Submission Format

### Activity Track

Submit a **`.parquet` or `.csv`** file with exactly **513 rows** and the following three columns:

| Column | Type | Description |
|--------|------|-------------|
| `SMILES` | string | SMILES string for the compound |
| `Molecule Name` | string | Compound identifier (e.g. `OADMET-00000`) |
| `pEC50` | float | Predicted pEC50 value — must not contain NaN or infinite values |

**Rules:**
- The file must contain predictions for all 513 compounds — no more, no fewer.
- Column names are case-sensitive and must match exactly.
- `pEC50` values must be finite floats (no `NaN`, `inf`, or `-inf`).
- `.parquet` is preferred; `.csv` is also accepted.

**Example (CSV):**
```
SMILES,Molecule Name,pEC50
CCO,OADMET-00000,6.23
c1ccccc1,OADMET-00001,5.87
...
```

---

### Structure Track

Submit a **`.zip`** archive containing exactly **184 `.pdb` files**, one per compound.

**File naming:** Each PDB file must be named after its compound identifier, e.g.:
```
structures.zip
├── x00011-1.pdb
├── x00011-2.pdb
├── ...
└── x00999-3.pdb
```

**PDB file requirements:**
- Each file must be a **full protein–ligand complex** — the complete receptor structure with the docked ligand included.
- The ligand must use the residue name **`LIG`** exactly. This is how the evaluator identifies the ligand within the structure. Files where no `LIG` residue can be found will receive worst-case penalty scores (LDDT-PLI = 0, BiSyRMSD = 20 Å).
- Standard PDB format is expected. Hydrogens may be included or omitted.
- The zip must contain exactly 184 `.pdb` files — no subdirectories, no extra files.

**Common mistakes to avoid:**
- Using a residue name other than `LIG` for the ligand (e.g. `UNK`, `MOL`, the compound name).
- Submitting ligand-only PDB files instead of full protein–ligand complexes.
- Including extra files or nested folders inside the zip.
- Mismatched filenames (e.g. wrong zero-padding or wrong prefix).
            """)
        if CURRENT_PHASE == 0:
            with gr.TabItem("🥊 Head-to-head Comparisons"):
                render_comparisons()


if __name__ == "__main__":
    logger.info("Starting PXR Challenge Gradio app...")
    demo.launch()
