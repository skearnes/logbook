import gradio as gr
import pandas as pd
from pathlib import Path
from typing import Optional
from about import (
    ENDPOINTS, API, METRICS,
    submissions_repo, 
    results_repo_test, 
    results_repo_validation, 
    test_repo,
    THROTTLE_MINUTES
)
from utils import (
    bootstrap_metrics, 
    clip_and_log_transform, 
    fetch_dataset_df,
    check_page_exists,
)
from huggingface_hub import hf_hub_download
import datetime
import io
import json, tempfile
import re
from pydantic import (
    BaseModel, 
    Field, 
    model_validator, 
    field_validator, 
    ValidationError
)
from loguru import logger

HF_USERNAME_RE = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-_]{1,38})$")
def _safeify_username(username: str) -> str:
    return str(username.strip()).replace("/", "_").replace(" ", "_")

def _unsafify_username(username: str) -> str:
    return str(username.strip()).replace("/", "_").replace(" ", "_")

def _check_required_columns(df: pd.DataFrame, name: str, cols: list[str]):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")

class ParticipantRecord(BaseModel):
    hf_username: str = Field(description="Hugging Face username")
    display_name: Optional[str] = Field(description="Name to display on leaderboard")
    participant_name: Optional[str] = Field(default=None, description="Participant's real name")
    discord_username: Optional[str] = Field(default=None, description="Discord username")
    email: Optional[str] = Field(default=None, description="Email address")
    affiliation: Optional[str] = Field(default=None, description="Affiliation")
    model_tag: Optional[str] = Field(default=None, description="Link to model description")
    anonymous: bool = Field(default=False, description="Whether to display username as 'anonymous'")
    consent_publication: bool = Field(default=False, description="Consent to be included in publications")

    @field_validator("hf_username")
    @classmethod
    def validate_hf_username(cls, v: str) -> str:
        """Validate hugging face username exists by checking the url

        Parameters
        ----------
        username : str
            Hugging Face username

        Returns
        -------
        bool
            Whether Hugging Face username is real
        """
        v = v.strip()
        if not HF_USERNAME_RE.match(v):
            raise gr.Error("Invalid Hugging Face username (letters, numbers, -, _; min 2, max ~39).")
    
        hf_url = f"https://huggingface.co/{v}"
        if not check_page_exists(hf_url, delay=2):
            raise gr.Error("The Hugging Face username used doesn't exist. An account is required to participate.")
        return v

    @field_validator("display_name")
    @classmethod
    def validate_display_name(cls, v: Optional[str]) -> Optional[str]:
        if v is None:
            return None
        v = v.strip()
        if not v:
            return None
        if len(v) > 20:
            raise ValueError("Display name is too long (max 20 chars).")
        return v
    
    @field_validator("model_tag", mode="before")
    @classmethod
    def normalize_url(cls, v):
        if v is None:
            return v
        s = str(v).strip()
        if not s:
            return None
        if "://" not in s:
            s = "https://" + s
        return s

    @model_validator(mode="after")
    def require_display_name_if_anonymous(self) -> "ParticipantRecord":
        if self.anonymous and not self.display_name:
            raise ValueError("Alias is required when anonymous box is checked.")
        return self

class SubmissionMetadata(BaseModel):
    submission_time_utc: str
    user: str
    original_filename: str
    evaluated: bool
    participant: ParticipantRecord


def submit_data(predictions_file: str, 
                user_state,
                participant_name: str = "",
                discord_username: str = "",
                email: str = "",
                affiliation: str = "",
                model_tag: str = "",
                user_display: str = "",
                anon_checkbox: bool = False,
                paper_checkbox: bool = False
):
    
    if user_state is None:
        raise gr.Error("Username or alias is required for submission.")
    

    # check the last time the user submitted
    data, __ = fetch_dataset_df()
    if not data[data['user'] == user_state].empty:
        last_time = data[data['user'] == user_state]['submission time'].max()
        delta = datetime.datetime.now(datetime.timezone.utc) - last_time.to_pydatetime()
        if delta < datetime.timedelta(minutes=THROTTLE_MINUTES):
            raise gr.Error(f"You have submitted within the last {THROTTLE_MINUTES} minutes. Please wait {THROTTLE_MINUTES - int(delta.total_seconds() // 60)} minutes before submitting again.")


    file_path = Path(predictions_file).resolve()
    if not file_path.exists():
        raise gr.Error("Uploaded file object does not have a valid file path.")

    # Read results file 
    try:
        results_df = pd.read_csv(file_path)            
    except Exception as e:
        raise gr.Error(f"❌ Error reading results file: {str(e)}")

    if results_df.empty:
        raise gr.Error("The uploaded file is empty.")
    
    missing = set(ENDPOINTS) - set(results_df.columns)
    if missing:
        raise gr.Error(f"The uploaded file must contain all endpoint predictions {ENDPOINTS} as columns, missing: {missing}")

    # Save participant record
    try:
        participant_record = ParticipantRecord(
            hf_username=user_state,
            participant_name=participant_name,
            discord_username=discord_username,
            email=email,
            affiliation=affiliation,
            model_tag=model_tag,
            display_name=user_display,
            anonymous=anon_checkbox,
            consent_publication=paper_checkbox
        )
    except ValidationError as e:
        raise gr.Error(f"❌ Error in participant information: {str(e)}")

    # Build destination filename in the dataset
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds") # should keep default time so can be deserialized correctly
    try:
        meta = SubmissionMetadata(
            submission_time_utc=ts,
            user=user_state,
            original_filename=file_path.name,
            evaluated=False,
            participant=participant_record
        )
    except ValidationError as e:
        raise gr.Error(f"❌ Error in metadata information: {str(e)}")

    safe_user = _safeify_username(user_state)
    destination_csv = f"submissions/{safe_user}_{ts}.csv"
    destination_json = destination_csv.replace(".csv", ".json")

    # Upload the CSV file
    API.upload_file(
        path_or_fileobj=str(file_path),
        path_in_repo=destination_csv,
        repo_id=submissions_repo,
        repo_type="dataset",
        commit_message=f"Add submission for user at {ts}"
    )
    # Upload the metadata JSON file
    meta_bytes = io.BytesIO(json.dumps(meta.model_dump(), indent=2).encode("utf-8"))
    API.upload_file(
        path_or_fileobj=meta_bytes,
        path_in_repo=destination_json,
        repo_id=submissions_repo,
        repo_type="dataset",
        commit_message=f"Add metadata for user submission at {ts}"
    )

    return "✅ Your submission has been received! Your scores will appear on the leaderboard shortly.", destination_csv

def evaluate_data(filename: str) -> None:
    # do test set first as a more stringent check of the submission w.r.t matching molecules
    logger.info(f"Evaluating submission file {filename}")
    # evaluate on the test set
    _evaluate_data(filename, test_repo=test_repo, split_filename="data/expansion_data_test.csv", results_repo=results_repo_test)
    # evaluate on the validation set
    _evaluate_data(filename, test_repo=test_repo, split_filename="data/expansion_data_test_validation.csv", results_repo=results_repo_validation)
    logger.info(f"Finished evaluating submission file {filename}")

def _evaluate_data(filename: str, test_repo: str, split_filename: str, results_repo: str) -> None:

    # Load the submission csv
    try:
        local_path = hf_hub_download(
            repo_id=submissions_repo,
            repo_type="dataset",
            filename=filename,
        )
    except Exception as e:
        raise gr.Error(f"Failed to download submission file: {e}")
    
    # Load the test set
    try: 
        test_path = hf_hub_download(
            repo_id=test_repo,
            repo_type="dataset",
            filename=split_filename
        )
    except Exception as e:
        raise gr.Error(f"Failed to download test file: {e}")
    
    data_df = pd.read_csv(local_path)
    test_df = pd.read_csv(test_path)
    try:
        results_df, results_raw_df = calculate_metrics(data_df, test_df)
        if not isinstance(results_df, pd.DataFrame) or results_df.empty:
            raise gr.Error("Evaluation produced no results.")
    except Exception as e:
        raise gr.Error(f'Evaluation failed: {e}. No results written to results dataset.')
    
    # Load metadata file
    meta_filename = filename.replace(".csv", ".json")
    try:
        meta_path = hf_hub_download(
                repo_id=submissions_repo,
                repo_type="dataset",
                filename=meta_filename,
            )
        with open(meta_path, "r", encoding="utf-8") as f:
            _meta = json.load(f)
        meta = SubmissionMetadata(**_meta)
        username = meta.participant.hf_username
        timestamp = meta.submission_time_utc
        report = meta.participant.model_tag
        if meta.participant.anonymous:
            display_name = meta.participant.display_name
        else:
            display_name = username
    except Exception as e:
        raise gr.Error(f"Failed to load metadata file: {e}. No results written to results dataset.")

    # Write results to results dataset
    results_df['user'] = display_name
    results_df['submission_time'] = timestamp
    results_df['model_report'] = report
    results_df['anonymous'] = meta.participant.anonymous
    results_df['hf_username'] = username

    results_raw_df = results_raw_df[results_raw_df['Endpoint']=='Average'] # Save ONLY for average endpoint, otherwise file is too large
    results_raw_df['user'] = display_name
    results_raw_df['submission_time'] = timestamp
    results_raw_df['model_report'] = report
    results_raw_df['anonymous'] = meta.participant.anonymous
    results_raw_df['hf_username'] = username   

    safe_user = _unsafify_username(username)
    destination_path = f"results/{safe_user}_{timestamp}_results.csv"
    tmp_name = None
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as tmp:
        results_df.to_csv(tmp, index=False)
        tmp.flush()
        tmp_name = tmp.name
    
    API.upload_file(
            path_or_fileobj=tmp_name, 
            path_in_repo=destination_path,
            repo_id=results_repo,
            repo_type="dataset",
            commit_message=f"Add result data for {username}"
        )
    Path(tmp_name).unlink()

    # Same for raw file
    destination_path_raw = f"results/{safe_user}_{timestamp}_results_raw.csv"
    tmp_name = None
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as tmp:
        results_raw_df.to_csv(tmp, index=False)
        tmp.flush()
        tmp_name = tmp.name
    
    API.upload_file(
            path_or_fileobj=tmp_name, 
            path_in_repo=destination_path_raw,
            repo_id=results_repo,
            repo_type="dataset",
            commit_message=f"Add raw result data for {username}"
        )
    Path(tmp_name).unlink()

def calculate_metrics(
        results_dataframe: pd.DataFrame,
        test_dataframe: pd.DataFrame
    ):
    import numpy as np
    
    # Do some checks

    # 1) Check all columns are present
    _check_required_columns(results_dataframe, "Results file", ["Molecule Name"] + ENDPOINTS)
    _check_required_columns(test_dataframe, "Test file", ["Molecule Name"] + ENDPOINTS)
    
    
       # 2) Check all Molecules in the test set are present in the predictions
    if not (test_dataframe['Molecule Name']).isin(results_dataframe['Molecule Name']).all():
        raise gr.Error("Some molecules in the test set are missing from the predictions file. Please ensure all molecules are included.")


    # 3) check no duplicated molecules in the predictions file
    if results_dataframe['Molecule Name'].duplicated().any():
        raise gr.Error("The predictions file contains duplicated molecules. Please ensure each molecule is only listed once.")
    
    # 4) Merge dataframes to ensure alignment
    merged_df = results_dataframe.merge(
        test_dataframe,
        on="Molecule Name",
        suffixes=('_pred', '_true'),
        how="inner"
    )
    merged_df = merged_df.sort_values("Molecule Name")

    # 5) loop over endpoints 

    final_cols = ["MAE", "RAE", "R2", "Spearman R", "Kendall's Tau"]
    all_endpoint_results = []
    all_endpoint_results_raw = []

    for ept in ENDPOINTS:
        pred_col = f"{ept}_pred"
        true_col = f"{ept}_true"

        # cast to numeric, coerce errors to NaN
        merged_df[pred_col] = pd.to_numeric(merged_df[pred_col], errors="coerce")
        merged_df[true_col] = pd.to_numeric(merged_df[true_col], errors="coerce")

        if merged_df[pred_col].isnull().all():
            raise gr.Error(f"All predictions are missing for endpoint {ept}. Please provide valid predictions.")
        
        # subset and drop NaNs
        subset = merged_df[[pred_col, true_col]].dropna()
        if subset.empty:
            raise gr.Error(f"No valid data available for endpoint {ept} after removing NaNs.")
        
        # extract numpy arrays
        y_pred = subset[pred_col].to_numpy()
        y_true = subset[true_col].to_numpy()

        # apply log10 + 1 transform except for logD
        if ept.lower() not in ['logd']:
            y_true_log = clip_and_log_transform(y_true)
            y_pred_log = clip_and_log_transform(y_pred)

        else:
            y_true_log = y_true
            y_pred_log = y_pred

        # calculate metrics with bootstrapping
        bootstrap_df = bootstrap_metrics(y_pred_log, y_true_log, ept, n_bootstrap_samples=1000)
        # Longer pivot alternative for the cases where all metric results are NaN, as pivot ignores those columns 
        grouped = bootstrap_df.groupby(["Endpoint", "Metric"])["Value"].agg(["mean", "std"])
        df_unstacked = grouped.unstack(level="Metric")
        df_reindexed = df_unstacked.reindex(columns=list(METRICS), level=1)

        df_reindexed.columns = [f"{agg}_{metric}" for agg, metric in df_reindexed.columns]
        df_endpoint = df_reindexed.reset_index()
        all_endpoint_results.append(df_endpoint)

        # Also save a raw dataframe with all the bootstrapping samples
        df_endpoint_raw = bootstrap_df.pivot_table(
            index=["Sample", "Endpoint"],  
            columns="Metric",
            values="Value"
        ).reset_index()
        df_endpoint_raw.columns.name = None
        df_endpoint_raw['Sample'] = df_endpoint_raw['Sample'].astype(int)
        all_endpoint_results_raw.append(df_endpoint_raw)

    df_results = pd.concat(all_endpoint_results, ignore_index=True)
    df_results_raw = pd.concat(all_endpoint_results_raw, ignore_index=True)

    # Average results
    mean_cols = [f'mean_{m}' for m in final_cols]
    std_cols = [f'std_{m}' for m in final_cols]
    macro_means = df_results[mean_cols].mean()
    macro_stds = df_results[std_cols].mean()
    avg_row = {"Endpoint": "Average"}
    avg_row.update(macro_means.to_dict())
    avg_row.update(macro_stds.to_dict())    
    df_with_average = pd.concat([df_results, pd.DataFrame([avg_row])], ignore_index=True)
    # Fix order of columns
    df_with_average = df_with_average[["Endpoint"]+mean_cols+std_cols]

    # Average results for raw dataframe
    macro_results_by_sample = df_results_raw.groupby('Sample')[final_cols].mean().reset_index()
    macro_results_by_sample["Endpoint"] = "Average"
    df_raw_with_average = pd.concat([df_results_raw, macro_results_by_sample], ignore_index=True)
    df_raw_with_average = df_raw_with_average[["Sample","Endpoint"] + final_cols]
    
    return df_with_average, df_raw_with_average