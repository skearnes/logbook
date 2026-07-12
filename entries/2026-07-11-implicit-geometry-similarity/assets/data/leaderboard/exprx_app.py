import gradio as gr
from gradio_leaderboard import Leaderboard
from gradio.themes.utils import sizes
import pandas as pd
import numpy as np
import os

from evaluate import submit_data, evaluate_data
from utils import (
    make_tag_clickable, 
    make_user_clickable, 
    fetch_dataset_df,
    map_metric_to_stats,
)
from cld import add_cld_to_leaderboard
from datasets import load_dataset
import tempfile
from loguru import logger
from about import ENDPOINTS, LB_COLS, LB_AVG, LB_DTYPES
import time
import threading



ALL_EPS = ['Average'] + ENDPOINTS
final_cols = ["rank", "user", "CLD", "MA-RAE", "R2", "Spearman R", "Kendall's Tau", "model details"]
ep_cols = ["rank", "user", "MAE", "R2", "Spearman R", "Kendall's Tau", "model details"]

def build_leaderboard(df_results, df_results_raw):
    logger.info("Rebuilding leaderboard data...")
    per_ep = {}
    for ep in ALL_EPS:
        df = df_results[df_results["Endpoint"] == ep].copy()
        if df is None:
            print(f"[refresh] {ep} returned None; using empty DF")
        if df.empty:
            per_ep[ep] = pd.DataFrame(columns=LB_COLS) # Empty df
            continue

        # Make model details clickable if it's a huggingface user
        df['model details'] = df['model_report'].apply(lambda x: make_tag_clickable(x)).astype(str)

        if ep == "Average":
            # MA-RAE is the average of the RAE per endpoint
            df = df.rename(columns={"mean_RAE": "mean_MA-RAE", 
                                    "std_RAE": "std_MA-RAE"}) 
            sorted_df = df.sort_values(by='mean_MA-RAE', ascending=True, kind="stable")
            sorted_df = map_metric_to_stats(sorted_df, average=True)
            # Add ranking column
            sorted_df['rank'] = np.arange(1, len(sorted_df) + 1)
            avg_leaderboard = sorted_df.copy()
            avg_cols = LB_AVG
            # Add CLD
            if df_results_raw is not None:
                df_raw = df_results_raw[df_results_raw["Endpoint"] == ep].copy()
                df_raw = df_raw.rename(columns={"RAE": "MA-RAE"})
                avg_leaderboard = add_cld_to_leaderboard(
                    sorted_df,
                    df_raw,
                    "MA-RAE",
                )
                avg_cols = ["rank", "user", "CLD", "MA-RAE", "R2", "Spearman R", "Kendall's Tau", "submission time", "model details"]

            # Make user and model details clickable if it's a huggingface user
            avg_leaderboard['user'] = avg_leaderboard.apply(
            lambda row: make_user_clickable(row['user']) if not row['anonymous'] else row['user'], 
            axis=1).astype(str)
            per_ep[ep] = avg_leaderboard[avg_cols]
        else:
            sorted_df = df.sort_values(by="mean_MAE", ascending=True, kind="stable")
            sorted_df = map_metric_to_stats(sorted_df)
            sorted_df['user'] = sorted_df.apply(
                lambda row: make_user_clickable(row['user']) if not row['anonymous'] else row['user'], 
                axis=1).astype(str)
            per_ep[ep] = sorted_df[LB_COLS]
    logger.info("Finished rebuilding leaderboard data.")
    return per_ep

def build_leaderboard_static():
    """Auxiliary function to add static leaderboars:
    - Last intermeadiate leaderboard 
    - Final leaderboard for each endpoint
    - Final leaderboard to be on a separate tab
    """
    # Load csv files of avg leaderboards
    inter_leaderboard = pd.read_csv("all_leaderboards/current_leaderboard.csv")
    final_leaderboard = pd.read_csv("all_leaderboards/leaderboard_cld_results.csv")

    per_ep = {}
    inter_leaderboard['model details'] = inter_leaderboard['model details'].apply(lambda x: make_tag_clickable(x)).astype(str)
    per_ep["Average"] = inter_leaderboard[LB_AVG]
    df = final_leaderboard.rename(columns={"final rank": "rank",
                                            "RAE_display": "MA-RAE",
                                            "R2_display": "R2",
                                            "Spearman R_display":"Spearman R",
                                            "Kendall's Tau_display": "Kendall's Tau",
                                            })
    # df['model details'] = df['model details'].apply(lambda x: make_tag_clickable(x)).astype(str)
    per_ep["Final"] = df[final_cols]

    for ep in ENDPOINTS:
        endpoint_safe = ep.replace(' ', '_').replace(">", '')
        df = pd.read_csv(f"all_leaderboards/leaderboard_{endpoint_safe}_results.csv")
        df = df.rename(columns={"final rank": "rank",
                                "MAE_display": "MAE",
                                "R2_display": "R2",
                                "Spearman R_display":"Spearman R",
                                "Kendall's Tau_display": "Kendall's Tau",
                                })
        per_ep[ep] = df[ep_cols]

    return per_ep

def get_leaderboard_csv(df_results: pd.DataFrame) -> str:
    """Auxiliary function to save leaderboard to temp csv file
    """
    # Generate csv-friendly version of leaderboard
    df = df_results[df_results["Endpoint"] == "Average"].copy()

    if df.empty:
        current_lb = pd.DataFrame(columns=LB_COLS) # Empty df
    df = df.rename(columns={"mean_RAE": "mean_MA-RAE", 
                            "std_RAE": "std_MA-RAE",
                            "model_report": "model details"}) 
    sorted_df = df.sort_values(by='mean_MA-RAE', ascending=True, kind="stable")
    sorted_df = map_metric_to_stats(sorted_df, average=True)
    # Add ranking column
    sorted_df['rank'] = np.arange(1, len(sorted_df) + 1)
    avg_leaderboard = sorted_df.copy()
    current_lb = avg_leaderboard[LB_AVG]

    # Temporary file
    temp_dir = tempfile.gettempdir()
    temp_filepath = os.path.join(temp_dir, "current_leaderboard.csv")
    current_lb.to_csv(temp_filepath, index=False)
    return temp_filepath
    
def get_leaderboard_csv_aux(df_lb: pd.DataFrame, lb_type: str) -> str:
    """Aux function for pre-formatted LBs
    """
    # Temporary file
    temp_dir = tempfile.gettempdir()
    temp_filepath = os.path.join(temp_dir, f"{lb_type}_leaderboard.csv")
    df_lb.to_csv(temp_filepath, index=False)
    return temp_filepath
# Initialize global dataframe
# current_df, current_df_raw = fetch_dataset_df()

# # Initialize global counter
# data_version_counter = 0

def update_current_dataframe():
    global current_df # ugly but works
    while True:
        logger.info("Fetching latest dataset for leaderboard...")
        current_df, current_df_raw = fetch_dataset_df()
        logger.debug(f"Dataset version updated")
        time.sleep(60)  # Check for updates every 60 sec

# threading.Thread(target=update_current_dataframe, daemon=True).start()




with gr.Blocks(title="OpenADMET ADMET Challenge", fill_height=False,
               theme=gr.themes.Default(text_size=sizes.text_lg)) as demo:
    # timer = gr.Timer(30)  # Run every 30 seconds
    data_version = gr.State(value=0)
    def increment_data_version(current_version):
        logger.debug("Incrementing data version counter... to " + str(current_version + 1))
        return current_version + 1
    
    # timer.tick(fn=increment_data_version, inputs=[data_version], outputs=data_version)
    
    ### Header
    with gr.Row():
        with gr.Column(scale=7):  # bigger text area
            gr.Markdown("""
                ## Welcome to the OpenADMET + ExpansionRx Blind Challenge!
                Your task is to develop and submit predictive models for key ADMET properties on a blinded test set of real world drug discovery data 🧑‍🔬
                        
                Go to the **Leaderboard** to check out how the challenge is going. 
                To participate, head out to the **Submit** tab and upload your results as a `CSV` file.
                    
                We released an intermediate leaderboard on December 2nd, with submissions evaluated against the full blinded test set.
                    Check it out on the OpenADMET [blog](https://openadmet.ghost.io/openadmet-expansionrx-blind-challenge-2/)!
                        
                **This challenge closed on January 19, 2026**. The final leaderboard can be found in the **FINAL** tab.
                        
                        
                A new challenge will be announced soon. Stay tuned!
                """
                )
        with gr.Column(scale=2):  # smaller side column for logo
            gr.Image(
                value="./_static/challenge_logo.png",
                show_label=False,
                show_download_button=False,
                width="5vw",  # Take up the width of the column (2/8 = 1/4)
            )
    # --- Welcome markdown message ---
    welcome_md = """
    # 💊 OpenADMET + ExpansionRx  
    ## Computational Blind Challenge in ADMET
    This challenge is a community-driven initiative to benchmark predictive models for ADMET properties in drug discovery,
    hosted by **OpenADMET** in collaboration with **ExpansionRx**. 
    ## Why are ADMET properties important in drug discovery?
    Small molecules continue to be the bricks and mortar of drug discovery globally, accounting for ~75% of FDA approvals over the last decade. 
    Oral bioavailability, easily tunable properties, modulation of a wide range of mechanisms, and ease of manufacturing make small molecules highly attractive as therapeutic agents. 
    Moreover, emerging small molecule modalities such as degraders, expression modulators, molecular glues, and antibody-drug conjugates (to name a few) have vastly expanded what we thought small molecules were capable of. 
    It is fairly difficult to predict the lifetime and distribution of small molecules within the body. Additionally, 
    interaction with off-targets can cause safety issues and toxicity. Collectively these *Absorption*, *Distribution*, *Metabolism*, *Excretion*, *Toxicology*--or **ADMET**--properties 
    sit in the middle of the assay cascade and can make or break preclinical candidate molecules. 
    **OpenADMET** aims to address these challenges through an open science effort to build predictive models of ADMET properties by characterizing the proteins and mechanisms 
    that give rise to these properties through integrated structural biology, high throughput experimentation and integrative computational models. 
    Read more about our strategy to transform drug discovery on our [website](https://openadmet.ghost.io/what-is-openadmet/). 
    Critical to our mission is developing open datasets and running community blind challenges to assess the current state of the art in ADMET modeling. 
    Building on the sucess of the recent [ASAP-Polaris-OpenADMET blind challenge](https://chemrxiv.org/engage/chemrxiv/article-details/68ac00d1728bf9025e22fe45) in computational methods for drug discovery,
    we bring you a brand new challenge in collaboration with **ExpansionRx**. During a recent series of drug discovery campaigns for RNA mediated diseases, 
    ExpansionRX collected a variety of ADMET data for off-targets and properties of interest, which they are generously sharing with the community for this challenge.

    ## 🧪 The Challenge
    Participants will be tasked with solving real-world ADMET prediction problems ExpansionRx faced during lead optimization. 
    Specifically, you will be asked to predict the ADMET properties of late-stage molecules based on earlier-stage data from the same campaigns.
    For this challenge we selected nine (9) crucial endpoints for the community to predict:
    - LogD
    - Kinetic Solubility **KSOL**: uM
    - Mouse Liver Microsomal (**MLM**) *CLint*: mL/min/kg
    - Human Liver Microsomal (**HLM**) *Clint*: mL/min/kg
    - Caco-2 Efflux Ratio
    - Caco-2 Papp A>B (10^-6 cm/s)
    - Mouse Plasma Protein Binding (**MPPB**): % Unbound
    - Mouse Brain Protein Binding (**MBPB**): % Unbound
    - Mouse Gastrocnemius Muscle Binding (**MGMB**): % Unbound
    
    Find more information about these endpoints on our [blog](https://openadmet.ghost.io/openadmet-expansionrx-blind-challenge/).
    
    **UPDATE:** The Challenge is now live! Data available at the following Hugging Face Datasets

    - Training: https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-train-data
    - Test: https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-test-data-blinded
    
    You can also watch a [Webinar](https://www.youtube.com/watch?v=9v0Ej_FL6k0) where we introduce the challenge, hosted by [Collaborative Drug Discovery (CDD)](https://www.collaborativedrug.com/).

    **UPDATE:** [OpenEye Cadence Molecular Sciences](https://www.eyesopen.com/) is generously providing access to their Toolkit for interested participants during the duration of the challenge. 
    Request access by filling out this [Google Form](https://forms.gle/2maSE1ne213TDmz87).
    
    ## ✅ How to Participate
    1. **Register**: Create an account with Hugging Face.
    2. **Walk through the tutorials**: We have prepared a [Tutorial](https://github.com/OpenADMET/ExpansionRx-Challenge-Tutorial/blob/main/expansion_tutorial.ipynb) showing how to train a model and submit to the leaderboard.
    3. **Download the Public Dataset**: Download the ExpansionRx [training](https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-train-data) and [blinded test](https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-test-data-blinded) sets from Hugging Face.
    4. **Train Your Model**: Use the provided training data for each ADMET property of your choice.
    5. **Submit Predictions**: Follow the instructions in the *Submit* tab to upload your predictions.
    6. Join the discussion on the [Challenge Discord](https://discord.gg/MY5cEFHH3D)!
    ## 📊 Data:
    The training set contains the following parameters:
    
    | Column                       | Unit        | Type      | Description                                   |
    |:---------------------------- |:----------: |:--------: |:----------------------------------------------|
    | Molecule Name                |             |    str    | Identifier for the molecule |
    | Smiles                       |             |    str    | Text representation of the 2D molecular structure |
    | LogD                         |             |   float   | LogD |
    | KSol                         |    uM       |   float   | Kinetic Solubility |
    | MLM CLint                    | mL/min/kg   |   float   | Mouse Liver Microsomal |
    | HLM CLint                    | mL/min/kg   |   float   | Human Liver Microsomal |
    | Caco-2 Permeability Efflux   |             |   float   | Caco-2 Permeability Efflux Ratio |
    | Caco-2 Permeability Papp A>B | 10^-6 cm/s  |   float   | Caco-2 Permeability Papp A>B |
    | MPPB                         | % Unbound   |   float   | Mouse Plasma Protein Binding |
    | MBPB                         | % Unbound   |   float   | Mouse Brain Protein Binding |
    | MGMB                         | % Unbound   |   float   | Mouse Gastrocnemius Muscle Binding |
 
    You can download the training data from the [Hugging Face dataset](https://huggingface.co/datasets/openadmet/openadmet-challenge-train-data).

    The test set will remained blinded until the challenge submission deadline. You will be tasked with predicting the same set of ADMET endpoints for the test set molecules.
    
    The training and blinded test set will also be made available on the [CDD Vault](https://www.collaborativedrug.com/). An account to access the CDD Vault can be requested by filling out this [form](https://forms.gle/KiviZ7AaGcuqtrwH8).
    Note that by joining the Vault, your account will be visible to other participants, so this option is **not recommended for those wishing to remain anonymous.**

    ### Experiment uncertainties
    For each endpoint, we present the experimental error statistics, measured with different control compounds. Special thanks to Jon Ainsley for providing this data!

    | Property | Compound | Count | Min | Max | SD | Mean abs. dev. | Median abs. dev. | Mean | Median |
    | :---- | :---- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
    | LogD | Progesterone | 456 | 3.7 | 4.02 | 0.07 | 0.06 | 0.05 | 3.85 | 3.85 |
    | KSOL | Progesterone | 336 | 12.33 | 25.69 | 2.28 | 1.79 | 1.48 | 16.88 | 16.8 |
    | MLM CLint | Verapamil | 276 | 237.67 | 514.78 | 76.03 | 69.34 | 36.78 | 334.41 | 292.73 |
    | HLM CLint | Verapamil | 359 | 104.68 | 304.07 | 61.28 | 57.58 | 61.9 | 197.84 | 210.25 |
    | Caco-2 Permeability Efflux | Atenolol | 113 | 0.85 | 1.71 | 0.2 | 0.17 | 0.16 | 1.28 | 1.27 |
    |  | Digoxin | 383 | 25.08 | 69.71 | 10.8 | 8.77 | 6.17 | 37.76 | 34.22 |
    |  | Minoxidil | 383 | 0.74 | 1.39 | 0.12 | 0.1 | 0.08 | 1.06 | 1.05 |
    |  | Propranolol | 270 | 0.66 | 1.29 | 0.13 | 0.11 | 0.09 | 0.96 | 0.94 |
    | Caco-2 Permeability Papp A\>B | Atenolol | 113 | 0.23 | 0.96 | 0.17 | 0.14 | 0.13 | 0.57 | 0.56 |
    |  | Digoxin | 383 | 0.22 | 0.82 | 0.13 | 0.11 | 0.1 | 0.54 | 0.55 |
    |  | Minoxidil | 383 | 4.99 | 9.98 | 1.17 | 0.98 | 0.87 | 7.8 | 7.81 |
    |  | Propranolol | 270 | 16.16 | 33.93 | 3.57 | 2.88 | 2.35 | 25.34 | 25.52 |
    | MPPB % Unbound | Ketoconazole | 265 | 98.96 | 99.62 | 0.13 | 0.1 | 0.08 | 99.35 | 99.35 |
    | MBPB % Unbound | Propranolol | 260 | 96.55 | 97.88 | 0.26 | 0.21 | 0.18 | 97.3 | 97.33 |
    | MGPB % Unbound | Propranolol | 95 | 92.28 | 98.53 | 1.06 | 0.76 | 0.54 | 94.18 | 93.94 |


    ## 📝 Evaluation
    The challenge will be judged based on the following criteria:
    - We welcome submissions of any kind, including machine learning and physics-based approaches. You can also employ pre-training approaches as you see fit, 
    as well as incorporate data from external sources into your models and submissions. 
    - In the spirit of open science and open source we would love to see code showing how you created your submission if possible, in the form of a Github Repository. 
    If not possible due to IP or other constraints you must at a minimum provide a short report written methodology based on the template [here](https://docs.google.com/document/d/1bttGiBQcLiSXFngmzUdEqVchzPhj-hcYLtYMszaOqP8/edit?usp=sharing).
    **Make sure your lat submission before the deadline includes a link to a report or to a Github repository.**
    - Each participant can submit as many times as they like, up to a limit of once per day. **Only your latest submission will be considered for the final leaderboard.**
    - The endpoints will be judged individually by mean absolute error (**MAE**), while an overall leaderboard will be judged by the macro-averaged relative absolute error (**MA-RAE**). 
    - For endpoints that are not already on a log scale (e.g LogD) they will be transformed to log scale to minimize the impact of outliers on evaluation.
    - We will estimate errors on the metrics using bootstrapping and use the statistical testing workflow outlined in [this paper](https://chemrxiv.org/engage/chemrxiv/article-details/672a91bd7be152b1d01a926b) to determine if model performance is statistically distinct.
    
    ## 📅 **Timeline**:  
    - **September 16:** Challenge announcement
    - **October 14:** Second announcement and sample data release
    - **October 27:** Challenge starts
    - **October-November:** Online Q&A sessions and support via the Discord channel
    - **December 1st:** Intermediate leaderboard release
    - **January 19, 2026:** Submission closes
    - **January 27, 2026:** Winners announced

    ## Acknowledgements
    We gratefully acknowledge Jon Ainsley, Andrew Good, Elyse Bourque, Lakshminarayana Vogeti, Renato Skerlj, Tiansheng Wang, and Mark Ledeboer for generously
    providing the Expansion Therapeutics dataset used in this challenge as an in-kind contribution.
    ---
    """
    # --- Gradio Interface ---
    gr.HTML("""
            <style>
            /* bold only the "Overall" tab label */
                #lb_subtabs [role="tab"][aria-controls="all_tab"] {
                    font-weight: 700 !important;
                }
            </style>
            <style>
            #welcome-md table {
                width: 80%;
                border-collapse: collapse;
                font-size: 0.95rem;     
                line-height: 1.2;    
            } 
            #welcome-md th, #welcome-md td {
                padding: 6px 10px;      
                border: 1px solid rgba(0,0,0,0.9);
                vertical-align: middle; 
            }
            #welcome-md thead th {
                background: var(--panel-background-fill, #f5f5f7);
                font-weight: 1000;
            }
            /* Header shading */
            #welcome-md thead th:nth-child(2),
            #welcome-md thead th:nth-child(3) {
                text-align: center;
            }
            /* Zebra striping */
            #welcome-md tbody tr:nth-child(odd)  { background: rgba(0,0,0,0.03); }
            #welcome-md tbody tr:hover          { background: rgba(0,0,0,0.06); }
            /* Align columns */
            #welcome-md td:nth-child(2),
            #welcome-md td:nth-child(3) { text-align: center; white-space: nowrap; }
            </style>
            """)
    with gr.Tabs(elem_classes="tab-buttons"):
        lboard_dict = {}
        with gr.TabItem("📖 About"):
            gr.Markdown(welcome_md, elem_id="welcome-md")
        with gr.TabItem("🚀 Leaderboard", elem_id="lb_subtabs"):
            gr.Markdown("""
                        View the leaderboard for each ADMET endpoint by selecting the appropiate tab.

                        Latest results from the live leaderboard (evaluated with the validation set), can be found in the "VALIDATION" tab, while the final leaderboard (evaluated on the full blinded test set) is in the "FINAL" tab.
                        **For the per-endpoint leaderboards, the results shown correspond to evaluating each endpoint with the full blinded set. Only entries in the final leaderboard are shown.**

                        """)
            # Make separate leaderboards in separate tabs  
            #per_ep = build_leaderboard()
            df_per_ep = build_leaderboard_static()
            # Aggregated leaderboard
            with gr.TabItem('VALIDATION', elem_id="all_tab"):
                current_df = pd.read_csv("all_leaderboards/current_leaderboard.csv")
                lboard_dict['Average'] = Leaderboard(
                    value=df_per_ep['Average'], #(current_df, current_df_raw)['Average'],
                    datatype=['number'] + LB_DTYPES,
                    select_columns=LB_AVG,
                    search_columns=["user"],
                    render=True,
                    # every=60,
                )
                # Set up button to download leaderboard as csv file 
                download_lb = gr.DownloadButton(
                        label="📥 Download current leaderboard as a csv file",
                        value=None,
                        variant="secondary",
                        )
                download_lb.click(
                    fn=get_leaderboard_csv_aux,
                    inputs=[gr.State(current_df), gr.State("validation")],
                    outputs=download_lb,
                    show_progress="hidden",
                )
            with gr.TabItem('FINAL', elem_id="all_tab"):
                df = df_per_ep['Final'].copy() # making a copy so that model details in csv version aren't html
                df['model details'] = df['model details'].apply(lambda x: make_tag_clickable(x)).astype(str)
                lboard_dict['Final'] = Leaderboard(
                    value=df,
                    datatype=['number', 'markdown', 'str', 'number', 'number', 'number', 'number', 'markdown'],
                    select_columns=final_cols,
                    search_columns=["user"],
                    render=True,
                )
                # Set up button to download leaderboard as csv file 
                download_lb = gr.DownloadButton(
                        label="📥 Download final leaderboard as a csv file",
                        value=None,
                        variant="secondary",
                        )
                download_lb.click(
                    fn=get_leaderboard_csv_aux,
                    inputs=[gr.State(df_per_ep['Final'].copy()),gr.State("final")], 
                    outputs=download_lb, 
                    show_progress="hidden", 
                )
            # per-endpoint leaderboard
            for endpoint in ENDPOINTS:
                with gr.TabItem(endpoint):
                    df = df_per_ep[endpoint].copy()
                    df['model details'] = df['model details'].apply(lambda x: make_tag_clickable(x)).astype(str)
                    lboard_dict[endpoint] = Leaderboard(
                        value=df, #(current_df, current_df_raw)[endpoint],
                        datatype=['number', 'markdown', 'number', 'number', 'number', 'number', 'markdown'],
                        select_columns=ep_cols,
                        search_columns=["user"],
                        render=True,
                        # every=60,
                    )
                    # Set up download button
                    download_lb = gr.DownloadButton(
                            label="📥 Download leaderboard as a csv file",
                            value=None,
                            variant="secondary",
                            )
                    download_lb.click(
                        fn=get_leaderboard_csv_aux,
                        inputs=[gr.State(df_per_ep[endpoint].copy()),gr.State(endpoint)],  
                        outputs=download_lb, 
                        show_progress="hidden", 
                    )
            # Auto-refresh 
            # def refresh_if_changed():
            #     logger.info("Refreshing on timer tick...")
            #     per_ep = build_leaderboard(current_df, current_df_raw)
            #     #return [gr.update(value=per_ep.get(ep, pd.DataFrame(columns=LB_COLS))) for ep in ALL_EPS]
            #     return [per_ep[ep] for ep in ALL_EPS]
            # data_version.change(fn=refresh_if_changed, outputs=[lboard_dict[ep] for ep in ALL_EPS])
        with gr.TabItem("✉️ Submit"):
            gr.Markdown(
            """
            # ADMET Endpoints Submission
            Upload your prediction files here as a csv file.
            """
            )
            filename = gr.State(value=None) 
            eval_state = gr.State(value=None) 
            user_state = gr.State(value=None)

            with gr.Row():
                
                with gr.Column():
                    gr.Markdown(
                        """
                        ## Participant Information
                        To participate, **we require a Hugging Face username**, which will be used to track multiple submissions.
                        Your username will be displayed on the leaderboard, unless you check the *anonymous* box. If you want to remain anonymous, please provide an alias to be used for the leaderboard (we'll keep the username hidden).
                        If you wish to be included in Challenge discussions, please provide your Discord username and email. 
                        If you wish to be included in a future publication with the Challenge results, please provide your name and affiliation (and check the box below).
                        We also ask you to provide a link to a report decribing your method. While not mandatory at the time of participation, 
                        you need to submit the link before the challenge deadline in order to be considered for the final leaderboard.
     
                        """
                        )
                    
                    username_input = gr.Textbox(
                        label="Username", 
                        placeholder="Enter your Hugging Face username",
                        # info="This will be displayed on the leaderboard."
                    )
                    user_alias = gr.Textbox(
                        label="Optional Alias", 
                        placeholder="Enter an identifying alias for the leaderboard if you wish to remain anonymous",
                        # info="This will be displayed on the leaderboard."
                    )
                    anon_checkbox = gr.Checkbox(
                        label="I want to submit anonymously",
                        info="If checked, your username will be replaced with the given *alias* on the leaderboard.",
                        value=False,
                    )
                with gr.Column():
                    # Info to track participant, that will not be displayed publicly
                    participant_name = gr.Textbox(
                        label="Participant Name",
                        placeholder="Enter your name (optional)",
                        info="This will not be displayed on the leaderboard but will be used for tracking participation."
                    )
                    discord_username= gr.Textbox(
                        label="Discord Username",
                        placeholder="Enter your Discord username (optional)",
                        info="Enter the username you will use for the Discord channel (if you are planning to engage in the discussion)."
                    )
                    email = gr.Textbox(
                        label="Email",
                        placeholder="Enter your email (optional)",
                    )
                    affiliation = gr.Textbox(
                        label="Affiliation",
                        placeholder="Enter your school/company affiliation (optional)",
                    )
                    model_tag = gr.Textbox(
                        label="Model Report",
                        placeholder="Link to a report describing your method (optional)",
                    )
                    paper_checkbox = gr.Checkbox(
                        label="I want to be included in a future publication detailing the Challenge results",
                        value=False,
                    )
            with gr.Row():
                with gr.Column():
                    gr.Markdown(
                        """
                        ## Submission Instructions
                        After training your model with the [ExpansionRx trainining set](https://huggingface.co/datasets/openadmet/openadmet-challenge-train-data),
                        please upload a single CSV file containing your predictions for all compounds in the test set.
                        Only your latest submission will be considered.
                        Download a CSV file with the compounds in the test set here:
                        **NOTE: Submission can sometimes take a few minutes to process**
                        **Please be patient and wait for the status message to update and your submission to reach the leaderboard.**
                        """
                    )
                    download_btn = gr.DownloadButton(
                        label="📥 Download Test Set Compounds",
                        value="./data/expansion_data_test_blinded.csv",
                        variant="secondary",
                        )
                with gr.Column():
                    predictions_file = gr.File(label="Single file with ADMET predictions (.csv)",
                                            file_types=[".csv"],
                                            file_count="single",)
            username_input.change(
                fn=lambda x: x if x.strip() else None,
                inputs=username_input,
                outputs=user_state
            )  
            submit_btn = gr.Button("📤 Submit Predictions (closed)", interactive=False)
            message = gr.Textbox(label="Status", lines=1, visible=False)
            #submit_btn.click(
            #    submit_data,
            #    inputs=[predictions_file, 
            #            user_state, 
            #            participant_name, 
            #            discord_username, 
            #            email, 
            #            affiliation, 
            #            model_tag, 
            #            user_alias,
            #            anon_checkbox,
            #            paper_checkbox],
            #    outputs=[message, filename],
            #).success(
            #    fn=lambda m: gr.update(value=m, visible=True),
            #    inputs=[message],
            #    outputs=[message],
            #).success(
            #    fn=evaluate_data,
            #    inputs=[filename],
            #    outputs=[eval_state]
            #)
        with gr.TabItem("🛠️ FAQ"):
            with gr.Column():
                gr.Markdown(
                    """
                    1. How long does it normally take for submissions to reach the leaderboard?

                    > At most, this should take 2 to 3 minutes. If it's taking longer, please ping us on Discord and let us know.

                    2. The leaderboard isn't updating.

                    > Early in the challenge we had a problem with the leaderboard not updating.  We believe this has been fixed. Please reach out on [Discord](https://discord.com/channels/1412827471488745545/1413320650281455657) if you have an issue.

                    2. My submission didn't upload, what's wrong.
                    
                    > Please check your submission and confirm that the column names are the same as those in the test set file. If you run into an issue, please reach out on [Discord](https://discord.com/channels/1412827471488745545/1413320650281455657), we'd be happy to help. 

                    3. I only want to submit for a one endpoint, what should I do.
                    
                    > Right now, you have to submit all columns. If you'd like to just submit LogD, or some other column, please include the other columns and put zeros for all values. This should still rank you on the LogD leaderboard. Your rankings on the other leaderboards will be low. 

                    4. What is the formula for macro-averaged relative absolute error (MA-RAE)? What is this relative to?
                    
                    > The code we use to compute all the metrics is available [here](https://huggingface.co/spaces/openadmet/OpenADMET-ExpansionRx-Challenge/tree/main)

                    5. Organizers, where are you?

                    > We try to answer questions on [Discord](https://discord.com/channels/1412827471488745545/1413320650281455657) as quickly as possible, but we occassionally need to sleep. 🙂

                    6. Do I have to upload my model to HuggingFace?
                    
                    > No, you only need to upload the results.

                    7. Do I have to make the code for my model public?
                    
                    > No, while we love open source, you don't have to make your model public.  We would appreciate a brief description of how you built your models.
                    > NOTE: You must provide a link to a report or github repository before the challenge deadline in order to be considered for the final leaderboard.

                    8. Can I use data beyond the training set to train my model.
                    
                    > Yes, absolutely, you're free to use any data you'd like to train your model.

                    9. How are you handling log transforms of zero values?
                    
                    > Please see the function function [clip_and_log_transform](https://huggingface.co/spaces/openadmet/OpenADMET-ExpansionRx-Challenge/blob/main/utils.py). This function adds 1 to both the training and test data before doing the log transform.

                    """)
if __name__ == "__main__":
    logger.info("Starting Gradio app...")
    demo.launch(ssr_mode=False)
    logger.info("Gradio app closed.")