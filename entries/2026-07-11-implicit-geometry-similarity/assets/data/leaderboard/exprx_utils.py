
import pandas as pd
import numpy as np
from typing import Tuple
from datasets import load_dataset, Features, Value
from about import results_repo_validation, results_repo_test
from about import METRICS, STANDARD_COLS
from loguru import logger
import time
import requests

import requests
import time

def check_page_exists(url: str, delay=0.2, max_retries=3, current_retries=0):
    """Checks if a web page exists at the given URL with a retry limit for 429 errors.

    Parameters
    ----------
    url : str
        Url of the page
    delay : float, optional
        Seconds to wait until submitting another request, by default 0.2
    max_retries : int, optional
        Maximum number of times to retry on a 429 error, by default 3
    current_retries : int, optional
        Current number of retries performed (internal counter), by default 0

    Returns
    -------
    bool
       If the page exists
    """
    safe_url = str(url).strip()
    
    # Attempt to fix url
    if not safe_url.startswith(('http://', 'https://')):
        safe_url = f"https://{safe_url}"
        
    try:
        response = requests.get(safe_url, timeout=5) 
        
        # Check for Rate Limit Error and retry if under the limit
        if response.status_code == 429:
            if current_retries < max_retries:
                # Make wait time exponential
                wait_time = 5 * (2 ** current_retries)
                print(f"Warning: Rate limit hit on {safe_url}. Attempt {current_retries + 1}/{max_retries}. Waiting for {wait_time} seconds...")
                time.sleep(wait_time)
                # Recurse with an incremented retry counter
                return check_page_exists(safe_url, delay=delay, max_retries=max_retries, current_retries=current_retries + 1) 
            else:
                print(f"Error: Max retries ({max_retries}) reached for rate limit on {safe_url}.")
                return False # Give up after max retries
        
        # Return True only for a successful status code (200)
        return response.status_code == 200
    
    except requests.exceptions.RequestException as e:
        print(f"Error checking URL {safe_url}: {e}")
        return False
        
    finally:
        # Sleep after every request to avoid HTTPS error
        time.sleep(delay)

def make_user_clickable(name: str):
    link =f'https://huggingface.co/{name}'
    return f'<a target="_blank" href="{link}" style="color: var(--link-text-color); text-decoration: underline;text-decoration-style: dotted;">{name}</a>'
def make_tag_clickable(tag: str):
    if tag is None:
        return "Not submitted"
    return f'<a target="_blank" href="{tag}" style="color: var(--link-text-color); text-decoration: underline;text-decoration-style: dotted;">link</a>'

def fetch_dataset_df(download_raw=False, test_repo=results_repo_validation): # Change download_raw to True for the final leaderboard
    logger.info("Fetching latest results dataset from Hugging Face Hub...")
    # Specify feature types to load results dataset
    metric_features = {
        f'mean_{m}': Value('float64') for m in METRICS
    }
    metric_features.update({
        f'std_{m}': Value('float64') for m in METRICS
    })
    other_features = {
        'user': Value('string'),
        'Endpoint': Value('string'),
        'submission_time': Value('string'),
        'model_report': Value('string'), 
        'anonymous': Value('bool'), 
        'hf_username': Value('string')
    }
    feature_schema = Features(metric_features | other_features)

    dset = load_dataset(test_repo, 
                        name='default',
                        split='train', 
                        features=feature_schema,
                        download_mode="force_redownload")
    full_df = dset.to_pandas()
    expected_mean_cols = [f"mean_{col}" for col in METRICS]
    expected_std_cols = [f"std_{col}" for col in METRICS]
    expected_all_cols = STANDARD_COLS + expected_mean_cols + expected_std_cols
    assert all(
        col in full_df.columns for col in expected_all_cols
    ), f"Expected columns not found in {full_df.columns}. Missing columns: {set(expected_all_cols) - set(full_df.columns)}"

    df = full_df.copy()
    df = df[df["user"] != "test"].copy()
    df["submission_time"] = pd.to_datetime(df["submission_time"], errors="coerce")
    df = df.dropna(subset=["submission_time"])

    # Get the most recent submission per user & endpoint
    latest = (
        df.sort_values("submission_time")
          .drop_duplicates(subset=["Endpoint", "hf_username"], keep="last") #IMPORTANT: unique on HF username not display name
          .sort_values(["Endpoint", "user"])
          .reset_index(drop=True)
    )
    latest.rename(columns={"submission_time": "submission time"}, inplace=True)

    # Also fetch raw dataset
    # We'll set download_raw to False for the live leaderboard, as it's too large to load
    latest_raw = None
    if download_raw:
        raw_metric_features = {
            m: Value('float64') for m in METRICS
        }
        other_features_raw = other_features.copy()
        other_features_raw.update({'Sample': Value("float32")})
        feature_schema = Features(raw_metric_features | other_features_raw)
        logger.info("Fetching raw bootstrapping dataset from Hugging Face Hub...")
        # Because the raw file is so long, we have to load it with delay and multiple retries
        max_retries = 10
        base_delay = 5
        for attempt in range(max_retries):
            try:
                logger.info("Attempting to download raw data")
                dset_raw = load_dataset(test_repo, 
                                    name='raw',
                                    split='train', 
                                    features=feature_schema,
                                    download_mode="force_redownload")
                raw_df = dset_raw.to_pandas()
                df_raw = raw_df.copy()
                df_raw["submission_time"] = pd.to_datetime(df_raw["submission_time"], errors="coerce")
                df_raw = df_raw.dropna(subset=["submission_time"])
                latest_raw = (
                    df_raw.sort_values("submission_time")
                    .drop_duplicates(subset=["Sample", "Endpoint", "hf_username"], keep="last") 
                    .sort_values(["Sample","Endpoint", "user"])
                    .reset_index(drop=True)
                )
                break # Exit try loop if successful
            except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, Exception) as e:
                logger.warning(f"Error fetching raw dataset on attempt {attempt + 1}: {e.__class__.__name__}. Retrying...")
                if attempt < max_retries - 1:
                    # Exponential backoff
                    sleep_time = base_delay * (2 ** attempt) 
                    logger.info(f"Waiting for {sleep_time} seconds before next attempt.")
                    time.sleep(sleep_time)
                else:
                    logger.error(f"Failed to fetch 'raw' dataset after {max_retries} retries.")


    return latest, latest_raw


def clip_and_log_transform(y: np.ndarray):
    """
    Clip to a detection limit and transform to log10 scale.

    Parameters
    ----------
    y : np.ndarray
        The array to be clipped and transformed.
    """
    y = np.clip(y, a_min=0, a_max=None)
    return np.log10(y + 1)


def bootstrap_sampling(size: int, n_samples: int) -> np.ndarray:
    """
    Generate bootstrap samples for a given size and number of samples.

    Parameters
    ----------
    size : int
        The size of the data.
    n_samples : int
        The number of samples to generate.

    Returns
    -------
    np.ndarray
        Returns a numpy array of the bootstrap samples.
    """
    rng = np.random.default_rng(0)
    return rng.choice(size, size=(n_samples, size), replace=True)


def metrics_per_ep(pred: np.ndarray, 
                   true: np.ndarray
    )->Tuple[float, float, float, float]:
    """Predict evaluation metrics for a single sample
    Parameters
    ----------
    pred : np.ndarray
        Array with predictions
    true : np.ndarray
        Array with actual values
    Returns
    -------
    Tuple[float, float, float, float]
        Resulting metrics: (MAE, RAE, R2, Spearman R, Kendall's Tau)
    """
    from scipy.stats import spearmanr, kendalltau
    from sklearn.metrics import mean_absolute_error, r2_score
    mae = mean_absolute_error(true, pred)
    rae = mae / np.mean(np.abs(true - np.mean(true)))
    if np.nanstd(true) == 0:
        r2=np.nan
    else:
        r2 = r2_score(true, pred)

    if np.nanstd(pred) < 0.0001:
        spr = np.nan
        ktau = np.nan
    else:
        spr = spearmanr(true, pred).statistic
        ktau = kendalltau(true, pred).statistic

    return mae, rae, r2, spr, ktau

def bootstrap_metrics(pred: np.ndarray, 
                      true: np.ndarray,
                      endpoint: str,
                      n_bootstrap_samples=1000
    )->pd.DataFrame:
    """Calculate bootstrap metrics given predicted and true values
    Parameters
    ----------
    pred : np.ndarray
        Predicted endpoints
    true : np.ndarray
        Actual endpoint values
    endpoint : str
        String with endpoint
    n_bootstrap_samples : int, optional
        Size of bootstrapsample, by default 1000
    Returns
    -------
    pd.DataFrame
        Dataframe with estimated metric per bootstrap sample for the given endpoint
    """
    cols = ["Sample", "Endpoint", "Metric", "Value"]
    bootstrap_results = pd.DataFrame(columns=cols) 
    for i, indx in enumerate(
        bootstrap_sampling(true.shape[0], n_bootstrap_samples)
    ):
        mae, rae, r2, spr, ktau = metrics_per_ep(pred[indx], true[indx])
        scores = pd.DataFrame(
            [
                [i, endpoint, "MAE", mae],
                [i, endpoint, "RAE", rae],
                [i, endpoint, "R2", r2],
                [i, endpoint, "Spearman R", spr],
                [i, endpoint, "Kendall's Tau", ktau]
            ],
            columns=cols
        )
        bootstrap_results = pd.concat([bootstrap_results, scores])
    return bootstrap_results

def map_metric_to_stats(df: pd.DataFrame, average=False) -> pd.DataFrame: 
    """Map mean and std to 'mean +/- std' string for each metric

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to modify
    average : bool, optional
        Whether the dataframe contains average info, by default False

    Returns
    -------
    pd.DataFrame
        Modified dataframe
    """
    metric_cols = METRICS[:] 
    if average:
        metric_cols[1] = "MA-RAE" 
    cols_drop = []
    for col in metric_cols:
        mean_col = f"mean_{col}"
        std_col = f"std_{col}"
        df[col] = df.apply(
            lambda row: f"{row[mean_col]:.2f} +/- {row[std_col]:.2f}", 
            axis=1
        )
        cols_drop.extend([mean_col, std_col])
    df = df.drop(columns=cols_drop)
    return df