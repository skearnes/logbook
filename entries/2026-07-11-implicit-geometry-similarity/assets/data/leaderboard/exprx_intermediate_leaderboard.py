from statsmodels.stats.multicomp import pairwise_tukeyhsd
from cld import cld
from utils import (
    check_page_exists,
    map_metric_to_stats,
    fetch_dataset_df,
)
from about import ENDPOINTS, LB_COLS, results_repo_test, METRICS
import pandas as pd



def validate_hf_username(username):
    username = str(username).strip()
    hf_url = f"https://huggingface.co/{username}"
    return check_page_exists(hf_url, delay=1, max_retries=10)
    # return True  # For testing purposes, assume all usernames are valid

def validate_model_details(tag):
    if tag is None:
        return "Not submitted"
    safe_tag = str(tag).strip()
    if not safe_tag.startswith("https://"):
        return "Invalid link"
    is_real_url = check_page_exists(safe_tag, delay=2)
    if not is_real_url:
        return "Invalid link"
    else:
        return safe_tag



def make_intermediate_lb():

    df_latest, df_latest_raw = fetch_dataset_df(
        download_raw=True, 
        test_repo=results_repo_test
    )

    # Make all usernames lowercase
    df_latest_raw["hf_username"] = df_latest_raw["hf_username"].str.lower()

    # HF username validation
    hf_usernames = df_latest_raw["hf_username"].unique()
    valid_hf_usernames = {username: validate_hf_username(username) for username in hf_usernames}

    # print all users and their validation status
    for username, is_valid in valid_hf_usernames.items():
        print(f"Username: {username}, Valid: {is_valid}")

    df_latest_raw["hf_user_valid"] = df_latest_raw["hf_username"].map(valid_hf_usernames)
    # drop invalid usernames
    df_latest_raw = df_latest_raw[df_latest_raw["hf_user_valid"]].reset_index(drop=True)

    # make sure to only keep the latest submission per user for the 'Average' endpoint
    df_latest_raw["submission_time"] = pd.to_datetime(df_latest_raw["submission_time"])
    df_latest_raw = df_latest_raw.query("Endpoint == 'Average'")
    df_latest_raw['latest_time_per_user'] = df_latest_raw.groupby('hf_username')['submission_time'].transform('max')
    latest_submissions_df = df_latest_raw[df_latest_raw['submission_time'] == df_latest_raw['latest_time_per_user']].copy()
    # Fix to order by the mean RAE and not the RAE of all samples (slight missmatch for some users)
    latest_submissions_df['mean_RAE'] = latest_submissions_df.groupby('hf_username')['RAE'].transform('mean')    
    latest_submissions_df = latest_submissions_df.sort_values(
        by=['mean_RAE', 'Sample'], ascending=True
    ).reset_index(drop=True)

    # Get the unique users in the order of their first appearance
    unique_users_ordered = latest_submissions_df['user'].unique()

    # Create a mapping dictionary: original_user -> prefixed_user
    user_mapping = {}
    for idx, user in enumerate(unique_users_ordered):
        # The prefix is the index starting from 0, formatted to be 3 digits (001, 002, etc.)
        # We use idx + 1 to start the sequence from 001 instead of 000
        prefix = f"{idx + 1:03d}"
        prefixed_user = f"{prefix}___{user}"
        user_mapping[user] = prefixed_user

    # Apply the mapping to create a new column with prefixed usernames
    latest_submissions_df['user'] = latest_submissions_df['user'].map(user_mapping)

    # Perform Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=latest_submissions_df['RAE'], groups=latest_submissions_df['user'], alpha=0.05)
    tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], 
                            columns=tukey._results_table.data[0])

    # add CLDs
    cld_dict = cld(tukey_df)

    cld_df = pd.DataFrame(cld_dict.items(),columns=["group","letter"]).sort_values("group")
    cld_df.letter = [",".join(x) for x in cld_df.letter]
    cld_df["user"] = cld_df.group

    # clean up CLD letters for extended alphabet (i.e with @ symbols)
    def clean_up(ser):
        ser = ser.split(",")
        # rejoin for late in alphabet
        if "@" in ser and len(ser) == 2:
            let = "@" + ser[1]
        elif "@" in ser and len(ser) == 4:
            let = "@" + ser[2] + "," + "@" + ser[3]
        else:
            let = ",".join(ser)
        return let
    
    cld_df["fixed_letter"] = cld_df["letter"].apply(lambda x: clean_up(x))
    report_cols = latest_submissions_df[['user', 'model_report']].drop_duplicates(keep='first')

    # gather means and stds for each metric for each user
    for metric in METRICS:
        metric_stats = latest_submissions_df.groupby('user')[metric].agg(['mean', 'std']).reset_index()
        metric_stats = metric_stats.rename(columns={'mean': f'{metric}_mean', 'std': f'{metric}_std'})
        metric_stats[f"{metric}_display"] = metric_stats.apply(
            lambda row: f"{row[f'{metric}_mean']:.4f} ± {row[f'{metric}_std']:.4f}", axis=1
        )
        cld_df = metric_stats[['user', f'{metric}_mean', f'{metric}_std', f'{metric}_display']].merge(cld_df, on='user', how='left')

    # re-sort by RAE mean, lowest is best
    cld_df = cld_df.sort_values(by='RAE_mean', ascending=True).reset_index(drop=True)
    cld_df = cld_df.merge(report_cols, on='user', how='inner')
    cld_df['user'] = cld_df['user'].str.split('___').str[1]

    cld_subset = cld_df[['user', 'fixed_letter'] + [f'{metric}_display' for metric in METRICS]+ ['model_report']]
    cld_subset = cld_subset.rename(columns={'user': 'user', 'fixed_letter': 'CLD', 'model_report': 'model details'})
    cld_subset['model details'] = cld_subset['model details'].apply(lambda x: validate_model_details(x)).astype(str)

    print(cld_subset.head())
    cld_subset.to_csv("leaderboard_cld_results.csv", index=False)

if __name__ == "__main__":
    make_intermediate_lb()