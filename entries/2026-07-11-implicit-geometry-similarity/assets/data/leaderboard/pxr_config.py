"""Configuration constants for the PXR Challenge HuggingFace Space."""

import os

# Phase can be set to 1, 2, or 0 for the end of the challenge
# (or -1 if the final leaderboard is not yet available)
CURRENT_PHASE = 0

ACTIVITY_DATASET_SIZE = 513
STRUCTURE_DATASET_SIZE = 184
REQUIRED_ACTIVITY_COLUMNS = {"SMILES", "Molecule Name", "pEC50"}
HOURS_BETWEEN_SUBMISSIONS = 0

S3_BUCKET: str = os.environ.get("S3_BUCKET", "")
ACTIVITY_LEADERBOARD_S3 = "leaderboard/interim/activity/leaderboard_latest.csv"
STRUCTURE_LEADERBOARD_S3 = "leaderboard/interim/structure/leaderboard_latest.csv"
PHASE_2_ACTIVITY_ENTRIES_S3 = (
    "submissions/manifest/phase_2_entries/activity/leaderboard_latest.csv"
)
