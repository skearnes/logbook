"""Submission persistence — uploads prediction files and metadata to S3.

Architecture
------------
Each submission is stored in S3 under a deterministic prefix::

    s3://{S3_BUCKET}/submissions/{track}/{username}/{submission_id}/
        metadata.json       <- Submission model serialised as JSON
        {original_filename} <- The uploaded prediction file

A scheduled AWS Lambda function (``backend/lambda_handler.py``) polls this
prefix for unscored submissions, runs evaluation, and writes scores back.

Environment variables
---------------------
S3_BUCKET
    Name of the S3 bucket (required at runtime, not at import time).
AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / AWS_DEFAULT_REGION
    Standard boto3 credentials — set via HuggingFace Space secrets.
"""

import os
from datetime import datetime, timezone
from pathlib import Path

import boto3
from loguru import logger
from models import Submission


def upload_submission(submission: Submission, file_path: Path) -> Submission:
    """Upload the prediction file and serialised metadata to S3.

    On success, returns the submission with ``s3_key`` populated.
    On failure, logs the error and returns the submission unchanged so the
    caller can still surface a user-facing message.

    Args:
        submission: Validated Submission instance (s3_key will be set here).
        file_path: Local path to the uploaded prediction file.

    Returns:
        The submission with s3_key set to the uploaded object key.

    Todo:
        - Set S3_BUCKET, AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY,
          AWS_DEFAULT_REGION as HuggingFace Space secrets before deploying.
        - Consider server-side encryption (SSE-S3 or SSE-KMS).
        - Add retry logic for transient S3 errors.

    """
    bucket = os.environ.get("S3_BUCKET")
    if not bucket:
        logger.warning(
            "S3_BUCKET not set — submission will not be persisted. "
            "Set S3_BUCKET and AWS credentials as Space secrets to enable storage."
        )
        return submission

    if submission.track == "Structure Prediction":
        canonical_filename = "structures.zip"
    else:
        canonical_filename = f"predictions{Path(file_path).suffix}"
    file_key = f"{submission.s3_prefix}/{canonical_filename}"
    metadata_key = f"{submission.s3_prefix}/metadata.json"
    submission = submission.model_copy(update={"s3_key": file_key})

    try:
        s3 = boto3.client("s3")

        # Upload prediction file
        logger.info(f"Uploading prediction file to s3://{bucket}/{file_key}")
        s3.upload_file(str(file_path), bucket, file_key)

        # Upload metadata JSON
        logger.info(f"Uploading metadata to s3://{bucket}/{metadata_key}")
        s3.put_object(
            Bucket=bucket,
            Key=metadata_key,
            Body=submission.model_dump_json(indent=2).encode(),
            ContentType="application/json",
        )

        logger.info(f"Submission {submission.submission_id!r} stored successfully.")

    except Exception as exc:
        logger.error(
            f"S3 upload failed for submission {submission.submission_id!r}: {exc}"
        )

    return submission


def _fetch_last_submission_date(track: str, user_id: str) -> datetime | None:
    """Fetch the submission date of the most recent submission for a track and user.

    Args:
        track (str): The track name (e.g., "activity" or "structure").
        user_id (str): The user ID to check for previous submissions.

    Returns:
        datetime | None: The submission date of the most recent submission, or None if
                         no previous submissions are found.

    """
    bucket = os.environ.get("S3_BUCKET")
    if not bucket:
        logger.warning(
            "S3_BUCKET not set — cannot fetch last submission date. "
            "Set S3_BUCKET and AWS credentials as Space secrets to enable this feature."
        )
        return None

    s3 = boto3.client("s3")
    prefix = f"submissions/{track}/{user_id}/"
    try:
        response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
        if response["IsTruncated"]:  # Unlikely to be > 1000 submissions per user
            logger.warning(
                f"ListObjectsV2 response truncated for prefix {prefix!r}. "
                "Only the first 1000 objects will be considered."
            )
        if "Contents" not in response:
            return None  # No submissions found

        logger.info(
            f"Found {len(response['Contents'])} objects under prefix {prefix!r}."
        )
        submission_dates = []
        for obj in response["Contents"]:
            if obj["Key"].endswith("metadata.json"):
                submission_dates.append(obj["LastModified"])

        if not submission_dates:  # Shouldn't be possible
            logger.warning(f"No metadata.json files found under prefix {prefix!r}.")
            return None

        return max(submission_dates).astimezone(timezone.utc)

    except Exception as exc:
        logger.error(f"Failed to fetch last submission date for {user_id!r}: {exc}")
        return None
