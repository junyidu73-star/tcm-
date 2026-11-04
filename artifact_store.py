"""Content-addressed local object store with an optional S3-compatible backend."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path


class LocalObjectStore:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def put_file(self, path: str | Path, object_key: str):
        source = Path(path)
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        target = self.root / digest[:2] / digest / Path(object_key).name
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(source.read_bytes())
        return {"uri": str(target.resolve()), "sha256": digest, "backend": "local_object_store"}


class S3ObjectStore:
    def __init__(self, bucket: str, endpoint_url: str | None = None):
        import boto3
        self.bucket = bucket
        self.client = boto3.client("s3", endpoint_url=endpoint_url or None)

    def put_file(self, path: str | Path, object_key: str):
        source = Path(path)
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        key = f"{digest[:2]}/{digest}/{Path(object_key).name}"
        self.client.upload_file(str(source), self.bucket, key, ExtraArgs={"Metadata": {"sha256": digest}})
        return {"uri": f"s3://{self.bucket}/{key}", "sha256": digest, "backend": "s3"}


def configured_object_store(local_root: str | Path):
    bucket = os.environ.get("ZHIHENG_OBJECT_BUCKET", "").strip()
    if bucket:
        try:
            return S3ObjectStore(bucket, os.environ.get("ZHIHENG_S3_ENDPOINT"))
        except ImportError as exc:
            raise RuntimeError("S3 storage is configured but boto3 is not installed") from exc
    return LocalObjectStore(local_root)
