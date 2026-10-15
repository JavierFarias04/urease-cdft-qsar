"""Input integrity. Step 01 calls this before reading anything.

The published results were computed from one specific snapshot of the three
input workbooks. If a reader's copies differ, every downstream number may
differ too, and a silent mismatch is worse than a stopped run.
"""
import hashlib

from . import config as cfg


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def verify_input_checksums(strict=True):
    """Compare the raw inputs against the recorded SHA-256 digests.

    Returns {filename: (status, digest)} with status one of 'ok', 'MISMATCH',
    'MISSING'. Raises on any problem when `strict`, so a pipeline run cannot
    proceed on the wrong data snapshot.
    """
    report, problems = {}, []
    for name, expected in sorted(cfg.INPUT_SHA256.items()):
        path = cfg.RAW / name
        if not path.exists():
            report[name] = ("MISSING", None)
            problems.append(f"{name}: not found at {path.relative_to(cfg.ROOT)}")
            continue
        got = sha256(path)
        if got == expected:
            report[name] = ("ok", got)
        else:
            report[name] = ("MISMATCH", got)
            problems.append(f"{name}: expected {expected[:16]}..., got {got[:16]}...")
    if problems and strict:
        raise SystemExit(
            "Input verification failed:\n  " + "\n  ".join(problems)
            + "\n\nSee data/raw/README.md for how to obtain the correct snapshot."
        )
    return report
