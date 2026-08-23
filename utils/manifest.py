"""Run manifest: what produced an output directory, and a resume guard.

The pipeline resumes from the newest COMPLETED per-round checkpoint plus the
incremental `round_NN/confirmed_families.tsv` files. Nothing checked that the
configuration being resumed under was the one that produced them, so changing
a parameter and resuming silently mixed two configurations inside one output
directory — and `summary.tsv` is written only at the very end, so afterwards
there was no record of which config produced which round. A silent mixture is
worse than a loud failure, hence the guard below.

The manifest (`run_manifest.json` in the output directory) is the record: run
id, git commit, resolved config, its hash, species-tree hash, input file
checksums, resolved tool versions, and one history entry per start/resume.

What the hash covers is the whole point of the guard, so it is a DENYLIST:
every Config field is hashed unless it is named in EXCLUDED_FROM_HASH. A
parameter added later is therefore covered by default — an allowlist would
leave exactly the newest knob unguarded. Excluded are the things that
legitimately differ between a run and its resume: how much machine is thrown
at the work (threads, workers, hmmsearch chunking and its sbatch extras),
where a tool lives (the *_bin paths — a rebuilt conda env moves them without
changing what they compute; the tools block records their resolved versions
instead), which device a model runs on, and the path of a reusable DIAMOND
cache. Command TEMPLATES (*_cmd) ARE hashed: they carry model names and
flags, not just locations.
"""

import hashlib
import json
import logging
import os
import subprocess
import uuid
from dataclasses import asdict, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from config import ToolNotFoundError, resolve_tool

logger = logging.getLogger("family_finder")

MANIFEST_FILE = "run_manifest.json"
MANIFEST_VERSION = 1

#: Config fields deliberately kept OUT of the config hash. See module docstring.
EXCLUDED_FROM_HASH = {
    # how much machine, not what is computed
    "orthofinder_threads",
    "n_workers",
    "hmmer_chunk_size",
    "hmmer_chunk_concurrent",
    "hmmer_chunk_sbatch_extra",
    # where a tool lives; the manifest's tools block records the version
    "orthofinder_bin",
    "mafft_bin",
    "fasttree_bin",
    "iqtree_bin",
    "codeml_bin",
    "pal2nal_bin",
    "hmmbuild_bin",
    "hmmsearch_bin",
    "hmmpress_bin",
    "deeploc_bin",
    "foldseek_bin",
    "epa_ng_bin",
    "raxml_ng_bin",
    "hmmalign_bin",
    # device selection and a cache location
    "deeploc_device",
    "orthofinder_reuse_blast_dir",
}


class ConfigMismatchError(RuntimeError):
    """Resume was refused: the configuration is not the one on record."""


# ---------------------------------------------------------------------------
# hashing
# ---------------------------------------------------------------------------

def _normalize(value):
    """Canonical form of a config value for hashing.

    Lists of strings are sorted: `clustering_species_exclude` and
    `codeml_models` are sets in spirit, and reordering one is not a
    configuration change.
    """
    if isinstance(value, list) and all(isinstance(v, str) for v in value):
        return sorted(value)
    return value


def hashed_config_fields(config) -> dict:
    """The subset of config values the hash covers, normalized."""
    return {
        f.name: _normalize(getattr(config, f.name))
        for f in fields(config)
        if f.name not in EXCLUDED_FROM_HASH
    }


def _canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def sha256_file(path) -> str:
    """Streaming sha256 of a file (proteomes are large)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def config_hash(config, species_tree_path: str) -> str:
    """Hash of the resolved config plus the species tree's CONTENTS.

    The tree file is hashed, not its path: two hosts store the same tree at
    different paths, and an edited tree at the same path is a different
    configuration.
    """
    payload = {
        "config": hashed_config_fields(config),
        "species_tree_sha256": sha256_file(species_tree_path),
    }
    return "sha256:" + hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def git_state(repo_dir: Optional[Path] = None) -> dict:
    """Current commit and whether the working tree was dirty.

    Never fatal: a run from an exported tarball has no git at all.
    """
    repo_dir = Path(repo_dir) if repo_dir else Path(__file__).resolve().parent.parent
    state = {"commit": None, "dirty": None}
    try:
        commit = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=15,
        )
        if commit.returncode != 0:
            return state
        state["commit"] = commit.stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repo_dir), "status", "--porcelain"],
            capture_output=True, text=True, timeout=30,
        )
        if status.returncode == 0:
            state["dirty"] = bool(status.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return state


def _tool_version(path: str) -> Optional[str]:
    """First line of `<tool> --version`, or None.

    Best effort by design: several tools here print usage and exit non-zero,
    and codeml prints no banner at all. stdin is closed and the call is
    capped so a tool that waits for input cannot stall a run's startup.
    """
    try:
        proc = subprocess.run(
            [path, "--version"],
            capture_output=True, text=True, timeout=10,
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    for stream in (proc.stdout, proc.stderr):
        for line in (stream or "").splitlines():
            if line.strip():
                return line.strip()[:200]
    return None


def _tool_keys(config) -> list:
    """The tools this configuration will actually invoke."""
    keys = ["mafft_bin", "pal2nal_bin"]
    if getattr(config, "clustering_method", "orthofinder") == "orthofinder":
        keys.append("orthofinder_bin")
    keys.append("iqtree_bin" if config.tree_builder == "iqtree" else "fasttree_bin")
    if config.hmmer_rescue or config.profile_assign_per_round:
        keys += ["hmmbuild_bin", "hmmsearch_bin", "hmmpress_bin"]
    if config.run_codeml:
        keys.append("codeml_bin")
    return keys


def tool_records(config) -> dict:
    """Resolved path + version per tool, or the failure that names the setting.

    Resolution goes through config.resolve_tool, and a failure is recorded
    rather than raised: the pipeline prepends a conda env to PATH at runtime
    and resolves tools at the point of use, so a bare name that looks
    unresolvable at startup may well be fine later.
    """
    records = {}
    for key in _tool_keys(config):
        configured = getattr(config, key, "")
        rec = {"configured": configured, "resolved": None, "version": None}
        try:
            resolved = resolve_tool(configured, key)
        except ToolNotFoundError as e:
            rec["error"] = str(e)
        else:
            rec["resolved"] = resolved
            rec["version"] = _tool_version(resolved)
            try:
                st = os.stat(resolved)
                rec["size"] = st.st_size
                rec["mtime"] = int(st.st_mtime)
            except OSError:
                pass
        records[key] = rec
    return records


def checksum_dir(path: str) -> dict:
    """sha256 of every regular file in an input directory, by file name."""
    d = Path(path)
    files = {}
    if d.is_dir():
        for p in sorted(d.iterdir()):
            if p.is_file():
                files[p.name] = sha256_file(p)
    return {"path": str(d), "files": files}


# ---------------------------------------------------------------------------
# the guard
# ---------------------------------------------------------------------------

def _short(value):
    if isinstance(value, str) and len(value) == 64:
        return value[:12] + "..."
    return repr(value)


def _config_differences(stored: dict, config, tree_sha: str) -> list:
    """Per-setting (setting, old, new) between a stored manifest and now."""
    old_cfg = stored.get("config") or {}
    new_cfg = asdict(config)
    diffs = []
    for name in sorted(set(hashed_config_fields(config))):
        new = _normalize(new_cfg.get(name))
        if name not in old_cfg:
            diffs.append({"setting": name, "old": None, "new": new,
                          "note": "absent from the stored manifest"})
            continue
        old = _normalize(old_cfg[name])
        if old != new:
            diffs.append({"setting": name, "old": old, "new": new})
    for name in sorted(set(old_cfg) - set(new_cfg)):
        if name not in EXCLUDED_FROM_HASH:
            diffs.append({"setting": name, "old": _normalize(old_cfg[name]),
                          "new": None, "note": "no longer a config field"})

    old_tree = (stored.get("species_tree") or {}).get("sha256")
    if old_tree != tree_sha:
        diffs.append({"setting": "species_tree (contents)",
                      "old": old_tree, "new": tree_sha})
    return diffs


def _input_differences(stored: dict, inputs: dict) -> list:
    """Per-file (role/name, old, new) between stored input checksums and now."""
    old_inputs = stored.get("inputs") or {}
    diffs = []
    for role in sorted(inputs):
        old_files = (old_inputs.get(role) or {}).get("files", {})
        new_files = inputs[role]["files"]
        for name in sorted(set(old_files) | set(new_files)):
            old = old_files.get(name)
            new = new_files.get(name)
            if old != new:
                diffs.append({"setting": f"{role}/{name}", "old": old, "new": new})
    return diffs


def _format_refusal(outdir: Path, diffs: list) -> str:
    lines = [
        f"Resume refused: the configuration differs from the one that produced "
        f"{outdir}.",
        "",
        "  setting: on record -> now",
    ]
    for d in diffs:
        note = f"   ({d['note']})" if d.get("note") else ""
        lines.append(
            f"  {d['setting']}: {_short(d['old'])} -> {_short(d['new'])}{note}"
        )
    lines += [
        "",
        "Resuming would mix results from two configurations in one output "
        "directory, and per-round outputs carry no record of which produced "
        "them. Start a fresh --outdir, or pass --allow-config-change to do it "
        "deliberately (the change is then recorded in the run manifest).",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------

def load_manifest(outdir) -> Optional[dict]:
    path = Path(outdir) / MANIFEST_FILE
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.warning(f"Could not read {path}: {e}")
        return None


def _write_manifest(outdir, data: dict):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = outdir / f"{MANIFEST_FILE}.tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    tmp.rename(outdir / MANIFEST_FILE)


def start_manifest(outdir, config, species_tree_path: str,
                   protein_dir: str, cds_dir: str,
                   resume: bool = False,
                   allow_config_change: bool = False) -> dict:
    """Write (or extend) the run manifest, refusing a mismatched resume.

    Raises ConfigMismatchError when resuming under a configuration, species
    tree, or input set other than the one on record, unless
    allow_config_change is set — in which case the override is recorded.
    """
    outdir = Path(outdir)
    tree_sha = sha256_file(species_tree_path)
    current_hash = config_hash(config, species_tree_path)
    inputs = {
        "protein_dir": checksum_dir(protein_dir),
        "cds_dir": checksum_dir(cds_dir),
    }
    git = git_state()
    now = _now()

    stored = load_manifest(outdir) if resume else None
    unverified = False
    diffs: list = []

    if resume and stored is None:
        unverified = True
        logger.warning(
            f"Resuming {outdir} with no run_manifest.json on record — the "
            f"configuration that produced the existing rounds cannot be "
            f"verified. This run's manifest starts from here."
        )
    elif stored is not None:
        diffs = (_config_differences(stored, config, tree_sha)
                 + _input_differences(stored, inputs))
        if diffs and not allow_config_change:
            raise ConfigMismatchError(_format_refusal(outdir, diffs))
        if diffs:
            logger.error(
                f"--allow-config-change: resuming {outdir} under a DIFFERENT "
                f"configuration ({len(diffs)} setting(s) changed). This output "
                f"directory now mixes configurations; the manifest records "
                f"which."
            )
        old_commit = (stored.get("git") or {}).get("commit")
        if old_commit and git.get("commit") and old_commit != git["commit"]:
            logger.warning(
                f"Resuming under a different code revision: "
                f"{old_commit[:12]} -> {git['commit'][:12]}"
            )

    entry = {
        "event": "resume" if stored is not None or (resume and unverified) else "start",
        "at": now,
        "git": git,
        "config_hash": current_hash,
    }
    if diffs:
        entry["config_changes"] = diffs
        entry["allow_config_change"] = True

    if stored is None:
        data = {
            "manifest_version": MANIFEST_VERSION,
            "run_id": f"{now.replace('-', '').replace(':', '')}-{uuid.uuid4().hex[:8]}",
            "created_at": now,
            "resume_count": 0,
            "config_override_used": False,
            "unverified_resume": unverified,
            "history": [],
        }
    else:
        data = dict(stored)
        data["resume_count"] = int(stored.get("resume_count", 0)) + 1

    data.update({
        "updated_at": now,
        "status": "running",
        "finished_at": None,
        "git": git,
        "config": asdict(config),
        "config_hash": current_hash,
        "hashed_settings": sorted(hashed_config_fields(config)),
        "species_tree": {"path": str(species_tree_path), "sha256": tree_sha},
        "inputs": inputs,
        "tools": tool_records(config),
    })
    data["history"] = list(data.get("history", [])) + [entry]
    if diffs:
        data["config_override_used"] = True
    data.setdefault("config_override_used", False)
    data.setdefault("unverified_resume", unverified)

    _write_manifest(outdir, data)
    logger.info(
        f"Run manifest: {outdir / MANIFEST_FILE} "
        f"(run {data['run_id']}, resume #{data['resume_count']}, "
        f"config {current_hash[7:19]})"
    )
    return data


def finish_manifest(outdir, status: str, error: Optional[str] = None) -> Optional[dict]:
    """Move the manifest to a terminal status. Never raises."""
    data = load_manifest(outdir)
    if data is None:
        return None
    data["status"] = status
    data["finished_at"] = _now()
    data["updated_at"] = data["finished_at"]
    if error is not None:
        data["error"] = error
    try:
        _write_manifest(outdir, data)
    except OSError as e:
        logger.warning(f"Could not update run manifest status: {e}")
    return data
