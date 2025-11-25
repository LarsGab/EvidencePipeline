#!/usr/bin/env python3
"""Convenience wrapper for running the EvidencePipeline Nextflow workflow.

The script validates that:
  * The provided params YAML file and Nextflow config exist.
  * Required input data referenced in the params file are available.
  * Required executables are on the PATH (unless skipped).

It then launches the pipeline with the given configuration.
"""
from __future__ import annotations

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover - PyYAML is an external dependency
    raise SystemExit("PyYAML is required to run this wrapper (pip install pyyaml).")


REPO_ROOT = Path(__file__).resolve().parent
PIPELINE_MAIN = REPO_ROOT / "EvidencePipeline" / "main.nf"
BASE_CONFIG = REPO_ROOT / "conf" / "base.config"

GLOB_CHARS = set("*?[]{}")

DEFAULT_TOOL_BINARIES: Dict[str, str] = {
    "hisat2": "hisat2",
    "hisat2_build": "hisat2-build",
    "minimap2": "minimap2",
    "stringtie": "stringtie",
    "samtools": "samtools",
    "transdecoder_longorfs": "TransDecoder.LongOrfs",
    "transdecoder_predict": "TransDecoder.Predict",
    "transdecoder_util_gtf2fa": "gtf_genome_to_cdna_fasta.pl",
    "transdecoder_util_orf2genome": "cdna_alignment_orf_to_genome_orf.pl",
    "transdecoder_gtf2gff": "gtf_to_alignment_gff3.pl",
    "diamond": "diamond",
    "bedtools": "bedtools",
    "miniprot": "miniprot",
    "miniprot_boundary_scorer": "miniprot_boundary_scorer",
    "miniprothint": "miniprothint.py",
    "bam2hints": "bam2hints",
}

TOOL_DESCRIPTIONS: Dict[str, str] = {
    "hisat2": "HISAT2 aligner (params.tools.hisat2)",
    "hisat2_build": "HISAT2 indexer (params.tools.hisat2_build)",
    "minimap2": "Minimap2 aligner (params.tools.minimap2)",
    "stringtie": "StringTie assembler (params.tools.stringtie)",
    "samtools": "Samtools (params.tools.samtools)",
    "transdecoder_longorfs": "TransDecoder.LongOrfs (params.tools.transdecoder_longorfs)",
    "transdecoder_predict": "TransDecoder.Predict (params.tools.transdecoder_predict)",
    "transdecoder_util_gtf2fa": "gtf_genome_to_cdna_fasta.pl (params.tools.transdecoder_util_gtf2fa)",
    "transdecoder_util_orf2genome": "cdna_alignment_orf_to_genome_orf.pl (params.tools.transdecoder_util_orf2genome)",
    "transdecoder_gtf2gff": "gtf_to_alignment_gff3.pl (params.tools.transdecoder_gtf2gff)",
    "diamond": "DIAMOND aligner (params.tools.diamond)",
    "bedtools": "BEDTools (params.tools.bedtools)",
    "miniprot": "MiniProt aligner (params.tools.miniprot)",
    "miniprot_boundary_scorer": "MiniProt boundary scorer (params.tools.miniprot_boundary_scorer)",
    "miniprothint": "MiniProtHint (params.tools.miniprothint)",
    "bam2hints": "BAM2HINTS (params.tools.bam2hints)",
}

GENERAL_COMMANDS = {
    "nextflow": "Nextflow executable used to launch the pipeline",
    "java": "Java runtime (version 11 or newer)",
    "singularity": "Singularity/Apptainer runtime",
    "python3": "System Python 3 interpreter",
    "perl": "Perl interpreter (needed for aln2hints.pl)",
    "fasterq-dump": "SRA Toolkit fasterq-dump utility",
}

OPTIONAL_COMMANDS = {
    "prefetch": "SRA Toolkit prefetch utility (optional but recommended for SRA downloads)",
}


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate inputs/tools and run the EvidencePipeline Nextflow workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-p",
        "--params_yaml",
        required=True,
        help="Path to the params YAML file that configures the run.",
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="Path to the Nextflow config to use (cluster/local settings).",
    )
    parser.add_argument(
        "--profile",
        help="Nextflow profile(s) to activate (comma-separated).",
    )
    parser.add_argument(
        "--nextflow_bin",
        default="nextflow",
        help="Path to the Nextflow executable to use.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Pass -resume to Nextflow to continue a previous execution.",
    )
    parser.add_argument(
        "--work_dir",
        help="Optional custom Nextflow work directory (passed to -work-dir).",
    )
    parser.add_argument(
        "--check_tools",
        action="store_true",
        help="Validate native tool executables (useful when not relying on Singularity).",
    )
    parser.add_argument(
        "--skip_singularity_check",
        action="store_true",
        help="Skip Singularity executable validation.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Only run validation; do not start Nextflow.",
    )
    parser.add_argument(
        "nextflow_args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded verbatim to Nextflow (prefix them with '--').",
    )
    return parser.parse_args(argv)


def load_params(params_path: Path) -> Dict:
    with params_path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise SystemExit(f"Expected a mapping at the top-level of {params_path}, got {type(data).__name__}.")
    return data


def has_glob_char(value: str) -> bool:
    return any(char in value for char in GLOB_CHARS)


def expand_braces(pattern: str) -> List[str]:
    """Expand basic {foo,bar} brace expressions recursively."""
    start = pattern.find("{")
    if start == -1:
        return [pattern]
    end = pattern.find("}", start)
    if end == -1:
        return [pattern]
    prefix = pattern[:start]
    suffix = pattern[end + 1 :]
    choices = pattern[start + 1 : end].split(",")
    expanded = []
    for choice in choices:
        expanded.extend(expand_braces(f"{prefix}{choice}{suffix}"))
    return expanded


def iter_strings(value) -> Iterable[str]:
    if isinstance(value, (list, tuple, set)):
        for sub in value:
            yield from iter_strings(sub)
    elif isinstance(value, (str, Path)):
        yield str(value)
    elif value is None:
        return
    else:
        raise TypeError(f"Unsupported value type in params file: {type(value)}")


def resolve_data_entries(value, base_dir: Path) -> Tuple[List[Path], List[str]]:
    """Resolve params entries into concrete file paths.

    Returns (paths, errors). Errors capture patterns that matched nothing.
    """
    resolved: List[Path] = []
    errors: List[str] = []
    for raw in iter_strings(value):
        cleaned = raw.strip()
        if not cleaned:
            continue
        expanded = os.path.expandvars(os.path.expanduser(cleaned))
        candidate = Path(expanded)
        if not candidate.is_absolute():
            candidate = (base_dir / candidate).resolve()
        if has_glob_char(cleaned):
            patterns = expand_braces(str(candidate)) if "{" in cleaned else [str(candidate)]
            globbed: List[str] = []
            for pattern in patterns:
                globbed.extend(glob.glob(pattern))
            if globbed:
                resolved.extend(Path(match) for match in globbed)
            else:
                errors.append(f"No files matched pattern '{cleaned}'.")
        else:
            resolved.append(candidate)
    return resolved, errors


def validate_input_data(params: Dict, params_path: Path) -> List[str]:
    errors: List[str] = []
    base_dir = params_path.parent

    def ensure_required(key: str, label: str) -> None:
        value = params.get(key)
        if not value:
            errors.append(f"{label} ('{key}') is not set in {params_path}")
            return
        check_entries(value, label)

    def check_entries(value, label: str) -> None:
        files, errs = resolve_data_entries(value, base_dir)
        errors.extend(errs)
        if not files:
            return
        for file_path in files:
            if not file_path.exists():
                errors.append(f"{label} missing: {file_path}")
            elif not file_path.is_file():
                errors.append(f"{label} is not a file: {file_path}")

    ensure_required("genome", "Genome FASTA")
    ensure_required("proteins", "Protein FASTA")

    optional_fields = {
        "rnaseq_single": "RNA-Seq single-end FASTQ",
        "rnaseq_paired": "RNA-Seq paired-end FASTQ",
        "isoseq": "Iso-Seq FASTQ",
        "scoring_matrix": "Scoring matrix",
    }
    for key, label in optional_fields.items():
        if key in params and params.get(key):
            check_entries(params[key], label)

    tiberius = params.get("tiberius") or {}
    if isinstance(tiberius, dict) and tiberius.get("run"):
        model_cfg = tiberius.get("model_cfg")
        if not model_cfg:
            errors.append("Tiberius is enabled but params.tiberius.model_cfg is missing.")
        else:
            check_entries(model_cfg, "Tiberius model_cfg")

    return errors


def resolve_executable(command: str, relative_to: Path) -> Path | None:
    expanded = os.path.expandvars(os.path.expanduser(command))
    has_sep = os.sep in expanded or (os.altsep and os.altsep in expanded)
    if has_sep:
        candidate = Path(expanded)
        if not candidate.is_absolute():
            candidate = (relative_to / candidate).resolve()
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
        return None
    found = shutil.which(expanded)
    return Path(found).resolve() if found else None


def check_java_version(java_path: Path) -> Tuple[bool, str | None]:
    try:
        proc = subprocess.run(
            [str(java_path), "-version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:  # pragma: no cover - depends on system java
        return False, str(exc)
    output = proc.stderr or proc.stdout
    version_line = output.splitlines()[0] if output else ""
    marker = '"'
    if marker in version_line:
        version = version_line.split(marker)[1]
        major = version.split(".")[0]
        try:
            if int(major) >= 11:
                return True, None
        except ValueError:
            pass
        return False, f"Java version {version} detected, but 11+ is required."
    return False, "Unable to parse Java version output."


def build_tool_command_map(params: Dict) -> Dict[str, str]:
    overrides = params.get("tools") or {}
    command_map = DEFAULT_TOOL_BINARIES.copy()
    if isinstance(overrides, dict):
        for key, value in overrides.items():
            if key in command_map and value:
                command_map[key] = str(value)
    return command_map


def validate_executables(
    params: Dict,
    nextflow_bin: str,
    check_tool_binaries: bool,
    skip_singularity_check: bool,
) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    def check_command(command: str, description: str, mandatory: bool = True, java: bool = False) -> None:
        path = resolve_executable(command, REPO_ROOT)
        if not path:
            msg = f"{description} not found on PATH (looked for '{command}')."
            if mandatory:
                errors.append(msg)
            else:
                warnings.append(msg)
            return
        if java:
            ok, problem = check_java_version(path)
            if not ok:
                errors.append(problem or f"Unable to validate Java executable at {path}.")

    for cmd, desc in GENERAL_COMMANDS.items():
        if cmd == "nextflow":
            check_command(nextflow_bin, desc)
        elif cmd == "singularity" and skip_singularity_check:
            continue
        elif cmd == "java":
            check_command(cmd, desc, java=True)
        else:
            check_command(cmd, desc)


    if check_tool_binaries:        
        for cmd, desc in OPTIONAL_COMMANDS.items():
            check_command(cmd, desc, mandatory=False)
        tools = build_tool_command_map(params)
        for key, cmd in tools.items():
            label = TOOL_DESCRIPTIONS.get(key, f"Tool '{key}'")
            check_command(cmd, label)

        tiberius = params.get("tiberius") or {}
        if isinstance(tiberius, dict) and tiberius.get("run"):
            check_command("tiberius.py", "Tiberius CLI (tiberius.py)")

    return errors, warnings


def run_nextflow(
    params_path: Path,
    config_path: Path,
    profile: str | None,
    nextflow_bin: str,
    resume: bool,
    work_dir: str | None,
    extra_args: Sequence[str],
) -> int:
    launch_cwd = Path.cwd()
    cmd = [
        nextflow_bin,
        "run",
        str(PIPELINE_MAIN),
        "-params-file",
        str(params_path),
        "-c",
        str(BASE_CONFIG),
        "-c",
        str(config_path),
    ]
    if profile:
        cmd.extend(["-profile", profile])
    if resume:
        cmd.append("-resume")
    if work_dir:
        cmd.extend(["-work_dir", work_dir])
    if extra_args:
        cmd.extend(extra_args)

    print("[INFO] Launching Nextflow with command:")
    print("       " + " ".join(shlex.quote(part) for part in cmd))

    completed = subprocess.run(cmd, cwd=launch_cwd)
    return completed.returncode


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    extra_args = list(args.nextflow_args)
    if extra_args and extra_args[0] == "--":
        extra_args = extra_args[1:]
    args.nextflow_args = extra_args
    params_path = Path(args.params_yaml).expanduser().resolve()
    config_path = Path(args.config).expanduser().resolve()
    wdir: Path | None = None
    if args.work_dir:
        wdir = Path(args.work_dir).expanduser().resolve()

    if not params_path.exists():
        raise SystemExit(f"Params YAML not found: {params_path}")
    if not config_path.exists():
        raise SystemExit(f"Nextflow config not found: {config_path}")
    if not PIPELINE_MAIN.exists():
        raise SystemExit(f"Pipeline entry point missing: {PIPELINE_MAIN}")
    if not BASE_CONFIG.exists():
        raise SystemExit(f"Base config not found: {BASE_CONFIG}")

    params = load_params(params_path)

    print("[INFO] Validating input files...")
    data_errors = validate_input_data(params, params_path)
    if data_errors:
        for error in data_errors:
            print(f"[ERROR] {error}")
        raise SystemExit("Input validation failed.")

    print("[INFO] Validating required executables...")
    exec_errors, exec_warnings = validate_executables(
        params=params,
        nextflow_bin=args.nextflow_bin,
        check_tool_binaries=args.check_tools,
        skip_singularity_check=args.skip_singularity_check,
    )
    if exec_warnings:
        for warning in exec_warnings:
            print(f"[WARN] {warning}")
    if exec_errors:
        for error in exec_errors:
            print(f"[ERROR] {error}")
        raise SystemExit("Executable validation failed.")

    if args.dry_run:
        print("[INFO] Dry run requested; skipping Nextflow execution.")
        return

    returncode = run_nextflow(
        params_path=params_path,
        config_path=config_path,
        profile=args.profile,
        nextflow_bin=args.nextflow_bin,
        resume=args.resume,
        work_dir=wdir,
        extra_args=args.nextflow_args,
    )
    if returncode != 0:
        raise SystemExit(returncode)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except KeyboardInterrupt:  # pragma: no cover - user interruption
        raise SystemExit("\nAborted by user.")
