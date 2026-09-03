#!/usr/bin/env python3
"""Validate the Phase 2 method, runner, and environment registries.

The environment registry keeps historical evidence, import smoke results, and
public-specification readiness as separate contracts.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFESTS = ROOT / "manifests"

EXPECTED_CONFIG_IDS = [f"C{i:02d}" for i in range(1, 21)]
EXPECTED_CONFIG_TO_OPERATIONAL = {
    "C01": "scvelo_deterministic",
    "C02": "scvelo_stochastic",
    "C03": "scvelo_dynamic",
    "C04": "multivelo",
    "C05": "unitvelo",
    "C06": "unitvelo_mode2",
    "C07": "dynamo",
    "C08": "TIvelo",
    "C09": "cell2fate",
    "C10": "TFvelo",
    "C11": "veloae",
    "C12": "veloVI",
    "C13": "DeepKINET",
    "C14": "TSvelo",
    "C15": "celldancer",
    "C16": "deepvelo",
    "C17": "SDEvelo",
    "C18": "deepvelo2022",
    "C19": "graphvelo_scv",
    "C20": "graphvelo_dynamo",
}
EXPECTED_OPERATIONAL_IDS = set(EXPECTED_CONFIG_TO_OPERATIONAL.values())

CONFIG_COLUMNS = {
    "config_id",
    "operational_id",
    "manuscript_label",
    "base_method",
    "model_type",
    "manuscript_version",
    "execution_device",
    "deviation_from_default",
    "implementation_source",
    "runner_id",
    "environment_id",
    "result_path_pattern",
    "velocity_layer",
    "source_embedding",
    "velocity_embedding",
    "release_strategy",
    "evidence_status",
    "notes",
}
RUNNER_COLUMNS = {
    "runner_id",
    "config_id",
    "candidate_runner",
    "supplemental_runners",
    "input_contract",
    "output_contract",
    "dataset_coverage",
    "result_path_pattern",
    "path_evidence",
    "evidence_status",
    "release_action",
    "notes",
}
ENVIRONMENT_COLUMNS = {
    "environment_id",
    "config_ids",
    "audited_conda_path",
    "assignment_status",
    "python_version",
    "scanpy_version",
    "scvelo_version",
    "scvi_tools_version",
    "method_package",
    "manuscript_version",
    "local_spec_source",
    "upstream_repository",
    "upstream_revision",
    "repository_status",
    "license_evidence",
    "release_action",
    "evidence_status",
    "notes",
    "runtime_status",
    "release_spec_status",
}

ALLOWED_CONFIG_RELEASE = {"curate_thin_runner", "recover_runner"}
ALLOWED_CONFIG_STATUS = {"confirmed", "conflicted"}
ALLOWED_RUNNER_RELEASE = {"curate_after_review", "recover_before_release"}
ALLOWED_RUNNER_STATUS = {"candidate", "partial", "ambiguous", "missing"}
ALLOWED_ENV_ASSIGNMENT = {
    "inferred",
    "ambiguous",
    "observed_name_only",
    "conflicted",
    "partial",
    "inferred_shared",
}
ALLOWED_ENV_RELEASE = {"create_reviewed_spec_later", "resolve_evidence_first"}
ALLOWED_ENV_STATUS = {"partial", "ambiguous", "conflicted", "missing"}
ALLOWED_ENV_RUNTIME = {
    "import_smoke_passed",
    "import_smoke_passed_with_warning",
    "import_smoke_failed",
}
ALLOWED_ENV_SPEC = {"not_yet_curated", "curation_required", "reviewed_clean"}


class ValidationError(RuntimeError):
    pass


def read_tsv(name: str, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    path = MANIFESTS / name
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        missing = required.difference(columns)
        extra = set(columns).difference(required)
        if missing or extra:
            raise ValidationError(
                f"{name}: schema mismatch; missing={sorted(missing)}, extra={sorted(extra)}"
            )
        rows = list(reader)
    for index, row in enumerate(rows, start=2):
        if None in row:
            raise ValidationError(f"{name}:{index}: non-rectangular row")
        empty_columns = [column for column in columns if row[column] == ""]
        if empty_columns:
            raise ValidationError(
                f"{name}:{index}: empty fields are not allowed: {empty_columns}"
            )
    return columns, rows


def require_unique(rows: list[dict[str, str]], field: str, name: str) -> None:
    values = [row[field] for row in rows]
    duplicates = sorted({value for value in values if values.count(value) > 1})
    if duplicates:
        raise ValidationError(f"{name}: duplicate {field}: {duplicates}")


def check_vocab(
    rows: list[dict[str, str]], field: str, allowed: set[str], name: str
) -> None:
    invalid = sorted({row[field] for row in rows}.difference(allowed))
    if invalid:
        raise ValidationError(f"{name}: invalid {field}: {invalid}")


def check_absolute_paths(
    name: str, rows: list[dict[str, str]], allowed_field: str | None = None
) -> None:
    windows_absolute = re.compile(r"^[A-Za-z]:[\\/]")
    for index, row in enumerate(rows, start=2):
        for field, value in row.items():
            if field == allowed_field:
                continue
            if value.startswith("/") or windows_absolute.match(value):
                raise ValidationError(
                    f"{name}:{index}: absolute path is not allowed in {field}: {value}"
                )
            if "/mnt/second19T/" in value or "/Users/" in value:
                raise ValidationError(
                    f"{name}:{index}: embedded private path in {field}: {value}"
                )


def main() -> int:
    _, configs = read_tsv("method_configurations.tsv", CONFIG_COLUMNS)
    _, runners = read_tsv("runner_inventory.tsv", RUNNER_COLUMNS)
    _, environments = read_tsv("environments.tsv", ENVIRONMENT_COLUMNS)

    require_unique(configs, "config_id", "method_configurations.tsv")
    require_unique(configs, "operational_id", "method_configurations.tsv")
    require_unique(configs, "runner_id", "method_configurations.tsv")
    require_unique(runners, "runner_id", "runner_inventory.tsv")
    require_unique(runners, "config_id", "runner_inventory.tsv")
    require_unique(environments, "environment_id", "environments.tsv")

    config_ids = [row["config_id"] for row in configs]
    if config_ids != EXPECTED_CONFIG_IDS:
        raise ValidationError(
            "method_configurations.tsv: config_id order/content must be exactly C01-C20"
        )
    operational_ids = {row["operational_id"] for row in configs}
    if operational_ids != EXPECTED_OPERATIONAL_IDS:
        raise ValidationError(
            "method_configurations.tsv: operational IDs differ from the 20-method result contract; "
            f"missing={sorted(EXPECTED_OPERATIONAL_IDS - operational_ids)}, "
            f"extra={sorted(operational_ids - EXPECTED_OPERATIONAL_IDS)}"
        )
    for row in configs:
        expected = EXPECTED_CONFIG_TO_OPERATIONAL[row["config_id"]]
        if row["operational_id"] != expected:
            raise ValidationError(
                f"{row['config_id']}: expected operational_id={expected}, "
                f"found {row['operational_id']}"
            )

    config_id_set = set(config_ids)
    runner_ids = {row["runner_id"] for row in runners}
    if {row["config_id"] for row in runners} != config_id_set:
        raise ValidationError("runner_inventory.tsv: expected exactly one row per config_id")
    missing_runner_refs = sorted(
        {row["runner_id"] for row in configs}.difference(runner_ids)
    )
    if missing_runner_refs:
        raise ValidationError(f"method configuration references unknown runners: {missing_runner_refs}")
    runner_by_config = {row["config_id"]: row["runner_id"] for row in runners}
    for row in configs:
        if runner_by_config[row["config_id"]] != row["runner_id"]:
            raise ValidationError(
                f"{row['config_id']}: runner_id differs between configuration and runner manifests"
            )
    runner_rows_by_config = {row["config_id"]: row for row in runners}
    for row in configs:
        runner = runner_rows_by_config[row["config_id"]]
        if runner["result_path_pattern"] != row["result_path_pattern"]:
            raise ValidationError(
                f"{row['config_id']}: result_path_pattern differs between manifests"
            )
        if runner["candidate_runner"] == "not_found":
            if runner["evidence_status"] != "missing" or row["release_strategy"] != "recover_runner":
                raise ValidationError(
                    f"{row['config_id']}: a missing runner must use missing/recover_runner states"
                )

    special_runner_checks = {
        "C10": {
            "candidate_prefix": "project:TFvelo-main/",
            "required_supplemental": "project:TFvelo-main/TFvelo_analysis_demo.py",
            "status": "partial",
            "release_action": "curate_after_review",
        },
        "C12": {
            "candidate": "velocity:methods/veloVI/veloVI_main_ori.py",
            "required_supplemental": "velocity:methods/veloVI/veloVI_main.py",
            "status": "partial",
            "release_action": "curate_after_review",
        },
        "C14": {
            "candidate": "velocity:methods/TSvelo/TSvelo_run.py",
            "status": "partial",
            "release_action": "curate_after_review",
        },
    }
    for config_id, expected in special_runner_checks.items():
        runner = runner_rows_by_config[config_id]
        if "candidate" in expected and runner["candidate_runner"] != expected["candidate"]:
            raise ValidationError(
                f"{config_id}: candidate_runner does not match the audited entry point"
            )
        if "candidate_prefix" in expected and not runner["candidate_runner"].startswith(
            expected["candidate_prefix"]
        ):
            raise ValidationError(
                f"{config_id}: candidate_runner must use the project: provenance locator"
            )
        if "required_supplemental" in expected and expected["required_supplemental"] not in runner[
            "supplemental_runners"
        ]:
            raise ValidationError(
                f"{config_id}: required supplemental runner is not recorded"
            )
        if runner["evidence_status"] != expected["status"]:
            raise ValidationError(
                f"{config_id}: runner evidence_status must be {expected['status']}"
            )
        if runner["release_action"] != expected["release_action"]:
            raise ValidationError(
                f"{config_id}: runner release_action must be {expected['release_action']}"
            )

    environment_ids = {row["environment_id"] for row in environments}
    require_unique(environments, "audited_conda_path", "environments.tsv")
    missing_environment_refs = sorted(
        {row["environment_id"] for row in configs}.difference(environment_ids)
    )
    if missing_environment_refs:
        raise ValidationError(
            f"method configuration references unknown environments: {missing_environment_refs}"
        )

    environment_membership: dict[str, str] = {}
    for row in environments:
        members = row["config_ids"].split(";")
        unknown = sorted(set(members).difference(config_id_set))
        if unknown:
            raise ValidationError(
                f"{row['environment_id']}: unknown config_ids: {unknown}"
            )
        for member in members:
            if member in environment_membership:
                raise ValidationError(
                    f"{member}: assigned by multiple environment rows"
                )
            environment_membership[member] = row["environment_id"]
    if set(environment_membership) != config_id_set:
        raise ValidationError("environments.tsv: every configuration must appear exactly once")
    for row in configs:
        if environment_membership[row["config_id"]] != row["environment_id"]:
            raise ValidationError(
                f"{row['config_id']}: environment_id differs between manifests"
            )
    expected_conda_prefix = "/mnt/second19T/zhaozelin/anaconda3/envs/"
    for row in environments:
        if not row["audited_conda_path"].startswith(expected_conda_prefix):
            raise ValidationError(
                f"{row['environment_id']}: audited_conda_path is outside the audited env root"
            )

    check_vocab(
        configs, "release_strategy", ALLOWED_CONFIG_RELEASE, "method_configurations.tsv"
    )
    check_vocab(
        configs, "evidence_status", ALLOWED_CONFIG_STATUS, "method_configurations.tsv"
    )
    check_vocab(
        runners, "release_action", ALLOWED_RUNNER_RELEASE, "runner_inventory.tsv"
    )
    check_vocab(
        runners, "evidence_status", ALLOWED_RUNNER_STATUS, "runner_inventory.tsv"
    )
    check_vocab(
        environments, "assignment_status", ALLOWED_ENV_ASSIGNMENT, "environments.tsv"
    )
    check_vocab(
        environments, "release_action", ALLOWED_ENV_RELEASE, "environments.tsv"
    )
    check_vocab(
        environments, "evidence_status", ALLOWED_ENV_STATUS, "environments.tsv"
    )
    check_vocab(
        environments, "runtime_status", ALLOWED_ENV_RUNTIME, "environments.tsv"
    )
    check_vocab(
        environments, "release_spec_status", ALLOWED_ENV_SPEC, "environments.tsv"
    )

    environment_by_id = {row["environment_id"]: row for row in environments}
    sdevelo = environment_by_id.get("env_sdevelo")
    if sdevelo is None:
        raise ValidationError("environments.tsv: env_sdevelo row is required")
    if not sdevelo["audited_conda_path"].endswith("/sdevelo_py310"):
        raise ValidationError(
            "environments.tsv: C17 must use the audited sdevelo_py310 candidate"
        )
    if sdevelo["runtime_status"] != "import_smoke_passed":
        raise ValidationError(
            "environments.tsv: sdevelo_py310 must be marked import_smoke_passed"
        )

    check_absolute_paths("method_configurations.tsv", configs)
    check_absolute_paths("runner_inventory.tsv", runners)
    check_absolute_paths(
        "environments.tsv", environments, allowed_field="audited_conda_path"
    )

    print(
        "Phase 2 manifests valid: "
        f"{len(configs)} configurations, {len(runners)} runner rows, "
        f"{len(environments)} audited environments."
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
