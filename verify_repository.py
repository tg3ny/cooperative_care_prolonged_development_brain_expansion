#!/usr/bin/env python3
"""Validate repository integrity and principal manuscript results."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent
DATA = ROOT


EXPECTED_SHAPES = {
    "Data_S1_Summary.csv": (200, 7),
    "Data_S2_Sensitivity.csv": (81, 5),
    "Data_S3_Trajectories.csv.gz": (200_000, 10),
    "Data_S4_Alternative_Models.csv": (5, 16),
    "Data_S5_Equilibrium_Statistics.csv": (8, 7),
    "Data_S6_Pure_Summary.csv": (200, 14),
    "Data_S7_Mixed_Summary.csv": (200, 14),
    "experiment2_control_trajectories.csv.gz": (200_000, 14),
    "experiment2_treatment_trajectories.csv.gz": (200_000, 14),
    "gradual_ceiling_control_trajectories.csv.gz": (200_000, 14),
    "gradual_ceiling_treatment_trajectories.csv.gz": (200_000, 14),
}


def fail(message: str) -> None:
    raise AssertionError(message)


def verify_hashes() -> None:
    manifest = ROOT / "checksums" / "SHA256SUMS"
    if not manifest.exists():
        fail("Missing checksums/SHA256SUMS")
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, rel = line.split("  ", 1)
        path = ROOT / rel
        if not path.is_file():
            fail(f"Manifest file missing: {rel}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected:
            fail(f"Checksum mismatch: {rel}")


def verify_shapes() -> None:
    for name, expected in EXPECTED_SHAPES.items():
        frame = pd.read_csv(DATA / name)
        if frame.shape != expected:
            fail(f"Unexpected shape for {name}: {frame.shape}, expected {expected}")


def verify_groups_and_metadata() -> None:
    for name in ("Data_S1_Summary.csv", "Data_S6_Pure_Summary.csv", "Data_S7_Mixed_Summary.csv"):
        frame = pd.read_csv(DATA / name)
        counts = frame.groupby("group")["replicate"].nunique().to_dict()
        if counts != {"control": 100, "treatment": 100}:
            fail(f"Unexpected replicate counts in {name}: {counts}")

    split_groups = {
        "experiment2_control_trajectories.csv.gz": "control",
        "experiment2_treatment_trajectories.csv.gz": "treatment",
        "gradual_ceiling_control_trajectories.csv.gz": "control",
        "gradual_ceiling_treatment_trajectories.csv.gz": "treatment",
    }
    for name, expected_group in split_groups.items():
        frame = pd.read_csv(DATA / name, usecols=["group", "replicate"])
        groups = set(frame["group"])
        replicates = frame["replicate"].nunique()
        if groups != {expected_group} or replicates != 100:
            fail(
                f"Unexpected split-file contents in {name}: "
                f"groups={groups}, replicates={replicates}"
            )

    for name in ("Data_S8_Pure_Results.json", "Data_S9_Mixed_Results.json", "Data_S10_Ramp_Results.json"):
        payload = json.loads((DATA / name).read_text(encoding="utf-8"))
        if payload.get("master_seed") != 42:
            fail(f"Unexpected master seed in {name}")
        if payload.get("n_replicates_per_group") != 100:
            fail(f"Unexpected replicate metadata in {name}")


def comparison(frame: pd.DataFrame, prefix: str, trait: str) -> tuple[float, float, float, float]:
    control = frame.loc[frame.group == "control", f"{prefix}_{trait}"].to_numpy()
    treatment = frame.loc[frame.group == "treatment", f"{prefix}_{trait}"].to_numpy()
    test = stats.ttest_ind(treatment, control, equal_var=True)
    return control.mean(), treatment.mean(), float(test.statistic), float(test.pvalue)


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not np.isclose(actual, expected, atol=tolerance, rtol=0):
        fail(f"{label}: {actual} does not match {expected} ± {tolerance}")


def verify_results() -> None:
    exp1 = pd.read_csv(DATA / "Data_S1_Summary.csv")
    c, t, stat, p = comparison(exp1, "phase2", "L")
    close(c, 6.4136, 0.0001, "Simulation 1 control childhood mean")
    close(t, 8.5012, 0.0001, "Simulation 1 treatment childhood mean")
    close(stat, 8.5479, 0.0001, "Simulation 1 childhood t statistic")
    close(p, 3.3682854890025566e-15, 1e-25, "Simulation 1 two-sided childhood p value")
    if p >= 0.001:
        fail("Simulation 1 childhood comparison is not p < 0.001")

    exp2 = pd.read_csv(DATA / "Data_S6_Pure_Summary.csv")
    c, t, stat, p = comparison(exp2, "phase3", "B")
    close(c, 41.9342, 0.0001, "Simulation 2 control brain mean")
    close(t, 55.6594, 0.0001, "Simulation 2 treatment brain mean")
    delta_control = exp2.loc[exp2.group == "control", "delta_B"].to_numpy()
    delta_treatment = exp2.loc[exp2.group == "treatment", "delta_B"].to_numpy()
    delta_test = stats.ttest_ind(delta_treatment, delta_control, equal_var=True)
    close(float(delta_test.pvalue), 4.908025082303922e-31, 1e-40, "Simulation 2 two-sided expansion p value")
    if p >= 0.001:
        fail("Simulation 2 brain comparison is not p < 0.001")

    exp3 = pd.read_csv(DATA / "Data_S7_Mixed_Summary.csv")
    c, t, stat, p = comparison(exp3, "phase3", "L")
    close(c, 7.7992, 0.0001, "Simulation 3 care-independent-route childhood mean")
    close(t, 9.8706, 0.0001, "Simulation 3 cooperative-care-route childhood mean")
    close(stat, 5.5236, 0.0001, "Simulation 3 childhood t statistic")
    close(p, 1.0360986557087542e-07, 1e-16, "Simulation 3 two-sided childhood p value")
    if p >= 0.001:
        fail("Simulation 3 childhood comparison is not p < 0.001")

    json_expectations = {
        "Data_S8_Pure_Results.json": 4.908025082303922e-31,
        "Data_S9_Mixed_Results.json": 0.8458526640924088,
        "Data_S10_Ramp_Results.json": 1.6990488241466804e-25,
    }
    for filename, expected in json_expectations.items():
        payload = json.loads((DATA / filename).read_text(encoding="utf-8"))
        close(payload["analysis"]["dB_p_value"], expected, max(abs(expected) * 1e-12, 1e-40), f"{filename} two-sided expansion p value")


def main() -> int:
    verify_hashes()
    verify_shapes()
    verify_groups_and_metadata()
    verify_results()
    print("Repository validation passed.")
    print(f"Validated {len(EXPECTED_SHAPES)} tabular datasets and 3 JSON result files.")
    print("Replicate counts, master seed, integrity hashes, and principal results are correct.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print(f"VALIDATION FAILED: {error}", file=sys.stderr)
        sys.exit(1)
