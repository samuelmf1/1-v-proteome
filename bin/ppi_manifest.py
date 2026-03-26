#!/usr/bin/env python3

import argparse
import csv
import json
import os
from pathlib import Path
import sys
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a manifest of AlphaFold PPI structures from summary_confidences.json "
            "and model_interfaces.csv files."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="ppi",
        help="Root directory to scan. Default: ppi",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output CSV path. Use '-' to write to stdout. Default: -",
    )
    return parser.parse_args()


def find_summary_files(root: Path) -> list[Path]:
    return sorted(
        path
        for path in root.glob("**/*_summary_confidences.json")
        if not any(part.startswith("seed-") for part in path.parts)
    )


def count_interaction_sites(interfaces_path: Path) -> int:
    if not interfaces_path.exists():
        return 0

    with interfaces_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        return sum(1 for _ in reader)


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def summarize_interface_metrics(interfaces_path: Path) -> dict[str, float | int | str]:
    if not interfaces_path.exists():
        return {}

    numeric_values: dict[str, list[float]] = {}

    with interfaces_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return {}

        for row in reader:
            for column in reader.fieldnames:
                parsed_value = parse_float(row.get(column, ""))
                if parsed_value is None:
                    continue
                numeric_values.setdefault(column, []).append(parsed_value)

    summary: dict[str, float | int | str] = {}
    for column in sorted(numeric_values):
        values = numeric_values[column]
        summary[f"interface_{column}_min"] = min(values)
        summary[f"interface_{column}_max"] = max(values)

    return summary


def list_labels(value: list[object]) -> list[str]:
    if len(value) == 2:
        return ["A", "B"]
    return [str(index) for index in range(1, len(value) + 1)]


def flatten_metric(prefix: str, value: object) -> dict[str, object]:
    if isinstance(value, list):
        flattened: dict[str, object] = {}
        for label, item in zip(list_labels(value), value):
            flattened.update(flatten_metric(f"{prefix}_{label}", item))
        return flattened

    if isinstance(value, dict):
        flattened = {}
        for key in sorted(value):
            flattened.update(flatten_metric(f"{prefix}_{key}", value[key]))
        return flattened

    return {prefix: value}


def split_structure_name(structure_name: str) -> tuple[str, str]:
    target_protein, separator, interacting_protein = structure_name.partition("__")
    if not separator:
        return structure_name, ""
    return target_protein, interacting_protein


def build_record(summary_path: Path) -> dict[str, object]:
    with summary_path.open(encoding="utf-8") as handle:
        metrics = json.load(handle)

    structure_name = summary_path.name.removesuffix("_summary_confidences.json")
    target_protein, interacting_protein = split_structure_name(structure_name)
    interfaces_path = summary_path.with_name(f"{structure_name}_model_interfaces.csv")
    interface_summary = summarize_interface_metrics(interfaces_path)

    record: dict[str, object] = {
        "structure_name": structure_name,
        "target_protein": target_protein,
        "interacting_protein": interacting_protein,
        "interaction_site_count": count_interaction_sites(interfaces_path),
        "structure_dir": str(summary_path.parent),
    }

    for key in sorted(metrics):
        record.update(flatten_metric(key, metrics[key]))

    record.update(interface_summary)

    return record


def add_other_target_proteins(records: list[dict[str, object]]) -> None:
    targets_by_interactor: dict[str, set[str]] = {}

    for record in records:
        interacting_protein = str(record.get("interacting_protein", ""))
        target_protein = str(record.get("target_protein", ""))
        if not interacting_protein or not target_protein:
            continue
        targets_by_interactor.setdefault(interacting_protein, set()).add(target_protein)

    for record in records:
        interacting_protein = str(record.get("interacting_protein", ""))
        target_protein = str(record.get("target_protein", ""))
        other_targets = sorted(targets_by_interactor.get(interacting_protein, set()) - {target_protein})
        record["other_target_proteins"] = ";".join(other_targets)


def collect_records(summary_files: Iterable[Path]) -> list[dict[str, object]]:
    records = [build_record(summary_path) for summary_path in summary_files]
    add_other_target_proteins(records)
    return records


def write_manifest(records: list[dict[str, object]], output_path: str) -> None:
    base_fields = [
        "structure_name",
        "target_protein",
        "interacting_protein",
        "other_target_proteins",
        "interaction_site_count",
    ]
    path_fields = [
        "structure_dir",
    ]

    metric_fields = sorted(
        {
            key
            for record in records
            for key in record
            if key not in set(base_fields + path_fields)
        }
    )
    fieldnames = base_fields + metric_fields + path_fields

    if output_path == "-":
        import sys

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
        return

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    summary_files = find_summary_files(root)

    if not summary_files:
        raise SystemExit(f"No summary_confidences.json files found under {root}")

    records = collect_records(summary_files)
    write_manifest(records, args.output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
        raise SystemExit(141)