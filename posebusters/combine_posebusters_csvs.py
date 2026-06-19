#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


NON_BOOL_COLUMNS = {"CCD_ID", "PDB_ID", "file", "molecule", "position"}

CSV_NAMES = {
    "comparison_report.csv",
    "pred_posebusters_report.csv",
    "posebusters_report.csv",
}


def extract_ids(csv_path: Path):
    pdb_id = csv_path.parent.name
    ccd_folder = csv_path.parent.parent.name

    if ccd_folder.endswith("_ok"):
        ccd_id = ccd_folder[:-3]
    else:
        ccd_id = ccd_folder

    return ccd_id, pdb_id


def normalize_bool_like_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if col in NON_BOOL_COLUMNS:
            continue

        non_null = df[col].dropna()
        if non_null.empty:
            continue

        normalized_values = {str(v).strip().lower() for v in non_null.unique()}

        if normalized_values <= {"true", "false", "1", "0"}:
            df[col] = df[col].map(
                lambda v: v if pd.isna(v)
                else ("TRUE" if str(v).strip().lower() in {"true", "1"} else "FALSE")
            )

    return df


def is_valid_csv(p: Path) -> bool:
    return (
        p.is_file()
        and p.name in CSV_NAMES
        and not p.name.startswith("._")
    )


def combine_csvs(base_dir: Path, output_file: Path):

    csv_files = [
        p for p in base_dir.rglob("*")
        if is_valid_csv(p)
    ]

    csv_files = sorted(csv_files)

    if not csv_files:
        print("[ERROR] No PoseBusters CSV files found.")
        return

    print(f"[INFO] Found {len(csv_files)} CSV files.")

    all_rows = []
    seen = {}

    for csv_path in csv_files:
        try:
            ccd_id, pdb_id = extract_ids(csv_path)

            key = (ccd_id, pdb_id)
            seen[key] = seen.get(key, 0) + 1

            df = pd.read_csv(csv_path)

            df.insert(0, "PDB_ID", pdb_id)
            df.insert(0, "CCD_ID", ccd_id)

            df = normalize_bool_like_columns(df)

            all_rows.append(df)

            print(f"[OK] Added: {ccd_id} / {pdb_id}")

        except Exception as e:
            print(f"[WARN] Skipped {csv_path}: {e}")

    # report duplicates (this is likely your +2 issue)
    print("\n[CHECK] duplicate CSVs per structure:")
    for k, v in seen.items():
        if v > 1:
            print(f"  DUPLICATE: {k} -> {v} files")

    combined = pd.concat(all_rows, ignore_index=True)

    combined = combined.sort_values(
        by=["CCD_ID", "PDB_ID"],
        ascending=[True, True]
    ).reset_index(drop=True)

    combined.to_csv(output_file, index=False)

    print("\n[DONE]")
    print(f"Total rows: {len(combined)}")
    print(f"Output: {output_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", nargs="?", default="posebusters_outputs")
    parser.add_argument("--output", default="combined_posebusters.csv")

    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    output_file = Path(args.output).resolve()

    print("[INFO] scanning:", base_dir)
    combine_csvs(base_dir, output_file)


if __name__ == "__main__":
    main()
