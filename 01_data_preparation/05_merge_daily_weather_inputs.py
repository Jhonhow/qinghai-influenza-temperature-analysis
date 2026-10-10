from pathlib import Path
import os
import json
import re
import pandas as pd
import numpy as np

ROOT = Path(os.environ.get("QH_PROJECT_ROOT", Path.cwd())).expanduser().resolve()
OLD = ROOT / "0718开始分析" / "青海省16-24年气象数据.xlsx"
TEMP_DIR = ROOT / "20260713之前" / "气温数据16-24"
PRES_DIR = ROOT / "20260713之前" / "气压数据16-24"
OUT_DIR = ROOT / "_codex_weather_merge"

ID_COLS = ["CountyC", "CityC", "ProvinceC", "County", "City", "Province"]


def load_wide(folder: Path, prefix: str, value_name: str):
    frames = []
    yearly = []
    for year in range(2016, 2025):
        file = folder / f"区县层级_{prefix}_{year}.xlsx"
        x = pd.read_excel(file, sheet_name=0)
        missing_ids = [c for c in ID_COLS if c not in x.columns]
        if missing_ids:
            raise ValueError(f"{file.name}: missing identifier columns {missing_ids}")
        date_cols = [c for c in x.columns if re.fullmatch(r"\d{4}_\d{2}_\d{2}", str(c))]
        non_date_extra = [c for c in x.columns if c not in ID_COLS and c not in date_cols]
        long = x.melt(id_vars=ID_COLS, value_vars=date_cols,
                      var_name="date_text", value_name=value_name)
        long["date"] = pd.to_datetime(long["date_text"], format="%Y_%m_%d", errors="raise")
        long = long.drop(columns="date_text")
        frames.append(long)
        yearly.append({
            "year": year,
            "file": file.name,
            "rows_wide": int(len(x)),
            "date_columns": int(len(date_cols)),
            "extra_columns": [str(c) for c in non_date_extra],
            "records_long": int(len(long)),
            "unique_counties": int(x["CountyC"].nunique()),
            "duplicate_county_codes": int(x["CountyC"].duplicated().sum()),
            "missing_values": int(long[value_name].isna().sum()),
            "min": None if long[value_name].dropna().empty else float(long[value_name].min()),
            "max": None if long[value_name].dropna().empty else float(long[value_name].max()),
        })
    out = pd.concat(frames, ignore_index=True)
    return out, yearly


old = pd.read_excel(OLD, sheet_name="Qinghai_PREC_Long")
old["date"] = pd.to_datetime(old["date"], errors="raise")

temp, temp_yearly = load_wide(TEMP_DIR, "tmpmean", "平均气温")
pres, pres_yearly = load_wide(PRES_DIR, "pres", "平均气压")

key = ["CountyC", "date"]
temp_dup = int(temp.duplicated(key).sum())
pres_dup = int(pres.duplicated(key).sum())
old_dup = int(old.duplicated(key).sum())
if temp_dup or pres_dup or old_dup:
    raise ValueError(f"Duplicate keys: old={old_dup}, temperature={temp_dup}, pressure={pres_dup}")

temp_names = temp[ID_COLS].drop_duplicates("CountyC")
pres_names = pres[ID_COLS].drop_duplicates("CountyC")
old_names = old[ID_COLS].drop_duplicates("CountyC")

merged = old.merge(temp[key + ["平均气温"]], on=key, how="left", validate="one_to_one")
merged = merged.merge(pres[key + ["平均气压"]], on=key, how="left", validate="one_to_one")

old_keys = set(map(tuple, old[key].itertuples(index=False, name=None)))
temp_keys = set(map(tuple, temp[key].itertuples(index=False, name=None)))
pres_keys = set(map(tuple, pres[key].itertuples(index=False, name=None)))

def name_mismatches(source_names):
    cmp = old_names.merge(source_names, on="CountyC", how="outer", suffixes=("_old", "_new"), indicator=True)
    mismatch_rows = []
    for _, row in cmp.iterrows():
        mismatched = row["_merge"] != "both"
        for c in ["County", "City", "Province", "CityC", "ProvinceC"]:
            if row.get(f"{c}_old") != row.get(f"{c}_new"):
                mismatched = True
        if mismatched:
            mismatch_rows.append({k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()})
    return mismatch_rows

audit = {
    "old_file": str(OLD),
    "old_rows": int(len(old)),
    "old_columns": old.columns.tolist(),
    "old_unique_counties": int(old["CountyC"].nunique()),
    "old_date_min": old["date"].min().strftime("%Y-%m-%d"),
    "old_date_max": old["date"].max().strftime("%Y-%m-%d"),
    "old_duplicate_keys": old_dup,
    "temperature_yearly": temp_yearly,
    "pressure_yearly": pres_yearly,
    "temperature_total_records": int(len(temp)),
    "pressure_total_records": int(len(pres)),
    "temperature_duplicate_keys": temp_dup,
    "pressure_duplicate_keys": pres_dup,
    "temperature_name_mismatches": name_mismatches(temp_names),
    "pressure_name_mismatches": name_mismatches(pres_names),
    "old_keys_missing_temperature": int(len(old_keys - temp_keys)),
    "old_keys_missing_pressure": int(len(old_keys - pres_keys)),
    "temperature_keys_not_in_old": int(len(temp_keys - old_keys)),
    "pressure_keys_not_in_old": int(len(pres_keys - old_keys)),
    "merged_missing_temperature": int(merged["平均气温"].isna().sum()),
    "merged_missing_pressure": int(merged["平均气压"].isna().sum()),
    "temperature_summary": {
        "min": float(merged["平均气温"].min()),
        "p01": float(merged["平均气温"].quantile(0.01)),
        "median": float(merged["平均气温"].median()),
        "p99": float(merged["平均气温"].quantile(0.99)),
        "max": float(merged["平均气温"].max()),
    },
    "pressure_summary": {
        "min": float(merged["平均气压"].min()),
        "p01": float(merged["平均气压"].quantile(0.01)),
        "median": float(merged["平均气压"].median()),
        "p99": float(merged["平均气压"].quantile(0.99)),
        "max": float(merged["平均气压"].max()),
    },
}

OUT_DIR.mkdir(parents=True, exist_ok=True)
with open(OUT_DIR / "merge_audit.json", "w", encoding="utf-8") as f:
    json.dump(audit, f, ensure_ascii=False, indent=2)

# Preserve exact original row order. JSON nulls are not expected after validation.
payload = {
    "headers": ["平均气温", "平均气压"],
    "values": merged[["平均气温", "平均气压"]].where(pd.notna(merged), None).values.tolist(),
}
with open(OUT_DIR / "merge_values.json", "w", encoding="utf-8") as f:
    json.dump(payload, f, ensure_ascii=False, separators=(",", ":"))

merged_export = merged.copy()
merged_export["date"] = (merged_export["date"] - pd.Timestamp("1899-12-30")).dt.days
merged_export.to_csv(
    OUT_DIR / "merged_full.csv",
    index=False,
    encoding="utf-8",
    date_format="%Y-%m-%d",
)

print(json.dumps(audit, ensure_ascii=False, indent=2))
