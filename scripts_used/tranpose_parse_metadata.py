import pandas as pd
import glob
import os

# Since you are running the script FROM the metadata folder:
metadata_dir = "." 

# Find all files ending in series_matrix.txt in the current folder
files = glob.glob(os.path.join(metadata_dir, "*series_matrix.txt"))

print(f"Found {len(files)} files: {files}")

all_samples = []

for filepath in files:
    # Use a DIFFERENT variable name for the ID so you don't overwrite 'filepath'
    dataset_id = os.path.basename(filepath).split("_")[0]
    
    print(f"Processing {dataset_id}...")

    rows = []
    # Use the full filepath variable here
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith("!Sample_"):
                parts = line.strip().split("\t")
                key = parts[0].replace("!Sample_", "")
                values = [v.strip('"') for v in parts[1:]]
                rows.append([key] + values)

    if len(rows) == 0:
        continue

    df = pd.DataFrame(rows)

    # first column = metadata field
    df = df.set_index(0)

    # transpose so rows = samples
    df = df.T

    df["dataset"] = dataset_id

    all_samples.append(df)


# merge all datasets
meta = pd.concat(all_samples, ignore_index=True)


# ------------------------------
# extract useful fields
# ------------------------------

def extract_field(series, keyword):
    return series.str.extract(f"{keyword}: ([^;]+)", expand=False)


# combine characteristics columns
char_cols = [c for c in meta.columns if "characteristics" in c]

meta["characteristics"] = meta[char_cols].astype(str).agg(";".join, axis=1)


# extract disease info
meta["disease"] = extract_field(meta["characteristics"], "disease state|disease")
meta["health_state"] = extract_field(meta["characteristics"], "health state|group")


# extract brain region
meta["brain_region"] = extract_field(meta["characteristics"], "brain")
meta["tissue"] = extract_field(meta["characteristics"], "tissue")


# extract sex
meta["sex"] = extract_field(meta["characteristics"], "sex|Sex")


# extract age
meta["age"] = extract_field(meta["characteristics"], "age|Age")


# rename GSM column if present
if "geo_accession" in meta.columns:
    meta = meta.rename(columns={"geo_accession": "sample_id"})


# keep useful columns
cols = [
    "sample_id",
    "dataset",
    "organism_ch1",
    "source_name_ch1",
    "disease",
    "health_state",
    "brain_region",
    "tissue",
    "sex",
    "age"
]

cols = [c for c in cols if c in meta.columns]

meta = meta[cols]


# save result
output = os.path.join(metadata_dir, "combined_metadata.csv")
meta.to_csv(output, index=False)


print("Metadata saved to:", output)
print(meta.head())