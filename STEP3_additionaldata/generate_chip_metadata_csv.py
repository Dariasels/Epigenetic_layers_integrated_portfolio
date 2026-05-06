#!/usr/bin/env python3
"""
Generate clean metadata CSV for GSE102538 from downloaded GSM files.
Writes: data/chip_GSE102538 /clean_metadata_GSE102538.csv
"""
import os
import csv

dir_path = os.path.join("data", "chip_GSE102538 ")
out_path = os.path.join(dir_path, "clean_metadata_GSE102538.csv")

files = sorted([f for f in os.listdir(dir_path) if f.endswith('.txt.gz')])
rows = []
for fn in files:
    sample_id = fn.split('_')[0]
    # guess condition from filename and normalise to project-wide labels
    cond = 'unknown'
    if '_AD' in fn.upper() or '_AD' in fn:
        cond = 'Alzheimer'
    elif '_CONTROL' in fn.upper() or '_Control' in fn:
        cond = 'Control'

    # tissue/organism
    dataset = 'GSE102538'
    organism = 'Homo sapiens'
    tissue = 'entorhinal cortex'

    # sex/age unknown here
    sex = ''
    age = ''

    rows.append({
        'sample_id': sample_id,
        'dataset': dataset,
        'organism_ch1': organism,
        'source_name_ch1': tissue,
        'sex': sex,
        'age': age,
        'condition': cond,
    })

# Write CSV with header expected by 02a_import_metadata.py (COLUMN_MAP keys)
fieldnames = ['sample_id','dataset','organism_ch1','source_name_ch1','sex','age','condition']
with open(out_path, 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    for r in rows:
        writer.writerow(r)

print(f"Wrote {len(rows)} rows to {out_path}")
