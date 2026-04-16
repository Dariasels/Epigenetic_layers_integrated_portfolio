import os
import requests

# destination directory
output_dir = "/home/daria/Documents/UGent/databases/multiomics/chatgpt_version/multiomics_data/metadata"
os.makedirs(output_dir, exist_ok=True)

# GEO datasets
datasets = [
    "GSE203206",
    "GSE129040",
    "GSE59685",
    "GSE105194",
    "GSE33000"
]

def build_geo_url(gse):
    """
    GEO stores series_matrix files in:
    ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEXXXnnn/GSEXXXXXX/matrix/
    """
    prefix = gse[:-3] + "nnn"
    filename = f"{gse}_series_matrix.txt.gz"

    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{filename}"
    return url, filename


for gse in datasets:
    url, filename = build_geo_url(gse)
    filepath = os.path.join(output_dir, filename)

    print(f"Downloading {gse} metadata...")

    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open(filepath, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Saved to {filepath}")
    else:
        print(f"Failed to download {gse}")

print("All downloads attempted.")

