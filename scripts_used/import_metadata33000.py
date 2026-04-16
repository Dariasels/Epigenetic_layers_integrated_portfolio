#this script worked!

import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
import re

input_file = "combined_metadata.csv"
output_file = "clean_metadata_full.csv"

df = pd.read_csv(input_file)

def get_gsm_metadata(gsm):

    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"

    for attempt in range(5):
        try:
            r = requests.get(url, timeout=20)
            soup = BeautifulSoup(r.text, "html.parser")
            text = soup.get_text("\n")

            disease = None
            age = None
            sex = None

            # disease status
            disease_match = re.search(r'disease status:\s*([^\n]+)', text, re.IGNORECASE)
            if disease_match:
                disease = disease_match.group(1).strip()

            # age
            age_match = re.search(r'age:\s*(\d+)', text, re.IGNORECASE)
            if age_match:
                age = age_match.group(1)

            # gender
            sex_match = re.search(r'gender:\s*([^\n]+)', text, re.IGNORECASE)
            if sex_match:
                sex = sex_match.group(1).strip()

            return disease, age, sex

        except requests.exceptions.RequestException:
            print(f"Retry for {gsm}")
            time.sleep(5)

    return None, None, None


disease_list = []
age_list = []
sex_list = []

for gsm in df["sample_id"]:

    print("Fetching", gsm)

    disease, age, sex = get_gsm_metadata(gsm)

    disease_list.append(disease)
    age_list.append(age)
    sex_list.append(sex)

    time.sleep(0.3)


df["disease"] = disease_list
df["age"] = age_list
df["sex"] = sex_list

df.to_csv(output_file, index=False)

print("Saved:", output_file)

