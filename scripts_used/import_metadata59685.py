#script worked
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

            condition = None
            age = None
            sex = None

            # -------------------------
            # GSE59685 format
            # ad.disease.status: C / AD
            # -------------------------
            match = re.search(r'ad\.disease\.status:\s*([A-Za-z]+)', text, re.IGNORECASE)
            if match:
                val = match.group(1).strip()

                if val == "C":
                    condition = "Control"
                elif val == "AD":
                    condition = "Alzheimer"
                else:
                    condition = val

            # -------------------------
            # Generic disease field
            # -------------------------
            if not condition:
                disease_match = re.search(r'disease status:\s*([^\n]+)', text, re.IGNORECASE)
                if disease_match:
                    condition = disease_match.group(1).strip()

            # -------------------------
            # Age
            # -------------------------
            age_match = re.search(r'age:\s*(\d+)', text, re.IGNORECASE)
            if age_match:
                age = age_match.group(1)

            # -------------------------
            # Sex / Gender
            # -------------------------
            sex_match = re.search(r'(gender|sex):\s*([^\n]+)', text, re.IGNORECASE)
            if sex_match:
                sex = sex_match.group(2).strip()

            return condition, age, sex

        except requests.exceptions.RequestException:
            print(f"Retry for {gsm}")
            time.sleep(5)

    return None, None, None


condition_list = []
age_list = []
sex_list = []

for gsm in df["sample_id"]:

    print("Fetching", gsm)

    condition, age, sex = get_gsm_metadata(gsm)

    condition_list.append(condition)
    age_list.append(age)
    sex_list.append(sex)

    time.sleep(0.3)

df["condition"] = condition_list
df["age"] = age_list
df["sex"] = sex_list

df.to_csv(output_file, index=False)

print("Saved:", output_file)