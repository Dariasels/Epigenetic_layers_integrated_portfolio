import pandas as pd

# 'skiprows' bypasses the text header
# We also use 'low_memory=False' because these files have mixed data types
df = pd.read_csv(
    'GSE59685_betas.csv.gz', 
    compression='gzip', 
    skiprows=5, 
    low_memory=False
)

print(df.head())