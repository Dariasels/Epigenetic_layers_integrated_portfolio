import pandas as pd
import mysql.connector
import os
import glob
import math

# --- DATABASE CONFIG ---
DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""), 
    "database": "brain_multiomics"
}

def import_atac_folder(folder_path):
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    
    # 1. Ensure table exists
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS atac_peaks (
        peak_id INT AUTO_INCREMENT PRIMARY KEY,
        sample_id VARCHAR(50),
        chrom VARCHAR(20),
        chrom_start INT,
        chrom_end INT,
        peak_name VARCHAR(100),
        signal_value FLOAT,
        INDEX (chrom, chrom_start),
        INDEX (sample_id)
    ) ENGINE=InnoDB;
    """)

    # 2. Find ALL broadPeak files (both .gz and regular)
    # Using a list comprehension to catch both extensions
    all_files = glob.glob(os.path.join(folder_path, "*broadPeak*"))
    # Filter to avoid double-processing the same file if both .gz and .broadPeak exist
    files = list(set([f if not f.endswith('.gz') else f for f in all_files]))
    
    # Sort files to process them in order
    files.sort()
    
    print(f"Found {len(files)} potential ATAC files.")

    sql_insert = """
        INSERT INTO atac_peaks (sample_id, chrom, chrom_start, chrom_end, peak_name, signal_value)
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    for f in files:
        filename = os.path.basename(f)
        
        # Logic to extract GSM ID (GSM3692184)
        sample_id = filename.split('_')[0]
        
        print(f"Processing {sample_id} ({filename})...")

        # Add sample to metadata table
        cursor.execute("INSERT IGNORE INTO samples (sample_id, dataset) VALUES (%s, 'GSE129040')", (sample_id,))

        try:
            # compression='infer' handles .gz and plain text automatically
            df = pd.read_csv(f, sep="\t", header=None, compression='infer', comment='#')
        except Exception as e:
            print(f"  ✗ Error reading {filename}: {e}")
            continue
        
        batch = []
        for _, row in df.iterrows():
            try:
                # broadPeak: 0=chrom, 1=start, 2=end, 3=name, 6=signalValue
                batch.append((
                    sample_id, 
                    str(row[0]), 
                    int(row[1]), 
                    int(row[2]), 
                    str(row[3]), 
                    float(row[6])
                ))
            except:
                continue

            if len(batch) >= 10000:
                cursor.executemany(sql_insert, batch)
                conn.commit()
                batch = []

        if batch:
            cursor.executemany(sql_insert, batch)
            conn.commit()
            
        print(f"  ✓ Imported {len(df):,} peaks.")

    cursor.close()
    conn.close()
    print("\n✓ ALL ATAC FILES IMPORTED.")

if __name__ == "__main__":
    folder = "/home/daria/Documents/UGent/databases/multiomics/chatgpt_version/multiomics_data/Data_files/GSE129040_RAW"
    import_atac_folder(folder)