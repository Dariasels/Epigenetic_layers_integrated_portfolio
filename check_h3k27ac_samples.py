#!/usr/bin/env python3
"""Check H3K27ac sample counts for weighted mean calculation."""
import pymysql

conn = pymysql.connect(host='localhost', user='daria', password='simba', database='brain_multiomics')
cursor = conn.cursor()

# Check what's in the enhancers table
cursor.execute("SELECT COUNT(*) FROM enhancers WHERE h3k27ac_ad_mean IS NOT NULL")
ad_count = cursor.fetchone()[0]
print(f"Enhancers with h3k27ac_ad_mean: {ad_count}")

cursor.execute("SELECT COUNT(*) FROM enhancers WHERE h3k27ac_control_mean IS NOT NULL")
ctrl_count = cursor.fetchone()[0]
print(f"Enhancers with h3k27ac_control_mean: {ctrl_count}")

# Check sample structure in chip_peak_counts
cursor.execute("SELECT COUNT(DISTINCT sample_id) FROM chip_peak_counts")
unique_samples = cursor.fetchone()[0]
print(f"\nUnique samples in chip_peak_counts: {unique_samples}")

# Count by condition
cursor.execute("""
SELECT samples.condition, COUNT(DISTINCT chip_peak_counts.sample_id) 
FROM chip_peak_counts 
JOIN samples ON chip_peak_counts.sample_id = samples.sample_id
GROUP BY samples.condition
""")
results = cursor.fetchall()
print("\nSamples by condition in chip_peak_counts:")
for row in results:
    print(f"  {row[0]}: {row[1]}")

# Now compute weighted average for a sample enhancer
cursor.execute("""
SELECT h3k27ac_ad_mean, h3k27ac_control_mean, h3k27ac_ad_n, h3k27ac_control_n 
FROM enhancers 
WHERE h3k27ac_ad_mean IS NOT NULL LIMIT 5
""")
results = cursor.fetchall()
if results:
    print("\nSample enhancers and their mean values:")
    for row in results:
        ad_mean, ctrl_mean, ad_n, ctrl_n = row
        if ad_n and ctrl_n:
            print(f"  AD_mean={ad_mean}, CTRL_mean={ctrl_mean}, AD_n={ad_n}, CTRL_n={ctrl_n}")
            weighted = (ad_mean * ad_n + ctrl_mean * ctrl_n) / (ad_n + ctrl_n)
            print(f"    -> Weighted avg: {weighted:.2f}")
        elif ad_n or ctrl_n:
            print(f"  AD_mean={ad_mean}, CTRL_mean={ctrl_mean}, AD_n={ad_n}, CTRL_n={ctrl_n}")

conn.close()
