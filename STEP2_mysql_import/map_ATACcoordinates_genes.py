#map coordinates to genes
import mysql.connector
import requests
import time

# --- DATABASE CONFIG ---
DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": "simba", 
    "database": "brain_multiomics"
}

def fetch_and_map_genes():
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    
    # 1. Get the list of unique genes you already imported from GSE33000
    print("Fetching unique gene symbols from database...")
    cursor.execute("SELECT DISTINCT gene_symbol FROM rna_expression")
    all_genes = [row[0] for row in cursor.fetchall()]
    print(f"Found {len(all_genes):,} unique genes.")

    # 2. Use MyGene.info API to get coordinates
    # We process in batches of 1000 to be efficient
    url = "https://mygene.info/v3/query"
    batch_size = 1000
    inserted = 0

    print("Fetching coordinates from MyGene.info API...")
    
    for i in range(0, len(all_genes), batch_size):
        batch = all_genes[i:i + batch_size]
        
        # Prepare the API request (batch query expects a JSON list in `q`)
        payload = {
            'q': batch,
            'scopes': 'symbol',
            'fields': 'genomic_pos_hg19,genomic_pos', # Try both hg19 and hg38
            'species': 'human'
        }
        
        try:
            response = requests.post(url, json=payload, timeout=30).json()
        except Exception as e:
            print(f"API Error: {e}")
            continue

        # MyGene returns a list for batch queries; a dict for single queries
        if isinstance(response, dict):
            response = [response]

        for item in response:
            symbol = item.get('query')
            
            # Extract position data (preferring hg19 as it's common for these studies)
            pos = item.get('genomic_pos_hg19') or item.get('genomic_pos')
            
            # If a gene has multiple locations, take the first one
            if isinstance(pos, list): pos = pos[0]
            
            if pos and 'chr' in pos:
                chrom = "chr" + str(pos['chr'])
                start = int(pos['start'])
                end = int(pos['end'])
                strand = "+" if pos['strand'] == 1 else "-"
                
                try:
                    cursor.execute("""
                        INSERT IGNORE INTO gene_coordinates (gene_symbol, chrom, start_pos, end_pos, strand)
                        VALUES (%s, %s, %s, %s, %s)
                    """, (symbol, chrom, start, end, strand))
                    inserted += 1
                except:
                    continue
        
        conn.commit()
        print(f"  Processed {i + len(batch):,} genes. Saved {inserted:,} coordinates...", end="\r")
        time.sleep(0.5) # Be nice to the API

    cursor.close()
    conn.close()
    print(f"\n✓ DONE! Successfully mapped {inserted:,} genes to genomic locations.")

if __name__ == "__main__":
    fetch_and_map_genes()
