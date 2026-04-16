import mysql.connector

# --- DATABASE CONFIG ---
DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": "simba", 
    "database": "brain_multiomics"
}

# --- THE FULL MANUAL LIST (Used for updating info) ---
PLASTICITY_GENES = [
    # (Symbol, Category, Source)
    ("BDNF", "neurotrophin signalling", "KEGG hsa04722"),
    ("NTRK2", "neurotrophin signalling", "KEGG hsa04722"),
    ("GRIA1", "LTP/LTD", "KEGG hsa04720"),
    ("GRIA2", "LTP/LTD", "KEGG hsa04720"),
    ("GRIN1", "LTP/LTD", "KEGG hsa04720"),
    ("GRIN2A", "LTP/LTD", "KEGG hsa04720"),
    ("GRIN2B", "LTP/LTD", "KEGG hsa04720"),
    ("CAMK2A", "LTP/LTD", "KEGG hsa04720"),
    ("CREB1", "LTP/LTD", "KEGG hsa04720"),
    ("ARC", "immediate early gene", "manual curation"),
    ("FOS", "immediate early gene", "manual curation"),
    ("DLG4", "synaptic scaffolding", "SynGO"),
    ("SHANK3", "synaptic scaffolding", "SynGO"),
    ("DNMT1", "epigenetic regulation", "Narayanan 2014"),
    ("APP", "AD-related plasticity", "manual curation"),
    ("APOE", "AD-related plasticity", "manual curation"),
    ("MTOR", "mTOR pathway", "manual curation"),
    ("WNT7A", "Wnt signalling", "manual curation"),
    ("GSK3B", "Wnt signalling", "manual curation")
    # ... (The script will work for any gene in this format)
]

def patch_metadata():
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()

    print("Updating categories for genes found in the manual list...")
    
    # We use UPDATE instead of INSERT. 
    # This keeps the gene in the table but fills in the Category and Source columns.
    sql_update = """
        UPDATE plasticity_genes 
        SET category = %s, source = %s 
        WHERE gene_symbol = %s
    """
    
    # Prepare data for UPDATE (Note the order: category, source, symbol)
    data = [(cat, src, sym) for sym, cat, src in PLASTICITY_GENES]
    
    cursor.executemany(sql_update, data)
    conn.commit()
    
    print(f"✓ Successfully updated categories for {cursor.rowcount} genes.")

    # Show the result
    print("\n--- Current Gene Distribution ---")
    cursor.execute("SELECT category, COUNT(*) FROM plasticity_genes GROUP BY category")
    for cat, count in cursor.fetchall():
        print(f"  {cat}: {count} genes")

    cursor.close()
    conn.close()

if __name__ == "__main__":
    patch_metadata()