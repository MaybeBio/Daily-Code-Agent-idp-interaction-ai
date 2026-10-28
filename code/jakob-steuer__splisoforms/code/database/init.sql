CREATE TABLE genes (
    id SERIAL PRIMARY KEY,
    symbol VARCHAR(50) UNIQUE NOT NULL,
    hgnc_symbol VARCHAR(50),
    uniprot_id VARCHAR(50) UNIQUE,
    description TEXT
);

CREATE TABLE isoforms (
    id SERIAL PRIMARY KEY,
    gene_id INTEGER REFERENCES genes(id),
    transcript_id VARCHAR(50) UNIQUE NOT NULL,
    sequence TEXT,
    fasta_sequence TEXT,
    delta_length INTEGER,
    junction_coords VARCHAR(255)
);

CREATE TABLE structures (
    id SERIAL PRIMARY KEY,
    isoform_id INTEGER REFERENCES isoforms(id),
    storage_uri VARCHAR(255) NOT NULL,
    af3_cif_url TEXT,
    af3_pae_url TEXT,
    confidence_score FLOAT,
    avg_plddt FLOAT
);

CREATE TABLE clinical_frequencies (
    id SERIAL PRIMARY KEY,
    isoform_id INTEGER REFERENCES isoforms(id),
    cancer_type VARCHAR(100),
    frequency FLOAT
);

CREATE TABLE ppi_networks (
    id SERIAL PRIMARY KEY,
    isoform_id_1 INTEGER REFERENCES isoforms(id),
    target_uniprot_id VARCHAR(50),
    interaction_score FLOAT,
    interaction_retained BOOLEAN
);
