from multi_padlock_design import config


def get_gene_synonyms(
    species: str,
    symbol: str,
):
    """Return the official gene symbol for a (species, gene) using Ensembl data.

    Looks up an Ensembl gene ID from a synonyms dictionary, then maps that to the
    official gene_symbol in the headers file. Returns None if not found.
    """
    ensembl_id = None

    ensembl_headers = config.BASE_DIR + f"/ensembl/{species}.allheaders.txt"
    ensembl_synonym_dict = config.BASE_DIR + f"/synonyms/{species}.synonyms.txt"

    # Look up synonym in dictionary
    with open(ensembl_synonym_dict) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            if symbol in parts:
                ensembl_id = parts[0]
                break

    if ensembl_id is None:
        print(f"Synonym {symbol} not found in synonym dict")
        return None

    print(f"Ensembl ID for {symbol} in {species}: {ensembl_id}")

    # Look up the official gene symbol in the headers
    with open(ensembl_headers) as f:
        for line in f:
            if f"gene:{ensembl_id}" in line:
                parts = line.strip().split(" ")
                for part in parts:
                    if part.startswith("gene_symbol:"):
                        return part.split(":")[1]

    return None
