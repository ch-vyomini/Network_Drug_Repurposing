
"""
01_build_networks.py
====================

PURPOSE
-------
Build human STRING physical protein-protein interaction networks
at three confidence thresholds: 400, 700 and 900.

For each threshold:
1. Read the STRING physical interaction data.
2. Remove self-loops and duplicate undirected interactions.
3. Build an undirected protein interaction network.
4. Extract the largest connected component (LCC).
5. Save the resulting network.

Additionally:
- Identify proteins shared across all three LCC networks.
- Attach preferred gene symbols as node attributes.
- Calculate network statistics.
- Estimate memory needed for a full distance matrix.

INPUT FILES
-----------
data/raw/string/
    9606.protein.info.v12.5.txt.gz
    9606.protein.physical.links.v12.5.min400.onlyAB.txt.gz

OUTPUT FILES
------------
data/networks/
    physical_conf400_lcc.pkl
    physical_conf700_lcc.pkl
    physical_conf900_lcc.pkl

results/
    physical_network_statistics.csv
    physical_shared_proteins.csv

IMPORTANT
---------
This script uses the PHYSICAL interaction network,
not the full functional association network.
"""

# ==============================================================
# 1. IMPORT LIBRARIES
# ==============================================================

import pickle
from pathlib import Path

import networkx as nx
import pandas as pd


# ==============================================================
# 2. PROJECT SETTINGS
# ==============================================================

# Main project directory.
BASE_DIR = Path(__file__).resolve().parents[1]

# Input files are inside data/raw/string.
STRING_DIR = BASE_DIR / "data" / "raw" / "string"

# Output folders.
NETWORK_DIR = BASE_DIR / "data" / "networks"
RESULTS_DIR = BASE_DIR / "results"

# Input file patterns.
LINKS_PATTERN = "9606.protein.physical.links.v12.5.txt.gz"
INFO_PATTERN = "9606.protein.info*"

# STRING confidence thresholds.
THRESHOLDS = [400, 700, 900]


# ==============================================================
# 3. CREATE OUTPUT FOLDERS
# ==============================================================

# Create these folders automatically if they do not exist.
NETWORK_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Check that the input folder exists.
if not STRING_DIR.exists():
    raise FileNotFoundError(
        f"STRING input folder not found:\n{STRING_DIR}"
    )


# ==============================================================
# 4. FIND THE INPUT FILES
# ==============================================================

def find_one_file(pattern):
    """
    Find exactly one actual file matching a pattern.

    Important:
    The downloaded STRING archives may have extracted folders
    with similar names. item.is_file() ignores those folders.
    """

    matches = sorted(
        item
        for item in STRING_DIR.glob(pattern)
        if item.is_file()
    )

    if len(matches) == 0:
        raise FileNotFoundError(
            f"No actual file matching '{pattern}' was found in:\n"
            f"{STRING_DIR}"
        )

    if len(matches) > 1:
        names = "\n  ".join(item.name for item in matches)

        raise RuntimeError(
            f"More than one actual file matches '{pattern}':\n"
            f"  {names}\n"
            "Please check the files in the input folder."
        )

    return matches[0]


# Locate the physical links and protein information files.
LINKS_FILE = find_one_file(LINKS_PATTERN)
INFO_FILE = find_one_file(INFO_PATTERN)

NETWORK_LABEL = "physical"

print("=" * 65)
print("STRING NETWORK BUILDING")
print("=" * 65)

print("Network type :", NETWORK_LABEL)
print("Links file   :", LINKS_FILE.name)
print("Info file    :", INFO_FILE.name)


# ==============================================================
# 5. READ AND VALIDATE THE INTERACTION DATA
# ==============================================================

print("\nReading STRING physical interaction data...")

# STRING files are whitespace-separated.
# compression='infer' automatically reads .gz files.

links = pd.read_csv(
    LINKS_FILE,
    sep=r"\s+",
    compression="infer"
)

# Check the required columns.
required_columns = {
    "protein1",
    "protein2",
    "combined_score"
}

if not required_columns.issubset(links.columns):
    raise ValueError(
        "Unexpected columns in the STRING links file:\n"
        f"{list(links.columns)}"
    )

# Convert confidence scores to numeric values.
links["combined_score"] = pd.to_numeric(
    links["combined_score"],
    errors="raise"
)

print("\nInteraction data successfully loaded.")

print(f"Raw interaction rows : {len(links):,}")

print(
    "Confidence score range:",
    links["combined_score"].min(),
    "-",
    links["combined_score"].max()
)


# ==============================================================
# 6. CLEAN THE INTERACTION DATA
# ==============================================================

print("\nCleaning interaction data...")

# STRING networks are undirected.
# Therefore, A-B and B-A represent the same interaction.

# Put each protein pair into a consistent alphabetical order.
swap = links["protein1"] > links["protein2"]

edges = pd.DataFrame({

    "a": links["protein1"].where(
        ~swap,
        links["protein2"]
    ),

    "b": links["protein2"].where(
        ~swap,
        links["protein1"]
    ),

    "score": links["combined_score"]

})

# Remove self-loops (a protein interacting with itself).
self_loops = int((edges["a"] == edges["b"]).sum())

edges = edges[edges["a"] != edges["b"]].copy()

# Remove duplicate undirected pairs.
# If duplicates have different scores, retain the highest score.
rows_before = len(edges)

edges = (
    edges
    .sort_values("score", ascending=False)
    .drop_duplicates(subset=["a", "b"])
    .copy()
)

duplicates_removed = rows_before - len(edges)

unique_proteins = pd.concat(
    [edges["a"], edges["b"]]
).nunique()

print(f"Self-loops removed       : {self_loops:,}")
print(f"Duplicate pairs removed  : {duplicates_removed:,}")
print(f"Unique interaction edges : {len(edges):,}")
print(f"Unique proteins          : {unique_proteins:,}")


# ==============================================================
# 7. READ PROTEIN INFORMATION
# ==============================================================

print("\nReading STRING protein information...")

info = pd.read_csv(
    INFO_FILE,
    sep="\t",
    compression="infer"
)

# The first column in the STRING file begins with '#'.
# Rename it to a convenient column name.
info = info.rename(
    columns={
        info.columns[0]: "string_protein_id"
    }
)

# Check that the expected columns exist.
if not {
    "string_protein_id",
    "preferred_name"
}.issubset(info.columns):

    raise ValueError(
        "Unexpected columns in the STRING protein-info file:\n"
        f"{list(info.columns)}"
    )

# Create a lookup dictionary:
# STRING protein ID -> preferred gene symbol.
symbol_lookup = dict(
    zip(
        info["string_protein_id"],
        info["preferred_name"]
    )
)

print(f"Protein information records: {len(info):,}")


# ==============================================================
# 8. BUILD NETWORKS AT EACH CONFIDENCE THRESHOLD
# ==============================================================

print("\nBuilding networks at the three confidence thresholds...")

# Store the LCC graphs and their statistics.
lcc_graphs = {}

stats = []

for threshold in THRESHOLDS:

    print("\n" + "-" * 55)
    print(f"Processing confidence threshold: {threshold}")
    print("-" * 55)

    # ----------------------------------------------------------
    # 8A. FILTER INTERACTIONS BY CONFIDENCE SCORE
    # ----------------------------------------------------------

    sub = edges[
        edges["score"] >= threshold
    ].copy()

    print(f"Interactions retained: {len(sub):,}")

    if sub.empty:
        raise ValueError(
            f"No interactions remain at threshold {threshold}."
        )

    # ----------------------------------------------------------
    # 8B. CREATE AN UNDIRECTED NETWORK
    # ----------------------------------------------------------

    G = nx.from_pandas_edgelist(
        sub,
        source="a",
        target="b",
        edge_attr="score",
        create_using=nx.Graph()
    )

    # Attach preferred gene symbols to protein nodes.
    nx.set_node_attributes(
        G,
        {
            protein: symbol_lookup.get(protein)
            for protein in G.nodes()
        },
        "preferred_name"
    )

    # ----------------------------------------------------------
    # 8C. IDENTIFY CONNECTED COMPONENTS
    # ----------------------------------------------------------

    components = list(
        nx.connected_components(G)
    )

    if not components:
        raise ValueError(
            f"No connected components found at threshold {threshold}."
        )

    # Find the largest connected component.
    biggest_component = max(
        components,
        key=len
    )

    # Extract the largest connected component.
    lcc = G.subgraph(
        biggest_component
    ).copy()

    # Save this threshold's LCC.
    lcc_graphs[threshold] = lcc

    # ----------------------------------------------------------
    # 8D. CALCULATE NETWORK STATISTICS
    # ----------------------------------------------------------

    nodes_all = G.number_of_nodes()
    edges_all = G.number_of_edges()

    lcc_nodes = lcc.number_of_nodes()
    lcc_edges = lcc.number_of_edges()

    degrees = [
        degree
        for node, degree in lcc.degree()
    ]

    # Estimate memory needed for a full float64 distance matrix.
    # This is an estimate only; no distance matrix is created here.
    matrix_gb = (
        lcc_nodes * lcc_nodes * 8
    ) / 1e9

    # Calculate the fraction of all network nodes retained in the LCC.
    lcc_fraction = (
        lcc_nodes / nodes_all
        if nodes_all > 0
        else 0
    )

    mean_degree = (
        sum(degrees) / lcc_nodes
        if lcc_nodes > 0
        else 0
    )

    max_degree = (
        max(degrees)
        if degrees
        else 0
    )

    # Store the statistics.
    stats.append({

        "network": NETWORK_LABEL,

        "threshold": threshold,

        "nodes_all": nodes_all,

        "edges_all": edges_all,

        "components": len(components),

        "lcc_nodes": lcc_nodes,

        "lcc_edges": lcc_edges,

        "lcc_fraction_of_nodes": round(
            lcc_fraction,
            4
        ),

        "lcc_mean_degree": round(
            mean_degree,
            2
        ),

        "lcc_max_degree": max_degree,

        "distance_matrix_GB_float64": round(
            matrix_gb,
            2
        )

    })

    # ----------------------------------------------------------
    # 8E. SAVE THE LCC GRAPH
    # ----------------------------------------------------------

    output_file = (
        NETWORK_DIR
        / f"{NETWORK_LABEL}_conf{threshold}_lcc.pkl"
    )

    with open(output_file, "wb") as file:

        pickle.dump(
            lcc,
            file,
            protocol=pickle.HIGHEST_PROTOCOL
        )

    print(f"LCC nodes       : {lcc_nodes:,}")
    print(f"LCC edges       : {lcc_edges:,}")
    print(f"Connected parts : {len(components):,}")
    print(f"Estimated matrix memory: {matrix_gb:.2f} GB")
    print(f"Saved graph     : {output_file}")


# ==============================================================
# 9. IDENTIFY PROTEINS SHARED ACROSS ALL THREE LCC NETWORKS
# ==============================================================

print("\n" + "=" * 65)
print("IDENTIFYING SHARED PROTEINS")
print("=" * 65)

# Convert the nodes of each LCC into a set.
node_sets = [
    set(graph.nodes())
    for graph in lcc_graphs.values()
]

# Find proteins present in every LCC.
shared_proteins = set.intersection(
    *node_sets
)

print(
    "Proteins shared across all three LCCs:",
    f"{len(shared_proteins):,}"
)

# Create a table containing the shared STRING protein IDs
# and their preferred gene symbols.
shared_df = pd.DataFrame({

    "string_protein_id": sorted(
        shared_proteins
    )

})

shared_df["preferred_name"] = (
    shared_df["string_protein_id"]
    .map(symbol_lookup)
)

# Save the shared protein table.
shared_file = (
    RESULTS_DIR
    / f"{NETWORK_LABEL}_shared_proteins.csv"
)

shared_df.to_csv(
    shared_file,
    index=False
)

print("Shared protein table saved:", shared_file)


# ==============================================================
# 10. IDENTIFIER QUALITY CHECKS
# ==============================================================

print("\nChecking protein identifier quality...")

# Use the LCC from the lowest confidence threshold.
loosest_lcc = lcc_graphs[min(THRESHOLDS)]

# Check how many proteins have no preferred gene symbol.
symbols = pd.Series({

    protein: symbol_lookup.get(protein)

    for protein in loosest_lcc.nodes()

})

missing_symbols = int(
    symbols.isna().sum()
)

# Check whether multiple STRING proteins share a preferred symbol.
duplicate_symbol_counts = (
    symbols
    .dropna()
    .value_counts()
)

repeated_symbols = int(
    (duplicate_symbol_counts > 1).sum()
)

print(
    "Proteins without preferred symbols:",
    f"{missing_symbols:,}"
)

print(
    "Preferred symbols shared by multiple STRING IDs:",
    f"{repeated_symbols:,}"
)

# Important:
# We retain STRING protein IDs as network nodes.
# We do NOT merge nodes just because they share a gene symbol.


# ==============================================================
# 11. SAVE NETWORK STATISTICS
# ==============================================================

stats_df = pd.DataFrame(stats)

stats_file = (
    RESULTS_DIR
    / f"{NETWORK_LABEL}_network_statistics.csv"
)

stats_df.to_csv(
    stats_file,
    index=False
)

print("\n" + "=" * 65)
print("FINAL NETWORK STATISTICS")
print("=" * 65)

print(
    stats_df.to_string(index=False)
)

print("\nStatistics saved:", stats_file)


# ==============================================================
# 12. FINAL SUMMARY
# ==============================================================

print("\n" + "=" * 65)
print("SCRIPT 01 COMPLETED SUCCESSFULLY")
print("=" * 65)

print("\nNetworks created:")

for threshold in THRESHOLDS:

    graph = lcc_graphs[threshold]

    print(
        f"Confidence {threshold}: "
        f"{graph.number_of_nodes():,} proteins, "
        f"{graph.number_of_edges():,} interactions"
    )

print(
    "\nAll network files are saved in:",
    NETWORK_DIR
)

print(
    "All result tables are saved in:",
    RESULTS_DIR
)
