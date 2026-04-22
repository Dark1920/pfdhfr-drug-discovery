"""
04_filter_lipinski.py
Filtrage Lipinski — adapté aux colonnes SwissADME réelles
"""

import pandas as pd
from pathlib import Path

# ─── FICHIERS ─────────────────────────────────────────────────
INPUT_CSV  = Path("data/ligands/african_plants/swissadme.csv")
OUTPUT_CSV = Path("data/ligands/prepared/molecules_filtered.csv")
REJECTED_CSV = Path("data/ligands/prepared/molecules_rejected.csv")

Path("data/ligands/prepared").mkdir(parents=True, exist_ok=True)

# ─── CHARGEMENT ───────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)
print(f"Molécules initiales : {len(df)}")

# ─── FILTRAGE LIPINSKI ────────────────────────────────────────
# Colonnes exactes de ton fichier SwissADME
mask = (
    (df["MW"] <= 500) &
    (df["Consensus Log P"] <= 5) &
    (df["#H-bond donors"] <= 5) &
    (df["#H-bond acceptors"] <= 10) &
    (df["Lipinski #violations"] == 0)
)

df_ok  = df[mask].copy()
df_rej = df[~mask].copy()

# ─── SAUVEGARDE ───────────────────────────────────────────────
df_ok.to_csv(OUTPUT_CSV, index=False)
df_rej.to_csv(REJECTED_CSV, index=False)

# ─── RAPPORT ──────────────────────────────────────────────────
print(f"\n✅ Retenues  : {len(df_ok)} molécules")
print(f"❌ Rejetées  : {len(df_rej)} molécules")

print("\n=== MOLÉCULES RETENUES ===")
for _, r in df_ok.iterrows():
    print(f"  {r['Molecule']:<12} | MW={r['MW']:>6.1f} "
          f"| LogP={r['Consensus Log P']:>5.2f} "
          f"| HBD={int(r['#H-bond donors'])} "
          f"| HBA={int(r['#H-bond acceptors'])} "
          f"| GI={r['GI absorption']}")

print("\n=== MOLÉCULES REJETÉES ===")
for _, r in df_rej.iterrows():
    print(f"  {r['Molecule']:<12} | MW={r['MW']:>6.1f} "
          f"| violations={int(r['Lipinski #violations'])}")