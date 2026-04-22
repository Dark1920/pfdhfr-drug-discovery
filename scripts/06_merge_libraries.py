"""
06_merge_libraries.py — VERSION CORRIGÉE
Fusion : african filtrées + COCONUT réduit à 500
"""

import pandas as pd
from rdkit import Chem
from pathlib import Path

# ─── FICHIERS ─────────────────────────────────────────────────
AFRICAN_FILTERED = Path("data/ligands/prepared/molecules_filtered.csv")
AFRICAN_ORIGINAL = Path("data/ligands/african_plants/African_molecule.csv")
COCONUT_CSV      = Path("data/ligands/prepared/coconut_filtered.csv")
OUTPUT_CSV       = Path("data/ligands/prepared/library_merged.csv")

Path("data/ligands/prepared").mkdir(parents=True, exist_ok=True)

# ══════════════════════════════════════════════════════════════
# PARTIE 1 — Molécules africaines avec vrais noms
# ══════════════════════════════════════════════════════════════
print("=" * 55)
print("  PARTIE 1 : Molécules africaines")
print("=" * 55)

df_filtered  = pd.read_csv(AFRICAN_FILTERED)
df_original = pd.read_csv(AFRICAN_ORIGINAL, sep=';')

# Réattribuer les vrais noms depuis le fichier original
# Les deux fichiers ont le même ordre de molécules
df_filtered["Name"] = df_original["Name"].values[:len(df_filtered)]

# Ajouter colonne Plant si elle existe dans l'original
if "Plant" in df_original.columns:
    df_filtered["Plant"] = df_original["Plant"].values[:len(df_filtered)]
else:
    df_filtered["Plant"] = "West African Plant"

df_filtered["Source"] = "African_Plants"

# Garder colonnes utiles
df_african = df_filtered[[
    "Name", "Canonical SMILES", "MW",
    "Consensus Log P", "#H-bond donors",
    "#H-bond acceptors", "TPSA",
    "#Rotatable bonds", "GI absorption",
    "Plant", "Source"
]].rename(columns={
    "Canonical SMILES"  : "SMILES",
    "Consensus Log P"   : "LogP",
    "#H-bond donors"    : "HBD",
    "#H-bond acceptors" : "HBA",
    "#Rotatable bonds"  : "RotB"
})

print(f"  Molécules africaines : {len(df_african)} ✅")
print(f"  Noms : {df_african['Name'].tolist()}")

# ══════════════════════════════════════════════════════════════
# PARTIE 2 — COCONUT réduit à 500 molécules diversifiées
# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 55)
print("  PARTIE 2 : COCONUT → sélection 500 molécules")
print("=" * 55)

df_coconut_full = pd.read_csv(COCONUT_CSV)
print(f"  COCONUT total disponible : {len(df_coconut_full):,}")

# Sélection de 500 molécules diversifiées par stratification sur MW
TARGET = 500

df_coconut_full["MW_bin"] = pd.cut(df_coconut_full["MW"], bins=5)

df_coconut = (
    df_coconut_full
    .groupby("MW_bin", observed=True)
    .apply(lambda x: x.sample(
        min(len(x), TARGET // 5),
        random_state=42
    ))
    .reset_index(drop=True)
)

# Compléter jusqu'à 500 si besoin
if len(df_coconut) < TARGET:
    already = df_coconut.index
    reste   = df_coconut_full[~df_coconut_full.index.isin(already)]
    extra   = reste.sample(
        min(TARGET - len(df_coconut), len(reste)),
        random_state=42
    )
    df_coconut = pd.concat([df_coconut, extra]).reset_index(drop=True)

df_coconut = df_coconut.drop(columns=["MW_bin"], errors="ignore")
df_coconut["Source"] = "COCONUT"
df_coconut["Plant"]  = "Natural Product"

# Harmoniser colonnes COCONUT
df_coconut = df_coconut[[
    "Name", "SMILES", "MW", "LogP",
    "HBD", "HBA", "TPSA", "RotB",
    "Plant", "Source"
]]

# Ajouter GI absorption si absente
if "GI absorption" not in df_coconut.columns:
    df_coconut["GI absorption"] = "Unknown"

print(f"  COCONUT sélectionné  : {len(df_coconut)} molécules ✅")

# ══════════════════════════════════════════════════════════════
# PARTIE 3 — Fusion et déduplication
# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 55)
print("  PARTIE 3 : Fusion & Déduplication")
print("=" * 55)

# Harmoniser les colonnes african pour le concat
if "GI absorption" not in df_african.columns:
    df_african["GI absorption"] = "Unknown"

COLS = ["Name", "SMILES", "MW", "LogP", "HBD",
        "HBA", "TPSA", "RotB", "GI absorption",
        "Plant", "Source"]

df_african = df_african[COLS]
df_coconut = df_coconut[COLS]

df_merged = pd.concat([df_african, df_coconut], ignore_index=True)
print(f"  Avant déduplication  : {len(df_merged)}")

# Calcul SMILES canoniques
def to_canonical(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

df_merged["SMILES_can"] = df_merged["SMILES"].apply(to_canonical)
df_merged = df_merged.dropna(subset=["SMILES_can"])
df_merged = df_merged.drop_duplicates(subset=["SMILES_can"], keep="first")
df_merged = df_merged.drop(columns=["SMILES_can"])
df_merged = df_merged.reset_index(drop=True)

# ID unique
df_merged.insert(0, "ID", [f"LIG_{i+1:04d}" for i in range(len(df_merged))])

# ══════════════════════════════════════════════════════════════
# PARTIE 4 — Sauvegarde + Rapport
# ══════════════════════════════════════════════════════════════
df_merged.to_csv(OUTPUT_CSV, index=False)

print(f"  Après déduplication  : {len(df_merged)}")
print("\n" + "=" * 55)
print("  RÉSULTAT FINAL")
print("=" * 55)
print(f"  Molécules africaines : {len(df_african)}")
print(f"  Molécules COCONUT    : {len(df_coconut)}")
print(f"  Total fusionné       : {len(df_merged)}")
print(f"  Fichier sauvegardé   → {OUTPUT_CSV}")

print("\n  Répartition par source :")
for source, count in df_merged["Source"].value_counts().items():
    pct = count / len(df_merged) * 100
    print(f"    {source:<20} : {count:>4} ({pct:.1f}%)")

print("\n  Statistiques :")
print(f"    MW   : {df_merged['MW'].mean():.1f} ± {df_merged['MW'].std():.1f} Da")
print(f"    LogP : {df_merged['LogP'].mean():.2f} ± {df_merged['LogP'].std():.2f}")
print(f"    TPSA : {df_merged['TPSA'].mean():.1f} ± {df_merged['TPSA'].std():.1f}")

print("\n  Aperçu des 5 premières molécules :")
print(df_merged[["ID","Name","MW","LogP","Source"]].head().to_string(index=False))

print("\n✅ Bibliothèque fusionnée prête pour la conversion PDBQT !")