
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from pathlib import Path
import random

# ─── CONFIGURATION ────────────────────────────────────────────
SDF_INPUT  = Path("data/ligands/african_plants/coconut_sdf_3d-03-2026.sdf")
OUTPUT_CSV = Path("data/ligands/prepared/coconut_filtered.csv")
FINAL_CSV  = Path("data/ligands/prepared/coconut_final_500.csv")

Path("data/ligands/prepared").mkdir(parents=True, exist_ok=True)

# Critères Lipinski stricts
MW_MAX    = 500
LOGP_MAX  = 5
HBD_MAX   = 5
HBA_MAX   = 10
TPSA_MAX  = 140   # perméabilité membranaire
ROTB_MAX  = 10    # flexibilité moléculaire
TARGET_N  = 500   # molécules finales à garder

# ─── FILTRAGE ─────────────────────────────────────────────────
print("Chargement COCONUT SDF (600k molécules)...")
print("Patience — cela peut prendre 5-10 minutes...\n")

suppl = Chem.SDMolSupplier(str(SDF_INPUT), removeHs=False)

results = []
total   = 0
passed  = 0

for mol in suppl:
    total += 1

    if total % 50000 == 0:
        print(f"  Traitement : {total:,} molécules | Retenues : {passed:,}")

    if mol is None:
        continue

    try:
        mw   = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

        # Filtre Lipinski + ADMET basique
        if (mw   <= MW_MAX  and
            logp <= LOGP_MAX and
            hbd  <= HBD_MAX  and
            hba  <= HBA_MAX  and
            tpsa <= TPSA_MAX and
            rotb <= ROTB_MAX):

            smiles = Chem.MolToSmiles(mol)
            name   = mol.GetProp("_Name") if mol.HasProp("_Name") else f"MOL_{passed}"

            results.append({
                "Name"   : name,
                "SMILES" : smiles,
                "MW"     : round(mw, 2),
                "LogP"   : round(logp, 2),
                "HBD"    : hbd,
                "HBA"    : hba,
                "TPSA"   : round(tpsa, 2),
                "RotB"   : rotb,
                "Source" : "COCONUT"
            })
            passed += 1

    except Exception:
        continue

# ─── SAUVEGARDE INTERMÉDIAIRE ─────────────────────────────────
print(f"\n✅ Filtrage terminé")
print(f"   Total traité  : {total:,}")
print(f"   Après filtre  : {passed:,}")

df = pd.DataFrame(results)
df.to_csv(OUTPUT_CSV, index=False)
print(f"   Sauvegardé   → {OUTPUT_CSV}")

# ─── SÉLECTION FINALE 500 ─────────────────────────────────────
print(f"\nSélection finale de {TARGET_N} molécules...")

if len(df) > TARGET_N:
    # Stratification par MW pour avoir de la diversité
    df["MW_bin"] = pd.cut(df["MW"], bins=5)
    
    df_final = (df.groupby("MW_bin", observed=True)
                  .apply(lambda x: x.sample(min(len(x), TARGET_N // 5),
                                             random_state=42))
                  .reset_index(drop=True))

    # drop MW_bin si elle existe encore
    if "MW_bin" in df_final.columns:
        df_final = df_final.drop(columns=["MW_bin"])

    # Compléter si moins de 500
    if len(df_final) < TARGET_N:
        df_no_bin = df.drop(columns=["MW_bin"], errors="ignore")
        remaining = df_no_bin[~df_no_bin.index.isin(df_final.index)]
        extra = remaining.sample(
            min(TARGET_N - len(df_final), len(remaining)),
            random_state=42
        )
        df_final = pd.concat([df_final, extra]).reset_index(drop=True)
else:
    df_final = df.drop(columns=["MW_bin"], errors="ignore")

# ─── RÉSUMÉ ───────────────────────────────────────────────────
print(f"\n{'='*50}")
print(f"  RÉSULTAT FINAL")
print(f"{'='*50}")
print(f"  COCONUT initial     : 600 000 molécules")
print(f"  Après filtre        : {passed:,} molécules")
print(f"  Sélection finale    : {len(df_final)} molécules")
print(f"  Fichier final       → {FINAL_CSV}")
print(f"\n  Stats des retenues :")
print(f"  MW   : {df_final['MW'].mean():.1f} ± {df_final['MW'].std():.1f} Da")
print(f"  LogP : {df_final['LogP'].mean():.2f} ± {df_final['LogP'].std():.2f}")
print(f"  TPSA : {df_final['TPSA'].mean():.1f} ± {df_final['TPSA'].std():.1f}")