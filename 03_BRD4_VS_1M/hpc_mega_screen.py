import pandas as pd
import numpy as np
import joblib
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

MODEL_PATH = "rf_optimized_fingerprints.joblib"
TRAIN_DATA_PATH = "DATASET_ecfp4_train.csv"
INPUT_FILE = "zinc22_1M_leadlike.csv"
OUTPUT_FILE = "BRD4_1M_Top2000.csv"
SMILES_COL = "SMILES"
NAME_COL = "zincid"
CHUNK_SIZE = 50000
SCORE_THRESHOLD = 6.0
TOP_N = 2000
BATCH_SIZE = 5000

def get_fp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
    except:
        pass
    return None

def process_chunk(chunk, model, X_train, train_sums):
    fps = [get_fp(s) for s in chunk[SMILES_COL]]
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    if not valid_idx:
        return None

    df_clean = chunk.iloc[valid_idx].copy()
    X = np.stack([fps[i] for i in valid_idx]).astype("float32")

    # Tanimoto similarity in batches
    max_sims = []
    for i in range(0, len(X), BATCH_SIZE):
        X_batch = X[i:i+BATCH_SIZE]
        intersection = np.dot(X_batch, X_train.T)
        zinc_sums = X_batch.sum(axis=1).reshape(-1, 1)
        union = zinc_sums + train_sums - intersection
        tanimoto = intersection / union
        knn = np.partition(tanimoto, -5, axis=1)[:, -5:].mean(axis=1)
        max_sims.extend(knn)

    df_clean["pIC50_pred"] = model.predict(X)
    df_clean["knn_similarity"] = np.array(max_sims)
    df_clean["reliability"] = np.where(df_clean["knn_similarity"] >= 0.3, "High", "Low")

    return df_clean[df_clean["pIC50_pred"] >= SCORE_THRESHOLD]

def main():
    start_time = time.time()

    print("Loading model...")
    model = joblib.load(MODEL_PATH)

    print("Loading training data...")
    train_df = pd.read_csv(TRAIN_DATA_PATH)
    feature_cols = [c for c in train_df.columns if c.startswith("ECFP4_")]
    X_train = train_df[feature_cols].values.astype("float32")
    train_sums = X_train.sum(axis=1).reshape(1, -1)

    hits = []
    total = 0

    for i, chunk in enumerate(pd.read_csv(INPUT_FILE, chunksize=CHUNK_SIZE), 1):
        try:
            result = process_chunk(chunk, model, X_train, train_sums)
            if result is not None and len(result) > 0:
                hits.append(result)
                print(f"Chunk {i}: {len(chunk)} processed, {len(result)} hits found.")
            else:
                print(f"Chunk {i}: {len(chunk)} processed, 0 hits found.")
        except Exception as e:
            print(f"Chunk {i} failed: {e}")
        total += len(chunk)

    elapsed = (time.time() - start_time) / 60
    print(f"\nDone. {total} molecules screened in {elapsed:.2f} min.")

    if hits:
        final = pd.concat(hits).sort_values("pIC50_pred", ascending=False).head(TOP_N)
        final.to_csv(OUTPUT_FILE, index=False)
        print(f"Top {len(final)} saved to {OUTPUT_FILE}")
        print(f"Best:   {final.iloc[0]['pIC50_pred']:.3f}")
        print(f"Cutoff: {final.iloc[-1]['pIC50_pred']:.3f}")
    else:
        print(f"No molecules exceeded pIC50 >= {SCORE_THRESHOLD}.")

if __name__ == "__main__":
    main()