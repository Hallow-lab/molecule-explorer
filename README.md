# Molecule Explorer (Streamlit)

## Features
- Search a starter molecule library and open each molecule.
- Search PubChem by compound name.
- Upload a CSV molecule library (`smiles` column required; `name` optional).
- Per-molecule 2D structure, interactive 3D view, canonical SMILES, and calculated properties.
- Molecular formula, molecular weight, exact mass, XLogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring count, heavy atoms, and formal charge.

## Run locally
Use Python 3.10 or 3.11 in a virtual environment:

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud
1. Upload `app.py`, `requirements.txt`, and this README to a GitHub repository.
2. In Streamlit Community Cloud, choose the repository and set the main file path to `app.py`.
3. Deploy. PubChem search needs internet access.

## Notes
- “All molecules” is not a finite set. This app provides a starter library, PubChem name search, and CSV import so the library can be expanded.
- Calculated values are descriptors, not experimental results.
- Some compounds may not have a PubChem 3D conformer. In those cases the app generates an illustrative conformer using RDKit.
