import io
import json
import requests
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem import Crippen
import py3Dmol
def get_2d_image(smiles):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        + requests.utils.quote(smiles, safe="")
        + "/PNG"
    )
    response = requests.get(url, timeout=20)
    response.raise_for_status()
    return response.content
from stmol import showmol

st.set_page_config(page_title="Molecule Explorer", page_icon="🧬", layout="wide")
st.title("🧬 Molecule Explorer")
st.caption("Explore molecular structures, 3D conformations, SMILES, and physicochemical properties.")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# A small starter library drawn from the chemistry reports/reference examples.
# Add more compounds using search or by uploading a CSV with name, smiles columns.
STARTER = [
    {"name":"Ethanol", "smiles":"CCO", "category":"Alcohol"},
    {"name":"Aspirin", "smiles":"CC(=O)Oc1ccccc1C(=O)O", "category":"Drug / analgesic"},
    {"name":"Benzene", "smiles":"c1ccccc1", "category":"Aromatic hydrocarbon"},
    {"name":"Acetic acid", "smiles":"CC(=O)O", "category":"Carboxylic acid"},
    {"name":"Phenol", "smiles":"Oc1ccccc1", "category":"Phenol"},
    {"name":"Caffeine", "smiles":"Cn1c(=O)c2c(ncn2C)n(C)c1=O", "category":"Alkaloid"},
    {"name":"Glucose", "smiles":"C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "category":"Carbohydrate"},
    {"name":"Paracetamol", "smiles":"CC(=O)NC1=CC=C(C=C1)O", "category":"Drug / analgesic"},
    {"name":"Ibuprofen", "smiles":"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "category":"Drug / NSAID"},
    {"name":"Nicotine", "smiles":"CN1CCC[C@H]1c2cccnc2", "category":"Alkaloid"},
    {"name":"Lactic acid", "smiles":"C[C@H](O)C(=O)O", "category":"Organic acid"},
    {"name":"Urea", "smiles":"NC(=O)N", "category":"Organic compound"},
]

@st.cache_data(show_spinner=False, ttl=3600)
def pubchem_search(query):
    url = f"{PUBCHEM}/compound/name/{requests.utils.quote(query)}/property/Title,CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight,InChIKey/JSON"
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        return None
    data = r.json()
    props = data.get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None

@st.cache_data(show_spinner=False, ttl=3600)
def pubchem_3d_sdf(query):
    # PubChem may not have a 3D conformer for every compound.
    url = f"{PUBCHEM}/compound/name/{requests.utils.quote(query)}/SDF?record_type=3d"
    r = requests.get(url, timeout=20)
    return r.text if r.status_code == 200 and "V2000" in r.text or (r.status_code == 200 and "V3000" in r.text) else None

def get_mol(smiles):
    return Chem.MolFromSmiles(str(smiles)) if smiles else None

def prop_table(mol):
    exact = Descriptors.ExactMolWt(mol)
    return pd.DataFrame([
        ("Molecular formula", rdMolDescriptors.CalcMolFormula(mol)),
        ("Molecular weight (g/mol)", round(Descriptors.MolWt(mol), 4)),
        ("Exact molecular weight (Da)", round(exact, 4)),
        ("XLogP", round(Crippen.MolLogP(mol), 3)),
        ("Topological polar surface area (Å²)", round(rdMolDescriptors.CalcTPSA(mol), 2)),
        ("Hydrogen-bond donors", rdMolDescriptors.CalcNumHBD(mol)),
        ("Hydrogen-bond acceptors", rdMolDescriptors.CalcNumHBA(mol)),
        ("Rotatable bonds", rdMolDescriptors.CalcNumRotatableBonds(mol)),
        ("Ring count", rdMolDescriptors.CalcNumRings(mol)),
        ("Heavy atoms", mol.GetNumHeavyAtoms()),
        ("Formal charge", Chem.GetFormalCharge(mol)),
    ], columns=["Property", "Value"])

def render_3d(mol, sdf_text=None):
    view = py3Dmol.view(width=700, height=430)
    if sdf_text:
        view.addModel(sdf_text, "sdf")
    else:
        m = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(m, params)
        if status != 0:
            st.info("A 3D conformer could not be generated for this structure.")
            return
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass
        block = Chem.MolToMolBlock(m)
        view.addModel(block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
    view.setBackgroundColor("#ffffff")
    view.zoomTo()
    showmol(view, height=450, width=700)

with st.sidebar:
    st.header("Find molecules")
    mode = st.radio("Browse mode", ["Molecule library", "Search PubChem", "Upload molecule CSV"])
    st.caption("The library is expandable. PubChem search requires an internet connection.")

if mode == "Molecule library":
    df = pd.DataFrame(STARTER)
    q = st.text_input("Filter library", placeholder="Search by molecule name, SMILES, or category")
    if q:
        mask = df.astype(str).apply(lambda col: col.str.contains(q, case=False, na=False)).any(axis=1)
        df = df[mask]
    st.subheader(f"Molecule library · {len(df)} shown")
    st.dataframe(df[["name", "category", "smiles"]], use_container_width=True, hide_index=True)
    choices = df["name"].tolist()
    if choices:
        selected = st.selectbox("Open a molecule to view details", choices)
        record = df[df["name"] == selected].iloc[0].to_dict()
        name, smiles = record["name"], record["smiles"]
        source_sdf = None
    else:
        st.info("No molecules match this filter.")
        st.stop()

elif mode == "Search PubChem":
    query = st.text_input("Molecule name", placeholder="e.g., caffeine, aspirin, glucose")
    if not query:
        st.info("Enter a compound name to search PubChem.")
        st.stop()
    with st.spinner("Searching PubChem…"):
        result = pubchem_search(query)
    if not result:
        st.error("No PubChem result found. Try a common name or the exact compound name.")
        st.stop()
    name = result.get("Title") or query
    smiles = result.get("ConnectivitySMILES") or result.get("CanonicalSMILES") or result.get("SMILES")
    if not smiles:
        st.error("PubChem did not return a usable SMILES string for this result.")
        st.stop()
    source_sdf = pubchem_3d_sdf(query)
    st.caption(f"PubChem result · CID {result.get('CID', 'not supplied')}")
    st.session_state["last_pubchem"] = {"name": name, "smiles": smiles}
    st.success(f"Loaded: {name}")

else:
    uploaded = st.file_uploader("Upload CSV", type=["csv"])
    st.caption("Required column: smiles. Optional column: name. Example: name,smiles")
    if uploaded is None:
        st.stop()
    up = pd.read_csv(uploaded)
    lower = {c.lower(): c for c in up.columns}
    if "smiles" not in lower:
        st.error("The CSV must include a column named 'smiles'.")
        st.stop()
    smiles_col = lower["smiles"]
    name_col = lower.get("name")
    up["__name"] = up[name_col].fillna("Unnamed molecule").astype(str) if name_col else [f"Molecule {i+1}" for i in range(len(up))]
    up["__smiles"] = up[smiles_col].astype(str)
    up = up[up["__smiles"].str.strip().ne("")]
    if up.empty:
        st.warning("No SMILES records found.")
        st.stop()
    st.dataframe(up[["__name", "__smiles"]].rename(columns={"__name":"name","__smiles":"smiles"}), use_container_width=True, hide_index=True)
    selected = st.selectbox("Open a molecule to view details", up["__name"].tolist())
    record = up[up["__name"] == selected].iloc[0]
    name, smiles, source_sdf = selected, record["__smiles"], None

mol = get_mol(smiles)
if mol is None:
    st.error("This SMILES string is invalid or unsupported by RDKit.")
    st.code(smiles)
    st.stop()

canonical = Chem.MolToSmiles(mol, canonical=True)
st.divider()
st.header(name)
left, right = st.columns([1, 1.15], gap="large")
with left:
    st.subheader("2D structure")
    img = get_2d_image(mol, size=(520, 360), kekulize=True)
    st.image(img, use_container_width=True)
    st.subheader("SMILES")
    st.text_input("Canonical SMILES", value=canonical, key=f"smiles_{name}", disabled=True)
    if st.button("Copy-ready SMILES details"):
        st.code(canonical)
with right:
    st.subheader("Interactive 3D view")
    render_3d(mol, source_sdf)
st.subheader("Molecular properties")
st.dataframe(prop_table(mol), use_container_width=True, hide_index=True)
st.caption("Calculated properties are cheminformatics descriptors, not experimental measurements. 3D conformers are illustrative predictions unless a deposited PubChem conformer is available.")
st.markdown("**Data and tools:** PubChem PUG REST for name-based retrieval; RDKit for SMILES parsing, 2D depiction, and calculated descriptors; py3Dmol/stmol for interactive 3D viewing.")
