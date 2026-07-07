import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import io
import sqlite3
import json
from datetime import datetime

# ── Page config ───────────────────────────────────────────────
st.set_page_config(
    page_title="PyMutScan – DNA Mutation Detector",
    page_icon="🧬",
    layout="wide"
)

# ══════════════════════════════════════════════════════════════
#  DATABASE SETUP & FUNCTIONS
# ══════════════════════════════════════════════════════════════

DB_PATH = "pymutscane_patients.db"

def init_db():
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS patients (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            name             TEXT NOT NULL,
            age              INTEGER,
            gender           TEXT,
            blood_group      TEXT,
            medical_history  TEXT,
            diagnosis        TEXT,
            gene             TEXT,
            sequence         TEXT,
            created_at       TEXT
        )
    ''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS analyses (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            patient_id      INTEGER,
            patient_name    TEXT,
            gene            TEXT,
            num_mutations   INTEGER,
            identity_pct    REAL,
            pathogenic      INTEGER,
            uncertain       INTEGER,
            benign          INTEGER,
            mutations_json  TEXT,
            created_at      TEXT,
            FOREIGN KEY(patient_id) REFERENCES patients(id)
        )
    ''')
    conn.commit()
    conn.close()

def save_patient_to_db(name, age, gender, blood_group, medical_history, diagnosis, gene, sequence):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("""
        INSERT INTO patients
            (name, age, gender, blood_group, medical_history, diagnosis, gene, sequence, created_at)
        VALUES (?,?,?,?,?,?,?,?,?)
    """, (name, age, gender, blood_group, medical_history, diagnosis, gene, sequence,
          datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    pid = c.lastrowid
    conn.commit()
    conn.close()
    return pid

def update_patient_sequence(patient_id, gene, sequence):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("UPDATE patients SET gene=?, sequence=? WHERE id=?", (gene, sequence, patient_id))
    conn.commit()
    conn.close()

def save_analysis_to_db(patient_id, patient_name, gene, mutations, identity_pct):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("""
        INSERT INTO analyses
            (patient_id, patient_name, gene, num_mutations, identity_pct,
             pathogenic, uncertain, benign, mutations_json, created_at)
        VALUES (?,?,?,?,?,?,?,?,?,?)
    """, (patient_id, patient_name, gene,
          len(mutations), round(identity_pct, 2),
          sum(1 for m in mutations if m["Clinical"] == "Pathogenic"),
          sum(1 for m in mutations if m["Clinical"] == "Uncertain"),
          sum(1 for m in mutations if m["Clinical"] == "Benign"),
          json.dumps(mutations),
          datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    conn.commit()
    conn.close()

def get_all_patients():
    conn = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql("SELECT * FROM patients ORDER BY created_at DESC", conn)
    except Exception:
        df = pd.DataFrame()
    conn.close()
    return df

def get_all_analyses():
    conn = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql("SELECT * FROM analyses ORDER BY created_at DESC", conn)
    except Exception:
        df = pd.DataFrame()
    conn.close()
    return df

def get_patient_analyses(patient_id):
    conn = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql(
            "SELECT * FROM analyses WHERE patient_id=? ORDER BY created_at DESC",
            conn, params=(patient_id,))
    except Exception:
        df = pd.DataFrame()
    conn.close()
    return df

def get_patient_by_id(patient_id):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * FROM patients WHERE id=?", (patient_id,))
    row = c.fetchone()
    cols = [desc[0] for desc in c.description] if c.description else []
    conn.close()
    return dict(zip(cols, row)) if row else None

def delete_patient_from_db(pid):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("DELETE FROM analyses WHERE patient_id=?", (pid,))
    c.execute("DELETE FROM patients WHERE id=?", (pid,))
    conn.commit()
    conn.close()

def delete_analysis_from_db(aid):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("DELETE FROM analyses WHERE id=?", (aid,))
    conn.commit()
    conn.close()

# Initialise database on every run
init_db()

# ══════════════════════════════════════════════════════════════
#  STATIC DATA
# ══════════════════════════════════════════════════════════════

REFERENCE_SEQUENCES = {
    "BRCA1": "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAG",
    "EGFR":  "ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTC",
    "TP53":  "ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATCCTGACTGTACCACCATCCACTACAACTACATG",
}

SAMPLE_MUTATIONS = {
    "BRCA1": [
        {"pos": 12, "ref": "G", "alt": "A", "type": "SNP", "effect": "Missense",  "clinical": "Pathogenic"},
        {"pos": 27, "ref": "C", "alt": "T", "type": "SNP", "effect": "Silent",    "clinical": "Benign"},
        {"pos": 45, "ref": "A", "alt": "G", "type": "SNP", "effect": "Missense",  "clinical": "Uncertain"},
    ],
    "EGFR": [
        {"pos": 8,  "ref": "C", "alt": "T", "type": "SNP", "effect": "Missense",  "clinical": "Pathogenic"},
        {"pos": 33, "ref": "G", "alt": "A", "type": "SNP", "effect": "Nonsense",  "clinical": "Pathogenic"},
    ],
    "TP53": [
        {"pos": 5,  "ref": "G", "alt": "C", "type": "SNP", "effect": "Missense",  "clinical": "Pathogenic"},
        {"pos": 19, "ref": "A", "alt": "T", "type": "SNP", "effect": "Silent",    "clinical": "Benign"},
        {"pos": 52, "ref": "C", "alt": "G", "type": "SNP", "effect": "Missense",  "clinical": "Uncertain"},
        {"pos": 61, "ref": "T", "alt": "A", "type": "SNP", "effect": "Missense",  "clinical": "Pathogenic"},
    ],
}

CODON_TABLE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "TAT":"Tyr","TAC":"Tyr","TAA":"Stop","TAG":"Stop",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "TGT":"Cys","TGC":"Cys","TGA":"Stop","TGG":"Trp",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

RESTRICTION_ENZYMES = {
    "EcoRI":  {"site": "GAATTC",   "cut": "G^AATTC"},
    "BamHI":  {"site": "GGATCC",   "cut": "G^GATCC"},
    "HindIII":{"site": "AAGCTT",   "cut": "A^AGCTT"},
    "NotI":   {"site": "GCGGCCGC", "cut": "GC^GGCCGC"},
    "XhoI":   {"site": "CTCGAG",   "cut": "C^TCGAG"},
    "NcoI":   {"site": "CCATGG",   "cut": "C^CATGG"},
    "SalI":   {"site": "GTCGAC",   "cut": "G^TCGAC"},
    "KpnI":   {"site": "GGTACC",   "cut": "GGTAC^C"},
    "SmaI":   {"site": "CCCGGG",   "cut": "CCC^GGG"},
    "ClaI":   {"site": "ATCGAT",   "cut": "AT^CGAT"},
    "EcoRV":  {"site": "GATATC",   "cut": "GAT^ATC"},
    "PstI":   {"site": "CTGCAG",   "cut": "CTGCA^G"},
}

DRUG_DATABASE = {
    "BRCA1": [
        {"mutation_pos": 12, "drug": "Olaparib (PARP inhibitor)", "effect": "Reduced sensitivity", "evidence": "Strong",   "class": "Targeted therapy"},
        {"mutation_pos": 12, "drug": "Cisplatin",                  "effect": "Increased sensitivity","evidence": "Moderate","class": "Chemotherapy"},
        {"mutation_pos": 45, "drug": "Talazoparib",                "effect": "Reduced sensitivity", "evidence": "Moderate","class": "Targeted therapy"},
    ],
    "EGFR": [
        {"mutation_pos": 8,  "drug": "Gefitinib (Iressa)",         "effect": "Resistance",          "evidence": "Strong",   "class": "TKI"},
        {"mutation_pos": 8,  "drug": "Erlotinib (Tarceva)",        "effect": "Resistance",          "evidence": "Strong",   "class": "TKI"},
        {"mutation_pos": 33, "drug": "Osimertinib (Tagrisso)",     "effect": "Sensitivity",         "evidence": "Strong",   "class": "3rd Gen TKI"},
    ],
    "TP53": [
        {"mutation_pos": 5,  "drug": "Doxorubicin",                "effect": "Reduced response",    "evidence": "Moderate", "class": "Chemotherapy"},
        {"mutation_pos": 5,  "drug": "APR-246 (Eprenetapopt)",     "effect": "Restored sensitivity","evidence": "Clinical trial","class": "p53 reactivator"},
        {"mutation_pos": 61, "drug": "5-Fluorouracil",             "effect": "Reduced response",    "evidence": "Moderate", "class": "Chemotherapy"},
    ],
}

# ══════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════

def parse_fasta(text):
    lines = text.strip().splitlines()
    seq = ""
    for line in lines:
        if not line.startswith(">"):
            seq += line.strip().upper()
    return "".join(c for c in seq if c in "ATCG")

def inject_mutations(ref_seq, gene):
    seq = list(ref_seq)
    for m in SAMPLE_MUTATIONS[gene]:
        seq[m["pos"] - 1] = m["alt"]
    return "".join(seq)

def detect_mutations(ref, query, gene):
    results = []
    for i in range(min(len(ref), len(query))):
        if ref[i] != query[i]:
            known = next((m for m in SAMPLE_MUTATIONS[gene] if m["pos"] == i + 1), None)
            results.append({
                "Position": i + 1,
                "Codon#":   (i // 3) + 1,
                "Ref":      ref[i],
                "Alt":      query[i],
                "Type":     known["type"]     if known else "SNP",
                "Effect":   known["effect"]   if known else "Missense",
                "Clinical": known["clinical"] if known else "Uncertain",
            })
    return results

def gc_content(seq):
    if not seq: return 0.0
    return round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)

def nucleotide_counts(seq):
    return {b: seq.count(b) for b in "ATCG"}

def translate_dna(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "?")
        protein.append(aa)
        if aa == "Stop": break
    return protein

def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G'}
    return ''.join(comp.get(b,'N') for b in reversed(seq))

def calc_tm(seq):
    a,t,g,c = seq.count('A'),seq.count('T'),seq.count('G'),seq.count('C')
    if len(seq) < 14: return 2*(a+t) + 4*(g+c)
    return 64.9 + 41*(g+c-16.4)/len(seq)

def calc_gc_primer(seq):
    return round((seq.count('G')+seq.count('C'))/len(seq)*100,1)

def color_clinical(val):
    colors = {
        "Pathogenic": "background-color:#ff4d6d33;color:#cc0000;font-weight:bold",
        "Benign":     "background-color:#00d97e33;color:#007a45;font-weight:bold",
        "Uncertain":  "background-color:#ffd16633;color:#a07800;font-weight:bold",
    }
    return colors.get(val, "")

def hamming_distance(s1, s2):
    length = min(len(s1), len(s2))
    return sum(c1 != c2 for c1, c2 in zip(s1[:length], s2[:length]))

def similarity_pct(s1, s2):
    length = min(len(s1), len(s2))
    matches = sum(c1 == c2 for c1, c2 in zip(s1[:length], s2[:length]))
    return round(matches / length * 100, 2) if length > 0 else 0

# ══════════════════════════════════════════════════════════════
#  NCBI LINKS
# ══════════════════════════════════════════════════════════════

NCBI_LINKS = {
    "BRCA1": {
        "gene_id":   "672",
        "refseq":    "NM_007294.4",
        "ncbi_url":  "https://www.ncbi.nlm.nih.gov/gene/672",
        "fasta_url": "https://www.ncbi.nlm.nih.gov/nuccore/NM_007294.4?report=fasta",
        "clinvar":   "https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1",
        "omim":      "https://omim.org/entry/113705",
        "full_length": "81,189 bp",
        "location":  "Chr 17q21.31",
        "function":  "DNA damage repair, tumor suppression",
    },
    "EGFR": {
        "gene_id":   "1956",
        "refseq":    "NM_005228.5",
        "ncbi_url":  "https://www.ncbi.nlm.nih.gov/gene/1956",
        "fasta_url": "https://www.ncbi.nlm.nih.gov/nuccore/NM_005228.5?report=fasta",
        "clinvar":   "https://www.ncbi.nlm.nih.gov/clinvar/?term=EGFR",
        "omim":      "https://omim.org/entry/131550",
        "full_length": "188,309 bp",
        "location":  "Chr 7p11.2",
        "function":  "Cell growth signaling, tyrosine kinase",
    },
    "TP53": {
        "gene_id":   "7157",
        "refseq":    "NM_000546.6",
        "ncbi_url":  "https://www.ncbi.nlm.nih.gov/gene/7157",
        "fasta_url": "https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.6?report=fasta",
        "clinvar":   "https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53",
        "omim":      "https://omim.org/entry/191170",
        "full_length": "19,149 bp",
        "location":  "Chr 17p13.1",
        "function":  "Guardian of the genome, apoptosis",
    },
}

# ══════════════════════════════════════════════════════════════
#  SIDEBAR
# ══════════════════════════════════════════════════════════════

GENE_OPTIONS = ["BRCA1", "EGFR", "TP53"]

with st.sidebar:
    st.title("🧬 PyMutScan")
    st.caption("Advanced DNA Mutation Detector")
    st.divider()

    # Gene selector — respects loaded patient gene
    gene_idx = st.session_state.get("loaded_gene_idx", 0)
    gene = st.selectbox("🔖 Select Reference Gene", GENE_OPTIONS, index=gene_idx)

    # ── Current Patient Badge ─────────────────────────────────
    if st.session_state.get("current_patient"):
        cp = st.session_state["current_patient"]
        st.divider()
        st.markdown("**👤 Active Patient**")
        st.success(
            f"**{cp['name']}**  \n"
            f"Age: {cp['age']} · {cp['gender']}  \n"
            f"Blood: {cp.get('blood_group','—')}  \n"
            f"Gene: {cp.get('gene', gene)}"
        )
        if st.button("❌ Clear Patient", use_container_width=True):
            st.session_state.pop("current_patient", None)
            st.session_state.pop("loaded_gene_idx", None)
            st.rerun()

    # ── NCBI Source Links ─────────────────────────────────────
    st.divider()
    st.markdown("**🔗 NCBI Reference Sources**")
    info = NCBI_LINKS[gene]
    st.markdown(f"**Gene:** {gene} | ID: `{info['gene_id']}`")
    st.markdown(f"**RefSeq:** `{info['refseq']}`")
    st.markdown(f"**Location:** {info['location']}")
    st.markdown(f"**Full Length:** {info['full_length']}")
    st.markdown(f"**Function:** {info['function']}")
    st.markdown(f"[🌐 NCBI Gene Page]({info['ncbi_url']})")
    st.markdown(f"[📄 FASTA Sequence]({info['fasta_url']})")
    st.markdown(f"[🔬 ClinVar Variants]({info['clinvar']})")
    st.markdown(f"[📚 OMIM Entry]({info['omim']})")

    st.divider()
    st.markdown("**📂 Upload FASTA File**")
    uploaded_file = st.file_uploader("Upload .fasta or .txt", type=["fasta","fa","txt"])
    st.divider()
    st.markdown("**ℹ️ About**")
    st.markdown(
        "PyMutScan aligns patient DNA against a reference gene, detects SNPs, "
        "translates codons, designs primers, predicts drug resistance, visualizes "
        "mutation hotspots, and stores patient records in a local SQLite database."
    )

# ══════════════════════════════════════════════════════════════
#  MAIN HEADER
# ══════════════════════════════════════════════════════════════

st.title("🧬 PyMutScan — Advanced DNA Mutation Detector")
st.caption("Sequence alignment · Variant calling · Codon translation · Primer design · Drug resistance · Mutation heatmap · Patient records")
st.divider()

ref_seq = REFERENCE_SEQUENCES[gene]

# ══════════════════════════════════════════════════════════════
#  👤  VIRTUAL PATIENT MANAGER  (always visible)
# ══════════════════════════════════════════════════════════════

with st.expander("👤 Virtual Patient Manager — Create or Load a Patient Profile", expanded=False):
    vp_new, vp_load = st.tabs(["➕ Create New Patient", "📂 Load Existing Patient"])

    # ─── CREATE NEW PATIENT ────────────────────────────────────
    with vp_new:
        st.markdown("Fill in patient details below. **Gene** and **Sequence** will be pulled from the current sidebar selection and query box when you click Save.")

        vp_c1, vp_c2, vp_c3, vp_c4 = st.columns(4)
        vp_name   = vp_c1.text_input("Patient Name *", placeholder="e.g. Ravi Kumar")
        vp_age    = vp_c2.number_input("Age", min_value=1, max_value=120, value=35)
        vp_gender = vp_c3.selectbox("Gender", ["Male", "Female", "Other", "Prefer not to say"])
        vp_blood  = vp_c4.selectbox("Blood Group", ["A+","A−","B+","B−","AB+","AB−","O+","O−","Unknown"])

        vp_c5, vp_c6 = st.columns(2)
        vp_history   = vp_c5.text_area("Medical History / Symptoms", placeholder="e.g. Family history of breast cancer, BRCA1 mutation suspected", height=100)
        vp_diagnosis = vp_c6.text_area("Preliminary Diagnosis / Notes", placeholder="e.g. Referred for pharmacogenomics screening", height=100)

        st.caption("⚠️ The current Gene (sidebar) and Query Sequence (below) will be saved with this patient profile.")

        vp_save_btn = st.button("💾 Save Patient Profile", type="primary", use_container_width=False)

        if vp_save_btn:
            if not vp_name.strip():
                st.error("Patient name is required.")
            else:
                # Sequence from session state / query (we grab it at save time — empty string is fine)
                current_seq = st.session_state.get("query", "")
                pid = save_patient_to_db(
                    vp_name.strip(), vp_age, vp_gender, vp_blood,
                    vp_history, vp_diagnosis, gene, current_seq
                )
                st.session_state["current_patient"] = {
                    "id": pid, "name": vp_name.strip(),
                    "age": vp_age, "gender": vp_gender,
                    "blood_group": vp_blood, "gene": gene
                }
                st.success(f"✅ Patient **{vp_name.strip()}** saved! (ID: {pid}). They are now the active patient.")
                st.rerun()

    # ─── LOAD EXISTING PATIENT ─────────────────────────────────
    with vp_load:
        patients_df = get_all_patients()

        if patients_df.empty:
            st.info("No saved patients yet. Create one in the 'Create New Patient' tab.")
        else:
            # Build display list
            patient_options = {
                f"{row['name']} (ID: {row['id']}) — {row['gene']} — {row['created_at'][:10]}": row['id']
                for _, row in patients_df.iterrows()
            }
            selected_label = st.selectbox("Select a patient to load", list(patient_options.keys()))
            selected_id    = patient_options[selected_label]
            patient_data   = get_patient_by_id(selected_id)

            if patient_data:
                # Show patient card
                lc1, lc2, lc3, lc4 = st.columns(4)
                lc1.metric("Name",        patient_data["name"])
                lc2.metric("Age",         patient_data["age"])
                lc3.metric("Gender",      patient_data["gender"])
                lc4.metric("Blood Group", patient_data.get("blood_group", "—"))

                lc5, lc6 = st.columns(2)
                lc5.markdown(f"**Medical History:** {patient_data.get('medical_history','—')}")
                lc6.markdown(f"**Diagnosis:** {patient_data.get('diagnosis','—')}")

                st.markdown(f"**Saved Gene:** `{patient_data.get('gene','—')}`")
                saved_seq = patient_data.get("sequence","")
                if saved_seq:
                    st.markdown(f"**Saved Sequence ({len(saved_seq)} bp):**")
                    st.code(saved_seq[:80] + ("..." if len(saved_seq)>80 else ""), language=None)
                else:
                    st.caption("No sequence saved with this patient.")

                # Past analyses for this patient
                past = get_patient_analyses(selected_id)
                if not past.empty:
                    st.markdown(f"**📊 Past Analyses ({len(past)}):**")
                    display_cols = ["gene","num_mutations","identity_pct","pathogenic","uncertain","benign","created_at"]
                    st.dataframe(past[[c for c in display_cols if c in past.columns]], use_container_width=True)

                btn_col1, btn_col2 = st.columns(2)
                with btn_col1:
                    if st.button("🔄 Load Patient into Analysis", type="primary", use_container_width=True):
                        # Set active patient
                        st.session_state["current_patient"] = {
                            "id":         patient_data["id"],
                            "name":       patient_data["name"],
                            "age":        patient_data["age"],
                            "gender":     patient_data["gender"],
                            "blood_group":patient_data.get("blood_group",""),
                            "gene":       patient_data.get("gene","BRCA1"),
                        }
                        # Auto-fill gene and sequence
                        gene_name = patient_data.get("gene","BRCA1")
                        if gene_name in GENE_OPTIONS:
                            st.session_state["loaded_gene_idx"] = GENE_OPTIONS.index(gene_name)
                        if saved_seq:
                            st.session_state["query"] = saved_seq
                        st.success(f"✅ Patient **{patient_data['name']}** loaded! Gene and sequence pre-filled.")
                        st.rerun()

                with btn_col2:
                    if st.button("🗑️ Delete Patient", type="secondary", use_container_width=True):
                        delete_patient_from_db(selected_id)
                        if st.session_state.get("current_patient", {}).get("id") == selected_id:
                            st.session_state.pop("current_patient", None)
                        st.warning(f"Patient ID {selected_id} deleted.")
                        st.rerun()

# ══════════════════════════════════════════════════════════════
#  REFERENCE & QUERY INPUT
# ══════════════════════════════════════════════════════════════

col1, col2 = st.columns(2)

with col1:
    st.subheader(f"📄 Reference — {gene}")
    st.code(ref_seq, language=None)
    st.caption(f"{len(ref_seq)} bp demo segment · Source: NCBI GenBank `{NCBI_LINKS[gene]['refseq']}`")
    with st.expander("📖 View Full Gene Info & Sources"):
        info = NCBI_LINKS[gene]
        c1, c2 = st.columns(2)
        c1.metric("Full Gene Length", info["full_length"])
        c2.metric("Chromosomal Location", info["location"])
        st.markdown(f"**Function:** {info['function']}")
        st.markdown(f"**RefSeq Accession:** `{info['refseq']}`")
        st.markdown(f"**NCBI Gene ID:** `{info['gene_id']}`")
        st.markdown("**⚠️ Note:** PyMutScan uses a 67 bp demo segment for educational purposes. Full clinical sequences available at NCBI.")
        col_a, col_b, col_c, col_d = st.columns(4)
        col_a.markdown(f"[🌐 NCBI Gene]({info['ncbi_url']})")
        col_b.markdown(f"[📄 FASTA]({info['fasta_url']})")
        col_c.markdown(f"[🔬 ClinVar]({info['clinvar']})")
        col_d.markdown(f"[📚 OMIM]({info['omim']})")

with col2:
    st.subheader("🔬 Query Sequence (Patient)")
    prefill = ""
    if uploaded_file:
        raw = uploaded_file.read().decode("utf-8")
        prefill = parse_fasta(raw)
        st.success(f"✅ FASTA loaded — {len(prefill)} bp detected")
    if st.button("⚡ Load Sample with Mutations"):
        st.session_state["query"] = inject_mutations(ref_seq, gene)
    query_raw = st.text_area(
        "Paste DNA sequence or upload FASTA above (A, T, C, G only)",
        value=prefill if prefill else st.session_state.get("query", ""),
        height=120,
    )
    query_seq = "".join(c for c in query_raw.upper() if c in "ATCG")
    st.caption(f"{len(query_seq)} bp entered")

    # Show active patient info inline
    if st.session_state.get("current_patient"):
        cp = st.session_state["current_patient"]
        st.info(f"👤 Active patient: **{cp['name']}** · Analysis will be auto-saved to database")

st.divider()

# ── Update patient sequence in DB if sequence changed and patient loaded ──
if st.session_state.get("current_patient") and query_seq:
    cp = st.session_state["current_patient"]
    # Only update if sequence is not empty — silently keep DB fresh
    if query_seq != st.session_state.get("_last_saved_seq", ""):
        update_patient_sequence(cp["id"], gene, query_seq)
        st.session_state["_last_saved_seq"] = query_seq

# ══════════════════════════════════════════════════════════════
#  RUN ANALYSIS BUTTON
# ══════════════════════════════════════════════════════════════

run = st.button("🔍 Run Mutation Analysis", type="primary", disabled=len(query_seq)==0)

if run:
    with st.spinner("Running alignment pipeline..."):
        mutations = detect_mutations(ref_seq, query_seq, gene)

    identity = ((len(query_seq)-len(mutations))/len(query_seq)*100) if query_seq else 0
    st.success(f"✅ Analysis complete — {len(mutations)} variant(s) detected")

    # ── Auto-save to DB if a patient is active ────────────────
    if st.session_state.get("current_patient"):
        cp = st.session_state["current_patient"]
        save_analysis_to_db(cp["id"], cp["name"], gene, mutations, identity)
        st.toast(f"💾 Analysis saved to database for patient {cp['name']}", icon="🗄️")

    m1,m2,m3,m4,m5 = st.columns(5)
    m1.metric("Total Variants",    len(mutations))
    m2.metric("Sequence Identity", f"{identity:.1f}%")
    m3.metric("🔴 Pathogenic",     sum(1 for m in mutations if m["Clinical"]=="Pathogenic"))
    m4.metric("🟡 Uncertain",      sum(1 for m in mutations if m["Clinical"]=="Uncertain"))
    m5.metric("🟢 Benign",         sum(1 for m in mutations if m["Clinical"]=="Benign"))

    st.divider()

    # ── All 9 Analysis Tabs ───────────────────────────────────
    tab1,tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9 = st.tabs([
        "🔬 Alignment",
        "📋 Mutations Table",
        "🧪 Sequence Stats & GC",
        "🔤 Codon Translation",
        "🗺️ Mutation Heatmap",
        "🧲 Primer Design",
        "✂️ Restriction Enzymes",
        "💊 Drug Resistance",
        "🌿 Phylogenetic Tree",
    ])

    # ── Tab 1: Alignment ──────────────────────────────────────
    with tab1:
        st.subheader("Pairwise Alignment")
        align_ref   = ref_seq[:len(query_seq)]
        align_query = query_seq[:len(ref_seq)]
        align_bar   = "".join("|" if align_ref[i]==align_query[i] else "X" for i in range(min(len(align_ref),len(align_query))))
        st.code(f"REF [{gene}]   {align_ref}\n               {align_bar}\nQRY [Patient]  {align_query}", language=None)
        st.caption("| = match    X = mismatch (mutation)")

    # ── Tab 2: Mutations Table ────────────────────────────────
    with tab2:
        st.subheader(f"Detected Variants — {gene}")
        if mutations:
            df = pd.DataFrame(mutations)
            st.dataframe(df.style.applymap(color_clinical, subset=["Clinical"]), use_container_width=True)
            st.download_button("⬇️ Download Results as CSV", data=df.to_csv(index=False), file_name=f"PyMutScan_{gene}_results.csv", mime="text/csv")
        else:
            st.info("No mutations detected — sequences are identical.")

    # ── Tab 3: Sequence Stats & GC ────────────────────────────
    with tab3:
        st.subheader("Sequence Statistics")
        c1,c2 = st.columns(2)
        with c1:
            st.markdown("**Reference Sequence**")
            ref_gc = gc_content(ref_seq)
            st.metric("GC Content", f"{ref_gc}%")
            st.metric("AT Content", f"{round(100-ref_gc,2)}%")
            st.metric("Length", f"{len(ref_seq)} bp")
            ref_nt_df = pd.DataFrame(nucleotide_counts(ref_seq).items(), columns=["Nucleotide","Count"])
            fig_ref = px.bar(ref_nt_df, x="Nucleotide", y="Count", color="Nucleotide",
                color_discrete_map={"A":"#00b4d8","T":"#ff6b6b","G":"#00d97e","C":"#ffd166"},
                title="Nucleotide Distribution — Reference")
            st.plotly_chart(fig_ref, use_container_width=True)
        with c2:
            st.markdown("**Query (Patient) Sequence**")
            qry_gc = gc_content(query_seq)
            st.metric("GC Content", f"{qry_gc}%")
            st.metric("AT Content", f"{round(100-qry_gc,2)}%")
            st.metric("Length", f"{len(query_seq)} bp")
            qry_nt_df = pd.DataFrame(nucleotide_counts(query_seq).items(), columns=["Nucleotide","Count"])
            fig_qry = px.bar(qry_nt_df, x="Nucleotide", y="Count", color="Nucleotide",
                color_discrete_map={"A":"#00b4d8","T":"#ff6b6b","G":"#00d97e","C":"#ffd166"},
                title="Nucleotide Distribution — Query")
            st.plotly_chart(fig_qry, use_container_width=True)
        st.markdown("**GC Content Sliding Window (10 bp)**")
        window = 10
        gc_vals = [gc_content(query_seq[i:i+window]) for i in range(0, len(query_seq)-window+1)]
        gc_df = pd.DataFrame({"Position": range(1,len(gc_vals)+1), "GC%": gc_vals})
        fig_gc = px.line(gc_df, x="Position", y="GC%", title="GC% Along Query Sequence")
        fig_gc.add_hline(y=50, line_dash="dash", line_color="gray", annotation_text="50%")
        st.plotly_chart(fig_gc, use_container_width=True)

    # ── Tab 4: Codon Translation ──────────────────────────────
    with tab4:
        st.subheader("DNA → Protein Translation")
        c1,c2 = st.columns(2)
        with c1:
            st.markdown(f"**Reference Protein — {gene}**")
            ref_protein = translate_dna(ref_seq)
            st.code(" – ".join(ref_protein), language=None)
            st.caption(f"{len(ref_protein)} amino acids")
        with c2:
            st.markdown("**Query Protein (Patient)**")
            qry_protein = translate_dna(query_seq)
            st.code(" – ".join(qry_protein), language=None)
            st.caption(f"{len(qry_protein)} amino acids")
        st.markdown("**Amino Acid Changes Due to Mutations**")
        aa_changes = []
        for i in range(min(len(ref_protein),len(qry_protein))):
            if ref_protein[i] != qry_protein[i]:
                aa_changes.append({"Codon #":i+1,"Ref AA":ref_protein[i],"Alt AA":qry_protein[i],"Change":f"{ref_protein[i]}{i+1}{qry_protein[i]}"})
        if aa_changes:
            st.dataframe(pd.DataFrame(aa_changes), use_container_width=True)
        else:
            st.info("No amino acid changes detected.")

    # ── Tab 5: Mutation Heatmap ───────────────────────────────
    with tab5:
        st.subheader("Mutation Position Heatmap")
        if mutations:
            seq_len  = min(len(ref_seq), len(query_seq))
            mut_pos  = {m["Position"]: m["Clinical"] for m in mutations}
            clinical_map = {"Pathogenic":3,"Uncertain":2,"Benign":1}
            heatmap_vals,hover_texts = [],[]
            for i in range(1, seq_len+1):
                if i in mut_pos:
                    heatmap_vals.append(clinical_map.get(mut_pos[i],2))
                    hover_texts.append(f"Pos {i}: {mut_pos[i]}")
                else:
                    heatmap_vals.append(0)
                    hover_texts.append(f"Pos {i}: No mutation")
            cols_n = 10
            rows_n = (seq_len+cols_n-1)//cols_n
            grid,htxt = [],[]
            for r in range(rows_n):
                row_vals = heatmap_vals[r*cols_n:(r+1)*cols_n]
                row_txt  = hover_texts[r*cols_n:(r+1)*cols_n]
                while len(row_vals) < cols_n: row_vals.append(-1); row_txt.append("")
                grid.append(row_vals); htxt.append(row_txt)
            fig_hm = go.Figure(go.Heatmap(z=grid,text=htxt,hoverinfo="text",
                colorscale=[[0.0,"#1e293b"],[0.25,"#1e293b"],[0.35,"#00d97e"],[0.55,"#ffd166"],[0.75,"#ff4d6d"],[1.0,"#ff4d6d"]],
                showscale=False,xgap=2,ygap=2))
            fig_hm.update_layout(title=f"Mutation Heatmap — {gene}",height=350)
            st.plotly_chart(fig_hm, use_container_width=True)
            st.caption("🟢 Benign   🟡 Uncertain   🔴 Pathogenic   ⬛ No mutation")
            mut_df = pd.DataFrame(mutations)
            color_map = {"Pathogenic":"#ff4d6d","Benign":"#00d97e","Uncertain":"#ffd166"}
            fig_bar = px.scatter(mut_df, x="Position", y="Effect", color="Clinical",
                color_discrete_map=color_map, size=[15]*len(mut_df),
                hover_data=["Ref","Alt","Clinical"], title="Mutation Positions by Effect Type")
            st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.info("No mutations to visualize.")

    # ── Tab 6: Primer Design ──────────────────────────────────
    with tab6:
        st.subheader("🧲 Primer Design Tool")
        st.caption("Design forward and reverse primers from your query sequence")
        primer_len    = st.slider("Primer Length (bp)", 18, 25, 20)
        amplicon_size = st.slider("Amplicon Size (bp)", 100, 500, 200)
        if len(query_seq) >= amplicon_size:
            fwd_primer = query_seq[:primer_len]
            rev_primer = reverse_complement(query_seq[amplicon_size-primer_len:amplicon_size])
            c1,c2 = st.columns(2)
            with c1:
                st.markdown("**→ Forward Primer**")
                st.code(fwd_primer, language=None)
                st.metric("Tm (°C)",    f"{calc_tm(fwd_primer):.1f}")
                st.metric("GC Content", f"{calc_gc_primer(fwd_primer)}%")
                st.metric("Length",     f"{len(fwd_primer)} bp")
            with c2:
                st.markdown("**← Reverse Primer**")
                st.code(rev_primer, language=None)
                st.metric("Tm (°C)",    f"{calc_tm(rev_primer):.1f}")
                st.metric("GC Content", f"{calc_gc_primer(rev_primer)}%")
                st.metric("Length",     f"{len(rev_primer)} bp")
            st.info(f"Expected amplicon size: **{amplicon_size} bp**")
            st.markdown("**Primer Quality Check**")
            checks = {
                "GC content 40–60%":     40 <= calc_gc_primer(fwd_primer) <= 60,
                "Tm between 55–65°C":    55 <= calc_tm(fwd_primer) <= 65,
                "No poly-A/T runs (>4)": "AAAAA" not in fwd_primer and "TTTTT" not in fwd_primer,
                "No poly-G/C runs (>4)": "GGGGG" not in fwd_primer and "CCCCC" not in fwd_primer,
            }
            for check, passed in checks.items():
                st.write(f"{'✅' if passed else '❌'} {check}")
            primer_df = pd.DataFrame({
                "Primer":["Forward","Reverse"],
                "Sequence":[fwd_primer,rev_primer],
                "Length (bp)":[len(fwd_primer),len(rev_primer)],
                "Tm (°C)":[f"{calc_tm(fwd_primer):.1f}",f"{calc_tm(rev_primer):.1f}"],
                "GC%":[f"{calc_gc_primer(fwd_primer)}%",f"{calc_gc_primer(rev_primer)}%"],
            })
            st.download_button("⬇️ Download Primers as CSV", data=primer_df.to_csv(index=False), file_name=f"PyMutScan_{gene}_primers.csv", mime="text/csv")
        else:
            st.warning(f"Query sequence must be at least {amplicon_size} bp for primer design.")

    # ── Tab 7: Restriction Enzymes ────────────────────────────
    with tab7:
        st.subheader("✂️ Restriction Enzyme Analysis")
        st.caption("Find restriction enzyme cut sites in your DNA sequence")
        found_enzymes = []
        for enzyme, info in RESTRICTION_ENZYMES.items():
            site  = info["site"]
            count = query_seq.count(site)
            positions = [i+1 for i in range(len(query_seq)) if query_seq[i:i+len(site)]==site]
            if count > 0:
                found_enzymes.append({"Enzyme":enzyme,"Cut Site":info["cut"],"Site Seq":site,"# Cuts":count,"Positions":", ".join(map(str,positions)),"Fragments":count+1})
        if found_enzymes:
            st.success(f"✅ {len(found_enzymes)} restriction enzyme(s) found!")
            df_re = pd.DataFrame(found_enzymes)
            st.dataframe(df_re, use_container_width=True)
            fig_re = px.bar(df_re, x="Enzyme", y="# Cuts", color="# Cuts", color_continuous_scale="teal", title="Restriction Enzyme Cut Frequency")
            st.plotly_chart(fig_re, use_container_width=True)
            st.download_button("⬇️ Download Restriction Sites CSV", data=df_re.to_csv(index=False), file_name=f"PyMutScan_{gene}_restriction_sites.csv", mime="text/csv")
        else:
            st.info("No common restriction enzyme cut sites found in this sequence.")

    # ── Tab 8: Drug Resistance ────────────────────────────────
    with tab8:
        st.subheader("💊 Drug Resistance Predictor")
        st.caption("Predict drug response based on detected mutations")
        if mutations:
            mut_positions = {m["Position"] for m in mutations}
            drug_results  = [d for d in DRUG_DATABASE.get(gene,[]) if d["mutation_pos"] in mut_positions]
            if drug_results:
                st.warning(f"⚠️ {len(drug_results)} drug interaction(s) detected!")
                for drug in drug_results:
                    effect_color = "🔴" if "Resistance" in drug["effect"] or "Reduced" in drug["effect"] else "🟢"
                    with st.expander(f"{effect_color} {drug['drug']} — {drug['effect']}"):
                        c1,c2,c3 = st.columns(3)
                        c1.metric("Drug Class",     drug["class"])
                        c2.metric("Effect",         drug["effect"])
                        c3.metric("Evidence Level", drug["evidence"])
                        st.caption(f"Triggered by mutation at position {drug['mutation_pos']}")
                drug_df = pd.DataFrame(drug_results)
                st.dataframe(drug_df, use_container_width=True)
                st.download_button("⬇️ Download Drug Report CSV", data=drug_df.to_csv(index=False), file_name=f"PyMutScan_{gene}_drug_resistance.csv", mime="text/csv")
                st.error("⚠️ Disclaimer: For educational purposes only. Always consult a clinical oncologist.")
            else:
                st.success("✅ No known drug resistance mutations detected.")
        else:
            st.info("Run mutation analysis first to see drug resistance predictions.")

    # ── Tab 9: Phylogenetic Tree ──────────────────────────────
    with tab9:
        st.subheader("🌿 Phylogenetic Analysis")
        st.caption("Compare evolutionary distance between all reference sequences and your query")
        sequences = {
            "BRCA1 (Ref)":   REFERENCE_SEQUENCES["BRCA1"],
            "EGFR (Ref)":    REFERENCE_SEQUENCES["EGFR"],
            "TP53 (Ref)":    REFERENCE_SEQUENCES["TP53"],
            "Patient Query": query_seq,
        }
        labels = list(sequences.keys())
        seqs   = list(sequences.values())
        n      = len(labels)
        sim_matrix = [[similarity_pct(seqs[i],seqs[j]) for j in range(n)] for i in range(n)]
        matrix     = [[round(hamming_distance(seqs[i],seqs[j])/min(len(seqs[i]),len(seqs[j]))*100,2) if i!=j else 0.0 for j in range(n)] for i in range(n)]
        st.markdown("**Pairwise Distance Matrix (%)**")
        dist_df = pd.DataFrame(matrix, index=labels, columns=labels)
        st.dataframe(dist_df.style.background_gradient(cmap="RdYlGn_r"), use_container_width=True)
        fig_phy = go.Figure(go.Heatmap(
            z=sim_matrix, x=labels, y=labels, colorscale="Greens",
            text=[[f"{v}%" for v in row] for row in sim_matrix],
            texttemplate="%{text}", showscale=True))
        fig_phy.update_layout(title="Sequence Similarity Matrix (%)", height=400)
        st.plotly_chart(fig_phy, use_container_width=True)
        st.markdown("**Evolutionary Closest Match to Patient Query**")
        ref_distances = {g: hamming_distance(query_seq, REFERENCE_SEQUENCES[g]) for g in ["BRCA1","EGFR","TP53"]}
        closest = min(ref_distances, key=ref_distances.get)
        c1,c2,c3 = st.columns(3)
        for col,(g_name,dist) in zip([c1,c2,c3],ref_distances.items()):
            col.metric(g_name, f"{dist} mismatches", delta="Closest ✅" if g_name==closest else None)
        st.success(f"🌿 Patient sequence is evolutionarily closest to **{closest}**")


# ══════════════════════════════════════════════════════════════
#  🗄️  PATIENT RECORDS & DATABASE  (always visible)
# ══════════════════════════════════════════════════════════════

st.divider()
st.subheader("🗄️ Patient Records & Database")
st.caption("Persistent SQLite storage — all patient profiles and analysis results are saved locally in `pymutscane_patients.db`")

pr_tab1, pr_tab2, pr_tab3 = st.tabs(["👥 Patient Profiles", "📊 Analysis History", "📤 Export Data"])

# ─── TAB 1 : All Patient Profiles ─────────────────────────────
with pr_tab1:
    all_patients = get_all_patients()
    if all_patients.empty:
        st.info("No patients saved yet. Use the **Virtual Patient Manager** above to create one.")
    else:
        st.success(f"✅ {len(all_patients)} patient(s) in database")

        # Summary cards
        cols_per_row = 3
        rows = [all_patients.iloc[i:i+cols_per_row] for i in range(0, len(all_patients), cols_per_row)]
        for row_df in rows:
            cols = st.columns(cols_per_row)
            for col, (_, p) in zip(cols, row_df.iterrows()):
                with col:
                    past_count = len(get_patient_analyses(p["id"]))
                    col.markdown(
                        f"**{p['name']}**  \n"
                        f"🩸 {p.get('blood_group','—')} · {p['gender']} · Age {p['age']}  \n"
                        f"🧬 Gene: `{p.get('gene','—')}`  \n"
                        f"📊 Analyses: **{past_count}**  \n"
                        f"🕒 {str(p.get('created_at',''))[:10]}"
                    )

        st.divider()
        st.markdown("**Full Patient Table**")
        display_cols = ["id","name","age","gender","blood_group","diagnosis","gene","created_at"]
        st.dataframe(all_patients[[c for c in display_cols if c in all_patients.columns]],
                     use_container_width=True)

        # Delete individual patients
        st.markdown("**🗑️ Delete a Patient**")
        del_id = st.number_input("Enter Patient ID to delete", min_value=1, step=1, value=1)
        if st.button("Delete Patient & All Their Analyses", type="secondary"):
            delete_patient_from_db(del_id)
            st.warning(f"Patient ID {del_id} and all their analyses deleted.")
            st.rerun()

# ─── TAB 2 : Analysis History ─────────────────────────────────
with pr_tab2:
    all_analyses = get_all_analyses()
    if all_analyses.empty:
        st.info("No analyses saved yet. Run an analysis with an active patient to save results.")
    else:
        st.success(f"✅ {len(all_analyses)} analysis record(s) in database")

        # Summary metrics
        a1,a2,a3,a4 = st.columns(4)
        a1.metric("Total Analyses",       len(all_analyses))
        a2.metric("Unique Patients",       all_analyses["patient_name"].nunique() if "patient_name" in all_analyses else "—")
        a3.metric("Avg Mutations/Analysis",f"{all_analyses['num_mutations'].mean():.1f}" if "num_mutations" in all_analyses else "—")
        a4.metric("Total Pathogenic Found",int(all_analyses['pathogenic'].sum()) if "pathogenic" in all_analyses else "—")

        # Analyses table
        disp = ["id","patient_name","gene","num_mutations","identity_pct","pathogenic","uncertain","benign","created_at"]
        st.dataframe(all_analyses[[c for c in disp if c in all_analyses.columns]],
                     use_container_width=True)

        # Breakdown chart
        if "gene" in all_analyses.columns and len(all_analyses) > 0:
            gene_counts = all_analyses["gene"].value_counts().reset_index()
            gene_counts.columns = ["Gene","Analyses"]
            fig_g = px.pie(gene_counts, names="Gene", values="Analyses",
                           title="Analyses by Gene", color_discrete_sequence=["#00d97e","#ffd166","#ff4d6d"])
            st.plotly_chart(fig_g, use_container_width=True)

        # Detailed view per patient
        st.divider()
        st.markdown("**🔍 Filter by Patient**")
        if "patient_name" in all_analyses.columns:
            patient_filter = st.selectbox("Select patient", ["All"] + sorted(all_analyses["patient_name"].unique().tolist()))
            if patient_filter != "All":
                filtered = all_analyses[all_analyses["patient_name"] == patient_filter]
                st.dataframe(filtered[[c for c in disp if c in filtered.columns]], use_container_width=True)

                # Expand mutations detail for last analysis
                if not filtered.empty and "mutations_json" in filtered.columns:
                    last = filtered.iloc[0]
                    with st.expander(f"📋 Mutation Detail — Latest Analysis ({last['created_at'][:10]})"):
                        try:
                            muts = json.loads(last["mutations_json"])
                            if muts:
                                st.dataframe(pd.DataFrame(muts).style.applymap(
                                    color_clinical, subset=["Clinical"]), use_container_width=True)
                            else:
                                st.info("No mutations in this analysis.")
                        except Exception:
                            st.warning("Could not parse mutation details.")

        # Delete analysis by ID
        st.divider()
        st.markdown("**🗑️ Delete an Analysis Record**")
        del_aid = st.number_input("Enter Analysis ID to delete", min_value=1, step=1, value=1, key="del_analysis_id")
        if st.button("Delete Analysis Record", type="secondary", key="del_analysis_btn"):
            delete_analysis_from_db(del_aid)
            st.warning(f"Analysis ID {del_aid} deleted.")
            st.rerun()

# ─── TAB 3 : Export Data ──────────────────────────────────────
with pr_tab3:
    st.markdown("### 📤 Export All Data")

    all_patients_exp  = get_all_patients()
    all_analyses_exp  = get_all_analyses()

    ec1, ec2 = st.columns(2)

    with ec1:
        st.markdown("**👥 Patient Profiles**")
        if not all_patients_exp.empty:
            # CSV
            st.download_button(
                "⬇️ Download Patients CSV",
                data=all_patients_exp.to_csv(index=False),
                file_name="pymutscane_patients.csv",
                mime="text/csv",
                use_container_width=True,
            )
            # JSON
            st.download_button(
                "⬇️ Download Patients JSON",
                data=all_patients_exp.to_json(orient="records", indent=2),
                file_name="pymutscane_patients.json",
                mime="application/json",
                use_container_width=True,
            )
        else:
            st.info("No patient data to export.")

    with ec2:
        st.markdown("**📊 Analysis History**")
        if not all_analyses_exp.empty:
            # CSV (without mutations_json column for readability)
            export_cols = [c for c in all_analyses_exp.columns if c != "mutations_json"]
            st.download_button(
                "⬇️ Download Analyses CSV",
                data=all_analyses_exp[export_cols].to_csv(index=False),
                file_name="pymutscane_analyses.csv",
                mime="text/csv",
                use_container_width=True,
            )
            # JSON (full, with mutations)
            st.download_button(
                "⬇️ Download Analyses JSON (Full)",
                data=all_analyses_exp.to_json(orient="records", indent=2),
                file_name="pymutscane_analyses_full.json",
                mime="application/json",
                use_container_width=True,
            )
        else:
            st.info("No analysis data to export.")

    # ── Combined export ────────────────────────────────────────
    st.divider()
    st.markdown("**📦 Combined Full Export (JSON)**")
    if not all_patients_exp.empty or not all_analyses_exp.empty:
        combined = {
            "exported_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "patients":    json.loads(all_patients_exp.to_json(orient="records")) if not all_patients_exp.empty else [],
            "analyses":    json.loads(all_analyses_exp.to_json(orient="records")) if not all_analyses_exp.empty else [],
        }
        st.download_button(
            "⬇️ Download Complete Database as JSON",
            data=json.dumps(combined, indent=2),
            file_name="pymutscane_full_export.json",
            mime="application/json",
            use_container_width=False,
        )
    else:
        st.info("Database is empty. No data to export yet.")

    st.caption("⚠️ Database is stored locally as `pymutscane_patients.db` in the same directory as this script. Back it up regularly.")
