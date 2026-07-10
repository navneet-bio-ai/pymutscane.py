🧬** PyMutScan — DNA Mutation Detection & Variant Analysis Platform.

🔗 Live App:https://pymutscanepy-amd.streamlit.app/

A Streamlit-based platform for detecting and interpreting genetic variants from patient-derived DNA sequences. PyMutScan performs pairwise sequence alignment against a reference gene, calls SNPs, translates affected codons, classifies variant impact, and maintains persistent patient-level records — packaged as a single interactive application.


📑 Table of Contents


Overview
Core Features
Workflow
Tech Stack
Project Structure
Installation
Example Run
Data & Storage
Limitations
Roadmap
License



🔍 Overview

Genetic variant interpretation typically requires stitching together separate tools — one for alignment, another for translation, another for annotation. PyMutScan consolidates this into a single pipeline: a patient's query sequence is aligned against a reference gene (sourced from NCBI RefSeq), and any mismatches are automatically resolved into position-level variant calls, codon-level translations, and a preliminary pathogenicity classification.

The application is built for exploratory and educational use — demonstrating an end-to-end variant-calling workflow rather than replacing validated clinical pipelines (e.g., GATK, ANNOVAR).


✨ Core Features

Reference Data


Reference gene selection (e.g., BRCA1) with live metadata: NCBI Gene ID, RefSeq accession, chromosomal location, sequence length, and functional annotation
Direct links to NCBI Gene, FASTA record, ClinVar, and OMIM entries


Sequence Input


Manual paste or FASTA/TXT upload for the query sequence
Built-in mutated sample sequence for quick testing


Variant Analysis


Pairwise sequence alignment (reference vs. query) with match/mismatch visualization
Automated SNP detection — position, codon number, reference/alternate base, mutation type (e.g., SNP)
Codon-level translation showing resulting amino acid change (missense, silent, nonsense)
Variant classification: Pathogenic / Uncertain / Benign
Sequence identity scoring and summary metrics


Downstream Tools


GC content and basic sequence composition statistics
Mutation density heatmap across the reference gene
Primer design around detected variant regions


Patient Records


Patient profile management (demographics, clinical history, notes) linked to each analysis
Persistent local storage (SQLite)
Analysis history log per patient
CSV export of variant tables; JSON export of full patient + analysis records



🔄 Workflow


Select a reference gene from the sidebar; metadata and external database links load automatically.
Create or load a patient profile.
Provide a query sequence (paste, upload, or load the sample).
Run the analysis — the app aligns query against reference and calls variants.
Review results across dedicated views: Alignment, Mutations Table, Sequence Stats & GC, Codon Translation, Mutation Heatmap, and Primer Design.
Export results (CSV) or persist them to the patient's record for later retrieval.



🛠️ Tech Stack

ComponentTechnologyApplication frameworkStreamlitLanguagePython 3.10+Sequence handlingBiopythonDatabaseSQLiteReference data sourceNCBI GenBank / RefSeq

(Update to reflect the exact dependencies in requirements.txt.)


📂 Project Structure

PyMutScan/
├── app.py                     # Main Streamlit application
├── requirements.txt           # Python dependencies
├── pymutscane_patients.db     # Local SQLite database (auto-generated)
├── data/                      # Reference sequences / sample data
├── utils/                     # Alignment, variant calling, primer design modules
└── README.md


⚙️ Installation

bashgit clone https://github.com/navneet-bio-ai/pymutscane.git
cd pymutscane

python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate

pip install -r requirements.txt
streamlit run app.py

The app runs locally at http://localhost:8501.


🧪 Example Run


Select BRCA1 as the reference gene.
Load the sample mutated sequence (or paste a query sequence).
Run the analysis.
Output: 3 variants detected, 95.5% sequence identity — classified as 1 Pathogenic, 1 Uncertain, 1 Benign.
Inspect the mutations table and codon translation, then save the result to the patient's record.



🗄️ Data & Storage

All patient profiles and analysis results are stored locally in a SQLite database (pymutscane_patients.db) within the project directory. No data is transmitted externally.


⚠️ Disclaimer: PyMutScan is intended for educational and research demonstration purposes only. It has not been validated against clinical genomics pipelines and must not be used for diagnostic or treatment decisions.




⚠️ Limitations


Variant calling is based on direct pairwise alignment, not a probabilistic aligner — structural variants and indels beyond simple substitutions are not supported.
Pathogenicity classification is illustrative and not derived from a curated clinical variant database.
Currently scoped to a single reference gene per session.



🚀 Roadmap


 Extend reference gene panel beyond BRCA1
 Support indel detection and multi-sequence alignment
 Live ClinVar API integration for real-time variant annotation
 Expanded pharmacogenomic / drug-resistance prediction models
 Optional cloud-hosted database backend



📄 License

Licensed under the MIT License — see LICENSE for details.
