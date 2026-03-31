# Part 3: Evaluating Drug-Likeness with Lipinski Ro5 and Veber's Rules with RDKit and Python
### <sub>Implementing Lipinski's Rule of 5 and Veber's Rules to Filter Chemical Datasets for Drug Discovery</sub>

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![RDKit](https://img.shields.io/badge/RDKit-2022.03+-green.svg)
![Pandas](https://img.shields.io/badge/Pandas-Latest-orange.svg)

## ## 🧪 Project Overview
In this part of the series, we move beyond SMILES standardization to functional drug discovery analysis. We implement:
* **Lipinski’s Rule of 5 (Ro5):** Assessing oral bioavailability.
* **Veber’s Rules:** Evaluating molecular flexibility and polar surface area.
* **The Rule of 3 (Ro3):** Identifying lead-like fragments for optimization.
* **Modern Exceptions:** A brief look at "Beyond Rule of 5" (bRo5) molecules like PROTACs.

## 📂 Repository Structure
* `antibiotic_data.csv`: The raw dataset of ~2,300 compounds retrieved from PubChem.
* `drug_likeness_analysis.py`: The main Python script using RDKit for descriptor calculation and vectorized data filtering.
* `requirements.txt`: List of necessary Python libraries.

## 🚀 How to Run
1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YourUsername/YourRepoName.git](https://github.com/YourUsername/YourRepoName.git)
   cd YourRepoName
2. Install dependencies:
    pip install -r requirements.txt

3. Run the script:
    python drug_likeness_analysis.py

📊 Results
The script performs automated data curation, removing duplicates through SMILES standardization and identifying compounds that satisfy both Lipinski and Veber criteria.

🔗 Related Articles

Part 1: [A Chemist’s Guide to RDKit & Python](https://medium.com/@kparama82/from-bench-to-bytes-a-chemists-guide-to-rdkit-python-fe2d870da447?source=friends_link&sk=c1b08c4635b4a1bde4104f1e72355168)

Part 2:[ SMILES Standardization and Tautomers](https://medium.com/@kparama82/mastering-molecular-identity-handling-tautomers-in-heterocyclic-datasets-with-rdkit-cf0e7ba4f9dd?source=friends_link&sk=219be3cd7d0d81da947353b880792c02)

Part 3: Analysis of Drug-Likeness

