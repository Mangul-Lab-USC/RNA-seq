# Software Development Practices in 430 RNA-Seq Tools
We conducted a comprehensive assessment of 430 RNA-seq tools developed between 2008 and 2024, categorizing them based on the type of analysis they perform. Our evaluation encompassed their software development and distribution methodologies, as well as the attributes contributing to their widespread adoption and dependability within the biomedical community. Our findings establish the first documented positive correlation between rigorous software development practices - quantified by factors such as package manager availability, Docker containerization, multithreading support, documentation quality, and example dataset inclusion - and the scholarly impact of published RNA-seq tools, as measured by citations (Mann-Whitney U test, p-value = 4.9 × 10⁻-26). By identifying key characteristics of widely adopted software, our findings provide guidance for developing robust and user-friendly RNA-seq tools, thereby reinforcing the call for new community-wide standards. 

---

## 📑 Table of Contents
- [How to Cite This Study](#-how-to-cite-this-study)
- [Datasets](#-datasets)
- [Reproducing Results](#-reproducing-results)
- [License](#license)
- [Contact](#-contact)

---

## 📚 How to Cite This Study

Sharma, S., et al.(2025) Robust software development practices improve citations of RNA-seq tools. **Biopolymers and Cell**. DOI: [10.7124/bc.000AFE](http://dx.doi.org/10.7124/bc.000AFE)

---

## 📊 Datasets

We compiled publications describing novel RNA-seq tools from [Google Scholar](https://scholar.google.com), [PubMed](https://pubmed.ncbi.nlm.nih.gov), and [Oxford Academic](https://academic.oup.com).  
Our approach for extracting and verifying software links is described in the **Methods** section of the manuscript. Timeout links were manually verified.

The dataset is provided as a CSV file and contains the following fields:
- Name of the tool  
- Year of publication  
- Software interface utilized  
- Package manager availability  
- Docker containerization support  
- Multithreading capacity  
- User guide availability  
- Sample dataset availability  
- Archival stability  
- Number of releases or updates  
- Benchmarking practices  
- License type  
- Citation count  

---

## 🧪 Reproducing Results

To reproduce our figures and results, we provide a Google Colab Notebook:

- [Google Colab Notebook](https://github.com/Mangul-Lab-USC/RNA-seq/blob/master/notebook/RNA_Seq_Project..ipynb)

All figures and analyses can be reproduced using the accompanying code and data.

---

## License
This repository is under MIT license. 

---

## 📬 Contact

Please contact us with comments, suggestions, or questions:

- serghei.mangul@gmail.com    
- mealser@gmail.com  
- susharma@usc.edu
