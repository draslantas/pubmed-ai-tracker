import json
import os
import re
from Bio import Entrez
from datetime import datetime, timedelta

# --- SETTINGS ---
Entrez.email = "senin_emailin@example.com"
SEARCH_TERM = "artificial intelligence"
DATA_FILE = "papers.json"

# --- CATEGORY DICTIONARIES (EXPANDED) ---
KEYWORDS = {
    "Branch": {
        "Neurosurgery": [
            "neurosurgery", "brain", "spine", "spinal", "glioma", "aneurysm", "craniotomy", 
            "resection", "meningioma", "deep brain stimulation", "hydrocephalus", "neurosurgical"
        ],
        "Cardiology": [
            "cardiology", "heart", "cardiac", "ecg", "atrial", "coronary", "myocardial", 
            "infarction", "arrhythmia", "stent", "valve", "heart failure", "aortic"
        ],
        "Orthopedics": [
            "orthopedics", "bone", "fracture", "joint", "knee", "hip", "arthroplasty", 
            "trauma", "osteosynthesis", "acl", "meniscus", "spine fusion", "musculoskeletal"
        ],
        "Radiology": [
            "radiology", "mri", "ct scan", "x-ray", "ultrasound", "imaging", "radiomics", 
            "pet scan", "magnetic resonance", "computed tomography", "sonography"
        ],
        "Dermatology": [
            "dermatology", "skin", "lesion", "melanoma", "dermoscopy", "cutaneous", "psoriasis"
        ],
        "Ophthalmology": [
            "ophthalmology", "eye", "retina", "glaucoma", "cataract", "cornea", "macular", "visual"
        ],
        "Pathology": [
            "pathology", "histology", "biopsy", "slide", "microscopic", "cytology", "immunohistochemistry"
        ],
        "Oncology": [
            "oncology", "cancer", "tumor", "malignancy", "metastasis", "carcinoma", "sarcoma", 
            "chemotherapy", "radiotherapy"
        ],
        "Neurology": [
            "neurology", "stroke", "alzheimer", "parkinson", "epilepsy", "multiple sclerosis", 
            "dementia", "neurodegenerative", "seizure"
        ],
        "General Surgery": [
            "surgery", "surgical", "laparoscopy", "robotic surgery", "abdominal", "resection", 
            "operation", "postoperative"
        ],
        "Gastroenterology": [
            "gastroenterology", "liver", "stomach", "colon", "endoscopy", "hepatitis", "bowel"
        ],
        "Urology": [
            "urology", "kidney", "bladder", "prostate", "renal", "urinary"
        ]
    },
    "AI_Type": {
        "Imaging": ["segmentation", "detection", "classification", "image analysis", "reconstruction", "computer vision"],
        "Decision Support": ["decision support", "prediction", "prognosis", "risk assessment", "diagnosis", "predictive model"],
        "NLP/LLM App": ["chatbot", "summary", "generation", "question answering", "report generation", "natural language processing"],
        "Drug Discovery": ["drug discovery", "molecule", "docking", "protein folding", "virtual screening"],
        "Genomics": ["genomics", "gene", "sequence", "mutation", "bioinformatics"]
    },
    "Method": {
        "Deep Learning": ["deep learning", "cnn", "rnn", "convolutional", "neural network", "unet", "resnet"],
        "Machine Learning": ["machine learning", "svm", "random forest", "logistic regression", "gradient boosting", "xgboost"],
        "LLM/Transformer": ["llm", "transformer", "gpt", "bert", "llama", "large language model", "rag", "retrieval-augmented", "generative ai"],
        "Reinforcement Learning": ["reinforcement learning", "rl", "agent"]
    }
}

def get_stored_papers():
    """Loads stored papers."""
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r") as f:
            return json.load(f)
    return []

def save_papers(papers):
    """Saves papers to file."""
    with open(DATA_FILE, "w") as f:
        json.dump(papers, f, indent=4)

def analyze_paper(text):
    """Analyzes text and extracts tags."""
    text = text.lower()
    tags = {"Branch": [], "AI_Type": [], "Method": []}
    
    for category, subcategories in KEYWORDS.items():
        for tag, words in subcategories.items():
            if any(word in text for word in words):
                tags[category].append(tag)
    
    # Mark empty categories as "Other"
    for cat in tags:
        if not tags[cat]:
            tags[cat].append("General/Other")
            
    return tags

def fetch_new_papers(days=1):
    """Fetches papers from the last N days."""
    print(f"🔍 Scanning PubMed for the last {days} days...")
    
    handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, reldate=days, retmax=1000, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    if not id_list:
        print("No new papers found.")
        return

    # Fetch details
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    stored_papers = get_stored_papers()
    stored_ids = {p['id'] for p in stored_papers}
    
    new_count = 0
    
    for paper in records['PubmedArticle']:
        try:
            pmid = str(paper['MedlineCitation']['PMID'])
            
            if pmid in stored_ids:
                continue
                
            article = paper['MedlineCitation']['Article']
            title = article['ArticleTitle']
            
            # Journal
            journal = article.get('Journal', {}).get('Title', 'Unknown Journal')
            
            # Date Parsing (Robust)
            pub_date_data = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            if 'Year' in pub_date_data:
                year = pub_date_data.get('Year', '')
                month = pub_date_data.get('Month', '')
                day = pub_date_data.get('Day', '')
                pub_date = f"{year} {month} {day}".strip()
            elif 'MedlineDate' in pub_date_data:
                pub_date = pub_date_data['MedlineDate']
            else:
                pub_date = "Unknown Date"
            
            # Abstract
            abstract_list = article.get('Abstract', {}).get('AbstractText', [])
            abstract_text = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            
            # Authors
            authors = article.get('AuthorList', [])
            author_names = [f"{a.get('LastName', '')} {a.get('Initials', '')}" for a in authors[:3]]
            author_str = ", ".join(author_names)
            if len(authors) > 3:
                author_str += " et al."
            
            # Analyze
            tags = analyze_paper(title + " " + abstract_text)
            
            new_paper = {
                "id": pmid,
                "title": title,
                "authors": author_str,
                "journal": journal,
                "pub_date": pub_date,
                "abstract": abstract_text,
                "tags": tags,
                "date_added": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "status": "unread"
            }
            
            stored_papers.append(new_paper)
            new_count += 1
            
        except Exception as e:
            print(f"Error ({pmid}): {e}")
            
    save_papers(stored_papers)
    print(f"✅ {new_count} new papers saved.")

if __name__ == "__main__":
    fetch_new_papers(days=10)
