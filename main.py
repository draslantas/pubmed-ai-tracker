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

def parse_pub_date(pub_date_data):
    """Parses PubMed date into YYYY-MM-DD string for sorting."""
    try:
        year = pub_date_data.get('Year', '')
        month = pub_date_data.get('Month', '01') # Default Jan
        day = pub_date_data.get('Day', '01')     # Default 1st
        
        # Convert month name to number
        month_map = {
            "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06",
            "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
        }
        if month in month_map:
            month = month_map[month]
            
        # Handle MedlineDate (e.g. "2024 Oct-Dec")
        if not year and 'MedlineDate' in pub_date_data:
            medline_date = pub_date_data['MedlineDate']
            match = re.search(r'(\d{4})', medline_date)
            if match:
                year = match.group(1)
                
        if year:
            return f"{year}-{month}-{day}"
            
    except Exception:
        pass
        
    return "0000-00-00" # Fallback for unknown dates

def fetch_new_papers(days=1, is_full_scan=False):
    """
    Fetches papers from PubMed.
    
    Args:
        days: Number of days to look back (used when is_full_scan=False)
        is_full_scan: If True, uses explicit date range from Oct 8, 2025 to today
    """
    
    # Statistics tracking
    stats = {
        "pubmed_total": 0,
        "already_stored": 0,
        "filtered_by_date": 0,
        "date_parse_failed": 0,
        "new_added": 0
    }
    
    cutoff_date = "2025-10-08"
    
    if is_full_scan:
        # Use explicit date range for true "since Oct 8, 2025" behavior
        print(f"🔍 Full Scan: Fetching all papers since {cutoff_date}...")
        today = datetime.now().strftime("%Y/%m/%d")
        date_query = f'("{cutoff_date.replace("-", "/")}"[Date - Publication] : "{today}"[Date - Publication])'
        full_query = f"{SEARCH_TERM} AND {date_query}"
        
        handle = Entrez.esearch(db="pubmed", term=full_query, retmax=9999, sort="pub+date")
    else:
        # Daily scan: just last N days
        print(f"🔍 Daily Scan: Fetching papers from last {days} day(s)...")
        handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, reldate=days, retmax=9999, sort="pub+date")
    
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    stats["pubmed_total"] = len(id_list)
    
    if not id_list:
        print("❌ No papers found in PubMed for this query.")
        return

    # Fetch details in batches
    batch_size = 100
    stored_papers = get_stored_papers()
    stored_ids = {p['id'] for p in stored_papers}

    print(f"📥 PubMed returned {len(id_list)} papers. Processing...")

    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        ids = ",".join(batch_ids)
        
        try:
            handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            for paper in records['PubmedArticle']:
                try:
                    pmid = str(paper['MedlineCitation']['PMID'])
                    
                    # Check if already stored
                    if pmid in stored_ids:
                        stats["already_stored"] += 1
                        continue
                        
                    article = paper['MedlineCitation']['Article']
                    title = article['ArticleTitle']
                    journal = article.get('Journal', {}).get('Title', 'Unknown Journal')
                    
                    # Date Parsing
                    pub_date_data = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                    sort_date = parse_pub_date(pub_date_data)
                    
                    # Display Date
                    if 'Year' in pub_date_data:
                        pub_date_display = f"{pub_date_data.get('Year')} {pub_date_data.get('Month', '')} {pub_date_data.get('Day', '')}".strip()
                    elif 'MedlineDate' in pub_date_data:
                        pub_date_display = pub_date_data['MedlineDate']
                    else:
                        pub_date_display = "Unknown Date"

                    # Date validation
                    if sort_date == "0000-00-00":
                        stats["date_parse_failed"] += 1
                        # Still add it, but mark it
                        sort_date = "2025-12-01"  # Default to recent date if parse failed
                    
                    # Only apply cutoff for full scans (daily scans are already time-filtered by PubMed)
                    if is_full_scan and sort_date < cutoff_date:
                        stats["filtered_by_date"] += 1
                        continue

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
                        "pub_date": pub_date_display,
                        "sort_date": sort_date,
                        "abstract": abstract_text,
                        "tags": tags,
                        "date_added": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "status": "unread"
                    }
                    
                    stored_papers.append(new_paper)
                    stats["new_added"] += 1
                    
                except Exception as e:
                    print(f"⚠️  Error processing paper {pmid}: {e}")
            
            # Save after each batch (silent)
            save_papers(stored_papers)
            
        except Exception as e:
            print(f"❌ Error fetching batch {i//batch_size + 1}: {e}")
            continue

    # Print summary
    print(f"""
{'='*60}
📊 TARAMA ÖZETİ:
{'='*60}
PubMed'den gelen toplam:        {stats['pubmed_total']} makale
Zaten veritabanında:            {stats['already_stored']} makale
Tarih filtresiyle elenen:       {stats['filtered_by_date']} makale
Tarih parse edilemeyen:         {stats['date_parse_failed']} makale
✅ YENİ EKLENEN:                {stats['new_added']} makale
{'='*60}
    """)

if __name__ == "__main__":
    # Initial run: fetch last 365 days to catch everything in 2025
    fetch_new_papers(days=365)
