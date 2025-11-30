# Personal PubMed AI Assistant 🧬

This program tracks the PubMed database for new articles related to "Artificial Intelligence", categorizes them, and presents them in a user-friendly dashboard.

## 📦 Installation
Open a terminal in this folder and install the required libraries:
```bash
pip3 install streamlit schedule biopython
```

## 🚀 How to Run

### 1. Launch the Dashboard
To open the interface with Blue/Green cards:
```bash
streamlit run app.py
```

### 2. Start the Daily Scheduler (08:00 AM)
To run the background job that fetches new papers every morning:
```bash
python3 scheduler.py
```

## ✨ Features
- **Auto-Tracking:** Fetches papers published in the last 24 hours daily.
- **Smart Tagging:** Categorizes papers by Medical Branch (Neurosurgery, Cardiology, etc.) and AI Method (Deep Learning, LLM, etc.).
- **Reading List:** Unread papers appear **Blue**, read ones turn **Green**.
- **Keep Notebook:** Save inspiring papers to a special **Gold** list.
