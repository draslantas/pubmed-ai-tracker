import streamlit as st
import json
import os
import main
from datetime import datetime

DATA_FILE = "papers.json"

st.set_page_config(page_title="PubMed AI Tracker", page_icon="🧬", layout="wide")

def load_papers():
    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, "r") as f:
            return json.load(f)
    return []

def save_papers(papers):
    with open(DATA_FILE, "w") as f:
        json.dump(papers, f, indent=4)

def update_status(paper_id, new_status):
    papers = load_papers()
    for p in papers:
        if p['id'] == paper_id:
            p['status'] = new_status
            break
    save_papers(papers)
    st.rerun()

# --- UI ---
st.title("🧬 Personal PubMed AI Assistant")
st.markdown("---")

# Sidebar
with st.sidebar:
    st.header("Control Panel")
    if st.button("🔄 Full Scan (2025)"):
        with st.spinner("Scanning PubMed (This may take a while)..."):
            main.fetch_new_papers(days=365)
        st.success("Scan complete!")
        st.rerun()
    
    st.info("Auto-scan runs daily at 08:00.")

# Load and Sort Papers (Newest First by Publication Date)
papers = load_papers()
# Sort by 'sort_date' (YYYY-MM-DD) which we added in main.py
papers.sort(key=lambda x: x.get('sort_date', '0000-00-00'), reverse=True)

# Daily Count
today_str = datetime.now().strftime("%Y-%m-%d")
today_count = sum(1 for p in papers if p.get('date_added', '').startswith(today_str))
st.metric(label="Papers Added Today", value=today_count)

# Tabs
tab1, tab2, tab3 = st.tabs(["🟦 Unread", "🟩 Read / Archive", "💡 Keep Notebook"])

unread_papers = [p for p in papers if p.get('status') == 'unread']
read_papers = [p for p in papers if p.get('status') == 'read']
keep_papers = [p for p in papers if p.get('status') == 'keep']

def render_paper_card(paper, current_tab):
    # Card Border Color
    if current_tab == "unread":
        border_color = "blue"
    elif current_tab == "read":
        border_color = "green"
    else:
        border_color = "gold" # Keep
    
    with st.container():
        st.markdown(f"""
        <div style="border-left: 5px solid {border_color}; padding-left: 10px; margin-bottom: 20px;">
            <h3><a href="https://pubmed.ncbi.nlm.nih.gov/{paper['id']}/" target="_blank" style="text-decoration: none; color: inherit;">{paper['title']}</a></h3>
            <p><b>{paper.get('journal', 'Unknown Journal')}</b> | <i>{paper.get('pub_date', 'Unknown Date')}</i></p>
            <p><i>{paper['authors']}</i></p>
            <p><a href="https://pubmed.ncbi.nlm.nih.gov/{paper['id']}/" target="_blank">🔗 View on PubMed</a></p>
        </div>
        """, unsafe_allow_html=True)
        
        # Tags
        tags = paper.get('tags', {})
        cols = st.columns(3)
        with cols[0]:
            st.caption(f"🏷️ **Branch:** {', '.join(tags.get('Branch', []))}")
        with cols[1]:
            st.caption(f"🤖 **AI Type:** {', '.join(tags.get('AI_Type', []))}")
        with cols[2]:
            st.caption(f"⚙️ **Method:** {', '.join(tags.get('Method', []))}")
            
        with st.expander("📄 Read Abstract"):
            st.write(paper['abstract'])
            
            # Action Buttons
            col1, col2 = st.columns(2)
            
            with col1:
                if current_tab == "unread":
                    if st.button("✅ Mark as Read", key=f"read_{paper['id']}"):
                        update_status(paper['id'], "read")
                elif current_tab == "read":
                    if st.button("↩️ Mark as Unread", key=f"unread_{paper['id']}"):
                        update_status(paper['id'], "unread")
            
            with col2:
                if current_tab != "keep":
                    if st.button("⭐ Keep / Inspiring", key=f"keep_{paper['id']}"):
                        update_status(paper['id'], "keep")
                else:
                    if st.button("❌ Remove from Keep", key=f"unkeep_{paper['id']}"):
                        update_status(paper['id'], "read") # Move back to read

with tab1:
    st.subheader(f"Unread ({len(unread_papers)})")
    for p in unread_papers:
        render_paper_card(p, "unread")
        st.divider()

with tab2:
    st.subheader(f"Archive ({len(read_papers)})")
    for p in read_papers:
        render_paper_card(p, "read")
        st.divider()

with tab3:
    st.subheader(f"Inspiring Papers ({len(keep_papers)})")
    for p in keep_papers:
        render_paper_card(p, "keep")
        st.divider()
