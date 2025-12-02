from Bio import Entrez
import datetime

Entrez.email = "senin_emailin@example.com"
SEARCH_TERM = "artificial intelligence"

print(f"🔍 Debugging Daily Scan for: '{SEARCH_TERM}'")

# 1. Try reldate=1 (Last 1 day)
print("\n--- Test 1: reldate=1 ---")
try:
    handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, reldate=1, retmax=10, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()
    print(f"Found {record['Count']} papers.")
    if int(record['Count']) > 0:
        print("First 5 IDs:", record['IdList'][:5])
except Exception as e:
    print(f"Error: {e}")

# 2. Try explicit date range (Yesterday to Today)
print("\n--- Test 2: Explicit Date Range ---")
today = datetime.datetime.now().strftime("%Y/%m/%d")
yesterday = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime("%Y/%m/%d")
date_query = f'("{yesterday}"[Date - Publication] : "{today}"[Date - Publication])'
full_query = f"{SEARCH_TERM} AND {date_query}"

print(f"Query: {full_query}")
try:
    handle = Entrez.esearch(db="pubmed", term=full_query, retmax=10, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()
    print(f"Found {record['Count']} papers.")
    if int(record['Count']) > 0:
        print("First 5 IDs:", record['IdList'][:5])
except Exception as e:
    print(f"Error: {e}")
