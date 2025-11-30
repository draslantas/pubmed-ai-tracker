import schedule
import time
import main
from datetime import datetime

def job():
    print(f"⏰ Daily scan started: {datetime.now()}")
    # Scan last 24h
    main.fetch_new_papers(days=1)
    print("✅ Scan complete.")

# Run every day at 08:00
schedule.every().day.at("08:00").do(job)

print("⏳ Scheduler started. Will run daily at 08:00.")
print("Press CTRL+C to stop.")

# Run once on startup
job()

while True:
    schedule.run_pending()
    time.sleep(60)
