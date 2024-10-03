from Bio import Entrez, Medline  # BioPython modules
import pandas as pd
from tqdm import tqdm

# Set your email (this is required by NCBI for identification)
Entrez.email = "mostafadesoki86@gmail.com"

# Function to fetch data from PubMed API
def fetch_pubmed_data(query, max_records=1000):
    # Search PubMed for the query
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_records)
    record = Entrez.read(handle)
    handle.close()

    # Get list of PubMed IDs (PMIDs)
    pmids = record["IdList"]

    # Fetch the article details in MEDLINE format using the PMIDs
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
    articles = Medline.parse(handle)
    return list(articles)

# Define the search query (AI and Dentistry example)
search_query = '("artificial intelligence"[Title/Abstract] OR "AI"[Title/Abstract]) AND ("dentistry"[Title/Abstract] OR "dental"[Title/Abstract])'

# Fetch the articles based on the query
articles = fetch_pubmed_data(search_query, max_records=100)

# Initialize list to store the data
alldata = []

# Process the retrieved articles
for article in tqdm(articles):
    try:
        pmid = article.get("PMID", "")
        author = article.get("FAU", "")
        title = article.get("TI", "")
        journal = article.get("JT", "")
        year = article.get("DP", "")
        abstract = article.get("AB", "")
        identifier = article.get("AID", [""])[-1][:-6] if "AID" in article else ""

        dic = {
            "PMID": pmid,
            "Author": author,
            "Title": title,
            "Journal": journal,
            "Year": year,
            "Abstract": abstract,
            "DOI": f"https://doi.org/{identifier}" if identifier else ""
        }
        alldata.append(dic)

    except Exception as e:
        print(f"Error processing article: {e}")

# Convert the data into a pandas DataFrame
df = pd.DataFrame(alldata)

# Save the data to CSV file
df.to_csv('pubmed_articles.csv', index=False)

print("Finished processing and saved to CSV.")
