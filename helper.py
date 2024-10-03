from Bio import Medline, Entrez 



# fetch data from PubMed API
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
