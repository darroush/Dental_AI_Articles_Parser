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


# exclude Q3 and Q4 journals
journals_blacklist = ["AMIA ... Annual Symposium proceedings. AMIA Symposium",
"Bioinformation","Compendium of continuing education in dentistry (Jamesburg, N.J. : 1995)",
"Cureus",
"Journal of the Indian Society of Pedodontics and Preventive Dentistry",
"Journal of veterinary dentistry",
"Journal of visualized experiments : JoVE",
"JPMA. The Journal of the Pakistan Medical Association",
"L' Orthodontie francaise",
"Medecine sciences : M/S",
"Special care in dentistry : official publication of the American Association of Hospital Dentists, the Academy of Dentistry for the Handicapped, and the American Society for Geriatric Dentistry",
"Studies in health technology and informatics",
"Technology and health care : official journal of the European Society for Engineering and Medicine",
"The journal of contemporary dental practice",
"The Journal of forensic odonto-stomatology",
"Zhonghua kou qiang yi xue za zhi = Zhonghua kouqiang yixue zazhi = Chinese journal of stomatology"]

