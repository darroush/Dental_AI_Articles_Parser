from Bio import Medline, Entrez
import pandas as pd
from tqdm import tqdm 


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


# PubMed parser
def pubmed_parser(articles):
    alldata = []
    # Process the retrieved articles
    for article in tqdm(articles):
        if article["JT"] not in journals_blacklist:
            try:
                pmid = article["PMID"]
            except: 
                pmid = ""
            try:
                author = article["FAU"]
            except:
                author = ""
            try:
                title = article["TI"]
            except:
                title = ""
            try:
                journal = article["JT"]
            except:
                journal = ""
            try:
                year = article["DP"]
            except:
                year = ""
            try:
                identifier = article["AID"][-1][:-6]
            except:
                identifier = ""
            try:
                abstract = article["AB"]
            except:
                abstract = ""

            dic = {
                "PMID" : pmid,
                "Author" : author,
                "Title" : title,
                "Journal" : journal,
                "Year" : year,
                "Abstract" : abstract,
                "DOI" : "https://doi.org/" + identifier
            }
            alldata.append(dic)

    # save the data as excel file
    df = pd.DataFrame(alldata)
    df.to_excel(f'parser_output.xlsx', index=False)
    print('excel file extracted..')

