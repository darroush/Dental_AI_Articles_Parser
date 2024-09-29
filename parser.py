from Bio import Medline # type: ignore
import pandas as pd
from tqdm import tqdm

alldata = []

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

with open("pubmed-Artificial-set-september.txt", encoding= "utf-8") as file:
    articles = Medline.parse(file)
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

df = pd.DataFrame(alldata)
df.to_csv('parser_output_september.csv', index=False)
print('finished')




    