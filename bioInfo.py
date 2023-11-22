from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "user@gmail.com"
# ////////// A //////////

# setup
RECHERCHE = "nucleotide",'(SARS-CoV-2[organism] AND refseq []) OR SARSr-CoV RaTG13 genome OR MP789 MT121216'
NEW_FILE = "seq_covid"

def createSeq(search,newFileName):
    # on recupere les id des séquences que l'on veut
    request = Entrez.esearch(search[0],search[1])
    research = Entrez.read(request)
    request.close()
    # on récupere le fetch des sequences en gb
    request = Entrez.efetch(db='nucleotide',id=research['IdList'],rettype='gb',retmode='text')
    # on les converties en fichier gb 
    seqs = SeqIO.parse(request,'gb')
    SeqIO.write(seqs,f"{newFileName}.gb","gb")
    request.close()
#createSeq(RECHERCHE,NEW_FILE)

# ////////// B //////////

seqs = list(SeqIO.parse('seq_covid.gb','gb'))
with open("info_seq_covid.txt",'w') as fd:
    for seq in seqs:
        fd.write("\n///////////////////////////////////////////////////////\n\n")
        fd.write(f"organism name : {seq.annotations['organism']}\n")
        fd.write(f"taxonomy : {seq.annotations['taxonomy']}\n")
        fd.write(f"accession number : {seq.annotations['accessions'][0]}\n")
        fd.write(f"creation date : {seq.annotations['date']}\n")



        
