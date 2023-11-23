from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "user@gmail.com"
# ////////// A //////////

# setup
RECHERCHE = "nucleotide",'(SARS-CoV-2[organism] AND refseq []) OR SARSr-CoV RaTG13 genome OR MP789 MT121216'
NEW_FILE = "seq_covid"

def researchToGbSeq(search,newFileName):
    '''Ecris le fichier gb correspondant à une recherche dans nucleotide  '''
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
#createGbSeqFromResearch(RECHERCHE,NEW_FILE)

# ////////// B //////////

#setup
FILE='seq_covid.gb'
def gbToInfo():
    '''Creer un txt contenant les infos d'une séquence depuis un fichier GenBank'''
    # settings
    seqs = list(SeqIO.parse(FILE,'gb'))

    # Lecture, creation et Ecriture des infos
    with open("info_seq_covid.txt",'w') as fd:
        for seq in seqs:
            cptGene=0
            infoGenes=[]
            for feature in seq.features:
                if feature.type == 'gene':                                 # compte le nb de gene
                    cptGene+=1
                if feature.type =='CDS':
                    infoGenes.append([feature.qualifiers['gene'][0],       # infoGenes contient des listes avec :
                                    feature.location,                      # nom du gene, sa position et l'id de la proteine codée
                                    feature.qualifiers['protein_id'][0]])  # de chaque CDS  
                                
            fd.write("\n///////////////////////////////////////////////////////\n\n")
            fd.write(f"Organism name: {seq.annotations['organism']}\n")                 # nom de l'organism
            fd.write(f"Taxonomy: {seq.annotations['taxonomy']}\n")                      # taxonomie
            fd.write(f"Accession number: {seq.annotations['accessions'][0]}\n")         # num d'accession
            fd.write(f"Creation date: {seq.annotations['date']}\n")                     # date de creation
            fd.write(f"Number of genes: {cptGene}\n")                                   # nb gene 
            fd.write(f"% of GC : {round(((seq.seq.count('C')+seq.seq.count('G'))/len(seq.seq))*100,2)} % \n") # % de GC

            fd.write(f"Genes :\n")                                                      # info sur les genes
            for infoGene in infoGenes:
                fd.write(" "*3+f"- name: {infoGene[0]}")       # Remi vient recentrer ici jsp faire (regarde le .txt)
                fd.write(" "*6 + f"location: {infoGene[1]}")
                fd.write(" "*6 + f"protein_ID: {infoGene[2]}")
                fd.write("\n")
        
        






        
