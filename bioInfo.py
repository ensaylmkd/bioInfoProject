from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
from Bio.Align.Applications import MafftCommandline

Entrez.email = "user@gmail.com"
# ////////// A //////////

# setup
RECHERCHE = "nucleotide",'(SARS-CoV-2[organism] AND refseq []) OR SARSr-CoV RaTG13 genome OR MP789 MT121216'
NEW_FILE = "seq_covid"

def researchToGbSeq(search,newFileName):
    '''Ecris le fichier gb correspondant à une recherche dans nucleotide '''
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
#researchToGbSeq(RECHERCHE,NEW_FILE)

# ////////// B //////////

#setup
FILE='seq_covid.gb'
def gbToInfo(file):
    '''Creer un txt contenant les infos d'une séquence depuis un fichier GenBank'''
    # settings
    seqs = list(SeqIO.parse(file,'gb'))

    # Lecture, creation et Ecriture des infos
    with open("info_seq_covid.txt",'w') as fd:
        for seq in seqs:
            cptGene=0
            infoGenes=[]
            for feature in seq.features:
                if feature.type == 'gene': # compte le nb de gene
                    cptGene+=1
            if feature.type =='CDS':
                infoGenes.append([feature.qualifiers['gene'][0],      # infoGenes contient des listes avec :
                                feature.location,                     # nom du gene, sa position et l'id de la proteine codée
                                feature.qualifiers['protein_id'][0]]) # de chaque CDS 

            fd.write("\n///////////////////////////////////////////////////////\n\n")
            fd.write(f"Organism name: {seq.annotations['organism']}\n")                 # nom de l'organism
            fd.write(f"Taxonomy: {seq.annotations['taxonomy']}\n")                      # taxonomie
            fd.write(f"Accession number: {seq.annotations['accessions'][0]}\n")         # num d'accession
            fd.write(f"Creation date: {seq.annotations['date']}\n")                     # date de creation
            fd.write(f"Number of genes: {cptGene}\n")                                   # nb gene 
            fd.write(f"% of GC : {round(((seq.seq.count('C')+seq.seq.count('G'))/len(seq.seq))*100,2)} % \n") # % de GC

            fd.write(f"Genes :\n")                                                      # info sur les genes
            for infoGene in infoGenes:
                fd.write(f" - name: {str(infoGene[0]): <12}")
                fd.write(f"location: {str(infoGene[1]): <45}")
                fd.write(f"protein_ID: {infoGene[2]}")
                fd.write("\n")

gbToInfo(FILE)
# //////////////////// C ///////////////////////

def multiFastaS(gene: str="S"):
    seqs = list(SeqIO.parse('seq_covid.gb','gb'))
    with open("spike.fasta",'w') as fd:                                           # On ecrit dans le fichier spike.fasta :
        for seq in seqs:                                                          #     (pour chaque sequence dans seq_covid.gb
            for feature in seq.features:                                          #     pour chaque gene == au gene demandé (si le l'element est de type CDS et a une traduction proteique => c'est un gene à prelever ) )
                if feature.type == 'CDS' and 'translation' in feature.qualifiers: # les fasta proteique de tout les genes demandé .
                    if feature.qualifiers['gene'][0] == gene:
                        fd.write(f">{feature.qualifiers['protein_id'][0]} {feature.qualifiers['product'][0]} [{seq.annotations['organism']}]\n")
                        fd.write("\n".join([feature.qualifiers['translation'][0][i:i+70] for i in range(0, len(feature.qualifiers['translation'][0]),70)])+"\n")


# //////////////////// D ///////////////////////
def alignement():                                   # code d'allignement fournit
    commande = MafftCommandline(input="spike.fasta")
    myStdout, myStderr = commande()
    with open("aln-spike.fasta", 'w') as w:
        w.write(myStdout)



# /////////////////// E /////////////////////////

def comparaison(gene: str="S"):
    resultat: list[list[int, str, str, str]] = []
    seqs = list(SeqIO.parse("aln-spike.fasta", "fasta"))
    difference = [0, 0] # CHAUVE-SOURIS et PANGOLIN
    for i, seq in enumerate(zip(seqs[0].seq, seqs[1].seq, seqs[2].seq)):
        if not seq[0] == seq[1] == seq[2]:
            resultat.append([i, seq[0], seq[1], seq[2]])
    # /////////////////// Debut du F /////////////////////////
    if seq[1] != seq[0]: # Si chauve souris différent de homme
        difference[0] += 1
    if seq[2] != seq[0]: # Si pangolin différent de homme
        difference[1] += 1

    with open(f"resultatComparaison_gene{gene}.txt", "w") as f:
        f.write(f"{'POSITION':^15}{'HOMME':^15}{'CHAUVE-SOURIS':^15}{'PANGOLIN':^15}\n")
    for i in resultat:  
        f.write(f"{i[0]:^15}{i[1]:^15}{i[2]:^15}{i[3]:^15}\n")

 # /////////////////// Fin du F /////////////////////////

    taille = len(seqs[0].seq)
    print(f"Il y a : {(taille-difference[0])/taille*100:.2f}% de lettres identiques en l'Homme et la Chauve-souris pour le gene {gene}.")
    print(f"Il y a : {(taille-difference[1])/taille*100:.2f}% de lettres identiques en l'Homme et le Pangolin pour le gene {gene}.")


# /////////////////// G /////////////////////////

def C_F(gene: str|None=None):
    if not gene:
        gene = input("Quel gene choisissez vous : (S/M/N)")
    if gene not in ["S", "M", "N"]:
        gene = "S"
    multiFastaS(gene)
    alignement()
    comparaison(gene)

# /////////////////// H /////////////////////////

C_F("S")
C_F("M")
C_F("N")
"""
Il y a : 90.04% de lettres identiques en l'Homme et la Chauve-souris pour le gene S.
Il y a : 97.41% de lettres identiques en l'Homme et le Pangolin pour le gene S.
Il y a : 98.20% de lettres identiques en l'Homme et la Chauve-souris pour le gene M.
Il y a : 99.10% de lettres identiques en l'Homme et le Pangolin pour le gene M.
Il y a : 97.85% de lettres identiques en l'Homme et la Chauve-souris pour le gene N.
Il y a : 99.05% de lettres identiques en l'Homme et le Pangolin pour le gene N.

Voici ce qu'on obtient en comparant les genes S, M et N, elles se ressemblent au moins a 90%,
ce qui est trop élevé pour être le fruit du hasard.
Et les genes du pangolin sont plus proche de ceux de l'Homme que la Chauve-souris.
"""




        
