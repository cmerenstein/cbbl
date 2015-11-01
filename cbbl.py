# Carter Merenstein

''' Downloads all the sequences of the cbbl gene on the NCBI nucleotide database. Sequences that are only given in
whole genomes are excluded because of the time it takes to download them. All sequences are then given in FASTA
with the organism and accession number as the identifier. This will be used to design cbbl primers.
I download everything as text instead of using their xml parser because their parser is junk'''


from Bio import Entrez as E
from Bio.Seq import Seq
from urllib.error import HTTPError
import time
import re

E.email = "cmerenstein@gmail.com" ## NCBI requires an email

not_number = re.compile("[^0-9]")
not_base = re.compile("[^atcg]")

def getBounds(genbank):
    ''' returns the position of the cbbl gene and whether or not it's reverse complimented.
    This is necessary for entries that have more than ust the cbbl gene included.'''
    
    compliment = 0 # False
    bounds = [0,0]
    
    for line in genbank.split("FEATURES")[1].split('\n'):
        
        if line.strip()[0:4]  == "gene":

            if  "complement" in line:
                compliment = 1
            else:
                compliment = 0

            bound_line = line.split('..')
            try:
                assert (len(bound_line) ==2)
            except AssertionError:
                print(genbank)
            bounds[0] = int(re.sub(not_number, '', bound_line[0])) - 1 # gb sequences count from 1, not 0
            bounds[1] = int(re.sub(not_number, '', bound_line[1]))
            

        if line.strip()[0:5] == "/gene":
            gene = line.split('=')[-1]
            gene = gene.strip('\"').lower()
            if gene == "cbbl":
                bounds.append(compliment) # whether or not to reverse compliment is saved in bounds as 1 or 0
                return(bounds)
    
def getSequence(genbank):
    ''' returns the whole sequence from an entry, as a string'''

    sequence = genbank.split("ORIGIN")[-1]
    return re.sub(not_base, '', sequence)

def getOrganism(genbank):
    ''' returns the organism '''

    organism = genbank.split("REFERENCE")[0].split("ORGANISM")[1]
    organism = organism.split('\n')[0].strip()
    return organism

def getAccession (genbank):
    ''' returns the Accension '''

    accession = genbank.split("VERSION")[0].split("ACCESSION")[1]
    accession = accession.strip()
    return accession

def fetch(search):
    ''' Biopython guidelines suggests breaking up downloading files in batches in case they are large.
    Documentation also suggests excepting HTTPError. Search is a result from Entrez.esearch,
    with usehistory on.
    All entries get put into a fasta file, "cbbl.fasta" '''
    output = open("cbbL.fasta", 'w')
   
    batch_size = 1 # get one result at a time
    count = int(search["Count"])
    gb_list = search["IdList"]
    webenv = search["WebEnv"]
    query_key = search["QueryKey"]
    
    for start in range(0, count):
        attempt = 1
        while attempt<=5:
            try:
                handle = E.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
                record = handle.read()

                fa_header = "> "
                fa_header = fa_header + getAccession(record) + " -- " + getOrganism(record) + '\n'
                output.write(fa_header)
                
                bounds = getBounds(record)
                sequence = getSequence(record)
                cbbL_seq = Seq(sequence[bounds[0]:bounds[1]])
                if bounds[2] == 1:
                    cbbL_seq = cbbL_seq.reverse_complement()
                output.write(str(cbbL_seq) + '\n')
                break
            except HTTPError as err: # Sometimes there are connectivity issues
                print(attempt)
                attempt +=1
                time.sleep(5)

    output.close()
    

''' Main'''

search_handle = E.esearch(db="nucleotide", term="0:4000[Sequence Length] AND cbbl[gene]", usehistory="y")
search_results = E.read(search_handle)
search_handle.close()

fetch(search_results)


