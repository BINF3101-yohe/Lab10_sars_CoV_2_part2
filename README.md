# Lab10 SARS-CoV-2: Learn how to identify gene function through orthology

Let's review. In the Lab 9, we wrote used a python script to retrieve ORFs from the SARS-CoV-2 genome. Some retrieved ORFs may correspond to the actual genes, some were there by chance, and some because the start codon also encodes for methionine, so the retrieved ORF was a (short) part of the more extended coding sequence. We used protein hydrophobicity to find potential transmembrane proteins. But what about other types of proteins? How can we tell which ORFs correspond to the actual genes?

This is where evolution comes. If two species are morphologically different, their genetic material must also differ substantially. Right? Not exactly. Evolution conserved many genes crucial for survival. The degree of conservation is surprising. For example, while we would expect to share many genes with the chimpanzee, our closest relative, we would expect almost no similarity in our genetic material with bananas. The numbers tell us otherwise. While sharing 96% of our genes with the chimpanzee, we also share 60% of our genes with bananas! These similar genes in humans and bananas encode proteins that carry similar functions.

Two genes with similar sequences may perform similar functions. Genes like these are called homologous genes. Or, more precisely, in our case, these would be orthologous genes (see the figure below). A pair of orthologous genes are genes from two related species originating from a common ancestor (e.g., humans and chimpanzees). Due to the evolution, the two nucleotide or protein sequences (random accumulation of mutations due to drift or nonrandom accumulation of mutations due to selection) might differ slightly but not overwhelmingly.

<img width="1039" alt="Screenshot 2025-03-24 at 10 41 24 PM" src="https://github.com/user-attachments/assets/d638a709-1009-4a3d-a3f3-0e7a937f02b0" />

Let's first set up our lab environment. Log into the cluster.

```bash
mkdir lab_10
cd lab_10
#we also need the sars_cov_2 genome for this lab
cp /projects/class/binf3101_001/sars_cov_2.fasta ~/lab_10
cd lab_10
```

The human-banana example demonstrates that these two organisms don't have to be that closely related! If we tried to determine the gene function in humans using bananas as a reference, we could infer the role of over 60% of the genes. However, using two closely related organisms would undoubtedly give us more robust results. We would be better off choosing the genome of a mouse or a chimpanzee as our reference.

Scientists have studied hundreds of viruses and determined the functions of most of their genes. Using viruses similar in sequence to SARS-Cov-2, we can characterize the SARS-CoV-2 ORFs we found in our first homework and determine if they are the actual genes, and think about their function. Knowing what each gene does might help us find a way to fight the infection and potentially even find treatments for COVID-19.

This was not our first rodeo with coronaviruses. In fact, we have known coronaviruses circualted in bat populations for decades. This seems like a reasonable starting point.

In order to find similar genes in other genomes that we know more about, we have to figure out with genes match with what...sounds like a great task for an alignment algorithm!

```bash
#start python
module load python
python
```

Let's get this party started. While in your python environment, lets import some important packages and commands.

```python
import numpy as np
from Bio import Entrez
# In order to import from the python file without hassle, we add the current directory to the python path
import sys; sys.path.append(".")
email = "<your email>@charlotte.edu" #It's nice to tell useful, trustworthy resources you are using them
```

## Finding SARS-CoV-2's closest relatives

In the first part of this assignment, we will analyze several coronavirus genomes to determine which are most similar to SARS-CoV-2. We will then use the viruses with the most similar sequence as reference genomes to assign functions to SARS-CoV-2 ORFs.

To find the closest reference genomes, we will align the SARS-CoV-2 sequence to each candidate coronavirus genome. The virus genomes with the highest alignment scores will be our most closely-related viruses.

Below is a list of viruses we are going to pull from Entrez. It is list of coronaviruses and their corresponding NCBI accession codes. We have already copied over our SARS-CoV-2 genome.  

```python
accession_codes = {
    # 6 known human coronaviruses
    "Human-SARS": "NC_004718",
    "Human-MERS": "NC_019843",
    "Human-HCoV-OC43": "NC_006213",
    "Human-HCoV-229E": "NC_002645",
    "Human-HCoV-NL63": "NC_005831",
    "Human-HCoV-HKU1": "NC_006577",
    
    # Bat
    "Bat-CoV MOP1": "EU420138",
    "Bat-CoV HKU8": "NC_010438",
    "Bat-CoV HKU2": "NC_009988",
    "Bat-CoV HKU5": "NC_009020",
    "Bat-CoV RaTG13": "MN996532",
    "Bat-CoV-ENT": "NC_003045",
    
    # Other animals
    "Hedgehog-CoV 2012-174/GER/2012": "NC_039207",
    "Pangolin-CoV MP789": "MT121216",
    "Rabbit-CoV HKU14": "NC_017083",
    "Duck-CoV isolate DK/GD/27/2014": "NC_048214",
    "Feline infectious peritonitis virus": "NC_002306",  # cat
    "Giraffe-CoV US/OH3/2003": "EF424623",
    "Murine-CoV MHV/BHKR_lab/USA/icA59_L94P/2012": "KF268338",  # mouse
    "Equine-CoV Obihiro12-2": "LC061274",  # horse
}
```

### LQ10.1
What is the 'type' of object we just made storing all the names of these sequences? Hint: We used this function last lab.

You can easily read FASTA files using the SeqIO.read function from biopython. Make the following function in your terminal. 

```python
from Bio import Entrez
import time

def fetch_coronavirus_fasta(accession_codes, email="your.email@example.com"):
    """
    Fetch FASTA sequences from Entrez using given accession codes.
    
    Args:
        accession_codes (dict): Dictionary of {name: accession_number}
        email (str): Valid email for NCBI API usage
    
    Returns:
        dict: {name: fasta_string} dictionary containing FASTA-formatted sequences
    """
    Entrez.email = email
    sequences = {}
    
    for name, accession in accession_codes.items():
        try:
            with Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text"
            ) as handle:
                fasta_data = handle.read().strip()
                sequences[name] = fasta_data
                print(f"Retrieved {name} ({accession})")
                time.sleep(0.35)  # NCBI recommends <3 requests/sec
            
        except Exception as e:
            print(f"Error retrieving {name} ({accession}): {str(e)}")
            sequences[name] = None  # Store None for failed retrievals
    
    return sequences


```
Now that we defined the function, we are going to run it.

## LQ10.2
What are the two variables we defined above that we use to run as input into the function?

```python
fasta_sequences = fetch_coronavirus_fasta(accession_codes, email)

# Save all sequences to a combined FASTA file
with open("coronavirus_sequences.fasta", "w") as f:
    for name, fasta in fasta_sequences.items():
        if fasta:
            f.write(fasta + "\n\n")

# Access individual FASTA string
print(fasta_sequences["Human-SARS"][:200])  # Show first 200 characters
```

## LQ10.3

Show the first 200 characters of the bat coronavirus sequence that the first SARS-CoV-2 genome was first aligned to. Hint, this is the paper we read for Reading Quiz #12.


Here is a function to get all of the sequences to print their lenghts nicesly.
```python
from Bio import Entrez
from Bio import SeqIO
from io import StringIO
import time

def calculate_sequence_lengths(sequences, accession_codes):
    print("\n{:40} | {:15} | {}".format("Virus Name", "Accession", "Sequence Length"))
    print("-" * 70)
    
    for name, accession in accession_codes.items():
        fasta = sequences.get(name)
        if not fasta:
            print(f"{name[:40]:40} | {accession:15} | {'Retrieval failed':15}")
            continue
            
        try:
            record = SeqIO.read(StringIO(fasta), "fasta")
            print(f"{name[:40]:40} | {accession:15} | {len(record.seq):,} bp")
        except Exception as e:
            print(f"{name[:40]:40} | {accession:15} | {'Invalid format':15}")
```

Now we run the function.
```python
# Calculate and display lengths
calculate_sequence_lengths(fasta_sequences, accession_codes)
```
## LQ10.4
What is the name of the longest coronavirus sequence?

Okay, now we need to do some aligning! Notice that computing an alignment between two sequences of length N and M and requires the computation of a dynamic programming table with NxM entries. The size of this matrix is sufficiently small for short sequences. But even the short genomes, like viral genomes, are generally too long for this approach. Instead of computing similarities from the entire genomes, we will only focus on the spike protein sequence. In general, we would prefer to align the whole nucleotide sequences, but for this assignment, considering only a spike protein will be sufficient. The spike protein is also one of the essential parts of any coronavirus, as it is the one that grants the virus entry to host cells.


Let's find some similar regions in each viral sequence. We spent a lot of time talking about the spike protein. Let's find the spike protein sequence first in SARS-CoV-2. Then we need to find that same sequence in all of the other genomes we just downloaded.

```python
from Bio import SeqIO
from Bio.Seq import Seq

def extract_subsequence(file_path: str, strand: str, start: int, stop: int) -> str:
    """
    Extracts a subsequence from a local FASTA file with strand support.
    
    Args:
        file_path: Path to the local FASTA file
        strand: '+' for forward, '-' for reverse complement
        start: 1-based start position (inclusive)
        stop: 1-based end position (inclusive)
        
    Returns:
        DNA sequence as string
    """
    # Read the FASTA file
    with open(file_path, "r") as handle:
        record = next(SeqIO.parse(handle, "fasta"))
    
    genome = record.seq
    
    # Validate positions
    max_length = len(genome)
    if not (1 <= start <= stop <= max_length):
        raise ValueError(f"Positions must satisfy 1 ≤ start ≤ stop ≤ {max_length}")
    
    # Extract subsequence (0-based slicing)
    subseq = genome[start-1:stop]
    
    # Handle reverse complement
    if strand == "-":
        subseq = subseq.reverse_complement()
    elif strand != "+":
        raise ValueError("Strand must be '+' or '-'")
    
    return str(subseq)

```
Now it is your turn to write out the function. The spike protein starts 21562 and ends at 25384. It is along the forward strand.

```
# Example usage
file = "sars_cov_2.fasta"  # Replace with your actual file name
#write out and run the function here with the proper inputs
spike = extract_subsequence()
print(spike)
```
## L10.5
Paste your command to extract the spike protein.


We are now going to get all of the spike proteins from our genomes we downloaded.

```python
from Bio import SeqIO
from Bio import Entrez
import re

def get_spike_sequence(accession):
    Entrez.email = "your_email@example.com"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in record.features:
        if feature.type == "CDS":
            if "gene" in feature.qualifiers and feature.qualifiers["gene"][0] == "S":
                return feature.extract(record.seq)
            elif "product" in feature.qualifiers and re.search(r"spike", feature.qualifiers["product"][0], re.IGNORECASE):
                return feature.extract(record.seq)
    
    return None
```

Iterate through the accession codes and extract spike sequences:
```python
spike_sequences = {}

for virus, accession in accession_codes.items():
    spike_seq = get_spike_sequence(accession)
    if spike_seq:
        spike_sequences[virus] = spike_seq
    else:
        print(f"Spike sequence not found for {virus}")

# Handle SARS-CoV-2 separately (note this gives slightly different output than our spike
sars_cov2 = SeqIO.read("sars_cov_2.fasta", "fasta")
spike_sequences["SARS-CoV-2"] = sars_cov2.seq[21561:25384]
```
## L10.6
How does "spike" and "spike_sequences["SARS-CoV-2"]" differ? 


We now are going to use a Needleman-Wusnch algorithm to align each sequence to our SARS-CoV-2 spike protein.

```python
import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1))
    
    for i in range(m + 1):
        score_matrix[i][0] = i * gap
    for j in range(n + 1):
        score_matrix[0][j] = j * gap
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = score_matrix[i - 1][j - 1] + (match if seq1[j - 1] == seq2[i - 1] else mismatch)
            delete = score_matrix[i - 1][j] + gap
            insert = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(match_score, delete, insert)
    
    return score_matrix[m][n]
```

```python
# Assuming spike_sequences is a dictionary containing the spike protein sequences
sars_cov2_spike = spike_sequences["SARS-CoV-2"]

# Example usage with another spike sequences
virus1 = "Human-SARS"

alignment_score1 = needleman_wunsch(sars_cov2_spike, spike_sequences[virus1])

print(f"Alignment score between SARS-CoV-2 and {virus1}: {alignment_score1}")
```
## L10.7
What is the alignment score for Human-SARS?

Your job is to now compare the alignment score to each genome. You can either code the functions by hand or, for five bonus points, you can write a loop to print each score. Happy scripting!

