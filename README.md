# Lab10 SARS-CoV-2: Learn how to identify gene function through orthology

Let's review. In the Lab 9, we wrote used a python script to retrieve ORFs from the SARS-CoV-2 genome. Some retrieved ORFs may correspond to the actual genes, some were there by chance, and some because the start codon also encodes for methionine, so the retrieved ORF was a (short) part of the more extended coding sequence. We used protein hydrophobicity to find potential transmembrane proteins. But what about other types of proteins? How can we tell which ORFs correspond to the actual genes?

This is where evolution comes. If two species are morphologically different, their genetic material must also differ substantially. Right? Not exactly. Evolution conserved many genes crucial for survival. The degree of conservation is surprising. For example, while we would expect to share many genes with the chimpanzee, our closest relative, we would expect almost no similarity in our genetic material with bananas. The numbers tell us otherwise. While sharing 96% of our genes with the chimpanzee, we also share 60% of our genes with bananas! These similar genes in humans and bananas encode proteins that carry similar functions.

Two genes with similar sequences may perform similar functions. Genes like these are called homologous genes. Or, more precisely, in our case, these would be orthologous genes (see the figure below). A pair of orthologous genes are genes from two related species originating from a common ancestor (e.g., humans and chimpanzees). Due to the evolution, the two nucleotide or protein sequences (random accumulation of mutations due to drift or nonrandom accumulation of mutations due to selection) might differ slightly but not overwhelmingly.

<img width="1039" alt="Screenshot 2025-03-24 at 10 41 24 PM" src="https://github.com/user-attachments/assets/d638a709-1009-4a3d-a3f3-0e7a937f02b0" />

Let's first set up our lab environment. Log into the cluster.

```bash
mkdir lab_10
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
    "Bat-CoV_MOP1": "EU420138",
    "Bat-CoV_HKU8": "NC_010438",
    "Bat-CoV_HKU2": "NC_009988",
    "Bat-CoV_HKU5": "NC_009020",
    "Bat-CoV_RaTG13": "MN996532",
    "Bat-CoV-ENT": "NC_003045",
    
    # Other animals
    "Hedgehog-CoV_2012-174/GER/2012": "NC_039207",
    "Pangolin-CoV_MP789": "MT121216",
    "Rabbit-CoV_HKU14": "NC_017083",
    "Duck-CoV_isolate DK/GD/27/2014": "NC_048214",
    "Feline_infectious_peritonitis_virus": "NC_002306",  # cat
    "Giraffe-CoV_US/OH3/2003": "EF424623",
    "Murine-CoV_MHV/BHKR_lab/USA/icA59_L94P/2012": "KF268338",  # mouse
    "Equine-CoV_Obihiro12-2": "LC061274",  # horse
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
#Write out and run the function here with the proper inputs, you will need to adjust the code below to include the proper parameters.
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
    Entrez.email = "lyohe1@charlotte.edu"
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

Write this to .fasta file so we can access this file later on in Part 2 of the lab.

```python
from Bio.SeqRecord import SeqRecord  # Add this import at the top

# Create records list for FASTA writing
records = []
for virus_name, sequence in spike_sequences.items():
    record = SeqRecord(sequence, id=virus_name, description="")
    records.append(record)

# Write to FASTA file
SeqIO.write(records, "spike_proteins_nucleotides.fasta", "fasta")
```

# PART 2: Multiple Sequence Alignment and Phylogenetics

Multiple sequence alignment aligns multiple sequences, but its inner workings are a bit complicated. This type of alignment is used on a large number of more or less related sequences in order to infer homology and build evolutionary trees. My multiple sequence aligner of choice is mafft, which can be called from inside python.

Assuming you are starting a new session, navigate to your lab_10 folder and in your bash terminal, we need to load MAFFT, 

```bash
cd lab_10
module load python
module load mafft
```

We also need to prep our sequences for alignment. We are going to align the nucleotide sequencees of the spike protein genes that we pulled from each coronavirus. However, a very pesky problem often occurs when FASTA files get written to output. Let's have a look.

Say we want to inspect how long each of our spike sequences are. (This function is run in bash; it comes from Lab9)
```bash
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'    spike_proteins_nucleotides.fasta
```
##L10.8a 
Why is it printing the length of sequences many times over for each species?

Let's clean this up.
```bash
awk -v ORS= '/^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' spike_proteins_nucleotides.fasta > spike_proteins_nucleotides_chomp.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'    spike_proteins_nucleotides_chomp.fasta
```
##L10.8b 
What did the first awk function in the above two lines of code do to transform our file?

Notice they are not all the same length

##L10.9
The longest spike protein in our dataset is _____ bp and it belongs to _____.


In the first part of the lab, we performed a global alignment (Needleman-Wuncsh). We have aslo learned how to do Smith-Waterman by hand, which performs a local alignment, and can be implemented within the context of a multiple sequence alignment.

Key differences between multiple sequence alignment (MSA) and pairwise alignment:

|Feature    |Pairwise Alignment	 |Multiple Sequence Alignment |
| -------- | ------- | ------- |
|Scope    | 	Aligns two sequences	| Aligns ≥3 sequences simultaneously| 
|Complexity	| Less computationally intensive	| Requires heuristics due to NP-hard complexity| 
|Biological context| 	Identifies local/global matches between two sequences| 	Preserves evolutionary relationships across all sequences| 
|Common algorithms| 	Needleman-Wunsch (global), Smith-Waterman (local)| 	Progressive methods (e.g., FFT-NS-2), iterative refinement (e.g., G-INS-i)
|Output usage| 	Direct sequence comparison| 	Input for phylogenetic analysis, structural modeling| 


In MAFFT, pairwise Smith-Waterman alignments inform the iterative refinement process in strategies like H-INS-i, where local alignment information gets incorporated into the multiple alignment objective function. We are first going to implement the H-INS-i function within MAFFT. MAFFT's H-INS-i method combines local pairwise alignment information with iterative refinement to improve multiple sequence alignment accuracy. 

Core mechanism of H-INS-i:
*Pairwise local alignment generation: Uses FASTA's Smith-Waterman implementation (fasta34) to compute all possible local alignments between sequence pairs.
*Objective function construction: Creates a scoring system combining Weighted Sum-of-Pairs (WSP) score and an importance matrix derived from pairwise alignment consistency
*Iterative refinement process: repeatedly optimizes the alignment through dynamic programming realignment of non-conserved regions and anchoring of conserved blocks to maintain structural integrity
*Progressive incorporation of pairwise alignment evidence into the multiple alignment


Here we are calling an external program within python. You can also run mafft as a standalone program, like we have previously done using the cluster.
```python
python
import subprocess

#load our extracted spike protein sequences
in_file="spike_proteins_nucleotides_chomp.fasta"
out_file="aligned_spike_local.fasta"
subprocess.call(["mafft", "--maxiterate", "1000", "--localpair", "--out", out_file, in_file])
```
This is going to take a minute or so to run. Pay attention to the output and the steps it is doing as it runs.

Now we are going to align our sequences using a different algorithm, and happens to be the default within MAFFT, FFT-NS-2 (fast but rough). The key advantage of this approach is that it drastically reduces CPU time compared to traditional methods like CLUSTALW while maintaining comparable accuracy and is suitable for large datasets due to its computational speed. In contrast, H-INS-i focuses on refining alignments iteratively by incorporating pairwise local alignment information, making it better suited for datasets with complex sequence variability.

```python
import subprocess
in_file="spike_proteins_nucleotides_chomp.fasta"

out_file_2="aligned_spike_fast.fasta"
subprocess.call(["mafft", "--out", "aligned_spike_fast.fasta", in_file])
```
Notice how much faster this ran! Depending on your data set, you need to decide and defend which algorithm is best to use. Faster can come at the cost of lower accuracy. We are going to compare outcomes of using different alignemnt appraoches. Here is a summary of the differences. 

|Feature	|FFT-NS-2	|H-INS-i    |
| -------- | ------- | ------- |
|Alignment Type    |Progressive    |Iterative refinement|
|Algorithm Used|    Fast Fourier Transform (FFT)|	Smith-Waterman for local alignments|
|Accuracy|    Good for conserved sequences|    Higher accuracy for variable domains|
|Speed|    Faster due to FFT acceleration|	Slower due to iterative refinement|
|Best Application|    Large datasets with conserved regions|    Complex alignments with variable domains|

Exit out of python.
```python
exit()
```

We again need to clean up the new lines in our files. In your bash terminal, perform the following:
```bash
awk -v ORS= '/^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' aligned_spike_fast.fasta > aligned_spike_fast_chomp.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'   aligned_spike_fast_chomp.fasta

awk -v ORS= '/^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' aligned_spike_local.fasta > aligned_spike_local_chomp.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'   aligned_spike_local_chomp.fasta
```

## L10.10a

TRUE or FALSE: All sequences in the an alignment file have the same sequence length.

## L10.10b
The aligned sequences are (shorter/longer) than the unaligned sequences. The (FFT-NS-2/H-INS-i) algorithm resulted in a longer alignment file. This means this algorithm resulted in more (mutations/gaps/matches).

The first step is to import the (already aligned) sequences into our program. Using AlignIO to read in the data, and store it in a variable aln_fast. Print aln_fast to confirm this step worked correctly.
```bash
python
```

```python
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Read the alignment from a FASTA file
aln_fast=AlignIO.read("aligned_spike_fast_chomp.fasta", "fasta")

# Print the alignment to confirm it was read correctly
print(aln_fast)
```
Repeat this with your local alignment file "aligned_spike_local_chomp.fasta" and name your variable "aln_local".

Now that we have our data, we’re interested in how similar (or different) the sequences are from each other. The more similar, the more likely they are to be evolutionarily close. Therefore, we use a distance matrix to store information on how similar each DNA sequence is from every other sequence. We will creat a Distance Calculator object with the ‘identity’ model. Use it to calculate the distance matrix of the alignment, and store it in a variable "dm_fast".

```python
# Create a DistanceCalculator object with the 'identity' model
calculator = DistanceCalculator('identity')

# Calculate the distance matrix using the alignment
dm_fast = calculator.get_distance(aln_fast)

# Print the distance matrix to confirm it worked correctly
print(dm_fast)
```

Repeat this for your "aln_local" variable and store your distance matrix for it as "dm_local".

Now, the real magic… creating the phylogenetic tree! We can use our Distance Tree Constructor to make a tree from our distance matrix. We will use the UPGMA algorithm in Distance Tree Constructor to create a tree and store it in a variable "tree_fast".

```python
# Create a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()

# Use the UPGMA algorithm to construct a tree
tree_fast = constructor.upgma(dm_fast)

# Print the resulting tree to confirm it worked correctly
print(tree_fast)
```

Repeat this for your "dm_local" and store it as variable "tree_local".

Great, you’re done! Except… you’re probably interested in how your tree turned out! The simplest way to visualize your tree is to use the following methods to plot.

```python
# Print the tree as ASCII art
Phylo.draw_ascii(tree_fast)
Phylo.draw_ascii(tree_local)
```
## L10.11a
Paste your two resulting trees. 

## L10.11b
Do your trees agree with one another in terms of relationships? Describe two differences that you see.

## L10.11c
If we know that "Feline_infectious_peritonitis_virus" is an outgroup, which of the two trees (fast or local) demonstrates a more accurate representation of this.

## L10.11d
What are the three most closely related viruses to SARS-CoV-2 in the local tree? Are they the same in both trees? 

## L10.11e
Were these three species the same species with the "top scores" of your Needleman-Wunsch pairwise alignments (Part 1 of the lab)?
