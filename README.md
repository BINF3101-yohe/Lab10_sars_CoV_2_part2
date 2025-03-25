# Lab10 SARS-CoV-2: Learn how to identify gene function through orthology

Let's review. In the Lab 9, we wrote used a python script to retrieve ORFs from the SARS-CoV-2 genome. Some retrieved ORFs may correspond to the actual genes, some were there by chance, and some because the start codon also encodes for methionine, so the retrieved ORF was a (short) part of the more extended coding sequence. We used protein hydrophobicity to find potential transmembrane proteins. But what about other types of proteins? How can we tell which ORFs correspond to the actual genes?

This is where evolution comes. If two species are morphologically different, their genetic material must also differ substantially. Right? Not exactly. Evolution conserved many genes crucial for survival. The degree of conservation is surprising. For example, while we would expect to share many genes with the chimpanzee, our closest relative, we would expect almost no similarity in our genetic material with bananas. The numbers tell us otherwise. While sharing 96% of our genes with the chimpanzee, we also share 60% of our genes with bananas! These similar genes in humans and bananas encode proteins that carry similar functions.

Two genes with similar sequences may perform similar functions. Genes like these are called homologous genes. Or, more precisely, in our case, these would be orthologous genes (see the figure below). A pair of orthologous genes are genes from two related species originating from a common ancestor (e.g., humans and chimpanzees). Due to the evolution, the two nucleotide or protein sequences (random accumulation of mutations due to drift or nonrandom accumulation of mutations due to selection) might differ slightly but not overwhelmingly.

<img width="1039" alt="Screenshot 2025-03-24 at 10 41 24 PM" src="https://github.com/user-attachments/assets/d638a709-1009-4a3d-a3f3-0e7a937f02b0" />

Let's first set up our lab environment. Log into the cluster.

```bash
mkdir lab_10
cd lab_10
module load python
#we also need the sars_cov_2 genome for this lab
cp /projects/class/binf3101_001/sars_cov_2.fasta ~/lab_10
```

The human-banana example demonstrates that these two organisms don't have to be that closely related! If we tried to determine the gene function in humans using bananas as a reference, we could infer the role of over 60% of the genes. However, using two closely related organisms would undoubtedly give us more robust results. We would be better off choosing the genome of a mouse or a chimpanzee as our reference.

Scientists have studied hundreds of viruses and determined the functions of most of their genes. Using viruses similar in sequence to SARS-Cov-2, we can characterize the SARS-CoV-2 ORFs we found in our first homework and determine if they are the actual genes, and think about their function. Knowing what each gene does might help us find a way to fight the infection and potentially even find treatments for COVID-19.

This was not our first rodeo with coronaviruses. In fact, we have known coronaviruses circualted in bat populations for decades. This seems like a reasonable starting point.

In order to find similar genes in other genomes that we know more about, we have to figure out with genes match with what...sounds like a great task for an alignment algorithm!

```bash
#start python
python
```
