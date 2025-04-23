# Nanopore-LSHTM-Omics
Nanopore LSHTM 'Omics Course
## Manual Installations
```
conda create -n nanopore_test
conda activate nanopore_test
conda install porechop_abi
```
## Identifying drug-resistance from *Escherichia coli* WGS sequenced using Oxford Nanopore reads
In the data repository you have been provided with seven *Escherichia coli* whole genome sequences that have been sequenced using Oxford Nanopore (MinION). These samples were collected from a long-term health care facility in Japan by a previous study (ENA Project Accession: PRJDB9189). These sequences have already been basecalled, trimmed and we have provided some QC reports. Your task is to analyse the Nanopore fastq data to identify drug-resistance genes. 

As we have done previously, the first step will be to perform QC before going onto downstream analysis. To trim the nanopore reads we have used **Porechop** (*ab initio*) to identify any adapter sequences and remove them. This is an important step, especially when performing genome assembly. 

1️⃣ Step1: Let's first check the quality of the nanopore fastq data by mapping to an E.coli reference genome. To do this we use a 'sequencing_summary.txt' file generated using dorado. We check quality using **PycoQC**:
```
mkdir pycoqc
pycoQC –f sequencing_summary.txt –o pycoqc/Ecoli_Japan_1_PycoQC.html
```
Open the html file. What can you see from the output? How do you think the read lengths compare to Illumina? How does read quality vary over time?

We should also check to see if the sequence contains any contaminants. We have used **Kraken2** to do this. The Kraken reports have been generated for you as they require large databases. We can visualise the results using **recentrifuge (rcf)**:
```
mkdir kraken
rcf -n ./taxdump/ -k Ecoli_Japan_1.koutput.txt -o kraken/Ecoli_Japan_1_kraken.html
```
What can you determine from the Kraken2 outputs? Is the data clean? There are a few other genera included in the output- do you think they are contaminants (or not)?

Just in case we will filter our fastq files to remove any contaminants. We won't be too stringent with the filtering. To do this we will use **KrakenTools**:






1️⃣ Step: Let's first check the quality of the nanopore fastq data by mapping to an E.coli reference genome using **Minimap2**. Notice the '-ax map-ont' signalling we are using Nanopore data. Minimap2 creates a SAM file, a type of alignment file. We need to convert this to a BAM file to make the alignment compatible with other software. We do this using **samtools**.
```
cd nanopore/bacteria/data
minimap2 -ax map-ont Ecoli_reference.fasta Ecoli_Japan_1_trim.fastq.gz  | samtools sort -o Ecoli_Japan_1_aln.bam
```

### Advanced
Now that you have identified the drug-resistance genes for one sample, how about performing the steps on the remaining samples e.g. Ecoli_Japan_1_trim.fastq.gz, what are your result



### Tips and Tricks

You may wish in the future to trim your own reads and make Kraken2 reports by yourself. These are the commands we have used to generate the data for the practical. These commands take some time or require large databases which is why we have performed these steps in advance.

```
# Kraken2 
kraken2_client --host-ip XX.XX.XX.XX --sequence "Ecoli_Japan_1_trim.fastq.gz" --report "Ecoli_Japan_1.kreport.txt" > "Ecoli_Japan_1.koutput.txt"
# Trimming
porechop_abi -abi -i "Ecoli_Japan_1.fastq.gz" -o "Ecoli_Japan_1_trim.fastq.gz"

```
