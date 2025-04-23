# Nanopore-LSHTM-Omics
Nanopore LSHTM 'Omics Course
## Manual Installations
```
conda create -n nanopore_test
conda activate nanopore_test
conda install porechop_abi
```
## Identifying drug-resistance from *Escherichia coli* WGS sequenced using Oxford Nanopore
In the data repository you have been provided with seven *Escherichia coli* whole genome sequences that have been sequenced using Oxford Nanopore (MinION). These samples were collected from a long-term health care facility in Japan by a previous study (ENA Project Accession: PRJDB9189). These sequences have already been basecalled, trimmed and we have provided some QC reports. Your task is to analyse the Nanopore fastq data to identify drug-resistance genes. 

As we have done previously, the first step will be to perform QC before going onto downstream analysis. To trim the nanopore reads we have used **Porechop** (*ab initio*) to identify any adapter sequences and remove them. This is an important step, especially when performing genome assembly. 

1️⃣ Step1: Let's first check the quality of the nanopore fastq data using **fastqc**:
```
cd nanopore/bacteria/data

```

As we have done before, let's check the quality of our sample. What do you think? Have all of the adapters been removed?

```

```

We should also check to see if the sequence contains any contaminants. We will use **Kraken2** to do this. The Kraken reports have been generated for you as they require large databases. 

### Advanced
Now that you have identified the drug-resistance genes for one sample, how about performing the steps on the remaining samples e.g. Ecoli_Japan_1_trim.fastq.gz, what are your result



### Tips and Tricks

You may wish in the future to trim your own reads and make Kraken2 reports by yourself. These are the commands we have used to generate the data for the practical. These commands take some time or require large databases which is why we have performed these steps in advance.

```
# Kraken2 

# Trimming
porechop_abi -abi -i "Ecoli_Japan_1.fastq.gz" -o "Ecoli_Japan_1_trim.fastq.gz"

```
