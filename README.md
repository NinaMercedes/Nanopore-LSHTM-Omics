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

### Step 1. Quality Control and Contamination
**Step 1**: Let's first check the quality of the nanopore fastq data by mapping to an E.coli reference genome. To do this we use a 'sequencing_summary.txt' file generated using dorado. We check quality using **PycoQC**:
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

Just in case we will filter our fastq files to remove any contaminants. We won't be too stringent with the filtering- we will use the taxa id from NCBI for Enterobacteriaceae (543). To do this we will use **KrakenTools**:

```
python ./taxdump/extract_kraken_reads.py --include-children --fastq-output -t 543 -k "Ecoli_Japan_1.koutput.txt" -r "Ecoli_Japan_1.kreport.txt" -s "Ecoli_Japan_1_trim.fastq.gz" -o "./kraken/Ecoli_Japan_1.kraken_filtered.fq" 
```

### Step 2. Assembly
**Step 2**: Let's first now perform genome assembly using our high-quality, trimmed and filtered nanopore fastq. There are a few different tools we can use, but for Nanopore data **Flye** performs well. This may take 5 minutes or so to run.

```
flye --nano-raw "./kraken/Ecoli_Japan_1.kraken_filtered.fq"  --genome-size 4.6m --out-dir flye_output --threads 16 
```

### Step 3. Assess Quality of your Assembly
**Step 3**: Now the assembly is ready we will test its quality using **BUSCO**. BUSCO assesses the quality of genome assemblies by looking at the percentage of conserved genes. 
```
busco -i flye_output/assembly.fasta -l enterobacteriaceae_odb12 -c 16 -m genome -o assembly_QC
```
What do you think, is it a good assembly?

### Step 4. Annotate your assembly
**Step 4**: Ok so the assembly could be better, but how about we annotate some genes? We will do this using **Prokka**. This is a specialist software specific to bacteria. There are lots of other tools for annotating Eukaryotes and beyond, however, this tool is well-developed for bacteria species.
```
prokka --outdir assembly_annotation --prefix Ecoli_Japan_1 flye_output/assembly.fasta
```
Inside the assembly_annotation directory you will see "Ecoli_Japan_1.gff". This is the most common format for annotations. If you like, you can load this into **igv** to visualise the annotation (Load Genome: Ecoli_Japan_1.fna; Load files:Ecoli_Japan_1.gff). We can analyse the annotations further downstream (See Step 6).

### Step 5. Predicting Drug-resistance
**Step 5**: Great! We have an assembly and we have annotated it. How are the assemblies useful for downstream analysis? For bacteria specifically, genome assemblies can be input into specialist software to mass screen the contigs (bits of the assembly) for antimicrobial resistance genes, virulence genes or plasmids. **Abricate** is a handy suite of tools that can help us do this by scanning databases of these genes.
Some resistance genes:
```
mkdir abricate_results
abricate --db resfinder --quiet flye_output/assembly.fasta > abricate_results/resistance_results.txt
```
Some virulence genes:
```
abricate --db vfdb --quiet flye_output/assembly.fasta > abricate_results/virulence_results.txt

```
Some plasmids:
```
abricate --db plasmidfinder --quiet flye_output/assembly.fasta > abricate_results/plasmid_results.txt
```
Take a look at some of the resistance genes and plasmids- why might the blaCTX-M-27 gene be of concern?


### Step 6. Pan-genome Analysis
**Step 6**: With growth in the size of datasets, there is a need to understand key processes such as selection and evolution taking place in bacteria populations. For example, bacteria can transfer genes to one another, otherwise known as 'horizontal gene transfer' which can spread virulence and resistance genes. One way to identify some of these diffrences is to look *core* or *accessory* genes within a population: A Pangenome. **Roary** and **Pirate** can be used to compare our gff annotation files to one another and construct a pangenome. 
```

```

### Step 7. Mapping and Variant Calling

Mapping and variant calling tools are different that illumina. Nanopore fastq data by mapping to an E.coli reference genome using **Minimap2**. Notice the '-ax map-ont' signalling we are using Nanopore data. Minimap2 creates a SAM file, a type of alignment file. We need to convert this to a BAM file to make the alignment compatible with other software. We do this using **samtools**.
```
cd nanopore/bacteria/data
minimap2 -ax map-ont Ecoli_reference.fasta Ecoli_Japan_1_trim.fastq.gz  | samtools sort -o Ecoli_Japan_1_aln.bam
```
For variant calling we tend to use variant callers that have been designed specifically for Nanopore data. This includes **Clair3** and **Freebayes**, which can generate outputs that are seemingly compatible with GATK software. These tools can be fairly slow so feel free to trial them in your own time!

### Advanced
Now that you have identified the drug-resistance genes for one sample, how about performing the steps 2-5 on the remaining samples e.g. Ecoli_Japan_2_trim.fastq.gz, what are your results?

### Tips and Tricks

You may wish in the future to trim your own reads and make Kraken2 reports by yourself. These are the commands we have used to generate the data for the practical. These commands take some time or require large databases which is why we have performed these steps in advance.

```
# Kraken2 
kraken2_client --host-ip XX.XX.XX.XX --sequence "Ecoli_Japan_1_trim.fastq.gz" --report "Ecoli_Japan_1.kreport.txt" > "Ecoli_Japan_1.koutput.txt"
# Trimming
porechop_abi -abi -i "Ecoli_Japan_1.fastq.gz" -o "Ecoli_Japan_1_trim.fastq.gz"

```
