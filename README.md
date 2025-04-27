# Nanopore-LSHTM-Omics
Nanopore LSHTM 'Omics Course

## Identifying drug-resistance from *Escherichia coli* WGS sequenced using Oxford Nanopore
In the data repository you have been provided with seven *Escherichia coli* whole genome sequences that have been sequenced using Oxford Nanopore (MinION). These samples were collected from a long-term health care facility in Japan by a previous study (ENA Project Accession: PRJDB9189). These sequences have already been basecalled and we have provided some QC reports. Your task is to analyse the Nanopore fastq data to identify drug-resistance genes. 

As we have done previously, the first step will be to perform QC before going onto downstream analysis. 

### Step 1. Quality Control and Contamination
**Step 1**: Let's first check the quality of the nanopore fastq data. To do this we use a 'sequencing_summary.txt' file generated using dorado. We check quality using **PycoQC**:
```
conda activate nanopore
cd nanopore/data
mkdir pycoqc
pycoQC -f sequencing_summary.txt -o pycoqc/Ecoli_Japan_1_PycoQC.html
```
!!! question Open the html file in firefox. What can you see from the output? How do you think the read lengths compare to Illumina? How does read quality vary over time?

Next we should trim our files and filter for high-quality reads. To trim the nanopore reads we have used **Chopper** (You could also used Porechop_abi, see Tips and Tricks) to identify any adapter sequences and remove them. This is an important step, especially when performing genome assembly. 
```
gunzip -c "Ecoli_Japan_1.fastq.gz" | chopper -q 10 -l 500 --headcrop 50 --tailcrop 50 | gzip > "Ecoli_Japan_1_trim.fastq.gz"
```

We should also check to see if the sequence contains any contaminants. We have used **Kraken2** to do this. The Kraken reports have been generated for you as they require large databases. We can visualise the results using **recentrifuge (rcf)**, here there is a dependency clash in our conda environment, so let's quickly install using **pip**:
```
mkdir kraken
pip install recentrifuge
rcf -n ./taxdump/ -k Ecoli_Japan_1.koutput.txt -o kraken/Ecoli_Japan_1_kraken.html
```
!!! question  What can you determine from the Kraken2 outputs? Is the data clean?  There are a few other genera included in the output — do you think they might be contaminants?

Just in case we will filter our fastq files to remove any contaminants. We won't be too stringent with the filtering- we will use the taxa id from NCBI for Enterobacteriaceae (543). To do this we will use **KrakenTools**:

```
python ./taxdump/extract_kraken_reads.py --include-children --fastq-output -t 543 -k "Ecoli_Japan_1.koutput.txt" -r "Ecoli_Japan_1.kreport.txt" -s "Ecoli_Japan_1_trim.fastq.gz" -o "./kraken/Ecoli_Japan_1.kraken_filtered.fq" 
```

### Step 2. Assembly
**Step 2**: Let's now perform genome assembly using our high-quality, trimmed and filtered nanopore fastq. There are a few different tools we can use, but for Nanopore data **Flye** performs well. This may take 5 minutes or so to run. 

```
flye --nano-raw "./kraken/Ecoli_Japan_1.kraken_filtered.fq" --read-error 0.03 --genome-size 4.6m --out-dir flye_output --threads 16 
```

### Step 3. Assess the Quality of your Assembly
**Step 3**: Now the assembly is ready we will test its quality using **BUSCO**. BUSCO assesses the quality of genome assemblies by looking at the percentage of conserved genes. 
```
busco -i flye_output/assembly.fasta -l enterobacteriaceae_odb12 -c 4 -m genome -o assembly_QC
```
!!! question  What do you think, is it a good assembly?

### Step 4. Annotate your assembly
**Step 4**: Ok so the assembly could be better, but how about we annotate some genes? We will do this using **Prokka**. This is a specialist software specific to bacteria. There are lots of other tools for annotating Eukaryotes and beyond, however, this tool is well-developed for bacteria species.
```
prokka --outdir assembly_annotation --prefix Ecoli_Japan_1 flye_output/assembly.fasta
```
Inside the assembly_annotation directory you will see "Ecoli_Japan_1.gff". This is the most common format for annotations. If you like, you can load this into **igv** to visualise the annotation (Load Genome: Ecoli_Japan_1.fna; Load files:Ecoli_Japan_1.gff). We can analyse the annotations further downstream (See Step 6).

### Step 5. Predicting Drug-resistance
**Step 5**: Great! We have an assembly and we have annotated it. How are the assemblies useful for downstream analysis? For bacteria specifically, genome assemblies can be input into specialist software to mass screen the contigs (bits of the assembly) for antimicrobial resistance genes, virulence genes or plasmids. **Abricate** is a handy suite of tools that can help us do this by scanning databases of these genes.
To identify resistance genes:
```
mkdir abricate_results
abricate --db resfinder --quiet flye_output/assembly.fasta > abricate_results/resistance_results.txt
```
To identify virulence genes:
```
abricate --db vfdb --quiet flye_output/assembly.fasta > abricate_results/virulence_results.txt

```
To identify plasmids:
```
abricate --db plasmidfinder --quiet flye_output/assembly.fasta > abricate_results/plasmid_results.txt
```
!!! question  Take a look at some of the resistance genes and plasmids- why might the blaCTX-M-27 gene be of concern?


### Step 6. Pan-genome Analysis
**Step 6**: With growth in the size of datasets, there is a need to understand key processes such as selection and evolution taking place in bacteria populations. For example, bacteria can transfer genes to one another, otherwise known as 'horizontal gene transfer' which can spread virulence and resistance genes. One way to identify some of these diffrences is to look *core* or *accessory* genes within a population: A Pangenome. **Roary** and **Pirate** can be used to compare our gff annotation files to one another and construct a pangenome. We will use a different conda environment. This may take a moment or two!
```
conda deactivate
conda activate roary
cd all_gffs
roary -e --mafft -p 4 *.gff
conda deactivate
```
Open firefox and load https://jameshadfield.github.io/phandango/#/ to view the roary output files. Drag and drop the "assembly_annotations/accessory_binary_genes.fa.newick" and "assembly_annotations/gene_presence_absence.csv" files into phandango. You should get something that looks like this:
![image](https://github.com/user-attachments/assets/11f76e30-3125-43ae-9614-59e5cb984cb5)

!!! question  What do you think the blue blocks mean?

### Step 7. Mapping and Variant Calling

The mapping and variant calling tools for Nanopore are different than those typically used with Illumina data. Nanopore fastq data can be mapped to an E.coli reference genome using **Minimap2**. Notice the '-ax map-ont' signalling we are using Nanopore data. Minimap2 creates a SAM file, a type of alignment file. We need to convert this to a BAM file to make the alignment compatible with other software. We do this using **samtools**.
```
conda activate nanopore
cd ~/nanopore/data
minimap2 -ax map-ont Ecoli_reference.fasta Ecoli_Japan_1_trim.fastq.gz  | samtools sort -o Ecoli_Japan_1_aln.bam
```
For variant calling we tend to use variant callers that have been designed specifically for Nanopore data. This includes **Clair3** and **Freebayes**, which can generate outputs that are generally compatible with GATK-based workflows. These tools can be fairly slow so feel free to trial them in your own time!

### Can you do this analysis with Illumina data?
**Yes!** Anything from steps 3 to 6 can be done with Illumina data. You might just need to switch up your tool box. For bacteria genome assembly with short read data, there are some fantastic tools available such as **shovill**. The short read assemblies can be annotated using **prokka** and analysed using **abricate** as above. Here, we wanted to demonstrate the utility of nanopore data for real-time sequencing and help you gain some experience in the downstream analysis.

### Advanced
1. Now that you have identified the drug-resistance genes for one sample, how about performing the steps on the remaining samples e.g. Ecoli_Japan_2_trim.fastq.gz, what are your results? Please note you will not have to perform PycoQC on the remaining samples.
2. You can see in step 2 used a read error of **0.03**, if you have time you could change this to **0.05**. How does this affect your BUSCO and QUAST results?

### Tips and Tricks
You may wish in the future to trim your own reads and make Kraken2 reports by yourself. These are the commands we have used to generate the data for the practical. These commands take some time or require large databases which is why we have performed these steps in advance.

```
# Kraken2 
kraken2_client --host-ip XX.XX.XX.XX --sequence "Ecoli_Japan_1.fastq.gz" --report "Ecoli_Japan_1.kreport.txt" > "Ecoli_Japan_1.koutput.txt"
# Trimming
porechop_abi -abi -i "Ecoli_Japan_1.fastq.gz" -o "Ecoli_Japan_1_trim.fastq.gz"

```
