## Define the path of your raw data and your scripts
data_path=$"/mnt/scratch/hussien/Sig6/data"
script_path=$"/mnt/scratch/hussien/Sig6/scripts"
out_path=$"/mnt/scratch/hussien/Sig6/output"
genome_path=$"/mnt/scratch/hussien/Sig6/genome"
#####################################
#Running the fastQC again
 ### script name fastQC.sb
#!/bin/bash --login
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hussien@msu.edu

module load FastQC/0.11.7-Java-1.8.0_162

cd /mnt/scratch/hussien/Sig6/data
fastqc -f fastq -noextract ${sample}

# to run the script 
for S in *.fastq; do sbatch  --output=${out_path}/$S.out --export=sample=$S ${script_path}/fastqc.sb; done
# the results output was generated in the same as input data
# moving the fasrtQC output to a specific folder in the output 
mkdir ${out_path}/fastqc
cd ${out_path}/fastqc
mv *.html ${out_path}/fastqc
mv *.zip ${out_path}/fastqc
#####################
# there atre some over represetation of some adapter like sequences in some samples 
###################
# Running Trimmomatic
###################
#!/bin/bash --login
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hussien@msu.edu
 Trimmomatic/0.38-Java-1.8.0_162
 
 java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 ${sample} ${sample}.out ILLUMINACLIP:/mnt/research/common-data/Bio/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 
 ### to run the trimmomatic
 for S in *.fastq; do sbatch  --output=${out_path}/$S.out --export=sample=$S ${script_path}/trimmomatic.sb; done
 
 mkdir trimmomatic_output
 cd trimmomatic_output/
 mv ../*.out .
 
 
  # rerun the fastQC to make sure everytinig is ok (no adapters)
  
for S in *.out; do sbatch  --output=${out_path}/$S.output --export=sample=$S ${script_path}/fastqc.sb; done
mkdir ${out_path}/FastQC-after-trimmomatic
cd mkdir ${out_path}/FastQC-after-trimmomatic
mv ${data_path}/trimmomatic_output/*.html .
mv ${data_path}/trimmomatic_output/*.zip .

###### the data looks good (No adapters)
## Running the STAR fo mapping 
##Downloading the genome data:
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff

-running STAR aligner
#1-convert gff3 to gtf file

module load cufflinks
gffread TAIR10_GFF3_genes.gff  -T -o TAIR10_gtf_genes.gtf
sed -i 's/>1 Chr/> CHROMOSOME1/g' TAIR10_gtf_genes.gtf
sed -i 's/>2 Chr/> CHROMOSOME2/g' TAIR10_gtf_genes.gtf
sed -i 's/>3 Chr/> CHROMOSOME3/g' TAIR10_gtf_genes.gtf
sed -i 's/>4 Chr/> CHROMOSOME4/g' TAIR10_gtf_genes.gtf
sed -i 's/>5 Chr/> CHROMOSOME5/g' TAIR10_gtf_genes.gtf
# for indexing
mkdir ${genome_path}/STAR_index
module load STAR/2.6.0c

STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${genome_path}/STAR_index --genomeFastaFiles ${genome_path}/TAIR10_chr_all.fas --sjdbGTFfeatureExon exon --sjdbGTFfile ${genome_path}/TAIR10_gtf_genes.gtf genomeSAindexNbases 14

##Running the alignment 
mkdir ${output_path}/STAR_out
cd ${output_path}/STAR_out
cat ${script_path}/sample_list.txt | xargs mkdir  ## creat outpot files


#STAR --runThreadN 1 \ # can be increased if sufficient computational power is available
#--genomeDir ${genome_path}/STAR_index \
#--readFilesIn ${data_path}/trimmomatic_output/35S-SIG6-D1-1_S1_L003_R1_001.fastq.out
#--outFileNamePrefix alignment_STAR / ${out_path}/STAR_out/\ #/path/to/output/dir/prefix
#--outFilterMultimapNmax 1 \ # only reads with 1 match in the reference will be returned as aligned
#--outReadsUnmapped Fastx \ # will generate an extra output file with the unaligned reads
#--twopassMode Basic\ # STAR will perform mapping , then extract novel junctions which will be inserted into the genome index which will then be used to re - map all reads
#--outSAMtype BAM SortedByCoordinate Unsorted \ #Mapping reads and generating unsorted and coordinate-sorted BAM files, BAM conversion and coordinate-sorting while mapping
#--alignIntronMin 21 --alignIntronMax 6000 \
#--quantMode GeneCounts \ #STAR will count number reads per gene while mapping, A read is counted if it overlaps (1nt or more) one and only one gene, The counts coincide with those produced by htseq-count with default parameters, STAR outputs read counts per gene into ReadsPerGene.out.tab file
#--sjdbGTFfeatureExon ${genome_path}/TAIR10_gtf_genes.gtf

cd ${data_path}/trimmomatic_out

for S in *.out; do sbatch  --export=sample=$S ${script_path}/STAR_align.sh; done

### STAR_align.sh
#!/bin/bash --login
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hussien@msu.edu
data_path=$"/mnt/scratch/hussien/Sig6/data"
script_path=$"/mnt/scratch/hussien/Sig6/scripts"
out_path=$"/mnt/scratch/hussien/Sig6/output"
genome_path=$"/mnt/scratch/hussien/Sig6/genome"

module load STAR/2.6.0c

STAR --runThreadN 1 --genomeDir ${genome_path}/STAR_index --readFilesIn ${sample} --outFileNamePrefix ${out_path}/STAR_out/${sample}/${sample}-alignment_STAR --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --twopassMode Basic --outSAMtype BAM SortedByCoordinate Unsorted --alignIntronMin 21 --alignIntronMax 6000 --quantMode GeneCounts --sjdbGTFfeatureExon ${genome_path}/TAIR10_gtf_genes.gtf
######
#running the HTSeq fo counts
###preparing the output folders 
mkdir ${output_path}/HTseq_out
####collecting all the bam files to start the HTSeq 
mkdir ${out_path}/STAR_out/all_bam_files
cd ${out_path}/STAR_out/all_bam_files
cp ${out_path}/STAR_out/*/*-alignment_STARAligned.sortedByCoord.out.bam . 

#######HTSeq.sh####

#!/bin/bash --login
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hussien@msu.edu
data_path=$"/mnt/scratch/hussien/Sig6/data"
script_path=$"/mnt/scratch/hussien/Sig6/scripts"
out_path=$"/mnt/scratch/hussien/Sig6/output"
genome_path=$"/mnt/scratch/hussien/Sig6/genome"

htseq-count --format=bam --stranded=reverse ${out_path}/STAR_out/all_bam_files/${sample} ${genome_path}/TAIR10_gtf_genes.gtf > ${out_path}/HTseq_out/${sample}.htseqcount.txt

#####to call the script for all the samples in STAR_out folder


for S in ${out_path}/STAR_out/all_bam_files/* ; do sbatch  --export=sample=$S ${script_path}/HTSeq.sh; done
