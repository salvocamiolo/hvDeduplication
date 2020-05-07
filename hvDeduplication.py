import os,sys
import argparse

parser = argparse.ArgumentParser(description="Deduplicate reads for hypervariable genes")
parser.add_argument("-1", "--read1", type=str, required=True, help="Path to the fastq file 1")
parser.add_argument("-2", "--read2", type=str, required=True, help="Path to the fastq file 2")
parser.add_argument("-c","--conda_directory",type=str, required=True,help="Path to the conda 3 environment")
parser.add_argument("-t","--num_threads",type=str, required=False,help="number of threads")
args = vars(parser.parse_args())
condaDir = args['conda_directory']
read1 = args['read1']
read2 = args['read2']
if not args['num_threads']:
    threads = 1
else:
    threads = args['num_threads']

 
hvg = ['RL12','RL13','RL5A','RL6','UL11','UL120','UL139','UL146','UL1','UL20','UL73','UL74','UL9']

for gene in hvg: 
    print("Anlysing gene %s" %gene)
    print("Creating index....")
    os.system(condaDir+"/bin/bowtie2-build ./elongedCDS/"+gene+"_elongedCDS.fasta reference >null 2>&1")
    print("Performing alignment....")
    os.system(condaDir+"/bin/bowtie2 -1 "+read1+" -2 "+read2+" -x reference -S alignment.sam -p "+threads+" >null 2>&1")
    print("Converting to bam....")
    os.system(condaDir+"/bin/samtools view -h -bS -F 4 alignment.sam > alignment.bam")
    print("Sorting....")
    os.system(condaDir+"/bin/samtools sort -o mapped.bam alignment.bam")
    #print("Number of reads before deduplication:")
    os.system(condaDir+"/bin/samtools view mapped.bam | wc -l >readCount")
    infile = open("readCount")
    beforeDedupReads = int(infile.readline().rstrip())
    infile.close()
    
    print("Picard add groups....")
    os.system(condaDir+"/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus >null 2>&1")
    print("Picard mark duplicates....")
    os.system(condaDir+"/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
    print("Discarding reads with duplicates....")
    os.system(condaDir+"/bin/samtools view -h -b -F 1024 dedupped.bam > alignment_removedDup.bam")
    #print("Number of reads before deduplication:")
    os.system(condaDir+"/bin/samtools view alignment_removedDup.bam | wc -l >readCount")
    infile = open("readCount")
    afterDedupReads = int(infile.readline().rstrip())
    infile.close()
    print(gene,beforeDedupReads,afterDedupReads)
    print("Extracting reads....")
    os.system(condaDir+"/bin/bam2fastq -o "+gene+"_dedup#.fq alignment_removedDup.bam ")

print("Concatenating reads from first fastq file")
os.system("cat *dedup_1.fastq > all_Dedup_1.fastq")
print("Concatenating reads from second fastq file")
os.system("cat *dedup_2.fastq > all_Dedup_2.fastq")
os.system("mkdir -p dedupFastqFles")
os.system("mv *edup* ./dedupFastqFiles/")
os.system("rm *.bam *.sam *.bai readCount reference* null output.metrics")




    

