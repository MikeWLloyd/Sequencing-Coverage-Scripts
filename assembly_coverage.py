#https://github.com/BinPro/CONCOCT/blob/master/scripts/map-bowtie2-markduplicates.sh
#https://github.com/faircloth-lab/phyluce
import os
import subprocess
import shutil
import argparse
import sys
import glob
from phyluce.pth import get_user_path, get_user_param
from phyluce.helpers import FullPaths, is_dir, is_file

JAVA = get_user_param("java", "executable")
JAVA_PARAMS = get_user_param("java", "mem")
JAR_PATH = get_user_path("java", "jar")
FNULL = open(os.devnull, 'w')

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Calculate coverages for Trinity assemblies."""
    )
    parser.add_argument(
        "--rawreads",
        required=True,
        action=FullPaths,
        type=is_dir,
        default='',
        help="""Path to the uce-clean directory."""
    )
    parser.add_argument(
        "--input",
        required=True,
        type=is_dir,
        action=FullPaths,
        default=None,
        help="""A directory containing directories with Trinity assemblies. Directory names are taxon names. Directory contains 'contigs.fasta'"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory in which to store the summarized coverages."""
    )
    parser.add_argument(
        "--config",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""A configuration file containing names of taxa."""
    )
    parser.add_argument(
        "--threads",
        type=str,
        default=1,
        help="""The number of compute cores/threads to run."""
    )
    return parser.parse_args()


args = get_args()

f = open(args.config, 'r')
path_to_reads = args.rawreads

if not os.path.exists(args.output):
    os.makedirs(args.output)
else: 
    sys.exit("Output directory exists. Exiting to avoid overwriting files.")


for line in f:
    line = line.rstrip('\n')
    print(line)
    print('-------------')
    #directory='./coverage/'+line+'/'
    file_start=args.input+'/'+line+'/'
   
    # Index reference, Burrows-Wheeler Transform
    print('Building Bowtie2 reference index')
    COMMAND='bowtie2-build -f --quiet '+file_start+'contigs.fasta '+file_start+'contigs.fasta_index'
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    # Align Paired end and bam it
    print('Mapping reads using Bowtie2')
    os.popen('bowtie2 -p '+args.threads+' --quiet -x '+file_start+'contigs.fasta_index -1 '+path_to_reads+'/'+line+'/split-adapter-quality-trimmed/'+line+'-READ1.fastq.gz -2 '+path_to_reads+'/'+line+'/split-adapter-quality-trimmed/'+line+'-READ2.fastq.gz -U '+path_to_reads+'/'+line+'/split-adapter-quality-trimmed/'+line+'-READ-singleton.fastq.gz -S '+file_start+'contigs.fasta.sam').read()

    print('Indexing, sorting, and re-indexing the BAM file')
    COMMAND='samtools faidx '+file_start+'contigs.fasta'
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    COMMAND='samtools view -bt '+file_start+'contigs.fasta.fai '+file_start+'contigs.fasta.sam > '+file_start+'contigs.fasta.bam'
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    COMMAND='samtools sort '+file_start+'contigs.fasta.bam '+file_start+'contigs.fasta-s'
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    COMMAND='samtools index '+file_start+'contigs.fasta-s.bam'
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    # Mark duplicates and sort
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-XX:ParallelGCThreads="+args.threads+"",
        "-jar",
        os.path.join(JAR_PATH, "MarkDuplicates.jar"),
        "I="+file_start+"contigs.fasta-s.bam",
        "O="+file_start+"contigs.fasta-smd.bam",
        "METRICS_FILE="+file_start+"contigs.fasta-smd.metrics",
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000",
        "ASSUME_SORTED=TRUE",
        "VALIDATION_STRINGENCY=SILENT",
        "REMOVE_DUPLICATES=TRUE",
        "QUIET=TRUE",
    ]
    print 'Marking Duplicates for '+line
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

    #sort and index the duplicate mark-up
    os.popen('samtools sort '+file_start+'contigs.fasta-smd.bam '+file_start+'contigs.fasta-smds').read()
    os.popen('samtools index '+file_start+'contigs.fasta-smds.bam').read()

    #per base coverage
    print('Calculating coverages for '+line+' now.')
    os.popen('bedtools genomecov -d -ibam  '+file_start+'contigs.fasta-smds.bam > \''+args.output+'/'+line+'-smds.per.base.coverage\'').read()
    
    #per locus coverage
    os.popen('bedtools genomecov -ibam '+file_start+'contigs.fasta-smds.bam > \''+args.output+'/'+line+'-smds.coverage.per.contig\'').read()

    COMMAND = "awk 'BEGIN {pc=\"\"} {c=$1; if (c == pc) {cov=cov+$2*$5;} else {print pc,cov;cov=$2*$5;pc=c}} END {print pc,cov}' "+args.output+"/"+line+"-smds.coverage.per.contig | tail -n +2 > \'"+args.output+"/"+line+"-smds.coverage.per.contig.summary\'"
    subprocess.call(COMMAND, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    #cleaning up the directory, because we don't need it anymore. 
    print('Successfully calculated coverage for '+line+'. Cleaning up after myself, and moving on.')
    
    for CleanUp in glob.glob(file_start+'/*.*'):
        #print CleanUp  
        if not CleanUp.endswith('contigs.fasta'):
            #print CleanUp    
            os.remove(CleanUp)

    print '###########'
