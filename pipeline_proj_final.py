import os
import sys
import argparse
from Bio import SeqIO

#command line arguments parser 
def parse_cli_args(cli_args=None):
    parser = argparse.ArgumentParser(
        description="Wrapper script that runs a pipeline to analyze a virus genome"
    )
    parser.add_argument("-i", "--input", help="input files", nargs="+", required=True)
    return parser.parse_args(cli_args)

#count reads in FASTQ file 
#counts the number of records/read pairs in fastq file 
def fastq_record_count(fq_path):
    recs = 0
    for _ in SeqIO.parse(fq_path, 'fastq'):
        recs += 1
    return recs

#identifies paired end files based on name 
#identify pairs of files that belong together based on shared name
#returns a dictionary where the key is one file and the value is it paired file
def match_paired_files(file_list):
    pairing = {}
    for f1 in file_list:
        for f2 in file_list:
            parts = f1.split('_')
            if f2.startswith(parts[0]) and f1 != f2:
                if f2.endswith('1.fq'):
                    pairing[f2] = f1
                else:
                    pairing[f1] = f2
    return pairing

#run the spades assembler step 
#identifies all fq files, matched foward and reverse read files, executes the spades command, logs the command 
def run_spades_assembler():
    #gather all .fq files from current directory.
    fq_inputs = [f for f in os.listdir(os.getcwd()) if f.endswith('.fq')]
    paired_files = match_paired_files(fq_inputs)
    forwards = list(paired_files.keys())
    #assume exactly two pairs.
    fwd1 = forwards[0]
    rev1 = paired_files[fwd1]
    fwd2 = forwards[1]
    rev2 = paired_files[fwd2]
    #the spades command for assembly 
    spades_cmd = (
        f"spades.py -k 99 -t 2 --only-assembler "
        f"--pe1-1 {fwd1} --pe1-2 {rev1} "
        f"--pe2-1 {fwd2} --pe2-2 {rev2} -o SPADES_assembly/"
    )
    os.system(spades_cmd)
    #log the spades command that was used 
    with open("PipelineProject.log", "a") as log_handle:
        log_handle.write(f"SPAdes command used for assembly: {spades_cmd}\n")

#calculates the assembly statistics on contigs 
#parse the contig files from the spades assembly 
#count the number of contigs longer than 1000bp
#compute the total assembly length and log the results 
def calc_contig_statistics():
    contig_file = './SPADES_assembly/contigs.fasta'
    contig_count = 0
    total_length = 0
    with open(contig_file) as contigs:
        for rec in SeqIO.parse(contigs, 'fasta'):
            if len(rec.seq) > 1000:
                contig_count += 1
                total_length += len(rec.seq)
    with open("PipelineProject.log", "a") as log_handle:
        log_handle.write(f"There are {contig_count} contigs > 1000 bp in the assembly\n")
        log_handle.write(f"There are {total_length} bp in the assembly\n")

#performs blast on the longest contig found 
#identifies the longest contig (length >100bp) from the spades assembly 
#write it to FASTA file and download database of Betaherpesvirinae genome
#build the database and run blast against the longest contig and log results 
def perform_blast():
    seqs_over_thresh = []
    contig_file = './SPADES_assembly/contigs.fasta'
    #find all contigs above the 1000bp threshold
    with open(contig_file) as contigs:
        for rec in SeqIO.parse(contigs, 'fasta'):
            if len(rec.seq) > 1000:
                seqs_over_thresh.append(rec.seq)
    #sort contigs in descending order and get the longest 
    seqs_over_thresh.sort(key=len, reverse=True)
    longest_seq = seqs_over_thresh[0]
    #write the longest contig to file
    with open('longest_contig.fasta', 'w') as out_f:
        out_f.write(str(longest_seq))
    
    #download the genome dataset
    db_download = 'datasets download virus genome taxon betaherpesvirinae --include genome'
    os.system(db_download)
    os.system('unzip ncbi_dataset.zip')
    #path to the downloaded fasta files for the blast database 
    db_path = './ncbi_dataset/data/*genomic.fna'
    #create blast database from the downloaded sequences 
    mkdb_cmd = f"makeblastdb -in {db_path} -out betaherpesvirinae -title beatherpesvirinae -dbtype nucl"
    os.system(mkdb_cmd)
    
    #run blast against the longest contig 
    blast_cmd = (
        'blastn -query longest_contig.fasta -db betaherpesvirinae '
        '-max_target_seqs 10 -max_hsps 1 -out blast_output.tsv '
        '-outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    )
    os.system(blast_cmd)
    #log results 
    with open("PipelineProject.log", "a") as log_handle:
        log_handle.write("sacc   pident  length  qstart  qend    sstart  send    bitscore    evalue  stitle\n")
        with open('blast_output.tsv', 'r') as blast_out:
            for line in blast_out:
                log_handle.write(line + "\n")

#pipeline workflow 
def main():
    #parse command-line arguments.
    options = parse_cli_args(sys.argv[1:])
    file_args = options.input

    #create a project directory and move input files there.
    proj_dir = "PipelineProject_Evelyn_Umana"
    os.system(f"mkdir {proj_dir}")
    for item in os.listdir(os.getcwd()):
        if item in file_args:
            os.rename(os.path.join(os.getcwd(), item), os.path.join(proj_dir, item))
    os.chdir(proj_dir) #change directory to the project folder 
    open("PipelineProject.log", "x").close()  #create the log file

    #group input files 
    donor_group1 = []
    donor_group2 = []
    sample_group = []
    if len(file_args) == 4 and file_args[0].startswith('S'):
        for f in file_args:
            if f.startswith('SRR5660030'):
                donor_group1.append(f)
            else:
                donor_group2.append(f)
    else:
        sample_group = file_args

    #download reference genome data and build Bowtie2 index.
    #download the HCMV dataset
    ref_download = (
        'datasets download genome accession GCF_000845245.1 '
        '--include gff3,rna,cds,protein,genome,seq-report --filename HCMV_dataset.zip'
    )
    os.system(ref_download)
    os.system('unzip HCMV_dataset.zip')
    #path to ref FASTA file
    ref_fasta = 'ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna'
    #build the bowtie 2 index for ref genome 
    os.system(f"bowtie2-build {ref_fasta} HCMV")

    #mapping reads to the reference.
    if file_args[0].startswith('S'):
        #for data with 4 input files 
        pre_donor1 = fastq_record_count(donor_group1[0])
        pre_donor2 = fastq_record_count(donor_group2[0])
        fwd_d1 = ""
        rev_d1 = ""
        #indentify forward and reverse files for donor group 1
        for entry in donor_group1:
            if entry.endswith('1.fastq'):
                fwd_d1 = entry
            else:
                rev_d1 = entry
        map_cmd_d1 = f'bowtie2 -x HCMV -1 {fwd_d1} -2 {rev_d1} --al-conc-gz two_dpi_mapped_%.fq.gz -p 2'
        os.system(map_cmd_d1)
        fwd_d2 = ""
        rev_d2 = ""
        #identify foward and reverse files for group 2
        for entry in donor_group2:
            if entry.endswith('1.fastq'):
                fwd_d2 = entry
            else:
                rev_d2 = entry
        map_cmd_d2 = f'bowtie2 -x HCMV -1 {fwd_d2} -2 {rev_d2} --al-conc-gz six_dpi_mapped_%.fq.gz -p 2'
        os.system(map_cmd_d2)
        #unzip files from bowtie 2 mapping 
        os.system('gunzip *.gz')
        #count reads
        post_donor1 = fastq_record_count('two_dpi_mapped_1.fq')
        post_donor2 = fastq_record_count('six_dpi_mapped_1.fq')
        #log read count stats 
        with open("PipelineProject.log", "a") as log_handle:
            log_handle.write(
                f"Donor 1 (2dpi) had {pre_donor1} read pairs before Bowtie2 filtering and {post_donor1} after.\n"
            )
            log_handle.write(
                f"Donor 1 (6dpi) had {pre_donor2} read pairs before Bowtie2 filtering and {post_donor2} after.\n"
            )
    else:
        #for non donor grouping 
        samp_fw1 = ""
        samp_rv1 = ""
        samp_fw2 = ""
        samp_rv2 = ""
        for s in sample_group:
            if s.startswith('sampledata_'):
                if s.endswith('1.fastq'):
                    samp_fw1 = s
                else:
                    samp_rv1 = s
            else:
                if s.endswith('1.fastq'):
                    samp_fw2 = s
                else:
                    samp_rv2 = s
        pre_samp1 = fastq_record_count(samp_fw1)
        pre_samp2 = fastq_record_count(samp_fw2)
        map_cmd_s1 = f'bowtie2 -x HCMV -1 {samp_fw1} -2 {samp_rv1} --al-conc-gz sample_1_mapped_%.fq.gz -p 2'
        os.system(map_cmd_s1)
        map_cmd_s2 = f'bowtie2 -x HCMV -1 {samp_fw2} -2 {samp_rv2} --al-conc-gz sample2_mapped_%.fq.gz -p 2'
        os.system(map_cmd_s2)
        os.system('gunzip *gz')
        post_samp1 = fastq_record_count('sample_1_mapped_1.fq')
        post_samp2 = fastq_record_count('sample2_mapped_1.fq')
        #log the read count stats 
        with open("PipelineProject.log", "a") as log_handle:
            log_handle.write(
                f"Donor 1 (2dpi) had {pre_samp1} read pairs before Bowtie2 filtering and {post_samp1} after.\n"
            )
            log_handle.write(
                f"Donor 1 (6dpi) had {pre_samp2} read pairs before Bowtie2 filtering and {post_samp2} after.\n"
            )

    #clean up downloaded reference files.
    os.system('rm -r ncbi_dataset')
    os.system('rm README.md')
    os.system('rm md5sum.txt')

    #run assembly with SPAdes.
    run_spades_assembler()

    #filter contigs and report assembly stats.
    calc_contig_statistics()

    #BLAST analysis using the longest contig.
    perform_blast()

if __name__ == '__main__':
    main()
