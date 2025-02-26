import os
import sys
import argparse
from Bio import SeqIO

def parse_cli_args(cli_args=None):
    parser = argparse.ArgumentParser(
        description="Wrapper script that runs a pipeline to analyze a virus genome"
    )
    parser.add_argument("-i", "--input", help="input files", nargs="+", required=True)
    return parser.parse_args(cli_args)

def fastq_record_count(fq_path):
    recs = 0
    for _ in SeqIO.parse(fq_path, 'fastq'):
        recs += 1
    return recs

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

def run_spades_assembler():
    # Gather all .fq files from current directory.
    fq_inputs = [f for f in os.listdir(os.getcwd()) if f.endswith('.fq')]
    paired_files = match_paired_files(fq_inputs)
    forwards = list(paired_files.keys())
    # Assume exactly two pairs.
    fwd1 = forwards[0]
    rev1 = paired_files[fwd1]
    fwd2 = forwards[1]
    rev2 = paired_files[fwd2]
    spades_cmd = (
        f"spades.py -k 99 -t 2 --only-assembler "
        f"--pe1-1 {fwd1} --pe1-2 {rev1} "
        f"--pe2-1 {fwd2} --pe2-2 {rev2} -o SPADES_assembly/"
    )
    os.system(spades_cmd)
    with open("PipelineProject.log", "a") as log_handle:
        log_handle.write(f"SPAdes command used for assembly: {spades_cmd}\n")

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

def perform_blast():
    seqs_over_thresh = []
    contig_file = './SPADES_assembly/contigs.fasta'
    with open(contig_file) as contigs:
        for rec in SeqIO.parse(contigs, 'fasta'):
            if len(rec.seq) > 1000:
                seqs_over_thresh.append(rec.seq)
    # Sort in descending order of length.
    seqs_over_thresh.sort(key=len, reverse=True)
    longest_seq = seqs_over_thresh[0]
    with open('longest_contig.fasta', 'w') as out_f:
        out_f.write(str(longest_seq))
    
    # Download and set up the BLAST database.
    db_download = 'datasets download virus genome taxon betaherpesvirinae --include genome'
    os.system(db_download)
    os.system('unzip ncbi_dataset.zip')
    db_path = './ncbi_dataset/data/*genomic.fna'
    mkdb_cmd = f"makeblastdb -in {db_path} -out betaherpesvirinae -title beatherpesvirinae -dbtype nucl"
    os.system(mkdb_cmd)
    
    blast_cmd = (
        'blastn -query longest_contig.fasta -db betaherpesvirinae '
        '-max_target_seqs 10 -max_hsps 1 -out blast_output.tsv '
        '-outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    )
    os.system(blast_cmd)
    with open("PipelineProject.log", "a") as log_handle:
        log_handle.write("sacc   pident  length  qstart  qend    sstart  send    bitscore    evalue  stitle\n")
        with open('blast_output.tsv', 'r') as blast_out:
            for line in blast_out:
                log_handle.write(line + "\n")

def main():
    # Parse command-line arguments.
    options = parse_cli_args(sys.argv[1:])
    file_args = options.input

    # Create a dedicated project directory and move input files there.
    proj_dir = "PipelineProject_Evelyn_Umana"
    os.system(f"mkdir {proj_dir}")
    for item in os.listdir(os.getcwd()):
        if item in file_args:
            os.rename(os.path.join(os.getcwd(), item), os.path.join(proj_dir, item))
    os.chdir(proj_dir)
    open("PipelineProject.log", "x").close()  # Create the log file

    # Determine file grouping based on the provided input names.
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

    # Step 2: Download reference genome data and build Bowtie2 index.
    ref_download = (
        'datasets download genome accession GCF_000845245.1 '
        '--include gff3,rna,cds,protein,genome,seq-report --filename HCMV_dataset.zip'
    )
    os.system(ref_download)
    os.system('unzip HCMV_dataset.zip')
    ref_fasta = 'ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna'
    os.system(f"bowtie2-build {ref_fasta} HCMV")

    # Mapping reads to the reference.
    if file_args[0].startswith('S'):
        pre_donor1 = fastq_record_count(donor_group1[0])
        pre_donor2 = fastq_record_count(donor_group2[0])
        fwd_d1 = ""
        rev_d1 = ""
        for entry in donor_group1:
            if entry.endswith('1.fastq'):
                fwd_d1 = entry
            else:
                rev_d1 = entry
        map_cmd_d1 = f'bowtie2 -x HCMV -1 {fwd_d1} -2 {rev_d1} --al-conc-gz two_dpi_mapped_%.fq.gz -p 2'
        os.system(map_cmd_d1)
        fwd_d2 = ""
        rev_d2 = ""
        for entry in donor_group2:
            if entry.endswith('1.fastq'):
                fwd_d2 = entry
            else:
                rev_d2 = entry
        map_cmd_d2 = f'bowtie2 -x HCMV -1 {fwd_d2} -2 {rev_d2} --al-conc-gz six_dpi_mapped_%.fq.gz -p 2'
        os.system(map_cmd_d2)
        os.system('gunzip *.gz')
        post_donor1 = fastq_record_count('two_dpi_mapped_1.fq')
        post_donor2 = fastq_record_count('six_dpi_mapped_1.fq')
        with open("PipelineProject.log", "a") as log_handle:
            log_handle.write(
                f"Donor 1 (2dpi) had {pre_donor1} read pairs before Bowtie2 filtering and {post_donor1} after.\n"
            )
            log_handle.write(
                f"Donor 1 (6dpi) had {pre_donor2} read pairs before Bowtie2 filtering and {post_donor2} after.\n"
            )
    else:
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
        with open("PipelineProject.log", "a") as log_handle:
            log_handle.write(
                f"The first sample input had {pre_samp1} read pairs before Bowtie2 filtering and {post_samp1} after.\n"
            )
            log_handle.write(
                f"The second sample input had {pre_samp2} read pairs before Bowtie2 filtering and {post_samp2} after.\n"
            )

    # Clean up downloaded reference files.
    os.system('rm -r ncbi_dataset')
    os.system('rm README.md')
    os.system('rm md5sum.txt')

    # Step 3: Run assembly with SPAdes.
    run_spades_assembler()

    # Step 4: Filter contigs and report assembly stats.
    calc_contig_statistics()

    # Step 5: BLAST analysis using the longest contig.
    perform_blast()

if __name__ == '__main__':
    main()
