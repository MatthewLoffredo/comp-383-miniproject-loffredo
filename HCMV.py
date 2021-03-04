import os
import sys
from Bio import SeqIO
from Bio import Entrez
from utils import *
import argparse

######## Arg Parsing #########
# parser to check for flag to use test data
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--test', required = False, dest='test', action='store_true')
args = parser.parse_args()
if args.test:
  print("Using test data.")
else: 
  print("Using full data.")

######## Setup #########
# create directory for files if not already created
directory = os.popen('pwd').read().rstrip()
path = (directory + '/miniProject_Matt_Loffredo/')
if not os.path.exists(path):
  os.system('mkdir ' + path)
  os.chdir(path)
  os.system('touch miniProject.log')
else:
  os.chdir(path)
  
# Create list of sample numbers
sample_names = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

# Get reads from sra numbers if not using test data
if not args.test:
  os.system("cp ../getReads.sh getReads.sh")
  os.system('sh getReads.sh')
else:
  # Create subset of reads if not already created
  if not os.path.exists('SRR5660030._1.sub.fastq'):
    for sample in sample_names:
      os.system("cp ../%s_1.sub.fastq %s_1.sub.fastq" % (sample, sample))
      os.system("cp ../%s_2.sub.fastq %s_2.sub.fastq" % (sample, sample))
      #os.system('head -n 10000 ../%s_1.fastq > %s.sub_1.fastq' % (sample, sample))
      #os.system('head -n 10000 ../%s_2.fastq > %s.sub_2.fastq' % (sample, sample))

"""    
os.system('head -n 10000 ../SRR5660030_1.fastq > SRR5660030_1.sub.fastq')
os.system('head -n 10000 ../SRR5660030_2.fastq > SRR5660030_2.sub.fastq')
os.system('head -n 10000 ../SRR5660033_1.fastq > SRR5660033_1.sub.fastq')
os.system('head -n 10000 ../SRR5660033_2.fastq > SRR5660033_2.sub.fastq')
os.system('head -n 10000 ../SRR5660044_1.fastq > SRR5660044_1.sub.fastq')
os.system('head -n 10000 ../SRR5660044_2.fastq > SRR5660044_2.sub.fastq')
os.system('head -n 10000 ../SRR5660045_1.fastq > SRR5660045_1.sub.fastq')
os.system('head -n 10000 ../SRR5660045_2.fastq > SRR5660045_2.sub.fastq')
"""

######## Create ref genome from cds #########
# retrieve ref genome cds
Entrez.email = "mloffredo@luc.edu"
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb")
record = SeqIO.read(handle, "genbank")
#record.features = [f for f in record.features if f.type == "CDS"]
cds_num = 0
with open("HCMV.fa", "w") as ref_cds:
  for feature in record.features:
    if feature.type == "CDS":
      cds_num += 1
      ref_cds.write(">%s.%s CDS %s\n%s\n" % (
          record.id,
          cds_num,
          cds_num,
          feature.extract(record.seq)))
  ref_cds.close()

# Open log file and write to it
log = open("miniProject.log", "w")
log.write("The HCMV genome (EF99921) has " + str(cds_num) + " CDS." + "\n")
log.close()

######## Kallisto #########
# build index
os.system("time kallisto index -i index.idx HCMV.fa")

# Quantify TPM of cds
os.system("mkdir kallisto_output")
# use script for test/full data
if not args.test:
  os.system("cp ../runQuants.sh runQuants.sh")
  os.system("sh runQuants.sh")
else:
  os.system("cp ../runSubQuants.sh runSubQuants.sh")
  os.system("sh runSubQuants.sh")


######## Sleuth #########
# Make sample table
os.system("touch sample_info.txt")
info_table = open("sample_info.txt", "w")

info_table.write("sample condition path\n")
info_table.write("SRR5660030.1 2dpi kallisto_output/SRR5660030.1\n")
info_table.write("SRR5660044.1 2dpi kallisto_output/SRR5660044.1\n")
info_table.write("SRR5660033.1 6dpi kallisto_output/SRR5660033.1\n")
info_table.write("SRR5660045.1 6dpi kallisto_output/SRR5660045.1\n")
info_table.close()

# run sleuth script
os.system("cp ../sleuth_script.R sleuth_script.R")
os.system("Rscript ./sleuth_script.R")

#quant_cmd = "time kallisto quant -i index.idx -b 30 -t 2 -o kallisto_output/"
#os.system(quant_cmd + "SRR5660030.1 SRR5660030_1.sub.fastq SRR5660030_2.sub.fastq")

######## Bowtie Mapping #########
# Retreive HCMV Genome
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="fasta")
record = SeqIO.read(handle, "fasta")
with open("HCMV_genome.fna", "w") as output_handle:
    SeqIO.write(record, output_handle, "fasta")
    output_handle.close()

# create index for genome
bowtie_build("../HCMV_genome.fna", "HCMV")

# align sequences to index
bowtie_align("HCMV","SRR5660030", "SRR5660030.1", args.test)
bowtie_align("HCMV","SRR5660033", "SRR5660033.1", args.test)
bowtie_align("HCMV","SRR5660044", "SRR5660044.1", args.test)
bowtie_align("HCMV","SRR5660045", "SRR5660045.1", args.test)

######## Spades Assembly #########
# make list of sample names
aligned_samples = ["SRR5660030.1", "SRR5660033.1", "SRR5660044.1", "SRR5660045.1"]

print("## Assembling "+str(aligned_samples)+" using spades")
# create spade command arguments for each sample
spade_inputs = ""
arg_count = 0
map_dir = "HCMV_mapping/"
for sample in aligned_samples:
  arg_count += 1
  spade_inputs += "--pe"+str(arg_count)+"-1 "+map_dir+sample+".mapped.1.fastq --pe"+str(arg_count)+"-2 "+map_dir+sample+".mapped.2.fastq "

# run assembly
assembly_name = "HCMV_assembly"
# output_exists = os.path.isdir(assembly_name)
# if not output_exists:
os.system("mkdir "+assembly_name)
spades_command = "spades -k 55,77,99,127 -t 2 --only-assembler "+spade_inputs+"-o "+assembly_name+"/"
print(spades_command)
os.system(spades_command)
# Open log file and write to it
log = open("miniProject.log", "a+")
log.write(spades_command+"\n")
log.close()
# else: print("Files have already been assembled: "+str(aligned_samples))

######## Contig Filtering #########
num_contigs, bp = filter_contigs(assembly_name, 1000)
log = open("miniProject.log", "a+")
log.write("There are "+str(num_contigs)+" contigs > 1000 in the assembly.\n")
log.write("There are "+str(bp)+" bp in the assembly.\n")
log.close()

######## BLAST+ #########
# this was for making the local blast db which is now stored in the repo
# input_file = 'Betahers.fasta'
# makeblast_command = "makeblastdb -in "+input_file+" -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl"
# os.system(makeblast_command)

# copy local blast db to directory
os.system("cp ../Betaherpesvirinae.nhr Betaherpesvirinae.nhr")
os.system("cp ../Betaherpesvirinae.nin Betaherpesvirinae.nin")
os.system("cp ../Betaherpesvirinae.nsq Betaherpesvirinae.nsq")

# query local blast db
get_longest_contig(assembly_name)
input_file = 'longest_contig.fasta'
output_file = "blast_results.csv"
headers = ["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"]
headers_str = " ".join(headers)
blastCMD = 'blastn -query '+input_file+' -db Betaherpesvirinae -out '+output_file+' -outfmt "10 '+headers_str+'"'
print(blastCMD)
os.system(blastCMD)
rows = parse_blast(output_file, headers)

log = open("miniProject.log", "a+")
# writes headers to log
log.write("\t".join(headers)+"\n")
for row in rows:
  # remove values that have a type of None, which occurs when stitle has commas in it
  values = list(row.values())
  for value in values:
    if type(value) == list:
      values.remove(value)
  # writes values to log
  log.write("\t".join(values)+"\n")
log.close()

sys.exit(0)
