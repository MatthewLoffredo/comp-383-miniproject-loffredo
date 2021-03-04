import os
from itertools import (takewhile, repeat)
from Bio import SeqIO
import csv

# builds bowtie index
def bowtie_build(input, out):
  # make assembly folder and go there
  path = "HCMV_mapping"
  if not os.path.exists(path):
    os.system("mkdir " + path)
  os.chdir(path)
  
  print("Building index for "+input+"...")
  bowtie_command = 'bowtie2'
  input_exists = os.path.isfile(input)
  if input_exists:
    bowtie_command += '-build --quiet -f '+input+' '+out
    print(bowtie_command)
  else:
    print("Error: could not locate "+input+"")
    return False
  os.system(bowtie_command)
  os.chdir("..")

# runs bowtie alignment
def bowtie_align(idx, input, out, use_test):
  os.chdir("./HCMV_mapping")
  
  # use test data accordingly
  if use_test:
    r1, r2 = "../"+input+"_1.sub.fastq", "../"+input+"_2.sub.fastq"
  else:
    r1, r2 = "../"+input+"_1.fastq", "../"+input+"_2.fastq"
    
  input_exists = os.path.isfile(r1) and os.path.isfile(r2)
  index_exists = os.path.isfile(idx+".1.bt2")
  # if no mapped fastq file exists, map input to index
  #if not os.path.isfile(out+".mapped.1.fastq"):
  print("Aligning "+input+" to "+idx+" with Bowtie2")
  if input_exists and index_exists:
    bowtie_command = "bowtie2 --al-conc "+out+".mapped.fastq -x "+idx+" -1 "+r1+" -2 "+r2+" -S "+out; print(bowtie_command)
  elif index_exists: print("Error: could not locate paired FASTQ files: "+r1+","+r2); return False
  elif input_exists: print("Error: could not locate index file: "+idx); return False
  else: print("Error: could not locate any files provided."); return False
  os.system(bowtie_command)
  print("Done!")
  #else: print(input+" has already been mapped to "+idx+" with Bowtie2. Skipping!")
  
  # compute # of read pairs before and after mapping
  before = str(line_count(r1))
  after = str(line_count(out+".mapped.1.fastq"))
  
  # Write to log file
  log = open("../miniProject.log", "a+")
  log_prefix = ""
  if input == "SRR5660030":
    log_prefix = "Donor 1 (2dpi)"
  if input == "SRR5660033":
    log_prefix = "Donor 1 (6dpi)"
  if input == "SRR5660044":
    log_prefix = "Donor 2 (2dpi)"
  if input == "SRR5660045":
    log_prefix = "Donor 2 (6dpi)"
  
  log_message = log_prefix+" had "+before+" read pairs before Bowtie2 filtering and "+after+" read pairs after."+"\n"
  print(log_message)
  log.write(log_message)
  log.close()
  os.chdir("..")
  
# gets number of read pairs in file
def line_count(fi):
  f = open(fi, 'rb')
  bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
  seq_count = int(sum( buf.count(b'\n') for buf in bufgen )/4)
  return seq_count

# filters for contigs with lengths greater than 1000
def filter_contigs(assembly, n):
  print("## Filtering contigs in "+assembly+" by size "+str(n))
  in_path = assembly+"/contigs.fasta"
  outfile = open("contigs."+str(n)+".fasta", "w+")
  num_contigs, size = 0,0
  # extract length from record ids in fasta
  for record in SeqIO.parse(in_path, "fasta"):
    # find index of length_
    s = record.id.find("length_")
    # go past that to next _ to get length num
    l = record.id[s+7:]; l = int(l[:l.find("_")])
    if l > n:
      outfile.write(">"+record.id+"\n"+str(record.seq)+"\n")
      num_contigs += 1
      size += l
  return num_contigs, size

# gets longest contig from bowtie assembly
def get_longest_contig(assembly):
  print("## Retreiving longest contig in "+assembly)
  in_path = assembly+"/contigs.fasta"
  outfile = open("longest_contig.fasta", "w+")
  size = 0
  longest_record = None
  # extract length from record ids in fasta
  for record in SeqIO.parse(in_path, "fasta"):
    # find index of length_
    s = record.id.find("length_")
    # go past that to next _ to get length num
    l = record.id[s+7:]; l = int(l[:l.find("_")])
    if l > size:
      size = l
      longest_record = record
  outfile.write(">"+longest_record.id+"\n"+str(longest_record.seq)+"\n")

# parse blast input
def parse_blast(filename, headers):
  x=[]
  blast_results=open(filename,'r')
  rows=csv.DictReader(blast_results,headers,delimiter=',')
  for row in rows:
    x.append(row)
  blast_results.close()
  return x
