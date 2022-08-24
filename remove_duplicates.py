from Bio import SeqIO
import time

start = time.time() 

seen = []
records = []

covid_fasta_input = "/project/exaptation/covid_loop/aligned_10covi_19-22_july.fasta"
dupli_covid_fasta_input = "/project/exaptation/covid_loop/no_dupli_aligned_10covi_19-22_july.fasta"

for record in SeqIO.parse(covid_fasta_input, "fasta"):  
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)


#writing to a fasta file
SeqIO.write(records, dupli_covid_fasta_input, "fasta")
end = time.time()

print(f"Run time is {(end- start)/60}")