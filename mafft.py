import subprocess as sub
import os
import sys

# do pairwise alignment using mafft
import os
import sys
if (os.path.exists("/project/exaptation/")):
    project_path = "/project/exaptation/"
elif (os.path.exists("/home/kalaghat/exaptation/")):
    project_path = "/home/kalaghat/exaptation/"
sys.path.insert(0,project_path+'Projects/MSTBasedForests/scripts/')
toolPath = project_path+"Tools/"

def ReadFasta(fileName):    
    fastaFile = open(fileName,'r')
    seq=''
    name=''
    sequenceAlignment={}
    for line in fastaFile:
        if line.startswith('>'):
            if seq != '':
                seq = seq.upper()
                sequenceAlignment[name] = seq
                seq = ''
            name = line.strip().split('>')[1]
        else:
            seq += line.strip()

    sequenceAlignment[name] = seq
    fastaFile.close()
    return (sequenceAlignment)

def ReadPhylip(fileName):
    phylipFile = open(fileName,'r')
    phylipFile.readline()
    sequenceAlignment={}
    for line in phylipFile:
        splitLine = line.strip().split(' ')
        if len(splitLine)>1:
            name = splitLine[0]
            seq = splitLine[len(splitLine)-1]
        else:
            name, seq = line.strip().split('\t')
        sequenceAlignment[name] = seq
    phylipFile.close()
    return (sequenceAlignment)



def ReadAlignment(fileName):
    alignmentFile = open(fileName,'r')
    firstLine = alignmentFile.readline()
    alignmentFile.close()
    if firstLine.startswith('>'):
        alignment = ReadFasta(fileName)
    else:
        alignment = ReadPhylip(fileName)
    return alignment
    

def WriteAlignment(alignment,fileName,fileFormat="fasta",convertToAlpha=False):
    alignmentFile = open(fileName,'w')
    if fileFormat == "fasta":
        nucList = ['A','C','T','G']
        if convertToAlpha:
            for seqId in sorted(alignment.keys()):
                alignmentFile.write('>'+str(seqId)+'\n')
                for char in alignment[seqId]:
                    alignmentFile.write(nucList[int(char)])
                alignmentFile.write('\n')
        else:
            for seqId in sorted(alignment.keys()):
                alignmentFile.write('>'+str(seqId)+'\n')
                alignmentFile.write(str(alignment[seqId])+'\n')
    elif fileFormat=="phylip":
        numberOfSequences = len(alignment.keys())
        sequenceLength = len(alignment.values()[0])
        alignmentFile.write(str(numberOfSequences)+'\t'+str(sequenceLength)+'\n')
        if convertToAlpha:
            nucList = ['A','C','T','G']
            for seqId in sorted(alignment.keys()):
                alignmentFile.write(str(seqId)+'\t')
                for char in alignment[seqId]:
                    alignmentFile.write(nucList[int(char)])
                alignmentFile.write('\n')
        else:            
            for sequenceName in alignment.keys():
                alignmentFile.write(sequenceName+'\t'+alignment[sequenceName]+'\n')
    alignmentFile.close()

referenceSequence = ReadAlignment(project_path+"covid_loop/EPI_ISL_402124_uppercase.fasta")

# print(referenceSequence)


unaliged_fasta_file_name = project_path + "covid_loop/10covi_19-22_july.fasta"
raw_alignment = ReadAlignment(unaliged_fasta_file_name)
# align sequences in batches of 200 sequences
noOfSeqsInBatch = 0
smallAlignment = {}
numberOfParts = 0
seqsRemaining = len(raw_alignment)
for seqId, seq in raw_alignment.items():
    #smallAlignment[seqId] = seq 
    noOfSeqsInBatch += 1
    if noOfSeqsInBatch == 1000:
        numberOfParts+= 1
        noOfSeqsInBatch = 0
        #smallAlignment.update(referenceSequence)
        #WriteAlignment(smallAlignment, project_path+'covid_loop/covid_unaligned/part_'+str(numberOfParts)+'.fasta')
        #smallAlignment = {}
        #seqsRemaining -= 1000
        #print ('seqs remaining is', seqsRemaining)

numberOfParts += 1
#smallAlignment.update(referenceSequence)
#WriteAlignment(smallAlignment, project_path+'covid_loop/covid_unaligned/part_'+str(numberOfParts)+'.fasta')


# Run mafft in batch
#toolPath = "/project/exaptation/Tools/"
#mafftFile = open(project_path+'covid_loop/scripts/batch_mafft.sh','w')

for part in range(1,numberOfParts+1):
    unaligned_sequences_batch_fileName = '/project/exaptation/covid_loop/covid_unaligned/part_'+str(part)+'.fasta'
    aligned_sequences_batch_fileName = '/project/exaptation/covid_loop/covid_aligned/part_'+str(part)+'.fasta'
    #mafftCommand = toolPath+'mafft --maxiterate 1000 --thread 5 ' + unaligned_sequences_batch_fileName + ' > ' + aligned_sequences_batch_fileName 
    #mafftFile.write(mafftCommand+'\n')

#mafftFile.close()

#command_for_mafft = sub.call(project_path+'covid_loop/scripts/batch_mafft.sh', shell= True)
# command_to_convert = sub.call("python3 " + diff_tool + "--reference H3N2_HA_ref_trimmed.fasta --fasta H3N2_156.fasta --output H3N2_156.diff", shell= True)
#print(command_for_mafft, "The dataset ran on mafft tool")

# process mafft alignments and create concatenated alignment
referenceSequence = ReadAlignment(project_path+"covid_loop/EPI_ISL_402124_uppercase.fasta")
# align sequences in batches of 200 sequences
multipleSequenceAlignment = {}

concatenated_alignment_fileName = '/project/exaptation/covid_loop/covid_aligned/part_'+str(part)+'.fasta'
for part in range(1,numberOfParts+1):
    aligned_sequences_batch_fileName = '/project/exaptation/covid_loop/covid_aligned/part_'+str(part)+'.fasta'
    alignedSequences = ReadAlignment(aligned_sequences_batch_fileName)
    print (len(alignedSequences))
    noOfPosInAlignment = len(list(alignedSequences.values())[0])
    posToKeep = [True]*noOfPosInAlignment
    alignedRefSeq = alignedSequences[list(referenceSequence.keys())[0]]
    for pos in range(noOfPosInAlignment):
        if alignedRefSeq[pos]=='-':
            posToKeep[pos] = False
    for seqid in alignedSequences.keys():
        if seqid != list(referenceSequence.keys())[0]:
            trimmedSeq = ""
            for pos in range(noOfPosInAlignment):
                if posToKeep[pos]:
                    trimmedSeq += alignedSequences[seqid][pos]
            multipleSequenceAlignment[seqid] = trimmedSeq.upper()
    print ('part', part, 'of', numberOfParts)
WriteAlignment(multipleSequenceAlignment, project_path+'covid_loop/aligned_10covi_19-22_july.fasta', 'fasta')