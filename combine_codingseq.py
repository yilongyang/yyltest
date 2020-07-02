import os
import re


aj="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/AJ-codingseq.fasta"
F1D="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/F1D.scafSeq-augustus.codingseq"
FRA202="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/FRA202.a.lines.augustus.codingseq"
FRA1358="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/FRA1358-a-lines.augustus.codingseq"
fni="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/Fragaria_nipponica_Genome_v1.0_FNI_r1.1_cds.fasta"
fnu="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/Fragaria_nubicola_Genome_v1.0_FNU_r1.1_cds.fasta"
fvesca="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/Fragaria_vesca_v4.0.a1_makerStandard_CDS.fasta"
fxa="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/Fxa_v1.2_makerStandard_transcripts_woTposases.fasta"
gs91="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/GS91.a.lines.fvesca-augustus.codingseq"
potentilla="/mnt/home/davis/ytj2/assemblies/orthofinder/codingseq/Potentilla_transcripts.fasta"

#infile_list=[aj, F1D, FRA1358, FRA202, fni, fnu, fvesca, fxa, gs91, potentilla]

def fasta2dict(infile, str):
    print infile
    outfile="all_codingseq.fasta"
    out=open(outfile, 'a')
    f = open(infile)
    lines = f.readlines()
    infilestr=''.join(lines)
    seq_list=infilestr.split(">")    
    seq_list.pop(0)
    for sequence in seq_list:
        one_seq_list=sequence.split("\n")
        header=one_seq_list.pop(0)
        seq=''.join(one_seq_list)
        find_space=re.search(r'(.*?)\s', header)
        if find_space:
            newheader=find_space.group(1)
        else:
            newheader=header
        newheader=newheader+str
        out.write(">"+newheader+"\n")
        out.write(seq+"\n")
    out.close()
    f.close()
    

fasta2dict(aj, '')
fasta2dict(F1D, "-F1D")
fasta2dict(FRA202, "-Fchinensis")
fasta2dict(FRA1358, '-Fnilgerrensis')
fasta2dict(fni, '')
fasta2dict(fnu, '')
fasta2dict(fvesca, '')
fasta2dict(fxa, '')
fasta2dict(gs91, '-Fviridis')
fasta2dict(potentilla, '')