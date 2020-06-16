import argparse
import os
import re
import glob
import sys

def usage():
    argp=argparse.ArgumentParser(description='auto analysis pipeline', formatter_class=argparse.RawTextHelpFormatter)
    argp.add_argument('-indir', dest='in_dir', help='path to the input raw data directory')
    argp.add_argument('-outdir', dest='out_dir', help='path to the output directory')
    args=argp.parse_args()
    
    if args.in_dir and args.out_dir:
        r1_dict={}
        r2_dict={}
        indir=os.path.join(os.path.abspath(args.in_dir), '')
        file_path=indir+'*.gz'
        files=glob.glob(file_path, recursive=True)
        for fq_file in files:
            find_sample1=re.search(r'.+\/(.*?)_R1.fastq.gz', fq_file)
            find_sample2=re.search(r'.+\/(.*?)_R2.fastq.gz', fq_file)
            print (fq_file)
            if find_sample1:
                sample_name=find_sample1.group(1)
                r1_dict[sample_name]=fq_file
            elif find_sample2:
                sample_name=find_sample2.group(1)
                r2_dict[sample_name]=fq_file
        print (r1_dict)
        if len(r1_dict)>0:
            outdir=os.path.join(os.path.abspath(args.out_dir), '')
            trimmed_dir=outdir+"trimmed/"
            mkdir_cmd="mkdir -p "+trimmed_dir
            print (mkdir_cmd)
            os.system(mkdir_cmd)
            
            #bam_dir=outdir+"bam/"
            #mkdir_cmd="mkdir -p "+bam_dir
            #os.system(mkdir_cmd)
            
            #down_dir=outdir+"downsampled/"
            #mkdir_cmd="mkdir -p "+down_dir
            #os.system(mkdir_cmd)
            stats_dir=outdir+"stats/"
            mkdir_cmd="mkdir -p "+stats_dir
            os.system(mkdir_cmd)
            #outsh=outdir+"run_"+sample+".sh"
            outsh=outdir+"soapdenovo.sh"
            out=open(outsh, 'w')
            for sample, r1_file in r1_dict.items():
                
                echo_cmd='echo " process sample '+sample+'"\n'
                out.write(echo_cmd)
                r2_file=r2_dict[sample]
                r1_pe=trimmed_dir+sample+"-PE_R1.fastq.gz"
                r2_pe=trimmed_dir+sample+"-PE_R2.fastq.gz"
                r1_se=trimmed_dir+sample+"-SE_R1.fastq.gz"
                r2_se=trimmed_dir+sample+"-SE_R2.fastq.gz"
                #trim_cmd=" java -jar -Xmx150g /home/thinkmate/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar PE "+r1_file+" "+r2_file+" "+r1_pe+" "+r1_se+" "+r2_pe+" "+r2_se+" ILLUMINACLIP:/home/thinkmate/Downloads/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
                trim_cmd=" java -jar -Xmx200g /home/thinkmate/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 "+r1_file+" "+r2_file+" "+r1_pe+" "+r1_se+" "+r2_pe+" "+r2_se+" ILLUMINACLIP:/home/thinkmate/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
                #outsam=bam_dir+sample+".sam"
               # out.write(trim_cmd+"\n")
                #r1_down=down_dir+sample+"-down_R1.fastq.gz"
                #r2_down=down_dir+sample+"-down_R1.fastq.gz"
                #down_cmd="/home/thinkmate/Yilong/tools/bbmap/reformat.sh in="+r1_pe+"  in2="+r2_pe+" out="+r1_down+" out2="r2_down+" reads=1000000 samplereadstarget=1000000 "
                #bwa_cmd="bwa mem -t 24 /home/thinkmate/Yilong/Published_assembly/Fragaria_vesca_v4.0/Fragaria_vesca_v4.0.a1.normalized.fasta "+ r1_pe+" "+r2_pe+ " > "+outsam
                #out.write(bwa_cmd+"\n")
                #samtools_fixmate_cmd="samtools fixmate -O bam "+outsam+" "+outsam+".bam"
                #out.write(samtools_fixmate_cmd+"\n")
                #samtools_sort_cmd="samtools sort -O bam -o "+outsam+"_sorted.bam "+outsam+".bam"
                #out.write(samtools_sort_cmd+"\n")
                merge_txt=stats_dir+sample+"-ihist_merge.txt"
                bbmerge_cmd="/home/thinkmate/Yilong/tools/bbmap/bbmerge.sh  in1="+r1_pe+"  in2="+r2_pe+" ihist="+merge_txt+" reads=400000"
                #out.write(bbmerge_cmd+"\n")
                assembly_dir=outdir+"assembly/"+sample+"/"
                mkdir_cmd="mkdir -p "+assembly_dir
                os.system(mkdir_cmd)
                insert=insertsize(merge_txt)
                cfg_file=assembly_dir+sample+"-lib.cfg"
                sample_cfg(r1_pe, r2_pe, insert, cfg_file)
                assembly_out=assembly_dir+sample
                log=assembly_dir+sample+"-assembly.log"
                out.write("cd "+assembly_dir+"\n")
                assembly_cmd="soapdenovo2-127mer all -s "+cfg_file+"  -K 61 -p 22 -o "+assembly_out+" >> "+log
                out.write(assembly_cmd+"\n")
            out.close()
            os.system("chmod +x "+outsh)

def insertsize(bbmerge_out):
    insert=200
    with open (bbmerge_out, 'r') as f:
        for line in f:
            line=line.strip()
            find_insert=re.search(r'#Mean\t(.*?)$', line)
            if find_insert:
                insert=find_insert.group(1)
                break
    f.close()
    return insert





def sample_cfg(r1_file, r2_file,insert, outfile):
    out=open(outfile, 'w')
    example_cfg="/home/thinkmate/Yilong/assembly_analysis/soapdenovo.cfg"
    with open(example_cfg, 'r') as f:
        for line in f:
            line=line.strip()
            if re.search(r'^q1=', line):
                newline="q1="+r1_file
            elif  re.search(r'^q2=', line):
                newline="q2="+r2_file
            elif re.search(r'avg_ins=', line):
                newline="avg_ins="+str(insert)
            else:
                newline=line
            out.write(newline+"\n")

    out.close()

        
        

if __name__ == "__main__":
    usage()    
