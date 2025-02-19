#dependencies
#enviroments prepare
import pandas as pd;import numpy as np;import os;import glob;import re;from datetime import datetime
import pickle
from Bio import AlignIO
from Bio import SeqIO
from Bio import Align

#read in the 11k reads of the bam input,make a subdir in the target_dir and 
#depends on My_Perl from wl
def from_bam_to_reads(name,bam_file,target_dir):
    #first read the reads covering 11k and then output them to two fq files in the target_dir
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    a=os.popen(f'echo "Supercontig_12.6 309013 320543" |  extract_bam_pairs.pl --input - --format fastq  --bam {bam_file}  --output {target_dir}/{name}')



#use abyss to aseemble the reads, final result should be a xxx_scaffold.fa
def from_fq_to_denovo(name,target_dir):
    #use abyss to aseemble reads into scafolds fasta file
    os.chdir(target_dir)
    os.popen(f"abyss-pe name={name} -C {target_dir} K=55 k=95 in='{target_dir}/{name}_1.fq {target_dir}/{name}_2.fq'")
    
#重构代码，每次输入一条scaffold，align，to frame，输出序列之后结束单次流程
def mutate_loci(sites,new_bases,old_seq):
    out_seq=old_seq
    for site,new_base in zip(sites,new_bases):
        out_seq=out_seq[:site]+new_base+out_seq[site+1:]
    return out_seq
def scaffold_align_to_frame(sample,AB_ref,parent_dict,working_parent_dir):
    """基本思路，输入sample名，ABref dict，parent dict，working_parent_dir
    读入scaffold，align，挑选最优结果，写突变frame，以frame为基础，输出假想fasta序列，输出突变表格
    """
    #1st. #do align and keep the results 
    target_seq=AB_ref['Aref']
    file=os.path.join(working_parent_dir,sample,f"{sample}-scaffolds.fa")
    if not os.path.exists(file):
        return "NoScaffold"
    query_seqs=SeqIO.to_dict(SeqIO.parse(os.path.join(working_parent_dir,sample,f"{sample}-scaffolds.fa"),'fasta'))
    aligner = Align.PairwiseAligner()
    aligner.mismatch_score=-1
    aligner.gap_score=-5
    aligner.mode='local'
    result_alignments=dict()
    for i,query_seq in query_seqs.items():
        target=target_seq.seq;query=query_seq.seq
        alignments=list(aligner.align(target,query))+list(aligner.align(target,query.reverse_complement()))
        result_alignments[i]=max(alignments,key=lambda x:x.score)
    #read in the xxx_scaffold.fa and align with 11k
    #store the realignments in a dict
    
    #3.from alignments to snp frame
    seq_parentA=parent_dict['Aref']
    seq_parentB=parent_dict['Bref']
    print(sample)
    sample_result_frame=list()
    working_dir=os.path.join(working_parent_dir,sample)
    for key,test_ali in result_alignments.items():
        target_slice=slice(test_ali.aligned[0][0][0],test_ali.aligned[0][0][1])
        query_slice=slice(test_ali.aligned[1][0][0],test_ali.aligned[1][0][1])
        seq_dict={"seqA":list(AB_ref['Aref'].seq[target_slice]),
                 "seqB":list(AB_ref['Bref'].seq[target_slice]),
                 "parentA":list(seq_parentA.seq[target_slice]),
                 "parentB":list(seq_parentB.seq[target_slice]),
                 "query":list(test_ali.query[query_slice])}
        seq_frame=pd.DataFrame(seq_dict,index=range(test_ali.aligned[0][0][0],test_ali.aligned[0][0][1]))
        seq_frame.loc[:,"site"]=seq_frame.index
        seq_frame.loc[:,'seqB_marker']=seq_frame['seqB']!=seq_frame['seqA']
        seq_frame.loc[:,'parentB_marker']=seq_frame['parentB']!=seq_frame['seqA']
        seq_frame.loc[:,'sample']=sample
        seq_frame.loc[:,'scaffold']=key
        #print(all(seq_frame[seq_frame['parentB_marker']]['parentB']==seq_frame[seq_frame['parentB_marker']]['query']))
        check=seq_frame[seq_frame['parentB_marker']]['parentB']==seq_frame[seq_frame['parentB_marker']]['query']
        if all(check) and len(check)>0:
            tag="B"
        elif len(check)<=1:
            tag="unknown"
        else:
            tag='A'
        seq_frame.loc[:,'template']=tag
        #print(sample,key,len(check),sum(check),tag)
        if tag=="A":
            seq_frame['Mut']=seq_frame['query']!=seq_frame['parentA']
            seq_frame['CT_type']=seq_frame['parentA']+seq_frame['query']
        else:
            seq_frame['Mut']=seq_frame['query']!=seq_frame['parentB']
            seq_frame["CT_type"]=seq_frame['parentB']+seq_frame['query']
            
        #set the series 
        query_aligned=(test_ali.aligned[1][0][1]-test_ali.aligned[1][0][0])/len(test_ali.query)
        #print(sample,key,tag,len(check),sum(check),test_ali.aligned,query_aligned,len(seq_frame.index),len(result_frame.index))
        if query_aligned >0.5:
            sample_result_frame.append(seq_frame)
    if sample_result_frame:
        sample_concat_frame=pd.concat(sample_result_frame)
        sample_concat_frame.to_csv(os.path.join(working_dir,"mut_result_frame.csv"))
    else:
        return "NoUsableScaffold"
    #读入parent dict
    #4.读入sample 对应之mut result frame
    #print(sample,'make frame')
    mut_result_frame=os.path.join(working_dir,"mut_result_frame.csv")
    present_result =pd.read_csv(mut_result_frame,index_col=0)
    records=list()
    for template in ["A","B"]:
        if template not in present_result['template'].values:
            continue
        name=template+"ref"
        input_frame=present_result[present_result['Mut'] & (present_result['template']==template)]
        input_seq=[seq for key,seq in parent_dict.items() if name in key][0]
        out_seq=mutate_loci(input_frame['site'],input_frame['query'],input_seq)
        out_seq.id=sample+'--'+name
        records.append(out_seq)
    out_fasta=os.path.join(working_dir,sample+"_AB.fasta")
    SeqIO.write(records,out_fasta,'fasta')
    #5.基于拼接fasta序列，生成突变表
    sample_AB=os.path.join(working_dir,sample+"_AB.fasta")
    print(sample,'make fasta')
    sample_dict=SeqIO.to_dict(SeqIO.parse(sample_AB,'fasta'))
    frame_ls=list()
    for name,s_seq in sample_dict.items():
        input_seq=[seq for key,seq in parent_dict.items() if name[-4:] in key][0]
        sample_mut=pd.DataFrame([list(input_seq),list(s_seq)],index=['parent','sample']).T
        sample_mut['Mut']=sample_mut['parent']!=sample_mut['sample']
        sample_mut['template']=name[-4:]
        sample_mut['site']=sample_mut.index
        frame_ls.append(sample_mut)
    pd.concat(frame_ls).to_csv(os.path.join(working_dir,"mut_from_fasta.csv"))
    return "FinishedMut"
def scaffold_align_to_frame_unknown(sample,AB_ref,parent_dict,working_parent_dir):
    """基本思路，输入sample名，ABref dict，parent dict，working_parent_dir
    读入scaffold，align，挑选最优结果，写突变frame，以frame为基础，输出假想fasta序列，输出突变表格
    """
    #1st. #do align and keep the results 
    target_seq=AB_ref['Aref']
    file=os.path.join(working_parent_dir,sample,f"{sample}-scaffolds.fa")
    if not os.path.exists(file):
        return "NoScaffold"
    query_seqs=SeqIO.to_dict(SeqIO.parse(os.path.join(working_parent_dir,sample,f"{sample}-scaffolds.fa"),'fasta'))
    aligner = Align.PairwiseAligner()
    aligner.mismatch_score=-1
    aligner.gap_score=-5
    aligner.mode='local'
    result_alignments=dict()
    for i,query_seq in query_seqs.items():
        target=target_seq.seq;query=query_seq.seq
        alignments=list(aligner.align(target,query))+list(aligner.align(target,query.reverse_complement()))
        result_alignments[i]=max(alignments,key=lambda x:x.score)
    #read in the xxx_scaffold.fa and align with 11k
    #store the realignments in a dict
    
    #3.from alignments to snp frame
    seq_parentA=parent_dict['Aref']
    seq_parentB=parent_dict['Bref']
    print(sample)
    sample_result_frame=list()
    working_dir=os.path.join(working_parent_dir,sample)
    for key,test_ali in result_alignments.items():
        target_slice=slice(test_ali.aligned[0][0][0],test_ali.aligned[0][0][1])
        query_slice=slice(test_ali.aligned[1][0][0],test_ali.aligned[1][0][1])
        seq_dict={"seqA":list(AB_ref['Aref'].seq[target_slice]),
                 "seqB":list(AB_ref['Bref'].seq[target_slice]),
                 "parentA":list(seq_parentA.seq[target_slice]),
                 "parentB":list(seq_parentB.seq[target_slice]),
                 "query":list(test_ali.query[query_slice])}
        seq_frame=pd.DataFrame(seq_dict,index=range(test_ali.aligned[0][0][0],test_ali.aligned[0][0][1]))
        seq_frame.loc[:,"site"]=seq_frame.index
        seq_frame.loc[:,'seqB_marker']=seq_frame['seqB']!=seq_frame['seqA']
        seq_frame.loc[:,'parentB_marker']=seq_frame['parentB']!=seq_frame['seqA']
        seq_frame.loc[:,'sample']=sample
        seq_frame.loc[:,'scaffold']=key
        #print(all(seq_frame[seq_frame['parentB_marker']]['parentB']==seq_frame[seq_frame['parentB_marker']]['query']))
        check=seq_frame[seq_frame['parentB_marker']]['parentB']==seq_frame[seq_frame['parentB_marker']]['query']
        if all(check) and len(check)>0:
            tag="B"
        elif len(check)<=1:
            tag="u"
        else:
            tag='A'
        seq_frame.loc[:,'template']=tag
        #print(sample,key,len(check),sum(check),tag)
        if tag=="A":
            seq_frame['Mut']=seq_frame['query']!=seq_frame['parentA']
            seq_frame['CT_type']=seq_frame['parentA']+seq_frame['query']
        else:
            seq_frame['Mut']=seq_frame['query']!=seq_frame['parentB']
            seq_frame["CT_type"]=seq_frame['parentB']+seq_frame['query']
            
        #set the series 
        query_aligned=(test_ali.aligned[1][0][1]-test_ali.aligned[1][0][0])/len(test_ali.query)
        #print(sample,key,tag,len(check),sum(check),test_ali.aligned,query_aligned,len(seq_frame.index),len(result_frame.index))
        if query_aligned >0.5:
            sample_result_frame.append(seq_frame)
    if sample_result_frame:
        sample_concat_frame=pd.concat(sample_result_frame)
        sample_concat_frame.to_csv(os.path.join(working_dir,"mut_result_frame.csv"))
    else:
        return "NoUsableScaffold"
    #读入parent dict
    #4.读入sample 对应之mut result frame
    #print(sample,'make frame')
    mut_result_frame=os.path.join(working_dir,"mut_result_frame.csv")
    present_result =pd.read_csv(mut_result_frame,index_col=0)
    records=list()
    for template in ["A","B",'u']:
        if template not in present_result['template'].values:
            continue
        name=template+"ref"
        input_frame=present_result[present_result['Mut'] & (present_result['template']==template)]
        input_seq=[seq for key,seq in parent_dict.items() if name in key][0]
        out_seq=mutate_loci(input_frame['site'],input_frame['query'],input_seq)[:]
        out_seq.id=sample+'--'+name
        print(out_seq.id,template)
        records.append(out_seq)
        print([rec.id for rec in records])
    out_fasta=os.path.join(working_dir,sample+"_AB.fasta")
    
    SeqIO.write(records,out_fasta,'fasta')
    #5.基于拼接fasta序列，生成突变表
    sample_AB=os.path.join(working_dir,sample+"_AB.fasta")
    print(sample,'make fasta')
    sample_dict=SeqIO.to_dict(SeqIO.parse(sample_AB,'fasta'))
    frame_ls=list()
    for name,s_seq in sample_dict.items():
        input_seq=[seq for key,seq in parent_dict.items() if name[-4:] in key][0]
        sample_mut=pd.DataFrame([list(input_seq),list(s_seq)],index=['parent','sample']).T
        sample_mut['Mut']=sample_mut['parent']!=sample_mut['sample']
        sample_mut['template']=name[-4:]
        sample_mut['site']=sample_mut.index
        frame_ls.append(sample_mut)
    pd.concat(frame_ls).to_csv(os.path.join(working_dir,"mut_from_fasta.csv"))
    return "FinishedMut"