#modify the DIR and the software parameters before usage

cd $WORKDIR
 find $FQDIR -name "*_1.fq.gz" |sed 's/_1.fq.gz//' | xargs -I FILE -P 1 sh -c '
NAME=`basename FILE`
echo $NAME
mkdir $NAME
cd /mnt/san1/usr/hj/Data/01.neurospora/08.methyl/01.ref.NC12v42contig20MtHph
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark --genome /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/ -o $NAME --basename $NAME --bowtie2 -p 4 -1 FILE_1.fq.gz -2 FILE_2.fq.gz &> $NAME/${NAME}.log 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/deduplicate_bismark --bam --output_dir $NAME $NAME/${NAME}_pe.bam &>> $NAME/${NAME}.log 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark_methylation_extractor --gzip --buffer_size 30G --comprehensive --bedGraph --CX --cytosine_report \
	--genome_folder /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/ -o $NAME $NAME/${NAME}_pe.deduplicated.bam &>> $NAME/${NAME}.log 
	cd $NAME &&	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2report 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2summary ${NAME}_pe.deduplicated.bam -o {NAME}.summary_report'
#End


#methyl genome batch of kuac19-3 F22FTSECKF12166_FUNvqenH

cd /mnt/san1/usr/hj/Data/01.neurospora/08.methyl/01.ref.NC12v42contig20MtHph
 find /mnt/nas_sy/usr/hj/neurospora/00.rawdata/methyl/F22FTSECKF12166_FUNvqenH -name "*_1.fq.gz" |sed 's/_1.fq.gz//' | xargs -I FILE -P 5 sh -c '
NAME=`basename FILE`
echo $NAME
mkdir $NAME
cd /mnt/san1/usr/hj/Data/01.neurospora/08.methyl/01.ref.NC12v42contig20MtHph
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark --genome /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/ -o $NAME --basename $NAME --bowtie2 -p 4 -1 FILE_1.fq.gz -2 FILE_2.fq.gz &> $NAME/${NAME}.log 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/deduplicate_bismark --bam --output_dir $NAME $NAME/${NAME}_pe.bam &>> $NAME/${NAME}.log 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark_methylation_extractor --gzip --buffer_size 30G --comprehensive --bedGraph --CX --cytosine_report \
	--genome_folder /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/ -o $NAME $NAME/${NAME}_pe.deduplicated.bam &>> $NAME/${NAME}.log 
	cd $NAME &&	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2report 
	/mnt/san1/usr/thw/tansoft/bismark_v0.22.1/bismark2summary ${NAME}_pe.deduplicated.bam -o {NAME}.summary_report'
#End
