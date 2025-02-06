#
#before usage, modify the line starts with find to map a new directory
FQ_DIR=""
find $FQ_DIR -maxdepth 5 -name "*_1.fq.gz"  | sed 's/_1.fq.gz$//' | \
    xargs -n 1 -P 4 -I PREFIX \
    sh -c '
		CPU=4
        lane_id=`basename PREFIX`
        sample=`basename PREFIX`
		work_dir="/mnt/san1/usr/hj/Data/01.neurospora/01.bams/DNA2489bam/batch_20230126"
        echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"

       ## Align reads with BWA-MEM algorithm
       /home/wl/Data/biosoft/bwa-0.7.10/bin/bwa mem -t $CPU -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
            /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa ${read1} ${read2} \
            > ${work_dir}/${sample}.bwa.sam \
            2> ${work_dir}/${sample}.bwa.log
	echo "[`date`]: Successfully aligned ${sample} using BWA-MEM ... "

        ## sort bam file
        java -Djava.io.tmpdir=/home/hj/Data/tmp -jar /mnt/san1/usr/thw/tansoft/picard-tools-1.124/picard.jar SortSam \
            INPUT=${work_dir}/${sample}.bwa.sam \
            OUTPUT=${work_dir}/${sample}.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> ${work_dir}/${sample}.bwa.log 2>&1 && \

        echo "[`date`]: Start marking duplicates ${sample} ... "

        ## mark duplicates
        java -jar /mnt/san1/usr/thw/tansoft/picard-tools-1.124/picard.jar MarkDuplicates \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=${work_dir}/${sample}.bwa.sort.bam \
            OUTPUT=${work_dir}/${sample}.bwa.sort.dedup.bam \
            METRICS_FILE=${work_dir}/${sample}.bwa.sort.dedup.metrics \
            >>${work_dir}/${sample}.bwa.log 2>&1 && \

        ## index bam file
        samtools-1.3.1 index ${work_dir}/${sample}.bwa.sort.dedup.bam

        echo "[`date`]: Start realigning ${sample} ... "

        ## realignment
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
            -R /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa \
            -T RealignerTargetCreator -nt $CPU \
            -o ${work_dir}/${sample}.bwa.sort.dedup.realn.intervals \
            -I ${work_dir}/${sample}.bwa.sort.dedup.bam \
            >> ${work_dir}/${sample}.bwa.log 2>&1

        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
            -R /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa \
            -T IndelRealigner \
            -targetIntervals ${work_dir}/${sample}.bwa.sort.dedup.realn.intervals \
            -o ${work_dir}/${sample}.realn.bam \
            -I ${work_dir}/${sample}.bwa.sort.dedup.bam \
            >> ${work_dir}.bwa.log 2>&1  \
'
