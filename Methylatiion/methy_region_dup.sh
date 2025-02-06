bismark_dir=/mnt/san3/usr/cjx/all-methyl-CX/methyl_bismark
cd $bismark_dir
for file in `find ${bismark_dir} -name "Nc.*"`
do
    sample=`basename ${file}`
    # echo ${sample##*.}
    ## merge window methylation from site

    less ${sample}/${sample##*.}_pe.deduplicated.CX_report.txt.gz \
    |awk -v OFS='\t' '$4+$5>5{print $1,$2-1,$2,$4/($4+$5)}' |LC_ALL=C sort -S 80% -k 1,1 -k 2n - \
    -o ${sample}/${sample##*.}.CX.bed

    LC_ALL=C sort -S 80% -k 1,1 -k 2n /home/cjx/Data/data/neurospora/Nc2489_200bin.bed \
    |bedtools map -a - \
    -b ${sample}/${sample##*.}.CX.bed \
    -c 4 -o mean | awk '$4!="."{print $0}' \
    > ${sample}/${sample##*.}.CX.200w.bed

    /home/cjx/Data/set_up_backage/UCSC_Kent/bedGraphToBigWig ${sample}/${sample##*.}.CX.200w.bed \
    /home/cjx/Data/data/neurospora/Nc2489.chrom.size \
    ${sample}/${sample##*.}.CX.200w.bw

    gzip ${sample}/${sample##*.}.CX.bed

    ## Calculation area methylation

    echo -n -e "${sample}\t" && bedtools intersect -a ${sample}/${sample##*.}.CX.200w.bed \
    -b /home/cjx/Data/data/neurospora/dup/duplication_3kinds.s.bed \
    -wo \
    |awk -v OFS='\t' '{print $8,$9,$4,$1,$2,$3}' \
    |sort \
    |awk '$2==200' \
    |groupBy -g 1 -c 3 -o mean \
    |xargs > ${sample}/${sample##*.}.dupkind.methylation.txt

done 