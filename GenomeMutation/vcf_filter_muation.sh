#1st stage. Prepare markers.
#prepare the sample_group file
mkdir readcounts
dos2unix bams
dos2unix samples.group.txt
vcf_process.pl --vcf samples.vcf.gz --quality 30 --rare-only 1 | \
      bgzip -c > samples.fq3.vcf.gz && tabix -p vcf samples.fq3.vcf.gz
  bgzip -dc samples.fq3.vcf.gz |sed 's/^##fileformat=VCFv4.2/##fileformat=VCFv4.1/'| vcf-annotate --fill-type |\
      perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1];
          print "$chrom\t",($pos-1),"\t$pos\n";' > samples.fq3.snp.bed
#Extract candidate regions for indels

  bgzip -dc samples.fq3.vcf.gz |sed 's/^##fileformat=VCFv4.2/##fileformat=VCFv4.1/'| vcf-annotate --fill-type | \
      perl -ne 'next if(/\#/); next unless(/ins/ || /del/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
          my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq | \
      bedtools slop -i - -b 10 -g /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa.fai > samples.fq3.indel.bed
#the right version
#Step2: Count accurate allele depths for each locus and each sample
grep -v "#" bams | xargs -n 1 -P 6 -I PREFIX \
      sh -c '
          sample=`basename PREFIX | cut -d"." -f1`
          
          echo "[`date`]: Start processing ${sample} ... "
          
          samtools mpileup -Ad100000 -q20 -f /mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa -l /mnt/san1/usr/hj/Data/01.neurospora/02.vcf/5AC-dim2_v2/samples.fq3.snp.bed  PREFIX | grep -vP "\t0\t" \
              >readcounts/samples.fq3.snp.${sample}.mpileup
         java -jar /home/wl/Data/biosoft/varscan-2.4.2/VarScan.v2.4.2.jar readcounts readcounts/samples.fq3.snp.${sample}.mpileup \
              --min-base-qual 20 --min-coverage 1 \
              --output-file readcounts/samples.fq3.snp.${sample}.readcounts
          echo "[`date`]: Finished processing ${sample}"
              
            '
#Get the file list of all counting results

  for f in `find readcounts/ -name "samples.fq3.snp.*.readcounts" | sort`;
  do
      library=`basename $f | cut -d"." -f4`
      sample=${library}
      
      echo "${sample} ${library} ${f}"
  done > samples.fq3.snp.readcounts.list
fillVcfDepth.pl --vcf samples.fq3.vcf.gz --list samples.fq3.snp.readcounts.list \
      --minimum-vcf --update-AD | bgzip -c \
      > samples.fq3.snp.vcf.gz && tabix -p vcf samples.fq3.snp.vcf.gz
    bgzip -dc samples.fq3.snp.vcf.gz |sed 's/^##fileformat=VCFv4.2/##fileformat=VCFv4.1/' |
detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 \
      --max-cmp-miss 0 --min-supp-depth 5 --min-supp-plus 1 --min-supp-minus 1 \
      --mask-only LowDepth -g samples.group.txt --controls  3-dim2AB7g 3-dim2a4 | \
      vcf-annotate -f c=3,150 --fill-type > samples.fq3.snp.mut.c2t3d5m5.vcf
