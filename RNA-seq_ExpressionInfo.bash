#Mapping by STAR

for i in `ls ./20210716cleandata |grep gz| perl -e 'while(<>){s/\d\.clean\.fq\.gz$//; print;}' | sort| uniq`; 
do STAR --genomeDir ./genomeandannotation/STARindexCM 
              --readFilesIn ./20210716cleandata/${i}1.clean.fq.gz ./20210716cleandata/${i}2.clean.fq.gz 
              --runThreadN 6 --outFileNamePrefix ${i} 
              --outSAMtype BAM SortedByCoordinate 
              --readFilesCommand gunzip -c 
              --alignIntronMax 4000 
              --alignIntronMin 1 
              --outSJfilterReads Unique 
              --outFilterMismatchNmax 2; 
done


#featureCounts
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_22_STAR_featureCounts_1 t12_22_STARsorted-1Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_22_STAR_featureCounts_2 t12_22_STARsorted-2Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_22_STAR_featureCounts_3 t12_22_STARsorted-3Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_28_STAR_featureCounts_1 t12_28_STARsorted-1Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_28_STAR_featureCounts_2 t12_28_STARsorted-2Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o t12_28_STAR_featureCounts_3 t12_28_STARsorted-3Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_22_STAR_featureCounts_1 WT_22_STARsorted-1Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_22_STAR_featureCounts_2 WT_22_STARsorted-2Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_22_STAR_featureCounts_3 WT_22_STARsorted-3Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_28_STAR_featureCounts_1 WT_28_STARsorted-1Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_28_STAR_featureCounts_2 WT_28_STARsorted-2Aligned.sortedByCoord.out.bam -T 6 &
nohup featureCounts -p -a /lustre/user/liclab/lisky/panyg/LanJQ/Reference/Ath_GENE_TE_Clean.nochr.gtf -o WT_28_STAR_featureCounts_3 WT_28_STARsorted-3Aligned.sortedByCoord.out.bam -T 6 &

