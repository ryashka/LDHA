featureCounts -M -T 32 -p -a \
~/Documents/Programming/reference_files/mm10_star_reference_files/gencode.vM25.primary_assembly.annotation.gtf  \
-t exon -g gene_id -o ~/Documents/Programming/Rahul_Sequencing/Project_10863_B/featureCounts_files/combinedFeatureCounts.txt \
~/Documents/Programming/Rahul_Sequencing/Project_10863_B/output_files/*.bam
