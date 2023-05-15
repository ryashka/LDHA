for sn in {1..6}
do

cat $(find . -name "GC_${sn}*R1*.fastq.gz" | sort) > input_files/"Sample_${sn}_R1.fq.gz"
cat $(find . -name "GC_${sn}*R2*.fastq.gz" | sort) > input_files/"Sample_${sn}_R2.fq.gz"

gatk FastqToSam \
-F1 input_files/"Sample_${sn}_R1.fq.gz" \
-F2 input_files/"Sample_${sn}_R2.fq.gz" \
-O samfiles/"${sn}_fastqtosam.bam" \
-SM $sn

gatk --java-options "-Xmx100G" MarkIlluminaAdapters \
-I samfiles/"${sn}_fastqtosam.bam" \
-O samfiles/"${sn}_markilluminaadapters.bam" \
-M samfiles/"${sn}_markilluminaadapters_metrics.txt" \

gatk --java-options "-Xmx100G" ClipReads \
-I samfiles/"${sn}_markilluminaadapters.bam" \
-O samfiles/"${sn}_clipped.bam" \
-QT 10

gatk --java-options "-Xmx100G" SamToFastq \
-I samfiles/"${sn}_clipped.bam" \
--FASTQ fastq_files/"${sn}_samtofastq_R1.fq" \
--SECOND_END_FASTQ fastq_files/"${sn}_samtofastq_R2.fq" \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE false --NON_PF true \

STAR \
--genomeDir ~/Documents/Programming/reference_files/mm10_star_reference_files/star_genome \
--readFilesIn fastq_files/"${sn}_samtofastq_R1.fq" fastq_files/"${sn}_samtofastq_R2.fq" \
--outFileNamePrefix samfiles/"${sn}_" \
--runThreadN 400

gatk --java-options "-Xmx100G" MergeBamAlignment \
--ALIGNED_BAM samfiles/"${sn}_Aligned.out.sam" \
--UNMAPPED_BAM samfiles/"${sn}_fastqtosam.bam" \
--OUTPUT samfiles/"${sn}_not_piped.bam" \
-R ~/Documents/Programming/reference_files/mm10_star_reference_files/GRCm38.primary_assembly.genome.fa --CREATE_INDEX true --ADD_MATE_CIGAR true \
--CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS \

gatk --java-options "-Xmx100G" MarkDuplicatesSpark \
-I samfiles/"${sn}_not_piped.bam" \
-O output_files/"${sn}_dms_with_clipping.bam" \
--remove-sequencing-duplicates

done

