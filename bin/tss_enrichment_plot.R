library(argparse)
library(readr)
#library(dplyr)
options(bitmapType='cairo')
options(stringsAsFactors = FALSE)

parser = argparse::ArgumentParser(description='Script to make plot of TSS enrichment.')
parser$add_argument('--bowtie_log', required=TRUE, help='Bowtie mapping log file.')
parser$add_argument('--picard_log', required=TRUE, help='Picard deduplication log file.')
parser$add_argument('--picard_inserts_file', required=TRUE, help='Picard insert size file.')
parser$add_argument('--bam_file', required=TRUE, help='Original BAM file.')
parser$add_argument('--dedup_file', required=TRUE, help='Deduplicated BAM file.')
parser$add_argument('--bed_file', required=TRUE, help='Deduplicated BED file.')
parser$add_argument('--peak_file', required=TRUE, help='Peak calls from MACS BED file.')
parser$add_argument('--per_base_tss_region_coverage_file', required=TRUE, help='File with per base coverage for TSS regions.')
parser$add_argument('--plot', required=TRUE, help='Plot summarizing results.')
args = parser$parse_args()

tots = as.numeric(system(paste0("head -n1 ",args$bowtie_log," | cut -f1 -d' '"),intern=T))*2
mito = as.numeric(system(paste0("samtools view ",args$bam_file," | grep chrM | wc -l"),intern=T))
uniquers = as.numeric(system(paste0("samtools view -c ",args$dedup_file),intern=T))
dedups = readr::read_delim(args$picard_log,skip=6,n_max=1,delim="\t")
inserts = readr::read_delim(args$picard_inserts_file,skip=10,delim="\t")
subnuc_enrich = max(inserts$All_Reads.fr_count[inserts$insert_size<150])/max(inserts$All_Reads.fr_count[inserts$insert_size>149])
#RiDHS = as.numeric(system(paste0("bedtools intersect -a ",args$bed_file," -b <(zcat ",args$peak_file," | bedtools merge -i stdin) | wc -l"),intern=T))
RiDHS = as.numeric(system(paste0("bedtools intersect -a ",args$bed_file," -b ",args$peak_file," | wc -l"),intern=T))
#peakcount = as.numeric(system(paste0("zcat ",args$peak_file," | bedtools merge -i stdin | wc -l"),intern=T))
peakcount = as.numeric(system(paste0("cat ",args$peak_file," | wc -l"),intern=T))
sample_tss_coverage = readr::read_delim(args$per_base_tss_region_coverage_file, delim='\t')
pdf(args$plot,height=8,width=16)
par(mfrow=c(1,2))
  # TSS enrichment
  
  ## compute the min at edge of window as mean of left and right for average of 100bp on each side
  left_count = mean(sample_tss_coverage$total[sample_tss_coverage$position < min(sample_tss_coverage$position) + 100])
  right_count = mean(sample_tss_coverage$total[sample_tss_coverage$position > max(sample_tss_coverage$position) - 100])
  min_count = mean(c(left_count, right_count))

  ## compute enrichment over min and then define tss enrichment as max in the 200bp preceeding the center of the window
  sample_tss_coverage$enrichment = sample_tss_coverage$total / min_count
  tss_enrichment = max(sample_tss_coverage$enrichment[sample_tss_coverage$position > -200 & sample_tss_coverage$position < 0])
  report_vec = format(c(tots,uniquers,dedups$ESTIMATED_LIBRARY_SIZE,peakcount),big.mark=",", trim=TRUE)
  report_vec = c(report_vec[1],round(mito/tots,4),report_vec[2:3],round(RiDHS/uniquers,4),report_vec[4])
  
  plot(inserts$insert_size,inserts$All_Reads.fr_count,xlim=c(0,1200),col="red",type="l",lwd=2,xlab="Fragment Size",ylab="Count", main="Fragment Size Distribution")
  abline(v=150,lty="dashed",col="lightgray")
  plot(sample_tss_coverage$position, sample_tss_coverage$enrichment, type='l', lwd=2, col='black', xlab='Position relative to TSS (bp)', ylab='Fold Enrichment', main='TSS Enrichment')
  #grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")
  abline(v=0,lty="dashed",col="lightgray")
  legend("topright",
    c(paste0("Total Reads: ",report_vec[1]),paste0("Fraction Mito: ",report_vec[2]),paste0("Total Unique Reads: ",report_vec[3]),
    paste0("Est. Complexity: ",report_vec[4]," fragments"),paste0("FRiP: ",report_vec[5]),paste0("Called Peaks: ",report_vec[6]),
    paste0("TSS Enrichment: ", round(tss_enrichment, 4)),paste0("Subnucleosomal Enrichment: ",round(subnuc_enrich,4))),bty="n", cex=0.75, pt.cex = 1, text.font=2)
dev.off()