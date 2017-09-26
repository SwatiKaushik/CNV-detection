setwd("/Users/swati/Desktop/cnv-enrichment-EGFRm/EGFR-mutant")
rm(list=ls())

file <- read.table("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt", header=T, sep="\t")
file <- file[,-c(2,3)]
t <- t(file)
write.table(t, file="LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt.t", sep="\t", col.names =FALSE)

cnv.file <- read.table("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt.t", header=T, sep="\t")
gene.sym <- cnv.file$Gene.Symbol
sub.str <- substr(gene.sym,1,15)
new.sample.name <- gsub(".",'-',sub.str, fixed=T)
cnv.file$Gene.Symbol <- new.sample.name

mutant.sample.file <- read.table("EGFR-m-nosplicesite-FSD.txt", header=T, sep="\t")
mutant.sample.file <- mutant.sample.file[c(1),]
mutant.sample.names <- gsub(".",'-',colnames(mutant.sample.file), fixed=T)
colnames(mutant.sample.file) <- mutant.sample.names
sampl <- t(mutant.sample.file)
write.table(sampl, file="EGFR-m-nosplicesite-FSD.T", sep="\t",col.names =FALSE)
mutant.sample.file <- read.table("EGFR-m-nosplicesite-FSD.T", header=T, sep="\t",)

merged.data <- merge(mutant.sample.file, cnv.file, by.x = 'A', by.y = 'Gene.Symbol')
write.table(merged.data, file="merged.file", sep="\t")

output <- vector()

for (i in 3:ncol(merged.data)){  
	
	amp <- subset(merged.data[,2], merged.data[,i] >1)
	del <- subset(merged.data[,2], merged.data[,i] <= 1) #& merged.data[,i] >= -1))
	amp.len <- length(amp)
	del.len <- length(del)
	cat (colnames(merged.data[i]),"\n")
	
	if ((amp.len >2) & (del.len >2)){

		p.value <- wilcox.test(amp,del)$p.value
		mean.amp <- mean(amp)
		mean.del <- mean(del)		
		#cat <- paste(colnames(merged.data[i]), p.value,amp.len, del.len,mean.amp, mean.del, sep="	")
		name.gene <- colnames(merged.data[i])
		output <- rbind(output, data.frame(name.gene, p.value,amp.len, del.len,mean.amp, mean.del))
	} else { 
		#cat (colnames(merged.data[i]),"NA\n")
		#no.data <- paste(colnames(merged.data[i]), NA, sep="	")
		name.gene <- colnames(merged.data[i])
		p.value <- NA
		amp.len <- NA
		del.len <- NA
		mean.amp <- NA
		mean.del <- NA
		output <- rbind(output, data.frame(name.gene,  p.value,amp.len, del.len,mean.amp, mean.del))
	}
			
}
output$padjust <- p.adjust(output$p.value, method="fdr")
write.table(output, file="output.amp.others", sep="\t", row.names= FALSE, col.names=FALSE, quote = FALSE)


