#visualize SNPs positions from vcf
args <- commandArgs(trailingOnly = T)
data <- read.table(args[1],stringsAsFactors=F, header=F)
fai <- read.table(args[2],stringsAsFactors=F, header=F)
png_name <- args[3]
chr_name <- fai[, 1]
len_chr <- fai[, 2]


png(png_name, height=1440, width=1440, res=200)
plot(-1, type="n", bty="o", ,xaxt="n", yaxt="n", xlab="", ylab="",
	main=paste0(nrow(data), " SNPs position"), 
	cex.main=2, xlim=c(0,length(len_chr)+1), ylim=c(max(len_chr),0))
axis(side=1, at=1:length(len_chr), las=2, cex.axis=1.7,
     labels=c(chr_name))
axis(side=2, at=seq(0, 1e+8, by = 1e+7), las=2, cex.axis=1.7,
		labels=c("0M", "10M", "20M", "30M", "40M",
		         "50M", "60M", "70M", "80M", "90M", "100M"))
for(i in 1:length(len_chr)){	
  segments(i, 0, i, len_chr[i])
  data_select <- data[data[1]==chr_name[i], ]
  pos <- data_select[, 2]
  if(length(pos) > 0){
    for(j in 1:length(pos)){
      segments(i-0.3, pos[j], i+0.3, pos[j])
    }
  }
}
dev.off()


