#!/usr/bin/env Rscript
#example to run
#./chrarm_deletions_plot.R -b /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1 -a /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490 -o out2.txt
#./chrarm_deletions_plot.R -b /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1 -a /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490 -r /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/Or114_covtot.txt -o out3.txt

library(argparser, quietly=TRUE,warn.conflicts=FALSE)

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-b", help="input file: before cut site")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-r", help="use external input files for coverage normalization")
argv <- add_argument(argv, "-w", help="window size", default=1000)
argv <- add_argument(argv, "-m", help="minimum position", default=1)
argv <- add_argument(argv, "-M", help="maximum position", default=10000)
args <- parse_args(argv)


inputfile.pre<-args$b
output.file<-args$o
cov_normalization<-args$r
window_size<-as.numeric(args$w)
minimum.pos<-as.numeric(args$m)
maximum.pos<-as.numeric(args$M)
if (is.na(cov_normalization)){ cov_normalization<-0 } else { cov_normalization<-as.character(cov_normalization) }
if (is.na(window_size)){ window_size<-10000 }
if (is.na(minimum.pos)){ minimum.pos<-1 }
if (is.na(maximum.pos)){ maximum.pos<-10000 }
#module load R/4.0.0

#inputfile.pre<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1"
#inputfile.post<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490"
#


library(data.table)
df.pre<-fread(inputfile.pre)
###TO PROCESS A SINGLE FILE
#df<-df[grep(paste0("^",mychromosome),df$Locus),]

print("loading files")
df2<-tstrsplit(df.pre$Locus,":",names=c("chr","pos"))
df.pre$pos<-as.numeric(df2$pos)

#rm(df2)
#df$pos>as.numeric(position)
#strsplit(df$Locus[1:10])
#i_pos<-which(df$Locus==paste0(mychromosome,":",position))
###TO PROCESS TWO FILES: BEFORE AND AFTER POSITION


df.pre$read<-rleid(df.pre$Total_Depth)

df.pre.backup<-df.pre

df.pre.read<-df.pre[!duplicated(df.pre$read),]

#--------------generate genome wide plot------------

df.pre.read$cuts.1<-cut(df.pre.read$pos,seq(from=min(df.pre.read$pos),to=max(df.pre.read$pos),by=window_size))
df.reads<-df.pre.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]

df.reads$where<-0

df.reads$initpos<-df.reads$initpos/1000000

#output.file<-"temp"

if (cov_normalization!=0){
cov_normalization<-read.table(cov_normalization,header=TRUE)
meancov.ref<-mean(cov_normalization$nreads)
} else { 
meancov.ref<-mean(df.reads$nreads)
}
df.reads$meancov<-df.reads$nreads/meancov.ref


df.reads<-df.reads[df.reads$initpos>(minimum.pos/1000000) & df.reads$initpos<(maximum.pos/1000000),]
write.table(df.reads,file=paste0(output.file,".table"))

plot.new()
pdf(paste0(output.file,"_wholegenome.coveragescan.pdf"))
plot(df.reads$initpos,df.reads$nreads/meancov.ref,type="n",xlab="genomic position (Mb)",ylab="coverage",yaxt="n")
axis(side=2,at=seq(0,5,0.5),labels=seq(0,5,0.5))
points(df.reads$initpos,df.reads$nreads/meancov.ref,pch=19,col="gray54")
segments(min(df.reads$initpos),1,max(df.reads$initpos),1,col="gray44",lty=2,lwd=4)
dev.off()

#df.reads$coverage<-df.reads$nreads/mean(df.pre.reads$nreads)
#library(ggplot2)
#data_summary <- function(x) {
#   m <- mean(x)
#   ymin <- m-sd(x)
#   ymax <- m+sd(x)
#   return(c(y=m,ymin=ymin,ymax=ymax))
#}
#p <- ggplot(df.reads, aes(x=where, y=coverage)) + geom_violin() + stat_summary(fun.data=data_summary)
#ggsave(p,file=paste0(output.file,"_genomewide.violincoverageplot.pdf"))


#df.reads$coverage<-df.reads$nreads/mean(df.pre.reads$nreads)
#library(ggplot2)
#data_summary <- function(x) {
#   m <- mean(x)
#   ymin <- m-sd(x)
#   ymax <- m+sd(x)
#   return(c(y=m,ymin=ymin,ymax=ymax))
#}
#p <- ggplot(df.reads, aes(x=where, y=coverage)) + geom_violin() + stat_summary(fun.data=data_summary)
#ggsave(p,file=paste0(output.file,"_zoomin.violincoverageplot.pdf"))

#------------------------------------------------

