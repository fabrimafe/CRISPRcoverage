#!/usr/bin/env Rscript
#example to run
#first, generate a txt with the list of files you want to analyze:
#ls Or114_R1R2.sorted.RG.SL4.*tab > list.files.Or114_R1R2.sorted.RG.SL4.txt
#./chrarm_deletions_combinechromosomes.R -b list.files.Or114_R1R2.sorted.RG.SL4.txt -o Or114_R1R2.sorted.RG.SL4.cov_allchr.txt  
library(argparser, quietly=TRUE,warn.conflicts=FALSE)
library(data.table)
argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-b", help="input file. tab separated")
#argv <- add_argument(argv, "-a", help="input file: after cut site")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-w", help="window size in base pairs", default=50000)

args <- parse_args(argv)

inputfile.pre<-args$b
output.file<-args$o
window_size<-args$w
window_size<-as.numeric(window_size)
#module load R/4.0.0

#inputfile.pre<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1"
#inputfile.post<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490"
#
	df.pre<-fread(inputfile.pre)
	###TO PROCESS A SINGLE FILE
	#df<-df[grep(paste0("^",mychromosome),df$Locus),]

	df2<-tstrsplit(df.pre$Locus,":",names=c("chr","pos"))
	df.pre$pos<-as.numeric(df2$pos)
	df.pre$read<-rleid(df.pre$Total_Depth)
	df.pre.backup<-df.pre
	df.pre.read<-df.pre[!duplicated(df.pre$read),]

	df.pre.read$cuts.1<-cut(df.pre.read$pos,seq(from=min(df.pre.read$pos),to=max(df.pre.read$pos),by=window_size))
	df.pre.reads<-df.pre.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]
#	df.pre.reads$where<-inputfile.pre
	df.pre.reads$initpos<-df.pre.reads$initpos/1000000
	write.table(df.pre.reads,file=output.file,quote = FALSE,row.names = FALSE,sep="\t")
