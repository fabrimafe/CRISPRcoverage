#!/usr/bin/env Rscript
#example to run
#./chrarm_deletions_plot.R -b /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1 -a /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490 -o out2.txt
#./chrarm_deletions_plot.R -b /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1 -a /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490 -r /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/Or114_covtot.txt -o out3.txt

library(argparser, quietly=TRUE,warn.conflicts=FALSE)

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-b", help="input file: before cut site")
argv <- add_argument(argv, "-a", help="input file: after cut site")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-r", help="use external input files for coverage normalization")

args <- parse_args(argv)


inputfile.pre<-args$b
inputfile.post<-args$a
output.file<-args$o
cov_normalization<-args$r
if (is.na(cov_normalization)){ cov_normalization<-0 } else { cov_normalization<-as.character(cov_normalization) }
#module load R/4.0.0

#inputfile.pre<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1"
#inputfile.post<-"/home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490"
#


library(data.table)
df.pre<-fread(inputfile.pre)
df.post<-fread(inputfile.post)
###TO PROCESS A SINGLE FILE
#df<-df[grep(paste0("^",mychromosome),df$Locus),]

print("loading files")
df2<-tstrsplit(df.pre$Locus,":",names=c("chr","pos"))
df.pre$pos<-as.numeric(df2$pos)
df2<-tstrsplit(df.post$Locus,":",names=c("chr","pos"))
df.post$pos<-as.numeric(df2$pos)

#rm(df2)
#df$pos>as.numeric(position)
#strsplit(df$Locus[1:10])
#i_pos<-which(df$Locus==paste0(mychromosome,":",position))
###TO PROCESS TWO FILES: BEFORE AND AFTER POSITION


df.pre$read<-rleid(df.pre$Total_Depth)
df.post$read<-rleid(df.post$Total_Depth)

df.pre.backup<-df.pre
df.post.backup<-df.post

df.pre.read<-df.pre[!duplicated(df.pre$read),]
df.post.read<-df.post[!duplicated(df.post$read),]

#--------------generate genome wide plot------------
window_size<-50000

df.pre.read$cuts.1<-cut(df.pre.read$pos,seq(from=min(df.pre.read$pos),to=max(df.pre.read$pos),by=window_size))
df.post.read$cuts.1<-cut(df.post.read$pos,seq(from=min(df.post.read$pos),to=max(df.post.read$pos),by=window_size))

df.post.reads<-df.post.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]
df.pre.reads<-df.pre.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]

df.pre.reads$where<-0
df.post.reads$where<-1

df.post.reads$initpos<-df.post.reads$initpos/1000000
df.pre.reads$initpos<-df.pre.reads$initpos/1000000

#output.file<-"temp"

if (cov_normalization!=0){
cov_normalization<-read.table(cov_normalization,header=TRUE)
meancov.ref<-mean(cov_normalization$nreads)
} else { 
meancov.ref<-mean(df.pre.reads$nreads)
}

plot.new()
pdf(paste0(output.file,"_wholegenome.coverageplot.pdf"))
df.reads<-rbind(df.pre.reads,df.post.reads)
plot(df.reads$initpos,df.reads$nreads/meancov.ref,type="n",xlab="genomic position (Mb)",ylab="coverage",yaxt="n")
axis(side=2,at=seq(0,5,0.5),labels=seq(0,5,0.5))
points(df.pre.reads$initpos,df.pre.reads$nreads/meancov.ref,pch=19,col="gray54")
points(df.post.reads$initpos,df.post.reads$nreads/meancov.ref,pch=19,col="gray75")
segments(min(df.pre.reads$initpos),1,max(df.pre.reads$initpos),1,col="gray44",lty=2,lwd=4)
segments(min(df.post.reads$initpos),mean(df.post.reads$nreads)/meancov.ref,max(df.post.reads$initpos),mean(df.post.reads$nreads)/meancov.ref,col="gray85",lty=2,lwd=4)
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

#--------------generate zoom in plots------------
window_size<-5000

meancov.ref<-meancov.ref/10

df.pre.read$cuts.1<-cut(df.pre.read$pos,seq(from=min(df.pre.read$pos),to=max(df.pre.read$pos),by=window_size))
df.post.read$cuts.1<-cut(df.post.read$pos,seq(from=min(df.post.read$pos),to=max(df.post.read$pos),by=window_size))

df.post.reads<-df.post.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]
df.pre.reads<-df.pre.read[,.(nreads=length(read),initpos=min(pos)),by="cuts.1"]

df.pre.reads$where<-"genome"
df.post.reads$where<-"cut site"

df.post.reads$initpos<-df.post.reads$initpos/1000000
df.pre.reads$initpos<-df.pre.reads$initpos/1000000

df.pre.reads<-tail(df.pre.reads,1000)
df.post.reads<-head(df.post.reads,1000)


plot.new()
pdf(paste0(output.file,"_zoomin.coverageplot.pdf"))
#plot.new()
df.reads<-rbind(df.pre.reads,df.post.reads)
plot(df.reads$initpos,df.reads$nreads/meancov.ref,type="n",xlab="genomic position (Mb)",ylab="coverage",yaxt="n")
axis(side=2,at=seq(0,5,0.5),labels=seq(0,5,0.5))
points(df.pre.reads$initpos,df.pre.reads$nreads/meancov.ref,pch=19,col="gray54")
points(df.post.reads$initpos,df.post.reads$nreads/meancov.ref,pch=19,col="gray75")
segments(min(df.pre.reads$initpos),1,max(df.pre.reads$initpos),1,col="gray44",lty=2,lwd=4)
segments(min(df.post.reads$initpos),mean(df.post.reads$nreads)/meancov.ref,max(df.post.reads$initpos),mean(df.post.reads$nreads)/meancov.ref,col="gray85",lty=2,lwd=4)
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
#ggsave(p,file=paste0(output.file,"_zoomin.violincoverageplot.pdf"))

#------------------------------------------------

print("calculating p.values")
counter<-1
for ( window_size in c(0,1000,1000000))
{
#do test also on smaller region to avoid repeated regions and so on.

window_size_t<-window_size
if (window_size==0){window_size_t<-999999999999}
df.pre<-df.pre.backup[df.pre.backup$pos>(max(df.pre.backup$pos)-window_size_t)]
df.post<-df.post.backup[df.post.backup$pos<(min(df.post.backup$pos)+window_size_t)]

df.pre.read<-df.pre[!duplicated(df.pre$read),]
df.post.read<-df.post[!duplicated(df.post$read),]



if ( window_size==0 )
{
    l.pre<-max(df.pre$pos)-min(df.pre$pos)
    l.post<-max(df.post$pos)-min(df.post$pos)
    } else 
    { 
    l.pre<-window_size
    l.post<-window_size
}

newreads<-df.pre.read$Total_Depth[-1]-df.pre.read$Total_Depth[-length(df.pre.read$Total_Depth)]
tot_reads.pre<-sum(newreads[newreads>0])+df.pre.read$Total_Depth[1]
newreads<-df.post.read$Total_Depth[-1]-df.post.read$Total_Depth[-length(df.post.read$Total_Depth)]
tot_reads.post<-sum(newreads[newreads>0])+df.post.read$Total_Depth[1]


ll.2<-dbinom(x=tot_reads.post,size=tot_reads.post+tot_reads.pre,prob=l.post/(l.post+l.pre),log=TRUE)
ll.1<-dbinom(x=tot_reads.post,size=tot_reads.post+tot_reads.pre,prob=l.post/2/(l.post+l.pre),log=TRUE)
ll.0<-dbinom(x=tot_reads.post,size=tot_reads.post+tot_reads.pre,prob=l.post/1000/(l.post+l.pre),log=TRUE)
ll_l<-c(ll.0,ll.1,ll.2)

namesmodels<-c("complete.deletion","1chr.deletion","no.deletion")
bestmodel<-which(max(ll_l)==ll_l)[1]
re.ll<-exp(ll_l-ll_l[bestmodel])

temp_res<-c(window_size,tot_reads.post/l.post/(tot_reads.pre/l.pre),re.ll,ll_l,namesmodels[bestmodel])
if (counter==1){ res<-temp_res } else { res<-rbind(res,temp_res) }
counter<-counter+1
}
rownames(res)<-NULL
colnames(res)<-c("window_size","coverage_pre/coverage_post","pvalue.complete.deletion","pvalue.1chr.deletion","pvalue.no.deletion","loglik.complete.deletion","loglik.1chr.deletion","loglik.no.deletion","best_model")
#output.file="out.txt"
write.table(res,file=output.file,quote=FALSE)

