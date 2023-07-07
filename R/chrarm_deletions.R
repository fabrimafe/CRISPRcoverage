#!/usr/bin/env Rscript
#example to run
#./chrarm_deletions.R -b /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr11_from1 -a /home/labs/alevy/Collaboration/Tomato_WGS/bam_sorted_Or/read_coverage/Or114/Or114_chr02_from_47124490 -o out.txt

library(argparser, quietly=TRUE,warn.conflicts=FALSE)

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-b", help="input file: before cut site")
argv <- add_argument(argv, "-a", help="input file: after cut site")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")

args <- parse_args(argv)

inputfile.pre<-args$b
inputfile.post<-args$a
output.file<-args$o
#module load R/4.0.0

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

print("calculating p.values")
counter<-1
for ( window_size in c(0,1000,20000,1000000))
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
