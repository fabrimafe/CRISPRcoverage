#!/usr/bin/env Rscript
#script to generate genome-wide coverage plots as in fig.6.
#for instructions run:
#Rscript allchr_coverageplot.R --help

library(argparser, quietly=TRUE,warn.conflicts=FALSE)

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-i", help="file identifier: string that has to appear in all input files", default="250k.cov.txt")
argv <- add_argument(argv, "-l", help="file with name of plants (one per row)")


args <- parse_args(argv)


inputfile.string<-args$i
output.file<-args$o
inputfile.plants<-args$l
if (is.na(inputfile.string)){ inputfile.string<-"250k.cov.txt"  }
myplants<-read.table(inputfile.plants,header=FALSE)
myplants<-unname(c(unlist(myplants[,1])))


pdf(output.file)
par(mfrow=c(3,1))
#for ( plant in c("Or114","Or541","Or9-1-2")){
for ( plant in myplants){
myfiles<-dir()[grep(plant,dir())]
myfiles<-myfiles[grep(inputfile.string,myfiles)]
for (i in 1:length(myfiles))
        {
        print(i)
        mydata_t<-read.table(myfiles[i],header=TRUE)
        mydata_t$chr<-i
        mydata_t$initpos_tot<-mydata_t$initpos
        if (i==1){mydata<-mydata_t} else {
        mydata_t$initpos_tot<-mydata_t$initpos_tot+0.1+max(mydata$initpos_tot)
        mydata<-rbind(mydata,mydata_t)}
        }
mydata<-mydata[mydata$chr!=1,]
plot(mydata$initpos_tot,mydata$nreads/mean(mydata$nreads),pch=19,col=mydata$chr,xlab="genomic position",ylab="coverage",xaxt="n",yaxt="n",ylim=c(0,2),main=plant)
axis(side=2,at=seq(0,5,0.5),labels=seq(0,5,0.5))
abline(h=1,lty=2)
abline(h=0.5,lty=2)
abline(h=1.5,lty=2)
write.table(mydata,file=paste0(plant,"_covtot.txt"),quote=FALSE,row.names = FALSE, sep="\t")
}
dev.off()

