library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)
x = as.data.frame(table(x))
colnames(x) = c("x", "Freq")
x$x = as.numeric(as.character(x$x))
print(summary(x))

txtFontSize=10
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.x=element_line(size=0.7, color="black"), axis.line.y=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))

# VAC
png("vac.png", height=600, width=800)
p1=ggplot(data=x, aes(x=x, y=Freq))+geom_line(colour="darkblue")+ xlab("Sample count") + ylab("#Non-ref k-mers") + scale_x_log10(labels=comma, breaks=c(1,10,100,1000,2500), limits=c(1, max(x$x))) + scale_y_log10(labels=comma, breaks=c(1,10,100,1000,10000,100000,1000000,10000000), limits=c(min(x$Freq), 10000000)) + scienceTheme
p1=p1 + geom_vline(xintercept=1, linetype="longdash") + geom_vline(xintercept=2, linetype="longdash") + geom_vline(xintercept=3, linetype="longdash") + geom_vline(xintercept=4, linetype="longdash") + geom_vline(xintercept=5, linetype="longdash")
p1
dev.off()
print(warnings())
