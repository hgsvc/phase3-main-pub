library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)

x = read.table(args[1], header=TRUE)

g1k = F
if ("FIN" %in% x$pop) {
 g1k = T
 x$pop=factor(x$pop, levels=c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"))
 popCol=c("#FFCD00", "#FFB900", "#CC9933", "#E1B919", "#FFB933", "#FF9900", "#FF6600", "#CC3333", "#E10033", "#FF0000", "#CC3300", "#339900", "#ADCD00", "#00FF00", "#008B00", "#00CC33", "#0000FF", "#00C5CD", "#00EBFF", "#6495ED", "#00008B", "#8B008B", "#9400D3", "#B03060", "#E11289", "#FF00FF")
 names(popCol) = levels(x$pop)
}

print(table(x$pop))

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.x=element_line(size=0.7, color="black"), axis.line.y=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))

pdf(paste0(args[1], ".pdf"), width=12, height=10)
#p1 = ggplot(data=x, aes(x=pop, y=kmer))
p1 = ggplot(data=x, aes(x=pop, y=singleton))
#p1 = p1 + geom_violin(aes(fill=pop), adjust=3) + geom_point(position=position_jitter(w=0.05, h=0.05), alpha=0.1, size=0.1)
p1 = p1 + geom_boxplot(aes(fill=pop))
p1 = p1 + scienceTheme
if (g1k) { p1=p1 + scale_fill_manual(values=popCol, guide=F); } else { p1=p1 + theme(legend.position = "none"); }
p1 = p1 + xlab("Populations") + coord_flip() + theme(axis.text.x = element_text(angle=45, hjust=1))
p1 = p1 + ylab("#non-ref k-mer")
p1 = p1 + scale_y_continuous(labels=comma)
p1
dev.off()
