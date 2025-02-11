#!/usr/bin/env Rscript

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 06/27/23

# Purpose: Plot collapsed and resolved duplications.
# Input: requires four strings as arguments:
#          workflowDirectory filepath;
#          tsv filepath of all results with cols gene_num, gene, and depth;
#          tsv filepath of resolved results with cols gene_num, gene, and depth;
#          tsv filepath of collapsed results with cols gene_num, gene, and depth;
#          dataOutputSuffix string
# Output: Saves four depth plots to "workflowDirectory/../results" directory.

library(ggplot2)
library(gridExtra)
library(grid)

# Parse input
params <- commandArgs(trailingOnly = TRUE)

workflowDir <- params[[1]]
inputPath_combo <- params[[2]]
inputPath_res <- params[[3]]
inputPath_col <- params[[4]]
dataOutputSuffix <- params[[5]]

# Plotting Function and Variables
themeBaseSize <- 12
plotWidthInInches <- 7.5
geneDepthPlotter <- function (depthsDF, tickDist, themeBaseSize=12) {
    axisTextSize <- 8
    xExpansion <- 20
    yExpansion <- 1
    pointShape <- 1
    pointSize <- 1
    pointColor <- "#696969" # grey

    xAxisSeq <- seq(0,(max(depthsDF$depth)*100*2),tickDist)
    xAxisSeqLabels <- list()
    even <- TRUE
    for (i in xAxisSeq) {
        if(even) { # IF EVEN
            xAxisSeqLabels <- append(xAxisSeqLabels, c(i))
            even <- FALSE
        } else {
            xAxisSeqLabels <- append(xAxisSeqLabels, "")
            even <- TRUE
        }
    }

    depths_plot <- ggplot() +
        geom_vline(xintercept=100, linewidth=0.25, linetype="dotted") +
        geom_point(data=depthsDF,
                aes(y=depthsDF$gene_num[length(depthsDF$gene_num)]-gene_num,x=depth*100,color='Asm'),
                size=pointSize,
                shape=pointShape) +
        theme_classic(base_size = themeBaseSize) +
        theme(axis.text.y.left=element_text(size=axisTextSize),
              axis.text=element_text(color='black'),
              axis.ticks=element_line(color='black')) +
        scale_y_continuous(expand = c(0,1),
                        name="Gene",
                        labels=unique(depthsDF$gene),
                        breaks=seq(length(unique(depthsDF$gene))-1,0),
                        limits=c(0,NA)) +
        scale_x_continuous(expand = expansion(add=c(0,xExpansion)),
                        limits=c(0,NA),
                        breaks=xAxisSeq,
                        labels=xAxisSeqLabels,
                        name="Mean Gene Coverage / Assembly Coverage (%)") +
        scale_color_manual(name='Assemblies',
                        breaks=c('Asm'),
                        values=c('Asm'=pointColor)) +
        guides(color = "none")

    return(depths_plot)
}

# Create Plots
dataframe_combo <- read.table(inputPath_combo, header=TRUE)
dataframe_res   <- read.table(inputPath_res,   header=TRUE)
dataframe_col   <- read.table(inputPath_col,   header=TRUE)

# Helper function to pick a tick distance as a multiple of 25,
# if none found pick a tick distance as a multiple of 5,
# if none found pick a tick distance of 1.
# Considers desire for 20 ticks.
pickTickDist <- function(maxDepth) {
    tickDist <- round(maxDepth/20/25)*25 # pick a tick distance as a multiple of 25
    if (tickDist==0) {
        tickDist <- round(maxDepth/20/5)*5
        if (tickDist==0) {
            tickDist <- 1
        }
    }
    return(tickDist)
}

tickDist_combo <- pickTickDist(max(dataframe_combo$depth)*100)
tickDist_res   <- pickTickDist(max(dataframe_res$depth)*100)
tickDist_col   <- pickTickDist(max(dataframe_col$depth)*100)

depthPlot_combo <- geneDepthPlotter(dataframe_combo, tickDist_combo, themeBaseSize=themeBaseSize)
depthPlot_res   <- geneDepthPlotter(dataframe_res,   tickDist_res,   themeBaseSize=themeBaseSize)
depthPlot_col   <- geneDepthPlotter(dataframe_col,   tickDist_col,   themeBaseSize=themeBaseSize)

ggsave(depthPlot_combo,
       filename=paste(workflowDir,"/../results/G05_depthPlot",dataOutputSuffix,"_combo.pdf",sep=""),
       width=plotWidthInInches,
       height=length(dataframe_combo$gene_num)/5+1,
       units="in",
       limitsize=FALSE)
ggsave(depthPlot_res,
       filename=paste(workflowDir,"/../results/G05_depthPlot",dataOutputSuffix,"_res.pdf",sep=""),
       width=plotWidthInInches,
       height=length(dataframe_res$gene_num)/5+1,
       units="in",
       limitsize=FALSE)
ggsave(depthPlot_col,
       filename=paste(workflowDir,"/../results/G05_depthPlot",dataOutputSuffix,"_col.pdf",sep=""),
       width=plotWidthInInches,
       height=length(dataframe_col$gene_num)/5+1,
       units="in",
       limitsize=FALSE)

# Create Merged Plot
depthPlot_res_forMerged <- depthPlot_res +
    theme(axis.title.y=element_blank())
depthPlot_col_forMerged <- depthPlot_col +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank())

# Align left y axes of both plots
 depthPlot_col_forMerged_grob <- ggplotGrob(depthPlot_col_forMerged)
 depthPlot_res_forMerged_grob <- ggplotGrob(depthPlot_res_forMerged)
 maxWidth = grid::unit.pmax(depthPlot_col_forMerged_grob$widths[2:5], depthPlot_res_forMerged_grob$widths[2:5])
 depthPlot_col_forMerged_grob$widths[2:5] <- as.list(maxWidth)
 depthPlot_res_forMerged_grob$widths[2:5] <- as.list(maxWidth)

# Combine plots
layoutList <- c(rep(1,length(dataframe_col$depth)+2),rep(2,length(dataframe_res$depth)+2))
layout <- matrix(layoutList,ncol=1,byrow=TRUE)
depthPlot_merged <- grid.arrange(depthPlot_col_forMerged_grob,depthPlot_res_forMerged_grob,
                                 nrow=2,ncol=1,
                                 left=textGrob("Gene", gp=gpar(fontsize=themeBaseSize, col="black"), rot=90),
                                 layout_matrix=layout)

ggsave(depthPlot_merged,
       filename=paste(workflowDir,"/../results/G05_depthPlot",dataOutputSuffix,"_merged.pdf",sep=""),
       width=plotWidthInInches,
       height=length(dataframe_combo$gene_num)/5+1,
       units="in",
       limitsize=FALSE)