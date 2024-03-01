#!/usr/bin/env Rscript

##### File names from command line #####
library(argparse)

parser <- ArgumentParser(description= 'Peak calling with CAGEfightR')
parser$add_argument('--forward', required=TRUE, 
                    help='Names of forward bigwig files separated by commas [Example: starts1.fwd.bw,starts2.fwd.bw,startsN.fwd.bw]')
parser$add_argument('--reverse', required=TRUE, 
                    help='Names of reverse bigwig files separated by commas [Example: starts1.rev.bw,starts2.rev.bw,startsN.rev.bw]')
parser$add_argument('--chromosome', required=TRUE, 
                    help='Name of the file with chromosome length')
parser$add_argument('--output', required=TRUE, 
                    help='Name of the output BED file')
args <- parser$parse_args()

forward.files <- unlist(strsplit(args$forward, split = ","))
reverse.files <- unlist(strsplit(args$reverse, split = ","))

##### Genome processing #####
setClass("Seqinfo", representation(seqnames = "character",
                                   seqlengths = "integer",
                                   is_circular = "logical",
                                   genome = "character"))

chr <- read.table(args$chromosome,
                         header=FALSE, sep="\t")

genomeInfo <- new("Seqinfo", seqnames = as.character(chr[,1]),
                  seqlengths = as.integer(chr[,2]),
                  is_circular = rep(FALSE, nrow(chr)),
                  genome = rep("our_genome", nrow(chr)))

##### Peak calling ##### 
library(CAGEfightR)

nano_plus <- BigWigFileList(forward.files)
nano_minus <- BigWigFileList(reverse.files)
rep.names <- paste("rep", seq(1, length(forward.files)), sep="")
names(nano_plus) <- rep.names
names(nano_minus) <- rep.names

nano.design <- data.frame(Name = rep.names,
                          BigWigPlus = c("nano_plus"),
                          BigWigMinus = c("nano_minus"))
rownames(nano.design) <- rep.names

nanoCTSSs <- quantifyCTSSs(plusStrand = nano_plus,
                           minusStrand = nano_minus,
                           genome = genomeInfo,
                           design = nano.design)

nanoCTSSs <- calcPooled(nanoCTSSs, inputAssay="counts")
filtnanoCTSSs <- calcSupport(nanoCTSSs, inputAssay="counts",
                             outputColumn="support",
                             unexpressed=0)
supportednanoCTSSs <- subset(filtnanoCTSSs, support > 1)
supportednanoCTSSs <- calcPooled(supportednanoCTSSs, inputAssay="counts")
prefiltered_nanoTCs <- clusterUnidirectionally(supportednanoCTSSs,
                                               pooledCutoff = 2,
                                               mergeDist = 10)
export(prefiltered_nanoTCs,
       args$output)