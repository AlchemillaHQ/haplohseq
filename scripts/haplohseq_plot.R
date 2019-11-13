# This plots haplohseq posterior probability data for hg19 runs.
# To support other genomes, the chromosome length definitions below
# need to be altered.

require("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="posterior output from haplohseq", metavar="character"),
  make_option(c("-c", "--chromosome"), type="character", default="all", 
              help="chromosome to plot; all chromosomes by default", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="output file prefix", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
posteriors = opt$file
outputDir = opt$out
outputPrefix = opt$prefix
chromosome = opt$chromosome

# hg19 chromosome lengths
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/index.shtml
chrLen = c()
chrLen["1"] <- 249250621
chrLen["2"] <- 243199373
chrLen["3"] <- 198022430
chrLen["4"] <- 191154276
chrLen["5"] <- 180915260
chrLen["6"] <- 171115067
chrLen["7"] <- 159138663
chrLen["8"] <- 146364022
chrLen["9"] <- 141213431
chrLen["10"] <- 135534747
chrLen["11"] <- 135006516
chrLen["12"] <- 133851895
chrLen["13"] <- 115169878
chrLen["14"] <- 107349540
chrLen["15"] <- 102531392
chrLen["16"] <- 90354753
chrLen["17"] <- 81195210
chrLen["18"] <- 78077248
chrLen["19"] <- 59128983
chrLen["20"] <- 63025520
chrLen["21"] <- 48129895
chrLen["22"] <- 51304566
chrLen["X"] <- 155270560
chrLen["Y"] <- 59373566

# For a given chromosome, this calculates a genomic start and end (for a plot with chromosomes appended to each other).
setPlotStartEnd <- function(chromosome) {
	plotStart <- c()
	plotEnd <- c()
	
	if (chromosome == "all") {
		plotStart["1"] <- 1
		plotEnd["1"] <- plotStart["1"] + chrLen["1"] - 1
		
		for (chr in 2:22) { 
			plotStart[as.character(chr)] <- plotEnd[as.character(chr-1)] + 1
			plotEnd[as.character(chr)] <- plotStart[as.character(chr)] + chrLen[as.character(chr)] - 1
		}
		plotStart["X"] <- plotEnd["22"] + 1
		plotEnd["X"] <- plotStart["X"] + chrLen["X"]
		plotStart["Y"] <- plotEnd["X"] + 1
		plotEnd["Y"] <- plotStart["Y"] + chrLen["Y"]
		
		maxPos <- plotEnd["22"]
		
	} 
	else {
		plotStart[chromosome] <- 1
		plotEnd[chromosome] <- chrLen[chromosome]
		maxPos <- plotEnd[chromosome]
	}
	
	list("maxPos" = maxPos, "plotStart" = plotStart, "plotEnd" = plotEnd)
}
 
chrNum <- function(chrName) {
	substring(as.character(chrName), if (nchar(as.character(chrName)) > 2 && substring(as.character(chrName),0,3) == 'chr') {4} else 1)
}

chrNums <- function(chrNames) {
	sapply(chrNames, chrNum)
}

# given a list of chromosomes, this returns a list of start positions (for an aggregate plot)
plotStarts <- function(chrNums) {
	starts <- rep(0, length(chrNums))
	for(i in 1:length(chrNums)) {
		starts[i] <- plotStart[chrNums[i]] 
		}
	return(starts)
}

# main plotting function
plotHaplohseq <- function(chromosome, posteriorFile, plotStart, plotEnd, maxPos, title, title2) {
	posteriorData <- read.table(posteriorFile, header=TRUE)	
	desc = tail(strsplit(posteriorFile,split="/")[[1]],n=1)		 
	if (chromosome != "all") {
		posteriorData <- posteriorData[chrNums(posteriorData$CHR)==chromosome,]
	}
	
	# plot RAFs
	hap1RefFreqs <- posteriorData[posteriorData$HAP=="1",]
	hap2RefFreqs <- posteriorData[posteriorData$HAP=="2",]
	plot((plotStarts(chrNums(hap1RefFreqs$CHR)) + hap1RefFreqs$POS), hap1RefFreqs$RAF, ylim=c(0,1), xlim=c(0,maxPos), col="#0000FF88", xlab="", ylab="event prob", xaxt="n" , main=title, cex=0.4)
	points((plotStarts(chrNums(hap2RefFreqs$CHR)) + hap2RefFreqs$POS), hap2RefFreqs$RAF, ylim=c(0,1), xlim=c(0,maxPos), col="#00FF0088", cex=0.4)

	# plot posterior probabilities for up to 2 event states
	lines(head((plotStarts(chrNums(posteriorData$CHR)) + posteriorData$POS), -1), head(posteriorData$S1, -1), col="#FF000088", lwd=3, cex=0.5)
	lines(head((plotStarts(chrNums(posteriorData$CHR)) + posteriorData$POS), -1), head(posteriorData$S2, -1), col="#0000FF", lwd=3)
	
	# add lines and text to demarcate chromosomes
	currentGenomePosition = 0
	abline(v=currentGenomePosition, lty=2, col="#AAAAAAFF")
	if (chromosome == "all") {
		for (chr in 1:22) {
			previousGenomePosition = currentGenomePosition
			abline(v=plotEnd[chr], lty=2, col="#AAAAAAFF")
			currentGenomePosition <- plotEnd[chr] + 1
			labelPosition = previousGenomePosition + (currentGenomePosition - previousGenomePosition)/2
			text(c(labelPosition),c(0.02),cex=1, labels=c(paste0(chr)), font=2)
		}
	}
}

# initialize supporting variables
params <- setPlotStartEnd(chromosome)
maxPos <- params$maxPos
plotStart <- params$plotStart
plotEnd <- params$plotEnd

# generate plot image
suffix = ""
if (chromosome != "all") {
	suffix = paste("_chr",chromosome,sep="")
}
png(paste(outputDir,"/",outputPrefix,suffix,".png",sep=""),width=12,height=2,units="in",res=400)
par(mar=c(0.5,1.75,1,1), oma=c(1.25,1,1,1), ps=8, mgp=c(1,0.2,0))
layout(matrix(c(1,1),1,1,byrow=TRUE),widths=c(1),heights=c(1))
plotHaplohseq(chromosome, posteriors, plotStart, plotEnd, maxPos, "", "")
dev.off()
cat(paste("haplohseq image generated: ", paste(outputDir,"/",outputPrefix,suffix,".png",sep=""), "\n", sep=""))


