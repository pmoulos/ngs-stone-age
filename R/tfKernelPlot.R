# This R script creates a combined kernel density plot that characterizes the
# binding profile of a transcription factor. At the same time, given a file with
# gene expression from a respective knockout/down experiment, it associates the
# transcription factor binding events with the up and down regulated genes and
# overalays the respective kernel density plots.
#
# Usage:
# kdp <- tfKernelPlot(
#           signalFiles,
#           geneFile,
#           peakFiles,
#           statIndex=c(9,12),
#           signalType="bed",
#           promoter=c(-10000,1000),
#           calc="rpm",
#           fragSize=200,
#           lowCount=5,
#			minGeneSize=1000,
#           fcut=c(-1,1),
#           pcut=0.05,
#			overlay=FALSE,
#           xLim=NULL,
#           output=NULL,
#           device="x11",
#			plot=TRUE
#       )
#
# Parameters:
# signalFiles  : A set of mapped read files in BAM/BED format, where their input 
#                order corresponds to the input order of peakFiles, if provided
#                (e.g. if peakFiles = c("A.peak","B.peak"), 
#                signalFiles = c("A.bam","B.bam"), where A and B represent the
#                order and not the naming of the files, e.g. signalFile = 
#                c("B.bam","A.bam") is WRONG). This parameter mandatory.
# geneFile     : One file containing the organism's annotation in BED format 
#                WITH a header line. The columns should be as follows:
#                (chromosome, start, end, gene_id, something, strand). All 
#                additional columns can contain whatever other information
#                (e.g. expression values, additional annotation). At least TWO
#                additional columns must be provided: a column with a
#                statistical score (e.g. a p-value) and a column with a fold
#                change regarding each gene and its control (e.g. KO/WT). If
#                statIndex is null, then these two columns are not necessary.
#				 This argument can also be a data frame instead of a file.
# peakFiles    : Optionally, a set of peak files for each TF, usually the result
#                of a peak calling algorithm, in text tab delimited format WITH 
#                a header line. The tab delimited format should be similar to 
#                the UCSC BED format, where the first four columns are 
#                (chromosome, start, end, peak_id). The rest of the columns can 
#                contain whatever additional information. If provided, then the
#				 kernel density calculations will be perofmed only using reads
#				 in the promoters including these peak regions.
# statIndex    : A vector of two integers: the first one is the column number
#                with the statistical score and the second is the column with 
#                the fold change cutoff in log2 (important!). For a metaseqR
#                typical output, statIndex = c(9,12). A simple check is
#                performed to check if the correct columns have been provided
#                based on the content. If statIndex is NULL (default), then
#                no association with up and down genes is made and only the
#                kernel density of the binding in the promoters is plotted.
# signalType:  : The file type of the signal files. It can be "bed" (default)
#                or "bam".
# promoter     : The promoter region that will be used for the calculation of 
#                binding signal for the genes. It defaults to -10kb upstream 
#                and +1kb downstream (promoter = c(-10000,1000)).
# calc         : How to calculate the occupancy in the promoter regions. It can
#                be "rpm" (default) for reads per million (normalized by the 
#                length of the promoter region) or "coverage" for using the
#                actual read coverage in the promoter region.
# fragSize     : The original estimated fragment size after sonication 
#                (defaults to 200bp). The reads from the signal files will be
#                extended to this size.
# lowCount     : A simple threshold to exlcude binding regions with very low
#                signal in terms of read counts. Defaults to 5.
# minGeneSize  : A simple threshold to exlcude noisy genes/transcrpts (usually
#				 very short ones) from the calculations. Defaults to 1000 to
#				 remove genes smaller than 1kb. Set to NULL for not using this
#				 filter.
#                signal in terms of read counts. Defaults to 5.
# fcut         : A vector of length two which contains the fold change cutoffs 
#                which together with the statistical score cutoff define the 
#                differentially expressed genes. It must be in log2 scale and 
#                the first element is the lower cutoff and the second element is
#                the upper cutoff. E.g. fcut = c(-1,1) (default).
# pcut         : A number between 0 and 1 which defines the statistical score
#                (e.g. p-value) which together with the fold change cutoff 
#                defines the differentially expressed genes. The default is
#                0.05.
# overlay	   : If TRUE, then one plot is created with all kernel density
#				 profiles (as many as signalFiles), overlayed. The default is
#				 FALSE, thus one plot per signal file is created.
# xLim		   : A desired x-axis limit for the final plot. Defaults to NULL to
#				 let ggplot do the work but useful for later figure refinement.
# output       : The output file names to name the figures. If they are not 
#                given (NA) or they do not meet criteria (e.g. they must be
#                equal to the number of input TFs), they will be autogenerated.
# device       : One of the supported R graphics devices for the output figure
#                file. The default is "x11" to pop-up the figure. If "pdf", the
#                figure will be exported to a pdf file (possibly given in 
#                output).
# plot		   : If TRUE (default), the output kernel density plots will be
#				 generated. If FALSE, only the data will be generated. Use the
#				 output list to create custom plots.
# rc           : Fraction of available cores to use for parallel calculations.
#                Defaults to NULL (one core).
#
# The kdp variable is a list with values required to reconstruct the kernel
# density profiles. They can be used e.g. for later clustering of the profiles.
#
# Version: 0.11

# TODO: Separate plotting from data generation
# TODO: Implement overlay option

tfKernelPlot <- function(signalFiles,geneFile,peakFiles=NULL,statIndex=NULL,
    fragSize=200,signalType="bed",promoter=c(-10000,1000),xLim=NULL,
    calc=c("rpm","coverage"),lowCount=5,minGeneSize=1000,fcut=c(-1,1),pcut=0.05,
    overlay=FALSE,output=NULL,device=c("x11","png","jpg","tiff","bmp","pdf",
    "ps"),plot=TRUE,rc=NULL) {
    if (missing(signalFiles) || missing(geneFile))
        stop("You must provide both signalFiles and geneFile arguments!")
    if (!is.data.frame(geneFile) && !file.exists(geneFile))
		stop("geneFile must be an existing file or a data frame!")
    if (any(!file.exists(signalFiles))) {
		m <- paste(signalFiles[which(!file.exists(signalFiles))],sep=", ")
		stop("Signal files ",m," are not found!")
	}
    if (!is.null(peakFiles) && length(peakFiles)!=length(signalFiles))
		stop("Every signal file must be accompanied by a peak file when ",
			"peakFiles argument is provided!")
	if (!is.null(peakFiles) && any(!file.exists(peakFiles))) {
		mm <- paste(peakFiles[which(!file.exists(peakFiles))],sep=", ")
		stop("Peak files ",mm," are not found!")
	}
    signalType = signalType[1]
    calc = calc[1]
    if (!signalType %in% c("bed","bam"))
        stop("signalType must be one of \"bed\" or \"bam\"!")
    device <- device[1]
    if (!device %in% c("x11","png","jpg","tiff","bmp","pdf","ps"))
        stop("The device argument must be one of R supported graphics types.")
    if (!is.null(statIndex) 
        && (any(is.na(suppressWarnings(as.integer(statIndex)))) 
        || length(statIndex) != 2))
        stop("The statIndex argument must be an integer vector of lemngth 2!")
    if (any(is.na(suppressWarnings(as.integer(promoter)))) 
        || length(promoter) != 2)
        stop("The promoter argument must be an integer vector of lemngth 2!")
    if (promoter[1]>promoter[2])
		stop("The promoter argument must be sorted (promoter[1]<=promoter[2])!")
    if (is.na(suppressWarnings(as.integer(fragSize))))
        stop("The fragSize argument must be a positive integer!")
    if (!calc %in% c("rpm","coverage"))
        stop("The calc argument must be one of \"rpm\" or \"coverage\".")
    if (calc=="rpm")
        if (is.na(suppressWarnings(as.integer(lowCount))))
            stop("The lowCount argument must be a positive integer!")
    else if (calc=="coverage")
        if (is.na(suppressWarnings(as.numeric(lowCount))))
            stop("The lowCount argument must be a positive number!")
    if (!is.null(minGeneSize) 
		&& is.na(suppressWarnings(as.integer(minGeneSize))))
            stop("The lowCount argument must be a positive integer or NULL!")
    if (any(is.na(suppressWarnings(as.numeric(fcut)))) 
        || length(fcut) != 2)
        stop("The fcut argument must be a numeric vector of lemngth 2!")
    if (is.na(suppressWarnings(as.numeric(fcut))) || pcut<0 || pcut>1)
        stop("The pcut argument must be a number between 0 and 1!")
    if (!is.null(xLim) && (!all(is.numeric(xLim)) || length(xLim)!=2))
		stop("The xLim argument must be a numeric vector of length 2 if given!")
    
    # Check packages
    if (!require(ggplot2))
        stop("R package ggplot2 is required!")
    if (!require(GenomicRanges))
        stop("Bioconductor package IRanges is required!")
    if (signalType=="bed" && !require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")
    if (signalType=="bam" && !require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required!")    
    
    # Check requested output files
    if (!is.null(peakFiles)) {
		nameFiles <- names(peakFiles)
		if (is.null(nameFiles))
			nameFiles <- peakFiles
	}
	else {
		nameFiles <- names(signalFiles)
		if (is.null(nameFiles))
			nameFiles <- signalFiles
	}
	if (is.null(output) && device != "x11") {
		output <- character(length(nameFiles))
		for (i in 1:length(nameFiles))
			output[i] <- file.path(dirname(nameFiles[i]),
				paste("kerden_out_",sub("^([^.]*).*","\\1",
				basename(nameFiles[i])),".",device,sep=""))
	}
    if (!is.null(output) & length(output)!=length(nameFiles)) {
        warning("The number of output figures is different than the number of ",
            "the input peak files! They will be autogenerated...",
            call.=FALSE,immediate=TRUE)
        output <- character(length(nameFiles))
        for (i in 1:length(nameFiles))
            output[i] <- file.path(dirname(nameFiles[i]),
                paste("kerden_out_",sub("^([^.]*).*","\\1",
                basename(nameFiles[i])),".",device,sep=""))
    }
    
    # Initiate an object to write the profiles that will be generated for
    # general binding as well as gene-specific binding profile
    out <- vector("list",length(nameFiles))
    names(out) <- sub("^([^.]*).*","\\1",basename(nameFiles))
    
    # Read gene file
    if (!is.data.frame(geneFile)) {
		message("\nReading genes file ",basename(geneFile),
			" to GenomicRanges object...")
		geneData <- read.delim(geneFile)
		rownames(geneData) <- as.character(geneData[,4])
	}
	else {
		geneData <- geneFile
		rownames(geneData) <- as.character(geneData[,4])
	}
    
    # Get and check stat scores
    p <- f <- NA
    if (!is.null(statIndex)) {
        p <- as.numeric(geneData[,statIndex[1]])
        f <- as.numeric(geneData[,statIndex[2]])
        if (any(p[!is.na(p)]>1) || any(p[!is.na(p)]<0)) # Might include NA
            stop("Negative or greater than one values found in the ",
                "statistical score column! Maybe you accidentally gave them ",
                "in wrong order?")
        if (all(f>0))
            warning("All fold change values were found positive! Are you sure ",
                "that all your genes are upregulated? Maybe the fold change ",
                "is not in log2 scale. Please review.",immediate.=TRUE)
        names(p) <- names(f) <- rownames(geneData)
    }
    
    # Continue with Genomic Ranges operations
    genesGr <- makeGRangesFromDataFrame(
        df=geneData,
        keep.extra.columns=TRUE,
        seqnames.field="chromosome"
    )
    names(genesGr) <- rownames(geneData)
    
    # Remove small genes if wished
    if (!is.null(minGeneSize)) {
		message("Excluding genes less than ",minGeneSize," base pairs")
		s <- which(width(genesGr)<minGeneSize)
		if (length(s)>0)
			genesGr <- genesGr[-s]
	}
    
    # Get promoters and fix broken intervals
    if (promoter[1]<=0 && promoter[2]>=0)
		ranges(genesGr) <- promoters(ranges(genesGr),upstream=abs(promoter[1]),
			downstream=abs(promoter[2]))
	else if (promoter[1]<=0 && promoter[2]<=0) {
		ranges(genesGr) <- promoters(ranges(genesGr),upstream=abs(promoter[1]),
			downstream=0)
		ranges(genesGr) <- resize(ranges(genesGr),
			width=abs(promoter[1])-abs(promoter[2]))
	}
	else if (promoter[1]>=0 && promoter[2]>=0) {
		ranges(genesGr) <- promoters(ranges(genesGr),upstream=0,
			downstream=promoter[1])
		ranges(genesGr) <- resize(ranges(genesGr),
			width=promoter[2]-promoter[1])
	}
	ranges(genesGr) <- restrict(ranges(genesGr),start=1)

    # Do job with peak files
    for (i in 1:length(signalFiles)) {
        message("Reading signal file ",basename(signalFiles[i]),"...")
        if (signalType=="bed")
            signalGr <- import.bed(signalFiles[i],trackLine=FALSE)
        else if (signalType=="bam")
            signalGr <- trim(as(readGAlignments(file=signalFiles[i]),"GRanges"))
        message("  Extending reads to ",fragSize,"bp...")
        signalGr <- trim(resize(signalGr,fragSize,fix="start"))
        
        if (!is.null(peakFiles)) {
			message("  Processing peaks file ",basename(peakFiles[i]),"...")
			message("    Reading file to GenomicRanges object...")
			peakData <- read.delim(peakFiles[i])
			peaksGr <- makeGRangesFromDataFrame(
				df=peakData,
				keep.extra.columns=TRUE,
				seqnames.field="chromosome"
			)
			names(peaksGr) <- as.character(peakData[,4])
        
			# Get the promoters with peaks to use for coverage calculation
			message("   Getting promoters containing actual peaks...")
			pinp <- findOverlaps(genesGr,peaksGr)
			promsWithPeaks <- GRanges(
				seqnames=seqnames(genesGr)[queryHits(pinp)],
				ranges=IRanges(
					start=start(genesGr)[queryHits(pinp)],
					end=end(genesGr)[queryHits(pinp)]
				),
				strand=strand(genesGr)[queryHits(pinp)],
				peak=names(peaksGr)[subjectHits(pinp)],
				gene=names(genesGr)[queryHits(pinp)]
			)
			names(promsWithPeaks) <- as.character(promsWithPeaks$gene)
			peakPresence <- rep("-",length(genesGr))
			names(peakPresence) <- names(genesGr)
			peakPresence[unique(names(promsWithPeaks))] <- "+"
			
			message("  Calculating binding strength...")
			if (calc=="rpm") {
				promReads <- summarizeOverlaps(promsWithPeaks,signalGr,
				    singleEnd=TRUE,ignore.strand=TRUE)
				promReads <- assays(promReads)$counts[,1]
				#names(promReads) <- names(promsWithPeaks)
				if (!is.na(lowCount) && lowCount>0) {
					#filt <- which(apply(promReads,1,function(x) 
					#	all(x<lowCount)))
					filt <- which(promReads<lowCount)
					if (length(filt)>0)
						#promReads <- promReads[-filt,]
						promReads <- promReads[-filt]
				}
				promReads <- (promReads/diff(promoter))*1e+6/length(signalGr)
			}
			else if (calc=="coverage") {
				promCov <- calcCoverage(signalGr,promsWithPeaks,rc=rc)
				promReads <- sapply(promCov,mean)
				if (!is.na(lowCount) && lowCount>0) {
					filt <- which(promReads<lowCount)
					if (length(filt)>0)
						promReads <- promReads[-filt]
				}
			}
        }
        else {
			if (calc=="rpm") {
				promReads <- summarizeOverlaps(genesGr,signalGr,singleEnd=TRUE,
					ignore.strand=TRUE)
				promReads <- assays(promReads)$counts[,1]
				#names(promReads) <- names(genesGr)
				if (!is.na(lowCount) && lowCount>0) {
					#filt <- which(apply(promReads,1,
					#	function(x) all(x<lowCount)))
					filt <- which(promReads<lowCount)
					if (length(filt)>0)
						#promReads <- promReads[-filt,]
						promReads <- promReads[-filt]
				}
				promReads <- (promReads/diff(promoter))*1e+6/length(signalGr)
			}
			else if (calc=="coverage") {
				promCov <- calcCoverage(signalGr,genesGr,rc=rc)
				promReads <- sapply(promCov,mean)
				if (!is.na(lowCount) && lowCount>0) {
					filt <- which(promReads<lowCount)
					if (length(filt)>0)
						promReads <- promReads[-filt]
				}
			}
		}
        
        # Data frame with rpkms
        allData <- data.frame(
            Density=log(promReads),
            Source=rep(nameFiles[i],length(promReads))
        )

        # Get up and down genes
        if (!is.null(statIndex)) {
			message("  Separating gene expression profiles...")
            upBound <- which(!is.na(p) & p<pcut & f>=fcut[2] 
                & peakPresence=="+")
            if (length(upBound) > 1) {
                upData <- allData[upBound,,drop=FALSE]
                upNa <- which(is.na(upData[,1]))
                if (length(upNa) > 0)
                    upData <- upData[-which(is.na(upData[,1])),,drop=FALSE]
            }
            else
                upData <- NULL
            downBound <- which(!is.na(p) & p<pcut & f<=fcut[1] 
                & peakPresence=="+")
            if (length(downBound) > 1) {
                downData <- allData[downBound,,drop=FALSE]
                downNa <- which(is.na(downData[,1]))
                if (length(downNa) > 0)
                    downData <- downData[-which(is.na(downData[,1])),,
                        drop=FALSE]
            }
            else
                downData <- NULL
        }
        else
            upData <- downData <- NULL
        
        # Now that we have used the f and p, we can remove NAs from allData
        # (if any)
        allNa <- which(is.na(allData[,1]))
        if (length(allNa) > 0)
            allData <- allData[-which(is.na(allData[,1])),]
        
        # Construct the final data frame
        theData <- allData
        theData <- 
            cbind(theData,rep(paste(nameFiles[i],"occupancy"),nrow(theData)))
        names(theData)[3] <- "Status"
        theData$Status <- as.character(theData$Status)
        if (!is.null(statIndex)) {
            theData[rownames(upData),3] <- paste(nameFiles[i],"+ Up",sep="")
            theData[rownames(downData),3] <- paste(nameFiles[i],"+ Down",sep="")
        }
        
        message("  Generating kernel density profiles...")
        if (!is.null(statIndex)) {
            if (!is.null(upData) && !is.null(downData)) {
				if (plot) {
					thePlot <-
						ggplot(theData,mapping=aes(x=Density)) + 
						geom_line(aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						geom_line(data=theData[which(theData$Status %in% c(
							paste(nameFiles[i],"+ Up",sep=""),
							paste(nameFiles[i],"+ Down",sep="")
						)),],aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						theme_bw() +
						ggtitle(paste(nameFiles[i],
							"normalized read density in promoters\n")) +
						xlab(paste("\n",nameFiles[i],
							"normalized read density (ln)")) +
						ylab("Kernel density\n") +
						theme(
							title=element_text(size=14),
							axis.title.x=element_text(size=12,face="bold"),
							axis.title.y=element_text(size=12,face="bold"),
							axis.text.x=element_text(size=12,face="bold"),
							axis.text.y=element_text(size=12,face="bold"),
							panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							legend.position="bottom",
							legend.key=element_blank()
						) + 
						scale_fill_manual(
							values=c("#00B400","#AEAEAE","#B40000")) +
						scale_color_manual(
							values=c("#00B400","#AEAEAE","#B40000"))
				}
				else
					thePlot <- NULL
                    
                # Store the output to out
                out[[i]]$total <- density(allData[,1])
                out[[i]]$up <- density(upData[,1])
                out[[i]]$down <- density(downData[,1])
                out[[i]]$plotData <- theData
                out[[i]]$plot <- thePlot
            }
            else if (!is.null(upData) && is.null(downData)) {
				if (plot) {
					thePlot <-
						ggplot(theData,mapping=aes(x=Density)) + 
						geom_line(aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						geom_line(data=theData[which(theData$Status %in% c(
							paste(nameFiles[i],"+ Up",sep="")
						)),],aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						theme_bw() +
						ggtitle(paste(nameFiles[i],
							"normalized read density in promoters\n")) +
						xlab(paste("\n",nameFiles[i],
							"normalized read density (ln)")) +
						ylab("Kernel density\n") +
						theme(
							title=element_text(size=14),
							axis.title.x=element_text(size=12,face="bold"),
							axis.title.y=element_text(size=12,face="bold"),
							axis.text.x=element_text(size=12,face="bold"),
							axis.text.y=element_text(size=12,face="bold"),
							panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							legend.position="bottom",
							legend.key=element_blank()
						) + 
						scale_fill_manual(values=c("#AEAEAE","#B40000")) +
						scale_color_manual(values=c("#AEAEAE","#B40000"))
				}
				else
					thePlot <- NULL
                    
                # Store the output to out
                out[[i]]$total <- density(allData[,1])
                out[[i]]$up <- density(upData[,1])
                out[[i]]$down <- NULL
                out[[i]]$plotData <- theData
                out[[i]]$plot <- thePlot
            }
            else if (is.null(upData) && !is.null(downData)) {
				if (plot) {
					thePlot <-
						ggplot(theData,mapping=aes(x=Density)) + 
						geom_line(aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						geom_line(data=theData[which(theData$Status %in% c(
							paste(nameFiles[i],"+ Down",sep="")
						)),],aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						theme_bw() +
						ggtitle(paste(nameFiles[i],
							"normalized read density in promoters\n")) +
						xlab(paste("\n",nameFiles[i],
							"normalized read density (ln)")) +
						ylab("Kernel density\n") +
						theme(
							title=element_text(size=14),
							axis.title.x=element_text(size=12,face="bold"),
							axis.title.y=element_text(size=12,face="bold"),
							axis.text.x=element_text(size=12,face="bold"),
							axis.text.y=element_text(size=12,face="bold"),
							panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							legend.position="bottom",
							legend.key=element_blank()
						) + 
						scale_fill_manual(values=c("#00B400","#AEAEAE")) +
						scale_color_manual(values=c("#00B400","#AEAEAE"))
				}
				else
					thePlot <- NULL
                    
                # Store the output to out
                out[[i]]$total <- density(allData[,1])
                out[[i]]$up <- NULL
                out[[i]]$down <- density(downData[,1])
                out[[i]]$plotData <- theData
                out[[i]]$plot <- thePlot
            }
            else if (is.null(upData) && is.null(downData)) {
				if (plot) {
					thePlot <-
						ggplot(theData,mapping=aes(x=Density)) + 
						geom_line(aes(colour=Status),stat="density",
							kernel="gaussian",size=1) +
						theme_bw() +
						ggtitle(paste(nameFiles[i],
							"normalized read density in promoters\n")) +
						xlab(paste("\n",nameFiles[i],
							"normalized read density (ln)")) +
						ylab("Kernel density\n") +
						theme(
							title=element_text(size=14),
							axis.title.x=element_text(size=12,face="bold"),
							axis.title.y=element_text(size=12,face="bold"),
							axis.text.x=element_text(size=12,face="bold"),
							axis.text.y=element_text(size=12,face="bold"),
							panel.grid.major=element_blank(),
							panel.grid.minor=element_blank(),
							legend.position="bottom",
							legend.key=element_blank()
						) + 
						scale_fill_manual(values=c("#AEAEAE")) +
						scale_color_manual(values=c("#AEAEAE"))
				}
				else
					thePlot <- NULL
                    
                out[[i]]$total <- density(allData[,1])
                out[[i]]$up <- NULL
                out[[i]]$down <- NULL
                out[[i]]$plotData <- theData
                out[[i]]$plot <- thePlot
            }
        }
        else {
			if (plot) {
				thePlot <-
					ggplot(theData,mapping=aes(x=Density)) + 
					geom_line(aes(colour=Status),stat="density",
						kernel="gaussian",size=1) +
					#geom_density(aes(colour=Status),kernel="gaussian",size=1) +
					#geom_hline(yintercept=0,colour="white",size=1) +
					theme_bw() +
					ggtitle(paste(nameFiles[i],
						"normalized read density in promoters\n")) +
					xlab(paste("\n",nameFiles[i],
						"normalized read density (ln)")) +
					ylab("Kernel density\n") +
					theme(
						title=element_text(size=14),
						axis.title.x=element_text(size=12,face="bold"),
						axis.title.y=element_text(size=12,face="bold"),
						axis.text.x=element_text(size=12,face="bold"),
						axis.text.y=element_text(size=12,face="bold"),
						panel.grid.major=element_blank(),
						panel.grid.minor=element_blank(),
						legend.position="bottom",
						legend.key=element_blank()
					) + 
					scale_fill_manual(values=c("#AEAEAE")) +
					scale_color_manual(values=c("#AEAEAE"))
			}
			else
				thePlot <- NULL
                    
            out[[i]]$total <- density(allData[,1])
            out[[i]]$up <- NULL
            out[[i]]$down <- NULL
            out[[i]]$plotData <- theData
            out[[i]]$plot <- thePlot
        }
        
        if (plot) {
			if (!is.null(xLim))
				thePlot <- thePlot + xlim(xLim[1],xLim[2])

			if (device != "x11")
				ggsave(filename=output[i],plot=thePlot,width=7,height=6)
			else
				print(thePlot)
		}
    }
    
    message("Done!")
    
    return(out)
}

calcCoverage <- function(input,mask,strand=NULL,ignoreStrand=TRUE,rc=NULL) {
    if (!is(input,"GRanges") && !is.list(input))
        stop("The input argument must be a GenomicRanges object or or a list ",
            "of GenomicRanges")
    if (!is(mask,"GRanges") && !is(mask,"GRangesList"))
        stop("The mask argument must be a GRanges or GRangesList object")
    if (!is.null(strand) && !is.list(strand)) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    input <- splitBySeqname(input)
    if (is(mask,"GRanges"))
        index <- split(1:length(mask),as.character(seqnames(mask)))
    cov <- unlist(lapply(index,coverageFromRanges,mask,input,
            ignoreStrand,rc=rc))
    theNames <- 
        unlist(lapply(strsplit(names(cov),"\\."),
            function(x) return(x[length(x)])))
    names(cov) <- theNames
    gc(verbose=FALSE)
    return(cov)
}

coverageFromRanges <- function(i,mask,input,ignore.strand,rc=NULL) {
    x <- mask[i]    
    y<-list(
        chromosome=as.character(seqnames(x))[1], 
        start=start(x),
        end=end(x),
        strand=as.character(strand(x)),
        reads=NULL,
        coverage=NULL
    )
    message("  processing ",y$chromosome)
    if (!is.null(input[[y$chromosome]])) {
        y$reads <- input[[y$chromosome]][
            unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
                ignore.strand=ignore.strand)))]
    }
    else {
        message(y$chromosome," not found!")
        y$reads <- NULL
    }
    if (length(y$reads)>0) {
        xCov <- coverage(x)
        cc <- as.character(seqnames(y$reads))[1]
        y$coverage <- coverage(y$reads)
        y$coverage <- y$coverage[[cc]]
        xCov <- xCov[[cc]]
        covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
            tryCatch({
                if (d[j] == "+")
                    return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
                else if (d[j] == "-")
                    return(rev(r[s[j]:e[j]]/xCov[s[j]:e[j]]))
                else
                    return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
            },
            error=function(e) {
                message("Caught invalid genomic area!")
                print(mask[i][j])
                message("Will return zero coverage")
                return(Rle(NA))
                #return(NA)
            },finally={})
        },y$start,y$end,y$strand,y$coverage,rc=rc)
        names(covs) <- names(x)
        return(covs)
    }
    else
        return(Rle(NA))
}

splitBySeqname <- function(gr,rc=NULL) {
    gr.list <- cmclapply(levels(seqnames(gr)),function(x,lib) {
        tmp <- lib[seqnames(lib)==x]
        if (length(tmp)>0) return(tmp) else return(NULL)
    },gr,rc=rc)
    names(gr.list) <- levels(seqnames(gr))
    null <- which(sapply(gr.list,is.null))
    if (length(null)>0)
        gr.list <- gr.list[-null]
    return(gr.list)
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}
