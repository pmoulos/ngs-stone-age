# General peak correlation script to host several methods

source("/media/HD4/Fleming/dev/diagplots.R")
# input must be be at least bed files (chromosome, start, end)... The rest will be done
# based on the mapfile... Does not need IDs...
# output can be the correlation matrix in text format
correlatePeaks <- function(mapfile,avg.win=100,fdr.cut=1,method=c("pearson","spearman"),
							type=c("circle","ellipse","heatmap","cgram.shade","cgram.pts","simple","image"),
							output="x11",fil=NULL,...)
{
	if (missing(mapfile))
		stop("You must provide a file mapping between peak files and their tracks!")
  
	map <- read.delim(mapfile,stringsAsFactors=FALSE)
	if (is.null(map$control))
		hasControl <- FALSE
	else
		hasControl <- TRUE
  
	# Make sure we have more than one file!
	if (length(map$filename)==1)
		stop("You must provide more than one file for correlation calculation!")

	# Required packages
	if (!require(IRanges))
		stop("Bioconductor package IRanges is required!")
	if (!require(GenomicRanges))
		stop("Bioconductor package GenomicRanges is required!")
	if (!require(rtracklayer))
		stop("Bioconductor package rtracklayer is required!")

	# Get also destination directory
	path <- paste(dirname(mapfile),.Platform$file.sep,sep="")
	
	# Prepare output figure filenames
	if (!is.null(fil)) # We will auto-generate some names
		real.fil <- prepareAutoGenName(fil,output)

	# First, import the peak files and store them, as they are small
	input.bed <- ovObj.treat <- ovObj.control <- list()
	input <- map$filename
	for (i in 1:length(input))
	{
		cat("Importing peak file ",basename(input[i])," as bed... Please wait...\n",sep="")
		flush.console()
		tmp <- tempfile()
		#tmp.frame <- read.delim(input[i],header=FALSE,stringsAsFactors=FALSE,comment.char="#")
		tmp.frame <- read.delim(input[i],stringsAsFactors=FALSE,comment.char="#")
		#tmp.frame <- tmp.frame[-1,]
		tmp.frame <- cbind(tmp.frame[,c(1:5,7,9)],matrix("*",nrow(tmp.frame),1))
		#write.table(tmp.frame,tmp,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
		#input.bed[[basename(input[[i]])]] <- as.data.frame(import.bed(tmp,trackLine=FALSE,asRangedData=FALSE))
		input.bed[[basename(input[[i]])]] <- tmp.frame
		unlink(tmp)
	}
  
	# Now we must combine all the peaks
	peaks.pooled.matrix <- do.call("rbind",input.bed)
	peaks.pooled.df <- data.frame(chr=as.character(peaks.pooled.matrix[,1]),
								  start=as.numeric(peaks.pooled.matrix[,2]),
								  end=as.numeric(peaks.pooled.matrix[,3]),
								  id=paste(paste(peaks.pooled.matrix[,1],
								  peaks.pooled.matrix[,2],sep=":"),
								  peaks.pooled.matrix[,3],sep="-"),
								  length=as.character(peaks.pooled.matrix[,3] - peaks.pooled.matrix[,2] + 1),
								  summit=peaks.pooled.matrix[,5] + peaks.pooled.matrix[,2],
								  significance=peaks.pooled.matrix[,6],
								  fdr=peaks.pooled.matrix[,7],
								  strand=as.character(peaks.pooled.matrix[,8]))
	peaks.pooled.df <- peaks.pooled.df[-which(peaks.pooled.df$fdr<fdr.cut),]
	# Stop if no peaks left after FDR control
	if (any(dim(peaks.pooled.df)==0))
		stop("---===--- No peaks left after FDR control! No plots will be generated! ---===---")
	peaks.pooled.df <- peaks.pooled.df[order(peaks.pooled.df$chr,peaks.pooled.df$start),]
	pooled.peaks <- GRanges(seqnames=Rle(peaks.pooled.df$chr),
							IRanges(start=peaks.pooled.df$start,end=peaks.pooled.df$end),
							strand=Rle(peaks.pooled.df$strand))

	# Then, import and process the big mother bed files...
	for (i in 1:length(map$treatment))
	{
		treat <- map$treatment[i]
		cat("Importing track ",basename(treat),"... Please wait...\n",sep="")
		flush.console()
		treat.bed <- import.bed(treat,trackLine=FALSE,asRangedData=FALSE)

		if (hasControl)
		{
			control <- map$control[i]
			cat("Importing track ",basename(control),"... Please wait...\n",sep="")
			flush.console()
			control.bed <- import.bed(control,trackLine=FALSE,asRangedData=FALSE)
		}

		nams <- map$filename[which(map$treatment==treat)]
		for (nam in basename(nams))
		{
			cat("  Counting treatment tags in the imported treatment track ",basename(treat)," over the pooled peak set for ",basename(nam),"...\n",sep="")
			flush.console()
			ovObj.treat[[nam]] <- countOverlaps(pooled.peaks,treat.bed)
			if (hasControl)
			{
				cat("  Counting control tags in the imported control track ",basename(control)," over the pooled peak set for ",basename(nam),"...\n",sep="")
				flush.console()
				ovObj.control[[nam]] <- countOverlaps(pooled.peaks,control.bed)
			}
		}
		
		## Destroy the memory eating objects every time... Does not work for the time being...
		#suicide("treat.bed")
		#if (hasControl)
		#	suicide("control.bed")
		gc(verbose=FALSE)
	}

	# Now that we have the counts, we have to put them in the peak files
	treat.counts <- control.counts <- count.matrix <- avg.counts.treat <-
		avg.counts.control <- rank.sums <- ratios <- list()
	diffs <- peaks.pooled.df$end - peaks.pooled.df$start + 1
	for (nam in basename(input))
	{
		treat.counts[[nam]] <- unlist(ovObj.treat[[nam]])
		if (hasControl)
		{
			control.counts[[nam]] <- unlist(ovObj.control[[nam]])
			count.matrix[[nam]] <- cbind(treat.counts[[nam]],control.counts[[nam]])
			# Ratios
			ratios[[nam]] <-
				log2(ifelse(count.matrix[[nam]][,1]==0,1,count.matrix[[nam]][,1])/
					 ifelse(count.matrix[[nam]][,2]==0,1,count.matrix[[nam]][,2]))
			# Average tags
			avg.counts.treat[[nam]] <- count.matrix[[nam]][,1]/(diffs/avg.win)
			avg.counts.control[[nam]] <- count.matrix[[nam]][,2]/(diffs/avg.win)
		}
		else # Ratio cannot be calculated
		{
			count.matrix[[nam]] <- as.matrix(treat.counts[[nam]])
			avg.counts.treat[[nam]] <- count.matrix[[nam]][,1]/(diffs/avg.win)
			avg.counts.control[[nam]] <- NULL
			ratios[[nam]] <- NULL
		}
	}
  
	# Construct the additional stats matrix
	stats.matrix <- list()
	for (nam in basename(input))
	{
		stats.matrix[[nam]] <-
			cbind(
				count.matrix[[nam]],
				avg.counts.treat[[nam]],
				avg.counts.control[[nam]],
				ratios[[nam]])
		if (hasControl)
			colnames(stats.matrix[[nam]]) <- c("treat","control","avg.counts.treat","avg.counts.control","fe")
		else
			colnames(stats.matrix[[nam]]) <- c("treat","avg.counts.treat")
	}

	# Now, the matrices that are going to be used for correlations
	cor.list <- list()
	cor.list$treat <- cor.list$avg.counts.treat <- matrix(0,nrow(peaks.pooled.df),length(input))
	if (hasControl)
		cor.list$fe <- matrix(0,nrow(peaks.pooled.df),length(input))
	colnames(cor.list$treat) <- colnames(cor.list$avg.counts.treat) <- basename(input)
	if (hasControl)
		colnames(cor.list$fe) <- basename(input)
	for (nam in basename(input))
	{
		cor.list$treat[,nam] <- stats.matrix[[nam]][,"treat"]
		cor.list$avg.counts.treat[,nam] <- stats.matrix[[nam]][,"avg.counts.treat"]
		if (hasControl)
			cor.list$fe[,nam] <- stats.matrix[[nam]][,"fe"]
	}
	
	if ("circle" %in% type)
	{
		real.fil[3] <- "_TREATMENT_CIRCLES"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$treat,method="spearman",type="circle",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_CIRCLES"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$avg.counts.treat,method="spearman",type="circle",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_CIRCLES"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.cor(cor.list$fe,method="spearman",type="circle",output=output,fil=fig.name,...)
		}
	}
	
	if ("ellipse" %in% type)
	{
		real.fil[3] <- "_TREATMENT_ELLIPSE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$treat,method="spearman",type="ellipse",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_ELLIPSE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$avg.counts.treat,method="spearman",type="ellipse",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_ELLIPSE"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.cor(cor.list$fe,method="spearman",type="ellipse",output=output,fil=fig.name,...)
		}
	}
	
	if ("heatmap" %in% type)
	{
		real.fil[3] <- "_TREATMENT_HEATMAP"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$treat,method="spearman",type="heatmap",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_HEATMAP"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$avg.counts.treat,method="spearman",type="heatmap",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_HEATMAP"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.cor(cor.list$fe,method="spearman",type="heatmap",output=output,fil=fig.name,...)
		}
	}
	
	if ("cgram.shade" %in% type)
	{
		real.fil[3] <- "_TREATMENT_CGRAMSHADE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$treat,method="spearman",type="cgram.shade",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_CGRAMSHADE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$avg.counts.treat,method="spearman",type="cgram.shade",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_CGRAMSHADE"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.cor(cor.list$fe,method="spearman",type="cgram.shade",output=output,fil=fig.name,...)
		}
	}
	
	if ("cgram.pts" %in% type)
	{
		real.fil[3] <- "_TREATMENT_CGRAMPTS"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$treat,method="spearman",type="cgram.pts",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_CGRAMPTS"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.cor(cor.list$avg.counts.treat,method="spearman",type="cgram.pts",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_CGRAMPTS"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.cor(cor.list$fe,method="spearman",type="cgram.pts",output=output,fil=fig.name,...)
		}
	}
	
	if ("simple" %in% type)
	{
		real.fil[3] <- "_TREATMENT_SIMPLE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.pairs(cor.list$treat,type="simple",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_SIMPLE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.pairs(cor.list$avg.counts.treat,type="simple",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_SIMPLE"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.pairs(cor.list$fe,type="simple",output=output,fil=fig.name,...)
		}
	}
	
	if ("image" %in% type)
	{
		real.fil[3] <- "_TREATMENT_IMAGE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.pairs(cor.list$treat,type="image",output=output,fil=fig.name,...)

		real.fil[3] <- "_AVGCOUNT_IMAGE"
		fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
							sep=.Platform$file.sep)
		plot.pairs(cor.list$avg.counts.treat,type="image",output=output,fil=fig.name,...)
		
		if (hasControl)
		{
			real.fil[3] <- "_FE_IMAGE"
			fig.name <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
								sep=.Platform$file.sep)
			plot.pairs(cor.list$fe,type="image",output=output,fil=fig.name,...)
		}
	}
}

prepareAutoGenName <- function(fil,output)
{
	real.fil <- character(4)
	real.fil[1] <- dirname(fil)
	real.fil[2] <- sub("^([^.]*).*","\\1",basename(fil))
	#real.fil[3] remains empty to be filled with the plot type
	real.fil[4] <- paste(".",output,sep="")
	return(real.fil)
}

suicide <- function(...) { rm(...,envir=.GlobalEnv) }
