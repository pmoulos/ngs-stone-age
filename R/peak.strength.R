peak.strength <- function(mapfile,type="boxplot")
{
	if (missing(mapfile))
        stop("You must provide a data mapping file for the script to work!")
	else
	{
		map <- read.delim(mapfile,stringsAsFactors=FALSE)
		map <- parse.map(map)
	}

	if (!require(GenomicRanges))
        stop("Bioconductor package IRanges is required!")
    if (!require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")

	# First, import the peak files and store them, as they are small
	region.bed <- vector("list",length(map$region))
	names(region.bed) <- names(map$region)
	for (n in names(region.bed))
	{
		cat("Importing peak file ",basename(map$region[n])," as bed... Please wait...\n",sep="")
		flush.console()
		tmp <- tempfile()
		tmp.reg <- read.delim(map$region[n],header=FALSE)
		tmp.frame <- cbind(tmp.reg[,1:3],paste(tmp.reg[,1],paste(tmp.reg[,2],tmp.reg[,3],sep="-"),sep=":"),rep(0,nrow(tmp.reg)),rep("*",nrow(tmp.reg)))
		write.table(tmp.frame,tmp,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
		region.bed[[n]] <- import.bed(tmp,trackLine=FALSE,asRangedData=FALSE)
		unlink(tmp)
	}

	count.mats <- init.count.matrix(map$data,region.bed)
	
    # Import tracks and process each region
    for (t in map$udata)
    {
		cat("Importing track ",basename(t),"... Please wait...\n",sep="")
        flush.console()
        bed <- import.bed(t,trackLine=FALSE,asRangedData=FALSE)

		# Which regions have reads in this file?
		for (n in names(map$data))
		{
			ix <- grep(t,map$data[[n]],fixed=TRUE)
			if (length(ix) != 0)
				count.mats[[n]][,basename(t)] <- (1e+6*countOverlaps(region.bed[[n]],bed,minoverlap=median(width(bed))-1))/(width(region.bed[[n]])*length(bed))
		}
	}

	z <- lapply(count.mats,function(x) apply(x,1,mean))
}

parse.map <- function(map)
{
	stru <- list()
	stru$region <- map$region
	stru$data <- vector("list",length(map$region))
	tmp <- map$data
	names(stru$region) <- names(stru$data) <- names(tmp) <- basename(stru$region)
	for (n in names(stru$region))
		stru$data[[n]] <- 	unlist(strsplit(tmp[[n]],split=","))
	stru$udata <- unique(unlist(stru$data))
	return(stru)
}

init.count.matrix <- function(x,r)
{
	# Prepare a list with names the regions and matrices of [nrow,ncol]=[length(region.bed[n]),length(map$data[[region]])]
	# They will be filled with the overlaps
	cmat <- vector("list",length(x))
	names(cmat) <- names(x)
	for (n in names(x))
	{
		cmat[[n]] <- matrix(0,length(r[[n]]),length(x[[n]]))
		colnames(cmat[[n]]) <- basename(x[[n]])
	}
	return(cmat)
}
