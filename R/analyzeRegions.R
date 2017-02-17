# Function to do some work with POLII... and not only... ;-)

# analyzeRegions is an R function that will accept as input one or multiple files that
# contain tag counts for specific genomic regions (e.g. genes or promoters) and perform
# certain operations on those regions (filtering, averaging, tag  count ratio calculation
# among different regions etc.) The user can have a detailed output with the original
# counts, the average counts and the counts ratios per given region. This script can be
# used as the next step of other scripts (e.g. getridofgalaxy.pl) to retrieve a list of
# genomic regions of interest. The script also collects certain statistics and the user
# can explore them inside R (e.g. indices of regions filtered etc.) so as to be able to
# overlap for example regions with a small number of counts and regions smaller than a
# predefined number of bps
#
# Usage: out <- analyzeRegions(files.in,
#                              file.out=NULL,
#                              ann.col=c(1:3),
#                              data.col=c(4,5),
#                              control.col=1,
#                              loc.col=c(2,3),
#                              avg.win=1000,
#                              filter.avg.count=5,
#                              filter.pct.count=1,
#                              filter.count=5,
#                              filter.std=1.5,
#                              filter.rat=0,
#                              calc.ratio=TRUE,
#                              rank.sum=FALSE,
#                              rank.sum.pctile=90,
#                              rank.sum.valid=1000,
#                              bin.rat.rank=NULL,
#                              norm.method="none",
#                              has.header=TRUE,
#                              plot.rvr=NULL,
#                              export.what="ratios",
#                              export.negative=FALSE)
#
# Input arguments:
#
# files.in         : One or more files that contain tag counts per region. The input files
#                    should have the same format (e.g. the 1st 3 columns should contain
#                    annotation elements, chromosome, start, end and the rest should
#                    contain tag counts per regions. The best way for the user is to gather
#                    different samples together in a file before submitting but the program
#                    the counts can also handle a set of files but with more careful use ;-)
#                    If not given a file chooser will open allowing the user to choose the
#                    files but only when having one file to submit, else they should be
#                    submitted manually.
# file.out         : A name for the output results file. It defaults to NULL and if not
#                    given it will be auto-generated.
# data.col         : The columns in the file that contain the data. For example if the
#                    first 3 columns contain annotation elements, the rest (c(4,...))
#                    should contain the data (tag counts).
# control.col      : The column(s) that correspond to the control counts and will be used
#                    for ratio calculation. If the user has one control sample and based
#                    on this ratios are calculated, then it should be a single number
#                    denoting the corresponding column. If the user has one control for
#                    each treated sample the controls should be a vector of equal size to
#                    the treated samples. If the user has multiple controls for multiple
#                    samples, still under development...
# loc.col          : The chromosome start and end columns. The file(s) submitted to this
#                    should containg OBLIGATORY this information.
# avg.win          : A window over which the tags should be averaged. For example if
#                    avg.win=500 then the average number of tags per 500bp windows will be
#                    calculated. This is useful for noise filtering processes. For example
#                    if a genomic is 15kb but has only 200 tags, while this might seem
#                    meaningful at the beginning, if one sees how many tags are per window,
#                    then realizes that this is probably artifacts. It defaults to 1000.
#                    It shouldn't be used for small regions. It is not used when set to 0.
# filter.avg.count : A threshold below which genomic regions with less counts than this
#                    number (per window size, given by the avg.win parameter) will be
#                    filtered out. For example, if filter.avg.count=5, then regions with
#                    less than 5counts/avg.win will be filtered out. Defaults to 5. Not
#                    used when set to 0. It should also be mentioned that in cases of
#                    multiple conditions, if one condition sattisfies the filter, then
#                    all of them will remain. To deal with this, the user should use the
#                    option at.least.in with the desired number of conditions that should
#                    have more than e.g. 5 counts
# at.least.in      : When having multiple conditions, this number defines in how many of
#                    them the number of counts should be larger than filter.avg.count.
#                    Defaults to 1.
# filter.pct.count : A threshold below which genomic regions with less counts than the
#                    percentage of the region length denoted by this number are filtered
#                    out. For example, if a region is 12kb in length, the 1% is 120 counts.
#                    If it contains less than 120 counts, it will be filtered out. Useful
#                    because it sets a threshold directly based on region length. Defaults
#                    to 1, not used when set to 0.
# filter.count     : A threshold of number of tags below which regions with less than this
#                    number are filtered out. Useful for small regions, dangerous for large
#                    as it might return noise. Defaults to 20. Not used when set to 0.
# filter.len       : A threshold below which smaller regions will be filtered. For example,
#                    the user might not want to include in final results regions less than
#                    1kb. Defaults to 1000. Not used when set to 0.
# filter.std       : The distance in standard deviations from the mean, for which the
#                    ratios will be filtered. For example, if filter.std=2, then regions
#                    with ratios (or counts if ratios not calculated) below +/- 2 stds
#                    from the mean will be filtered out. Defaults to 1.5. Not used when
#                    set to 0.
# filter.rat       : The fold change filter threshold. It defaults to 0. If the filter is
#                    not desired, it should take the value 0. The values are absolute and
#                    in log2 scale.
# calc.ratio       : Should log2 ratios be calculated? If not then the whole process will
#                    be based on the counts themselves. Defaults to TRUE
# rank.sum         : Calculate the RankSum metric to characterize genomic regions based on
#                    their content in tags. Useful when there is a control sample. RankSum
#                    introduced by Nagesh Rao and the formula is
#                    R = 1/(avg(#control tags per xbp)+(avg(#treated tags per xbp))
#                    It can be TRUE or FALSE (default).
# rank.sum.pctile  : The percentile of the RankSum distribution above which genomic regions
#                    are filtered. Defaults to 90 (if rank.sum=TRUE). Can also be NULL.
# rank.sum.val     : Validate the chosen threshold of RankSum by creating the bootstrap
#                    distribution of the chosen percentiles. This parameter is the number
#                    of bootstrap iterations that should be performed (default is 1000).
# has.header       : It should be set to TRUE if the file the user submits has a first line
#                    containing column names. Defaults to TRUE. These column names are also
#                    used in the output.
# bin.rat.rank     : Bin the regions according to their ratios to control samples and
#                    according to their ranksums. The bins aim to categorize the areas in
#                    3 categories for ratios and 3 categories for ranksums. The bins should
#                    be given in a list of 2 elements, the 1st being corresponding to ratios
#                    and the 2nd to ranksums. Each element of the list should be a vector
#                    of 2 elements corresponding to quantiles. For example, to bin ratios
#                    and ranksums according to their 1st, and 3rd quartiles, the list should
#                    be bin.rat.rank=list(c(25,75),c(25,75)). This will also create a 3x3
#                    grid in the ranksum plot. Defaults to NULL for no use. Note that
#                    categorization is performed ONLY in quality filtered regions, NOT all
#                    the original regions. If you wish to make a ranksum plot for ALL the
#                    regions, use bin.rat.rank=NULL.
# norm.method      : Normalize the log2 ratios in order to further standardize the variance
#                    among different conditions so that thresholds can be applied uniformly.
#                    It is recommended if the total vaiance of log2 ratio distributions vary
#                    much accross conditions. Available methods are "none" for no further
#                    normalization, "tukey" for 1-step Tukey's biweight normalization, "mad"
#                    for MAD centering normalization and "quantile" for Quantile normalization.
#                    It defaults to "none".
# plot.rvr         : Plot ranksum versus fold ratio (if calculated). Default is NULL (for
#                    not plotting. For plotting set this parameter to "x11" for displaying
#                    plots in R, "png" for creating a PNG image or "pdf" for a PDF file.
# export.what      : Choose what to export together with the annotation columns (filtered
#                    data of course) in the form of a vector (e.g. c("ratios","counts"))
#                    The following options are available:
#                    "ratios"                    : Export only ratios (if calculated, else counts)
#                    "counts"                    : Export counts (original ones)
#                    "avgcounts"                 : Export averaged counts (per avg.win)
#                    "ranksums"                  : Export ranksums
#                    "categories"                : Export strength categorized regions
#                    "all"                       : Export everything
#                    "nothing"                   : Exports nothing (returns only statistics in R)
#                    Defaults to export.what=c("ratios").
# export.negative  : Export also filtered regions, with the same settings as export.what.
#                    It defaults to FALSE.
#
# Output:
#
# Apart from the file exported, a list with the following slots is returned:
# Mean         : The mean ratios for all conditions
# StDev        : The standard deviations of the above
# Filter.Small : The filtered indices of small regions
# Filter.Avg   : The filtered indices of small averaged counts
# Filter.Pct   : The filtered indices of regions containing few counts compared to their
#                length
# Filter.Count : The filtered indices of regions containing fewer counts than a specified
#                threshold
# Filter.Std   : The filtered indices based on distance from the ratio means
# Filter.All   : All filtered indices
# Pass         : Remaining indices after all filters
#
# It should be noted, that in the case of multiple conditons (>2) if a region passes one
# filter (in counts mainly) for one condition, then it is maintained for all conditions.
# For example, if an average count filter is set to 5, regions with averaged counts<5 in
# all conditions will be filtered. If one passes, then it remains for all conditions.
#
# Usage examples:
#
# Analyze gene regions for transcription in the presence of POLII given a file with 4
# annotation columns, headers, 6 conditions from which the 1st is control:
# out <- analyzeRegions(ann.col=c(1:4), data.col=c(5:10), control.col=5, loc.col=c(2,3),
#                       avg.win=1000, filter.avg.count=5, filter.len=1000, filter.std=1.5,
#                       calc.ratio=TRUE, has.header=TRUE, export.what="all")
#
# Analyze small promoter regions for transcription in the presence of POLII given a file
# with 3 annotation columns, headers, 6 conditions from which the 1st is control:
# out <- analyzeRegions(ann.col=c(1:3), data.col=c(5:10), control.col=5, loc.col=c(2,3),
#                       avg.win=0, filter.avg.count=5, filter.len=0, filter.std=1.5,
#                       calc.ratio=TRUE, has.header=TRUE, export.what="all")
#
# Analyze gene regions for transcription in the presence of POLII given a file with 4
# annotation columns, no headers, 6 conditions from which the 1st 3 are controls and
# export only ratios:
# out <- analyzeRegions(ann.col=c(1:4), data.col=c(5:10), control.col=c(5:7),
#                       loc.col=c(2,3), avg.win=1000, filter.avg.count=5, filter.len=1000,
#                       filter.std=1.5, calc.ratio=TRUE, has.header=FALSE,
#                       export.what="ratios")
#
# Author: Panagiotis Moulos (pmoulos@eie.gr)
# Version: 1.21
# Creation date: 10-6-2008 (dd-mm-yyyy)
# Last update: 13-3-2009 (dd-mm-yyyy)
#
# TODO : Allow more categories in binning of ratios and ranksums
#

analyzeRegions <- function(files.in,file.out=NULL,ann.col=c(1:3),data.col=c(4,5),control.col=4,
                           loc.col=c(2,3),avg.win=1000,filter.avg.count=5,at.least.in=1,
                           filter.pct.count=1,filter.count=20,filter.len=1000,filter.std=1.5,
                           calc.ratio=TRUE,rank.sum=FALSE,rank.sum.pctile=90,rank.sum.val=1000,
                           filter.rat=0,bin.rat.rank=NULL,norm.method="none",has.header=TRUE,
                           plot.rvr=NULL,export.what="ratios",export.negative=FALSE)
{
    if (missing(files.in))
    {
        disp("Select a file from the lists directory to get the directory name\n")
        if (require(tcltk))
        {
            files.in <- tclvalue(tkgetOpenFile())
            if (!nchar(files.in))
            {
            	tkmessageBox(message="No file was selected!")
                return()
            }
        } else
            stop("You must provide a single or a list of files for the script to work!")
            
        # Fix problem with last file that becomes first
        if (length(files.in)>1)
            files.in <- files.in[c(2:length(files.in),1)]
    }
    
    # Get also destination directory
    v <- strsplit(files.in[1],"/")
    if (length(v[[1]])!=1)
    {
        des.dir <- v[[1]][1]
        for (i in 2:(length(v[[1]])-1))
            des.dir <- paste(des.dir,v[[1]][i],sep="/")
    } else
        des.dir=getwd();
    
    # Some further input checking...
    if (!is.numeric(ann.col))
        stop("The annotation columns vector must be numeric.")
    if (!is.numeric(loc.col))
        stop("The region start and end columns must be numeric.")
    if (!is.numeric(data.col) && !is.list(data.col))
        stop("The data columns must be a numeric vector or a list of numeric vectors.")
    if (!is.numeric(control.col) && !is.list(control.col))
        stop("The control columns must be a numeric vector or a list of numeric vectors.")
    if (!is.numeric(avg.win) || !is.numeric(filter.avg.count) || !is.numeric(filter.len) ||
        !is.numeric(filter.count) || !is.numeric(at.least.in) || !is.numeric(filter.std) ||
        !is.numeric(rank.sum.val) || !is.numeric(filter.rat))
         stop("avg.len, filter.avg.count, filter.std, rank.sum.val, filter.rat must be numeric.")
    if (!is.numeric(filter.pct.count) || filter.pct.count<0 || filter.pct.count>100)
        stop("filter.pct.count must be a numeric value between 0 and 100.")
    if (!is.list(bin.rat.rank) && !is.null(bin.rat.rank))
        stop("bin.rat.rank should be a list.")
    #if ((!is.numeric(rank.sum.pctile) && !is.null(rank.sum.pctile)) || rank.sum.pctile<0 ||
    #    rank.sum.pctile>100)
    #     stop("rank.sum.pctile must be a numeric value between 0 and 100.")
    if (length(norm.method)>1)
        stop("norm.method should be ONLY ONE of none, tukey, mad, quantile")
    if (!("none" %in% norm.method) && !("tukey" %in% norm.method) && !("mad" %in% norm.method)
        && !("quantile" %in% norm.method))
         stop("norm.method should be one of none, tukey, mad, quantile")
    
    # Check for the case of grouped data
    if (is.list(data.col) || is.list(control.col))
    {
        prop <- correctForGrouped(data.col,control.col,ann.col)
        control.cols <- prop[[1]]
        data.cols <- prop[[2]]
        control.col <- prop[[3]]
        data.col <- prop[[4]]
    } else
    {
        control.cols <- control.col
        data.cols <- data.col
    }

    disp("Reading input region files...")
    if (length(files.in)>1) # Multiple files, more work
    {
        # Determine annotation column from the first file
        first.file <- read.delim(files.in[1],header=has.header,quote="",check.names=FALSE)
        annot.cols <- first.file[,ann.col]
        all.num.data <- as.numeric(first.file[,data.cols[1]])
        for (i in 2:length(files.in))
        {
            temp.table <- read.delim(files.in[i],header=has.header,quote="",check.names=FALSE)
            all.num.data <- cbind(all.num.data,as.numeric(temp.table[,data.cols[i]]))
        }
        all.num.data <- as.matrix(all.num.data)
        control.cols <- control.cols-ncol(annot.cols)
    } else # One file with counts... better
    {
        temp.table <- read.delim(files.in,header=has.header,quote="",check.names=FALSE)
        annot.cols <- temp.table[,ann.col]
        all.num.data <- as.matrix(temp.table[,data.cols])
        control.cols <- control.cols-ncol(annot.cols)
    }
    data.names <- colnames(all.num.data)
    
    if (calc.ratio)
    {
        disp("Calculating ratios...")
        # Transform 0 -> 1 for log2
        all.num.data.corr <- all.num.data
        all.num.data.corr[which(all.num.data.corr==0)] <- 1
        if (length(control.cols)==1) # The case of only one control, fine
        {
            ctrl.num <- as.matrix(all.num.data.corr[,control.cols])
            rest.num <- as.matrix(all.num.data.corr[,-control.cols])
            for (i in 1:ncol(rest.num))
                rest.num[,i] <- log2(rest.num[,i]/ctrl.num)
        } else # The case of each data point and its correponding control
        if (ncol(all.num.data[,control.cols])==ncol(all.num.data[,-control.cols]))
        {
            ctrl.num <- as.matrix(all.num.data.corr[,control.cols])
            rest.num <- as.matrix(all.num.data.corr[,-control.cols])
            for (i in 1:ncol(rest.num))
                rest.num[,i] <- log2(rest.num[,i]/ctrl.num[,i])
        }
        #else if (is.list(control.col)) # Hardest part... in development...
        #{
        #    ctrl.num <- list()
        #    rest.num <- list()
        #    for (i in 1:length(control.col))
        #    {
        #        ctrl.num[[i]] <- all.num.data[,control.col[[i]]]
        #        rest.num[[i]] <- all.num.data[,-control.col[[i]]]
        #        for (j in 1:dim(rest.num[[i]])[2])
        #            rest.num[[i]][,j] <- log2(rest.num[[i]][,j]/ctrl.num[[i]])
        #    }
        #}
        else
            stop("Improper input configuration")
    } else rest.num <- all.num.data
    
    # Start quality filtering process
    disp("Filtering regions...")
    all.ind <- 1:nrow(all.num.data)
    starts <- as.numeric(annot.cols[,loc.col[1]])
    ends <- as.numeric(annot.cols[,loc.col[2]])
    diffs <- ends-starts
    
    # Less than filter.len regions
    #ifelse(filter.len>0,filt.less <- which(diffs<filter.len),filt.less <- NULL)
    if (filter.len>0)
        filt.less <- which(diffs<filter.len) else filt.less <- NULL

    # If average counts per window
    if (avg.win>0)
    {
        avg.factor <- diffs/avg.win
        avg.counts <- matrix(NA,nrow=nrow(all.num.data),ncol=ncol(all.num.data))
        for (i in 1:ncol(all.num.data))
            avg.counts[,i]<-all.num.data[,i]/avg.factor
    } else avg.counts <- all.num.data

    # Less than filter.avg.count counts per region
    if (filter.avg.count>0)
    {
        filt.avg.temp <- apply(avg.counts,1,comp.ext,filter.avg.count,at.least.in)
        filt.avg.count <- which(filt.avg.temp!=0)
    } else filt.avg.count <- NULL
    
    # Less than filter.pct.count counts per region
    if (filter.pct.count>0)
    {
        dyn.lens <- (filter.pct.count/100)*diffs
        filt.pct.temp <- numeric(length(diffs))
        for (i in 1:length(diffs))
            if(all(all.num.data[i,]<dyn.lens[i]))
                filt.pct.temp[i] <- 1
        filt.pct <- which(filt.avg.temp!=0)
    } else filt.pct <- NULL
    
    # Less than filter.count counts per region
    if (filter.count>0)
    {
        filt.count.temp <- apply(all.num.data,1,comp.ext,filter.count,at.least.in)
        filt.count <- which(filt.count.temp!=0)
    } else filt.count <- NULL
    
    # Rank sum filtering
    if (rank.sum)
    {
        if (avg.win==0)
        {
            warning("Averaging window not defined. Ignoring rank sum filtering...",
                    immediate.=TRUE)
            rank.filt <- NULL
            ranks <- all.num.data
        } else
        {
            if (length(control.cols)==1) # The case of only one control, fine
            {
                ctrl.ac <- as.matrix(avg.counts[,control.cols])
                rest.ac <- as.matrix(avg.counts[,-control.cols])
            } else # The case of each data point and its correponding control
            if (ncol(all.num.data[,control.cols])==ncol(all.num.data[,-control.cols]))
            {
                ctrl.ac <- as.matrix(avg.counts[,control.cols])
                rest.ac <- as.matrix(avg.counts[,-control.cols])
            }
            ranks <- matrix(NA,nrow=nrow(rest.ac),ncol=ncol(rest.ac))
            rank.pctile <- new.rank.pctile <- rep(0,ncol(rest.ac))
            rank.filt.list <- list()
            for (i in 1:ncol(rest.ac))
            {
                ranks[,i] <- 1/(ctrl.ac+rest.ac[,i])
                # Find requested percentile of ranks
                if (!is.null(rank.sum.pctile))
                    rank.pctile[i] <- quantile(ranks[,i],rank.sum.pctile/100,type=5)
                # Correct it through bootstraping
                if (rank.sum.val!=0)
                {
                    disp("Bootstrapping RankSum scores...")
                    new.rank.pctile[i] <- correctBoot(ranks[,i],rank.pctile,rank.sum.pctile,
                                                      rank.sum.val)
                } else
                    new.rank.pctile[i] <- rank.pctile[i]
                rank.filt.list[[i]] <- which(ranks[,i]>new.rank.pctile[i])
            }
            if (!is.null(rank.sum.pctile))
                rank.filt <- unique(unlist(rank.filt.list)) else rank.filt <- NULL
                
        }
    } else
    {
        rank.filt <- NULL
        ranks <- all.num.data
    }

    # Unify quality filters...
    all.filt <- union(filt.less,union(filt.avg.count,union(filt.pct,union(filt.count,rank.filt))))
    
    # Start filtering on standard deviations...
    if (filter.std>0)
    {
      means <- apply(rest.num,2,mean)
      stds <- apply(rest.num,2,sd)
      logical.rest <- matrix(1,nrow(rest.num),ncol(rest.num))
      inds.per.col <- list()
      for (i in 1:ncol(rest.num))
      {
          inds.per.col[[i]] <- which(rest.num[,i]<means[i]-filter.std*stds[i] |
                                     rest.num[,i]>means[i]+filter.std*stds[i])
          logical.rest[inds.per.col[[i]],i] <- 0
      }
      # Find the union of the above
      thresh.temp <- apply(logical.rest,1,greater.than,0)
      thresh.filt <- which(thresh.temp!=0)
    } else 
    {
      thresh.filt <- NULL
      means <- NULL
      stds <- NULL
    }
    
    # Filter on ratios...
    if (filter.rat>0)
    {
        logical.rest <- matrix(1,nrow(rest.num),ncol(rest.num))
        inds.per.col <- list()
        for (i in 1:ncol(rest.num))
        {
            inds.per.col[[i]] <- which(abs(rest.num[,i])>filter.rat)
            logical.rest[inds.per.col[[i]],i] <- 0
        }
        # Find the union of the above
        thresh.temp <- apply(logical.rest,1,greater.than,0)
        thresh.rat <- which(thresh.temp!=0)
    } else thresh.rat <- NULL
    
    # Final filters (quality and thresholds)
    final.filter.inds <- union(all.filt,union(thresh.filt,thresh.rat))
    # Good regions...
    ifelse(is.null(final.filter.inds),good.ind <- all.ind,
           good.ind <- all.ind[-final.filter.inds])
    if (length(good.ind)==0)
           good.ind <- all.ind
           
    # Categorize final distributions according to ratios and ranksums
    if (!is.null(bin.rat.rank) && calc.ratio && rank.sum)
    {
        good.rats <- as.matrix(rest.num[good.ind,])
        good.ranks <- as.matrix(ranks[good.ind,])
        categs.rat <- categs.rank <- categs.rat.rank <- matrix(NA,nrow(good.rats),ncol(good.rats))
        rat.pct <- c(0,bin.rat.rank[[1]]/100,1)
        rank.pct <- c(0,bin.rat.rank[[2]]/100,1)
        for (i in 1:ncol(good.rats))
        {
            rat.qnt <- quantile(good.rats[,i],rat.pct)
            rank.qnt <- quantile(good.ranks[,i],rank.pct)
            for (j in 2:length(rat.pct))
                categs.rat[which(good.rats[,i]>=rat.qnt[j-1] & good.rats[,i]<=rat.qnt[j]),i] <- j-1
            for (j in 2:length(rank.pct))
                categs.rank[which(good.ranks[,i]>=rank.qnt[j-1] & good.ranks[,i]<=rank.qnt[j]),i] <- length(rank.pct)-j+1
            categs.rat.rank[,i] <- t(paste(categs.rat[,i],categs.rank[,i],sep=""))
        }
    } else categs.rat.rank <- matrix(NA,length(good.ind),dim(all.num.data)[2])
    
    # Normalize data if desired
    if (norm.method!="none")
    {
        if (norm.method=="tukey")
        {
            if (!require(affy))
                stop("You need to install Bioconductor! See http://www.bioconductor.org") else
            rest.num <- apply(rest.num,2,tukey.biweight)
        }
        if (norm.method=="quantile")
        {
            if (!require(preprocessCore))
                stop("You need to install Bioconductor! See http://www.bioconductor.org") else
            rest.num <- normalize.quantiles(rest.num)
        }
        if (norm.method=="mad")
        {
            data.median <- apply(rest.num,2,median)
            data.mad <- apply(rest.num,2,mad)
            for (i in 1:ncol(rest.num))
                rest.num[,i] <- (rest.num[,i]-data.median[i])/data.mad[i]
        }
    }

    # Export data
    if (!any(export.what=="nothing"))
    {
        disp("Exporting analyzed and filtered regions...")
        if (is.null(file.out))
        {
            suffix <- gsub('[ |:]','-',date())
            file.out <- paste(des.dir,"/","sign_regions_",suffix,".txt",sep="")
        } else
            f.file.out <- paste(des.dir,file.out,sep="/")
        exportData(file.out,data.names,annot.cols,control.cols,good.ind,rest.num,
                   all.num.data,avg.counts,ranks,categs.rat.rank,rank.sum,bin.rat.rank,
                   export.what)
                   
        # Export negative regions
        if (export.negative)
        {
            disp("Exporting negative regions...")
            bad.ind <- all.ind[final.filter.inds]
            if (is.null(file.out))
                file.neg <- paste(des.dir,"/","negative_regions_",suffix,".txt",sep="")
            else
            {
                f<-unlist(strsplit(file.out,".",fixed=TRUE))
                l.f<-length(f)
                n.n<-character(0)
                if (l.f==1)
                    n.n<-f else
                {
                    t.n<-f[1:(l.f-1)]
                    n.n<-t.n[1]
                    if (l.f>2)
                        for (j in 2:(l.f-1))
                            n.n<-paste(n.n,t.n[j],sep=".")
                }
                file.neg <- paste(des.dir,"/",n.n,"_negative.txt",sep="")
            }
            exportData(file.neg,data.names,annot.cols,control.cols,bad.ind,rest.num,
                       all.num.data,avg.counts,ranks,categs.rat.rank,rank.sum,
                       bin.rat.rank,export.what)
        }
    }
    
    # Plot if required
    if (!is.null(plot.rvr))
    {
        if (calc.ratio && rank.sum)
        {
            if (!is.null(bin.rat.rank))
            {
                plotRankRatMultiFixed(as.matrix(ranks[good.ind,]),as.matrix(rest.num[good.ind,]),
                                      categs.rat.rank,rank.qnt,rat.qnt,
                                      main.title=data.names[-control.cols],out=plot.rvr)
            } else
                plotRankRat(ranks,rest.num,ran.thr=new.rank.pctile,rat.thr=filter.rat,
                            main.title=data.names[-control.cols],out=plot.rvr)

        }
    }


    # Return some statistics
    out <- list(means,stds,filt.less,filt.count,thresh.filt,final.filter.inds,good.ind)
    names(out) <- c("Mean","StDev","Filter.Small","Filter.Avg","Filter.Std",
                    "Filter.All","Pass")
    disp("Finished\n")
    return(out)
}


correctBoot <- function(data.vec,old.val,pct=90,n=1000)
{
    if (missing(old.val))
        old.val<-quantile(data.vec,pct/100,type=5)
    if (!require(boot))
        install.packages("boot",repos="http://cran-mirror.cs.uu.nl")
        
    myq <- function(x,ind,...) { quantile(x[ind],...) }
    boot.obj <- boot(data.vec,myq,n,stype="i",probs=pct/100,type=5)
    q.obj <- imp.quantile(boot.obj,alpha=pct)
    return(q.obj$rat)
}


# Helper exporting function
exportData <- function(file.out,data.names,annot.cols,control.cols,good.ind,
                       rest.num,all.num.data,avg.counts,ranks,categs.rat.rank,
                       rank.sum,bin.rat.rank,export.what="ratios")
{
    # Prepare to export
    ann.export <- annot.cols[good.ind,]
    if ("ratios" %in% export.what)
    {
        ratios <- as.matrix(rest.num[good.ind,])
        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")

    } else ratios <- NULL
    if ("counts" %in% export.what)
    {
        counts <- as.matrix(all.num.data[good.ind,])
        colnames(counts) <- data.names
    } else counts <- NULL
    if ("avgcounts" %in% export.what)
    {
        avgcounts <- as.matrix(avg.counts[good.ind,])
        colnames(avgcounts) <- paste(data.names,"_AVG_Counts",sep="")
    } else avgcounts <- NULL
    if ("ranksums" %in% export.what)
    {
        ranksums <- as.matrix(ranks[good.ind,])
        if (rank.sum)
            colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum",sep="")
        else
            colnames(ranksums) <- paste(data.names,"_RankSum",sep="")
    } else ranksums <- NULL
    if ("categories" %in% export.what)
    {
        categories <- as.matrix(categs.rat.rank)
        if (!is.null(bin.rat.rank))
            colnames(categories) <- paste(data.names[-control.cols],"_Flag",sep="")
        else
            colnames(categories) <- paste(data.names,"_Flag",sep="")
    } else categories <- NULL
    if ("all" %in% export.what)
    {
        ratios <- as.matrix(rest.num[good.ind,])
        counts <- as.matrix(all.num.data[good.ind,])
        avgcounts <- as.matrix(avg.counts[good.ind,])
        ranksums <- as.matrix(ranks[good.ind,])
        categories <- as.matrix(categs.rat.rank)
        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
        colnames(counts) <- data.names
        colnames(avgcounts) <- paste(data.names,"_AVG_Counts",sep="")
        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum",sep="")
        colnames(categories) <- paste(data.names[-control.cols],"_Flag",sep="")
    }
    # Export
    data.export <- as.data.frame(cbind(ratios,counts,avgcounts,ranksums,categories))
    all.export <- cbind(ann.export,data.export)
    write.table(all.export,file=file.out,quote=FALSE,sep="\t",row.names=FALSE)
}


# Plot RankSum against ratio
plotRankRat <- function(ranks,ratios,ran.thr=NULL,rat.thr=NULL,
                        main.title=rep("Data",ncol(ranks)),out="x11")
{
    if (out=="pdf")
    {
        pdf(file=paste(main.title[1],".pdf",sep=""),width=7,height=6)
        sex <- 0.4
        lsd <- 1
    } else
    {
        sex <- 0.7
        lsd <- 2
    }
        
    for (i in 1:ncol(ranks))
    {
        if (out=="png")
            png(filename=paste(main.title[i],".png",sep=""),width=640,height=480)
        if (out=="x11") x11()
        
        plot(ranks[,i],ratios[,i],pch=20,col="blue",main=paste("RankSum vs Ratio plot for",main.title),
             xlab="RankSum",ylab="log2(Ratio)",cex=sex)
        if (!is.null(ran.thr) && !all(ran.thr==0))
            abline(v=ran.thr[i],lty="dashed",col="green",lwd=lsd)
        if (!is.null(rat.thr) && rat.thr!=0)
        {
            abline(h=rat.thr,lty="dashed",col="red",lwd=lsd)
            abline(h=-rat.thr,lty="dashed",col="red",lwd=lsd)
        }
        
        if (out=="png") dev.off()
    }
    
    if (out=="pdf") dev.off()
}


# Plot RankSum against ratio
plotRankRatMultiFixed <- function(ranks,ratios,flags,rank.qnt,rat.qnt,
                                  main.title=rep("Data",ncol(ranks)),out="x11")
{
    if (out=="pdf")
    {
        pdf(file=paste(main.title[1],".pdf",sep=""),width=7,height=6)
        sex <- 0.5
        sex.axis <- 1.8
        sex.main <- 2.5
        sex.lab <- 1.8
        sex.leg <- 1.7
        lsd <- 1.7
    } else
    {
        sex <- 0.7
        sex.axis <- 1.3
        sex.main <- 1.5
        sex.lab <- 1.3
        sex.leg <- 1.3
        lsd <- 2
    }
    
    # We plot row-wise
    fixed.colors <- c("red","yellow","black","green","grey","magenta","blue","orange","brown")
    ab.rank <- rank.qnt[-c(1,length(rank.qnt))]
    ab.rat <- rat.qnt[-c(1,length(rat.qnt))]

    for (i in 1:ncol(ranks))
    {
        if (out=="png")
            png(filename=paste(main.title[i],".png",sep=""),width=640,height=480)
        if (out=="x11") x11()

        un.flags <- sort(unique(flags[,i]),decreasing=TRUE)
        par(cex=sex,cex.axis=sex.axis,cex.main=sex.main,cex.lab=sex.lab,font.lab=2)
        plot.new()
        plot.window(c(range(ranks[,i])[1],range(ranks[,i])[2]+0.025),range(ratios[,i]))
        axis(1,at=pretty(ranks[,i]))
        axis(2,at=pretty(ratios[,i]))
        title(main=paste("RankSum vs Ratio plot for",main.title),xlab="RankSum",ylab="log2(Ratio)")

        for (j in 1:length(un.flags))
            points(ranks[which(flags[,i]==un.flags[j]),i],ratios[which(flags[,i]==un.flags[j]),i],
                   pch=20,col=fixed.colors[j])

        for (j in 1:length(ab.rank))
            abline(v=ab.rank[j],lty="dashed",col="green",lwd=lsd)
        for (j in 1:length(ab.rat))
            #abline(h=ab.rat[j],lty="dashed",col="red",lwd=lsd)
            lines(seq(range(ranks[,i])[1],range(ranks[,i])[2],length.out=100),
                  rep(ab.rat[j],100),lty="dashed",col="red",lwd=lsd)
        
        # Legend stuff, put outside plotting region on the right
        old.par<-par(xpd=TRUE)
        y.pos<-mean(par("usr")[3:4])
        in2x.pos<-(par("usr")[2]-par("usr")[1])/par("pin")[1]
        x.pos<-par("usr")[2]-par("mai")[4]*in2x.pos
        legend(x.pos,y.pos,legend=un.flags,col=fixed.colors,text.col=fixed.colors,
               cex=sex.leg,pch=20,box.lwd=0.5,xjust=0.4,yjust=0.5)
        par(old.par)

        if (out=="png") dev.off()
    }

    if (out=="pdf") dev.off()
}


# Compare per row funcs... one per case instead of all-in-one for speed
less.than <- function(x,val) { ifelse(all(x<val),return(1),return(0)) }
greater.than <- function(x,val) { ifelse(all(x>val),return(1),return(0)) }
less.equal <- function(x,val) { ifelse(all(x<=val),return(1),return(0)) }
greater.equal <- function(x,val) { ifelse(all(x>=val),return(1),return(0)) }

comp.ext <- function(x,val,num)
{
    logic <- x>val
    ans <- length(which(logic==TRUE))
    ifelse(ans>=num,return(0),return(1))
}
    

# Case where we have e.g. 2 controls and 4 treated so we have to group controls and
# treated samples somehow, original structure remains
# Follows a function to perform a total correction in indices... since this feature is
# still in development (and probably will remain...) correct for these things in a wrapper
# function so as to avoid garbage code
correctForGrouped <- function(d.c,c.c,a.c)
{
    if ((is.list(d.c) && !is.list(c.c)) || (!is.list(d.c) && is.list(c.c)))
        stop(paste("If your data are grouped with single controls corresponding to multiple treatments",
                   "you should provide the columns grouped in lists in both cases"))
    if (is.list(c.c))
    {
        c.cs <- unlist(c.c)
        d.cs <- unlist(d.c)
        if (!is.numeric(c.cs) || !is.numeric(d.cs))
            stop("The vectors of columns must be numeric!")
        for (i in 1:length(c.c))
        {
            c.c[[i]] <- c.c[[i]]-length(a.c)
            d.c[[i]] <- d.c[[i]]-length(a.c)
        }
    }
    else
    {
        c.cs <- c.c
        d.cs <- d.c
    }
    return(list(c.cs,d.cs,c.c,d.c))
}


disp <- function(msg)
{
    cat(paste(msg,"\n",sep=""))
    flush.console()
}


## Old helper exporting function
#exportData <- function(file.out,data.names,annot.cols,control.cols,good.ind,
#                       rest.num,all.num.data,avg.counts,ranks,export.what="ratios")
#{
#    # Prepare to export
#    ann.export <- annot.cols[good.ind,]
#    if (export.what=="ratios")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        data.export <- as.data.frame(ratios)
#    } else
#    if (export.what=="counts")
#    {
#        counts <- as.matrix(all.num.data[good.ind,])
#        colnames(counts) <- data.names
#        data.export <- as.data.frame(counts)
#    } else
#    if (export.what=="avgcounts")
#    {
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        data.export <- as.data.frame(avgcounts)
#    } else
#    if (export.what=="ranksums")
#    {
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(ranksums)
#    } else
#    if (export.what=="ratios-counts")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        counts <- as.matrix(all.num.data[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(counts) <- data.names
#        data.export <- as.data.frame(cbind(ratios,counts))
#    } else
#    if (export.what=="ratios-avgcounts")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        as.data.frame(cbind(ratios,avgcounts))
#    } else
#    if (export.what=="counts-avgcounts")
#    {
#        counts <- as.matrix(all.num.data[good.ind,])
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        colnames(counts) <- data.names
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        as.data.frame(cbind(counts,avgcounts))
#    } else
#    if (export.what=="ratios-ranksums")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(cbind(ratios,ranksums))
#    } else
#    if (export.what=="counts-ranksums")
#    {
#        counts <- as.matrix(all.num.data[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(counts) <- data.names
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(cbind(counts,ranksums))
#    } else
#    if (export.what=="avgcounts-ranksums")
#    {
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(cbind(avgcounts,ranksums))
#    } else
#    if (export.what=="ratios-counts-ranksums")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        counts <- as.matrix(all.num.data[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(counts) <- data.names
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(cbind(ratios,counts,ranksums))
#    } else
#    if (export.what=="ratios-avgcounts-ranksums")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        as.data.frame(cbind(ratios,avgcounts,ranksums))
#    } else
#    if (export.what=="counts-avgcounts-ranksums")
#    {
#        counts <- as.matrix(all.num.data[good.ind,])
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(counts) <- data.names
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        as.data.frame(cbind(counts,avgcounts,ranksums))
#    } else
#    if (export.what=="all")
#    {
#        ratios <- as.matrix(rest.num[good.ind,])
#        counts <- as.matrix(all.num.data[good.ind,])
#        avgcounts <- as.matrix(avg.counts[good.ind,])
#        ranksums <- as.matrix(ranks[good.ind,])
#        colnames(ratios) <- paste(data.names[-control.cols],"_log2_Ratio",sep="")
#        colnames(counts) <- data.names
#        colnames(avgcounts) <- paste(data.names,"_AVG_Counts")
#        colnames(ranksums) <- paste(data.names[-control.cols],"_RankSum")
#        data.export <- as.data.frame(cbind(ratios,counts,avgcounts,ranksums))
#    }
#    # Export
#    all.export <- cbind(ann.export,data.export)
#    write.table(all.export,file=file.out,quote=FALSE,sep="\t",row.names=FALSE)
#}

# In the case of ever finding a way to read properly multiple files with tk...
#files.in <- strsplit(files.in,"} ")
#files.in <- files.in[[1]]
#files.in <- gsub('[{}]','',files.in)
