# Parse MACS output

#avg: average per xbps (100)
#bin.rr     : Bin the regions according to their ratios to control samples and
#             according to their ranksums. The bins aim to categorize the areas in
#             3 categories for ratios and 3 categories for ranksums. The bins should
#             be given in a list of 2 elements, named "ratio" and "rs" respectively.
#             Each element of the list should be a vector of 2 elements corresponding
#             to quantiles. For example, to bin ratios and ranksums according to their
#             1st and 3rd quartiles, the list should be bin.rr=list(ratio=c(25,75),rs=c(25,75)).
#             This will also create a 3x3 grid in the ranksum plot. Defaults to NULL for 
#             no use. Note that categorization is performed ONLY in quality filtered 
#             regions, NOT all the original regions.

# Rank sum becomes obsolette, we now work with average counts per avg.win
# bin.rr=list(ratio=c(25,75),rpw=c(25,75))
# bin.rr=list(ratio=c(5,50),rpw=c(5,25))

# TODO: When normalization != none but has not control, the control norm reads still try to export

processMACSOutput <- function(input,output=NA,mapfile=NA,fdr.cut=ifelse(ver==14,1,0.05),fc.cut=1,org="hg18",avg.win=100,win.size=10000,
                               normalize=c("none","linear","balance","rpkm","peakseq.original","peakseq.rlm"),chrom.info.file=NULL,ver=c(14,2),
                               eff.size=if (is.element(org,c("hg18","hg19","mm8","mm9","mm10","dm2","dm3","ce"))) org else NULL,sat=c(0.05,0.2),
                               pileup=5,write.output=TRUE,bin.rr=NULL,plot.rvr=FALSE,plot.svr=FALSE,image.format="x11",export.negative=FALSE)
{
    # If user does not provide files
    if (missing(input) && is.na(mapfile))
    {
        cat("Select your MACS peak file(s). Make sure they are all of the same format.\n\n")
        flush.console()
        if (require(tcltk))
        {
            input <- tclvalue(tkgetOpenFile())
            if (!nchar(input))
            {
                tkmessageBox(message="No file was selected!")
                return()
            }
        }
        else
            stop("You must provide MACS output file(s) for the script to work!")

        # Fix problem with last file that becomes first
        if (length(input)>1)
        input <- input[c(2:length(input),1)]
    }
    else if (missing(input) && !is.na(mapfile))
    {
        map <- read.delim(mapfile,stringsAsFactors=FALSE)
        input <- map$filename
    }
  
    # Get also destination directories
    path <- dirname(input)
    names(path) <- basename(input)
    
    # Some more input checking
    normalize <- tolower(normalize[1])
    if (!is.element(normalize,c("none","linear","balance","rpkm","peakseq.original","peakseq.rlm")))
        stop("normalize must be one of \"none\", \"linear\", \"balance\", \"rpkm\", \"peakseq.original\" or \"peakseq.rlm\"")
    ver <- ver[1]
    if (!is.element(ver,c(14,2)))
        stop("MACS version must be 14 or 2!")
    if (!is.numeric(fdr.cut))
        stop("The FDR cutoff must be numeric.")
    if (!is.numeric(fc.cut))
        stop("The fold change cutoff must be numeric.")
    if (!is.numeric(pileup))
        stop("The pileup cutoff must be numeric.")
    if (is.na(output))
        output="auto"
    if (!is.na(output) && length(input)!=length(output))
    {
        warning("Output file names do not have the same length as input or not given! They will be auto-generated...",
                call.=FALSE,immediate.=TRUE)
        output <- "auto"
    }
    if (!is.na(mapfile))
    {
        if (!require(IRanges))
            stop("Bioconductor package IRanges is required!")
        if (!require(GenomicRanges))
            stop("Bioconductor package GenomicRanges is required!")
        if (!require(rtracklayer))
            stop("Bioconductor package rtracklayer is required!")
        map <- read.delim(mapfile,stringsAsFactors=FALSE)
        if (is.null(map$control))
            hasControl <- FALSE
        else
            hasControl <- TRUE
    }
    if (normalize=="balance")
    {
        if (is.null(map$rawt))
            stop("You must provide also raw input bed files before external tag balance normalization!")
        if (hasControl && is.null(map$rawc))
            stop("You must provide also raw input bed files before external tag balance normalization for control samples!")
    }
    if (normalize=="peakseq")
    {
        if (!require(GenomicRanges))
            stop("Bioconductor package IRanges is required!")
        if (!is.null(chrom.info.file))
        {
            chrom.info <- read.delim(chr.info.file,header=FALSE)
            genome.size <- sum(as.numeric(chrom.info[,2]))
            chrom.info <- chrom.info[-grep("rand|chrM|hap",chrom.info[,1]),]
            attr(chrom.info,"genome.size") <- genome.size
        }
        else
            chrom.info <- getChromInfo(org)
        gen.win <- make.windows(chrlen=as.numeric(chrom.info[,2]),winsize=win.size,chrname=as.character(chrom.info[,1]))
        gen.win.names <- paste(seqnames(gen.win),paste(start(gen.win),end(gen.win),sep="-"),sep=":")
    }

    if (eff.size == org)
    {
        eff.gen.size <- switch(eff.size,
            hg18 = 2700000000,
            hg19 = 2700000000,
            mm8 = 1865500000,
            mm9 = 1865500000,
            mm10 = 1865500000,
            dm2 = 120000000,
            ce = 90000000)
    }
    else if (is.null(eff.size))
        stop("Please provide an effective genome size for your input organism!")
    else
        eff.gen.size <- eff.size

    n <- length(input)
    all.peaks <- fdr.fail <- fc.fail <- sat.fail <- filter.fail <- list()
    if (ver == 14) # MACS version 1.4
    {
        for (i in 1:n)
        {
            cat("Reading MACS 1.4 output ",basename(input[i])," ...\n",sep="")
            all.peaks[[basename(input[i])]] <- read.delim(input[i],skip=23,row.names=NULL)
            names(all.peaks[[basename(input[i])]]) <- 
                c("chromosome","start","end","length","summit","tags","significance","MACS_fe","fdr")
            # Correct the summit
            all.peaks[[basename(input[i])]]$summit <- 
                all.peaks[[basename(input[i])]]$summit+all.peaks[[basename(input[i])]]$start
            # FDR filtering
            fdr.fail[[basename(input[i])]] <- which(all.peaks[[basename(input[i])]]$fdr>=fdr.cut)
        }
    }
    else if (ver==2) # MACS version 2
    {
        for (i in 1:n)
        {
            cat("Reading MACS 2 output ",basename(input[i])," ...\n",sep="")
            all.peaks[[basename(input[i])]] <- tryCatch(
                read.delim(input[i],skip=26,row.names=NULL),
            error = function(e) {
                return(read.delim(input[i],skip=27,row.names=NULL))
            },finally=NULL)
            if (ncol(all.peaks[[basename(input[i])]])==10)
                names(all.peaks[[basename(input[i])]]) <- 
                    c("chromosome","start","end","length","summit","pileup",
                        "significance","MACS_fe","fdr","name")
            else
                names(all.peaks[[basename(input[i])]]) <- 
                    c("chromosome","start","end","length","pileup",
                        "significance","MACS_fe","fdr","name")
            # FDR filtering
            if (pileup != 0 && !is.na(pileup) && !is.null(pileup))
                fdr.fail[[basename(input[i])]] <- which(10^-all.peaks[[basename(input[i])]]$fdr>=fdr.cut | all.peaks[[basename(input[i])]]$pileup<pileup)
            else
                fdr.fail[[basename(input[i])]] <- which(10^-all.peaks[[basename(input[i])]]$fdr>=fdr.cut)
        }
    }

    # Now we must apply counting... IRanges and GRanges
    if (!is.na(mapfile))
    {
        # First, import the peak files and store them, as they are small
        input.bed <- ovObj.treat <- ovObj.control <- ovObj.rawt <- ovObj.rawc <-
        ovObj.nreads <- ovObj.peakseq <- ovObj.taglen <- list()
        for (i in 1:n)
        {
            cat("Importing peak file ",basename(input[i])," as bed... Please wait...\n",sep="")
            flush.console()
            tmp <- tempfile()
            tmp.frame <- cbind(all.peaks[[basename(input[[i]])]][,1:5],as.matrix(rep("*",nrow(all.peaks[[basename(input[[i]])]]))))
            write.table(tmp.frame,tmp,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
            input.bed[[basename(input[[i]])]] <- import.bed(tmp,trackLine=FALSE)
            unlink(tmp)
        }
        
        # Then, import and process the big mother bed files...
        for (i in 1:length(map$treatment))
        {
            treat <- map$treatment[i]
            cat("Importing track ",basename(treat),"... Please wait...\n",sep="")
            flush.console()
            treat.bed <- import.bed(treat,trackLine=FALSE)

            if (normalize=="balance")
            {
                rawt <- map$rawt[i]
                cat("Importing track ",basename(rawt),"... Please wait...\n",sep="")
                rawt.bed <- import.bed(rawt,trackLine=FALSE)
                flush.console()
            }
            
            if (hasControl)
            {
                control <- map$control[i]
                cat("Importing track ",basename(control),"... Please wait...\n",sep="")
                flush.console()
                control.bed <- import.bed(control,trackLine=FALSE)

                if (normalize=="balance")
                {
                    rawc <- map$rawc[i]
                    cat("Importing track ",basename(rawc),"... Please wait...\n",sep="")
                    rawc.bed <- import.bed(rawc,trackLine=FALSE)
                    flush.console()
                }
            }
        
            nams <- map$filename[which(map$treatment==treat)]
            for (nam in nams)
            {
                nam <- basename(nam)
                cat("  Counting tags in the treatment track ",basename(treat)," over the corresponding peak file ",nam,"... Please wait...\n",sep="")
                flush.console()
                ovObj.nreads[[nam]]$treatment = length(treat.bed)
                ovObj.taglen[[nam]]$treatment = median(width(treat.bed))
                ovObj.treat[[nam]] <- countOverlaps(input.bed[[nam]],treat.bed,minoverlap=ovObj.taglen[[nam]]$treatment-1)
                if (normalize=="balance")
                {
                    cat("  Counting tags in the raw treatment track ",basename(rawt)," over the corresponding peak file ",nam,"... Please wait...\n",sep="")
                    ovObj.rawt[[nam]] <- countOverlaps(input.bed[[nam]],rawt.bed,minoverlap=ovObj.taglen[[nam]]$treatment-1)
                }
                if (normalize=="peakseq.original" || normalize=="peakseq.rlm")
                {
                    cat("  Counting tags in the treatment track ",basename(treat)," over genomic windows for peakseq... Please wait...\n",sep="")
                    ovObj.peakseq[[nam]]$treatment <- countOverlaps(gen.win,treat.bed)
                    names(ovObj.peakseq[[nam]]$treatment) <- gen.win.names
                }
                if (hasControl)
                {
                    cat("  Counting tags in the control track ",basename(control)," over the corresponding peak file ",nam,"... Please wait...\n",sep="")
                    flush.console()
                    ovObj.nreads[[nam]]$control = length(control.bed)
                    ovObj.taglen[[nam]]$control = median(width(control.bed))
                    ovObj.control[[nam]] <- countOverlaps(input.bed[[nam]],control.bed,minoverlap=ovObj.taglen[[nam]]$control-1)
                    if (normalize=="balance")
                    {
                        cat("  Counting tags in the raw control track ",basename(rawc)," over the corresponding peak file ",nam,"... Please wait...\n",sep="")
                        ovObj.rawc[[nam]] <- countOverlaps(input.bed[[nam]],rawc.bed,minoverlap=ovObj.taglen[[nam]]$control-1)
                    }
                    if (normalize=="peakseq.original" || normalize=="peakseq.rlm")
                    {
                        cat("  Counting control tags in the control track ",control," over genomic windows for peakseq... Please wait...\n",sep="")
                        ovObj.peakseq[[nam]]$control <- countOverlaps(gen.win,control.bed)
                        names(ovObj.peakseq[[nam]]$control) <- gen.win.names
                    }
                }
            }
            
            gc(verbose=FALSE)
        }
            
        # Now that we have the counts, we have to put them in the peak files     
        treat.counts <- control.counts <- #rank.sums <- 
        count.matrix <- count.matrix.norm <-
        avg.counts.matrix <- avg.counts.matrix.norm <- 
        ratios <- ratios.norm <- saturation <- list()
        for (nam in basename(input))
        {
            treat.counts[[nam]] <- unlist(ovObj.treat[[nam]])
            saturation[[nam]]$treatment <- treat.counts[[nam]]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$treatment + 1)
            if (hasControl)
            {
                control.counts[[nam]] <- unlist(ovObj.control[[nam]])
                count.matrix[[nam]] <- cbind(treat.counts[[nam]],control.counts[[nam]])
                # Ratios
                ratios[[nam]] <- 
                    log2(count.matrix[[nam]][,1]/ifelse(count.matrix[[nam]][,2]==0,1,count.matrix[[nam]][,2]))
                # Average tags
                diffs <- all.peaks[[nam]]$end - all.peaks[[nam]]$start + 1
                avg.counts.matrix[[nam]] <- count.matrix[[nam]]/(diffs/avg.win)
                # Saturation index
                saturation[[nam]]$control <- control.counts[[nam]]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$control + 1)
                # Normalization
                if (normalize=="linear")
                {
                    if (ovObj.nreads[[nam]]$treatment > ovObj.nreads[[nam]]$control)
                    {
                        scale.factor <- ovObj.nreads[[nam]]$control/ovObj.nreads[[nam]]$treatment
                        count.matrix.norm[[nam]] <- cbind(round(scale.factor*count.matrix[[nam]][,1]),count.matrix[[nam]][,2])
                    }
                    else if (ovObj.nreads[[nam]]$control > ovObj.nreads[[nam]]$treatment)
                    {
                        scale.factor <- ovObj.nreads[[nam]]$treatment/ovObj.nreads[[nam]]$control
                        count.matrix.norm[[nam]] <- cbind(count.matrix[[nam]][,1],round(scale.factor*count.matrix[[nam]][,2]))
                    }
                    avg.counts.matrix.norm[[nam]] <- count.matrix.norm[[nam]]/(diffs/avg.win)
                    ratios.norm[[nam]] <- log2(count.matrix.norm[[nam]][,1]/ifelse(count.matrix.norm[[nam]][,2]==0,1,count.matrix.norm[[nam]][,2]))
                }
                else if (normalize=="balance") # Has been externally normalized, so swap...
                {
                    count.matrix.norm[[nam]] <- count.matrix[[nam]]
                    avg.counts.matrix.norm[[nam]] <- count.matrix.norm[[nam]]/(diffs/avg.win)
                    ratios.norm[[nam]] <- log2(count.matrix.norm[[nam]][,1]/ifelse(count.matrix.norm[[nam]][,2]==0,1,count.matrix.norm[[nam]][,2]))

                    count.matrix[[nam]] <- cbind(unlist(ovObj.rawt[[nam]]),unlist(ovObj.rawc[[nam]]))
                    avg.counts.matrix[[nam]] <- count.matrix[[nam]]/(diffs/avg.win)
                    ratios[[nam]] <- 
                        log2(count.matrix[[nam]][,1]/ifelse(count.matrix[[nam]][,2]==0,1,count.matrix[[nam]][,2]))

                    saturation[[nam]]$treatment <- count.matrix[[nam]][,1]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$treatment + 1)
                    saturation[[nam]]$control <- count.matrix[[nam]][,2]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$treatment + 1)
                }
                else if (normalize=="rpkm")
                {
                    count.matrix.norm[[nam]] <- cbind(count.matrix[[nam]][,1]*1e+6/ovObj.nreads[[nam]]$treatment,
                                                       count.matrix[[nam]][,2]*1e+6/ovObj.nreads[[nam]]$control)
                    avg.counts.matrix.norm[[nam]] <- count.matrix.norm[[nam]]*1e+6/(diffs/avg.win)
                    ratios.norm[[nam]] <- log2(count.matrix.norm[[nam]][,1]/ifelse(count.matrix.norm[[nam]][,2]==0,1,count.matrix.norm[[nam]][,2]))
                }
                else if (!is.na(pmatch("peakseq",normalize)))
                {
                    if (normalize=="peakseq.original")
                    {
                        # 1. Remove dead zones for each normalization pair
                        alive <- which(ovObj.peakseq[[nam]]$control!=0 | ovObj.peakseq[[nam]]$treatment!=0)
                        x <- ovObj.peakseq[[nam]]$control[alive]
                        y <- ovObj.peakseq[[nam]]$treatment[alive]
                        # 2. Calculate the vector of 10kb lamdas based on the alive treatment read distribution
                        #    Could be control, but now it's more conservative
                        lamda <- (ovObj.nreads[[nam]]$treatment*(y/win.size))/(eff.gen.size*win.size/attr(chrom.info,"genome.size"))
                        # 3. Exclude low quality windows using the above Poisson rate on the alive zones (p<1e-6)
                        exp.tag <- qpois(1-1e-6,lambda=lamda)
                        xx <- x[which(x > exp.tag & y > exp.tag)]
                        yy <- y[which(x > exp.tag & y > exp.tag)]
                        # TODO: Exclude low and high quartiles of the ratio to create unbiased curve
                        # TODO: Treat somehow the repeat areas. Look at notes for plotting those regions.
                        model <- lm(y ~ x)
                        #model$coefficients[2]*
                    }
                    else if (normalize=="peakseq.rlm")
                    {
                        x <- ovObj.peakseq[[nam]]$control
                        y <- ovObj.peakseq[[nam]]$treatment
                        rlm(x,y,maxit=100,...)
                    }
                }
                else if (normalize=="none")
                {
                    count.matrix.norm[[nam]] <- NULL
                    avg.counts.matrix.norm[[nam]] <- NULL
                    ratios.norm[[nam]] <- NULL
                }

                if (fc.cut!=0 && !is.null(fc.cut) && !is.na(fc.cut))
                {
                    if (normalize=="none")
                        fc.fail[[nam]] <- which(ratios[[nam]]<fc.cut)
                    else
                        fc.fail[[nam]] <- which(ratios.norm[[nam]]<fc.cut)
                }
                else fc.fail[[nam]] <- NULL

                ## Saturation filter
                #if (!is.null(sat))
                #{
                #   sat.fail[[nam]] <-
                #       which(saturation[[nam]]$treatment<sat[2] | saturation[[nam]]$control>sat[1])
                #       #which((saturation[[nam]]$treatment/saturation[[nam]]$max)<sat[2] & (saturation[[nam]]$control/saturation[[nam]]$max>sat[1]))
                #}
                #else sat.fail[[nam]] <- NULL
            }
            else # Ratio cannot be calculated and categorization cannot be done
            {
                count.matrix[[nam]] <- as.matrix(treat.counts[[nam]])
                diffs <- all.peaks[[nam]]$end - all.peaks[[nam]]$start + 1
                avg.counts.matrix[[nam]] <- as.matrix(count.matrix[[nam]][,1]/(diffs/avg.win))
                saturation[[nam]]$treatment <- count.matrix[[nam]][,1]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$treatment + 1)
                ratios[[nam]] <- NULL
                fc.fail[[nam]] <- NULL
                #rank.sums[[nam]] <- 1/ifelse(avg.counts.treat[[nam]])
                bin.rr <- NULL
                
                if (normalize=="balance") # Has been externally normalized, so swap...
                {
                    count.matrix.norm[[nam]] <- count.matrix[[nam]]
                    avg.counts.matrix.norm[[nam]] <- count.matrix.norm[[nam]]/(diffs/avg.win)
                    ratios.norm[[nam]] <- NULL

                    count.matrix[[nam]] <- as.matrix(unlist(ovObj.rawt[[nam]]))
                    avg.counts.matrix[[nam]] <- count.matrix[[nam]]/(diffs/avg.win)
                    saturation[[nam]]$treatment <- count.matrix[[nam]][,1]/(all.peaks[[nam]]$length - ovObj.taglen[[nam]]$treatment + 1)
                }
                if (normalize=="rpkm")
                {
                    count.matrix.norm[[nam]] <- count.matrix[[nam]][,1]*1e+6/ovObj.nreads[[nam]]$treatment
                    avg.counts.matrix.norm[[nam]] <- count.matrix.norm[[nam]][,1]*1e+6/(diffs/avg.win)
                    ratios.norm[[nam]] <- NULL
                }
                else if (normalize=="peakseq")
                {
                    # TODO
                }
                else if (normalize=="none")
                {
                    count.matrix.norm[[nam]] <- NULL
                    avg.counts.matrix.norm[[nam]] <- NULL
                    ratios.norm[[nam]] <- NULL
                }               
            }
        }

        # Saturation filter
        p.index <- 1:nrow(all.peaks[[nam]])
        if (!is.null(sat))
        {
            if (hasControl)
            {
                good <- which(saturation[[nam]]$treatment>sat[2] & saturation[[nam]]$control<sat[1])
                #good <- which((saturation[[nam]]$treatment/saturation[[nam]]$max)>sat[2] & (saturation[[nam]]$control/saturation[[nam]]$max)<sat[1])
            }
            else
            {
                good <- which(saturation[[nam]]$treatment>sat[2])
            }
            if (length(good)!=0)
                sat.fail[[nam]] <- p.index[-good]
            else sat.fail[[nam]] <- NULL
        }
        else sat.fail[[nam]] <- NULL

        # Unify all the filters
        for (nam in basename(input))
            filter.fail[[nam]] <- base:::union(base:::union(fdr.fail[[nam]],fc.fail[[nam]]),sat.fail[[nam]]) # Masked by IRanges... Argh!
        
        # Categorize final distributions according to ratios and ranksums
        rat.flags <- rank.flags <- iqc.flags <- rat.qnt <- rank.qnt <- list()
        if (!is.null(bin.rr) && hasControl)
        {
            for (nam in basename(input))
            {
                if (length(filter.fail[[nam]]) != 0)
                {
                    good.rats <- ratios[[nam]][-filter.fail[[nam]]]
                    #good.ranks <- rank.sums[[nam]][-filter.fail[[nam]]]
                    good.ranks <- avg.counts.matrix[[nam]][-filter.fail[[nam]],1]
                    rat.flags[[nam]] <- rank.flags[[nam]] <- iqc.flags[[nam]] <- numeric(length(good.rats))
                    rat.pct <- c(0,bin.rr$ratio/100,1)
                    rank.pct <- c(0,bin.rr$rpw/100,1)
                    rat.qnt[[nam]] <- quantile(good.rats,rat.pct)
                    rank.qnt[[nam]] <- quantile(good.ranks,rank.pct)
                    for (i in 2:length(rat.pct))
                        rat.flags[[nam]][which(good.rats>=rat.qnt[[nam]][i-1] & good.rats<=rat.qnt[[nam]][i])] <- i-1
                    for (i in 2:length(rank.pct))
                        #rank.flags[[nam]][which(good.ranks>=rank.qnt[[nam]][i-1] & good.ranks<=rank.qnt[[nam]][i])] <- length(rank.pct)-i+1
                        rank.flags[[nam]][which(good.ranks>=rank.qnt[[nam]][i-1] & good.ranks<=rank.qnt[[nam]][i])] <- i-1
                    iqc.flags[[nam]] <- as.matrix(paste(rat.flags[[nam]],rank.flags[[nam]],sep=""))
                }
                else
                {
                    good.rats <- ratios[[nam]]
                    #good.ranks <- rank.sums[[nam]]
                    good.ranks <- avg.counts.matrix[[nam]][,1]
                    rat.flags[[nam]] <- rank.flags[[nam]] <- iqc.flags[[nam]] <- numeric(length(good.rats))
                    rat.pct <- c(0,bin.rr$ratio/100,1)
                    rank.pct <- c(0,bin.rr$rpw/100,1)
                    rat.qnt[[nam]] <- quantile(good.rats,rat.pct)
                    rank.qnt[[nam]] <- quantile(good.ranks,rank.pct)
                    for (i in 2:length(rat.pct))
                        rat.flags[[nam]][which(good.rats>=rat.qnt[[nam]][i-1] & good.rats<=rat.qnt[[nam]][i])] <- i-1
                    for (i in 2:length(rank.pct))
                        #rank.flags[[nam]][which(good.ranks>=rank.qnt[[nam]][i-1] & good.ranks<=rank.qnt[[nam]][i])] <- length(rank.pct)-i+1
                        rank.flags[[nam]][which(good.ranks>=rank.qnt[[nam]][i-1] & good.ranks<=rank.qnt[[nam]][i])] <- i-1
                    iqc.flags[[nam]] <- as.matrix(paste(rat.flags[[nam]],rank.flags[[nam]],sep=""))
                }
            }
        }
        else 
        {
            warning("Binning areas not defined or ratios not calculated! Skipping categorization...",
                    immediate.=TRUE,call.=FALSE)
            iqc.flags[[nam]] <- NULL
            plot.rvr <- FALSE
        }
        
        # Construct the additional stats matrix
        add.matrix <- list()
        for (nam in basename(input))
        {
            add.matrix[[nam]] <- 
            cbind(
                avg.counts.matrix[[nam]],
                ratios[[nam]]
                #rank.sums[[nam]]
            )
        }
        
        #...and construct the output
        finalObj <- negObj <- list()
        for (nam in basename(input))
        {    
            if (length(filter.fail[[nam]]) != 0)
            {
                finalObj[[nam]] <- 
                data.frame(
                    #paste("chr",all.peaks[[nam]]$chromosome[-filter.fail[[nam]]],sep=""),
                    all.peaks[[nam]]$chromosome[-filter.fail[[nam]]],
                    all.peaks[[nam]]$start[-filter.fail[[nam]]],
                    all.peaks[[nam]]$end[-filter.fail[[nam]]],
                    paste(paste(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]],
                                all.peaks[[nam]]$start[-filter.fail[[nam]]],sep=":"),
                                all.peaks[[nam]]$end[-filter.fail[[nam]]],sep="-"),
                    all.peaks[[nam]]$length[-filter.fail[[nam]]],
                    if (!is.null(all.peaks[[nam]]$summit))
                        all.peaks[[nam]]$summit[-filter.fail[[nam]]]
                    else
                        rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])),
                    count.matrix[[nam]][-filter.fail[[nam]],],
                    if (!is.null(count.matrix.norm[[nam]])) count.matrix.norm[[nam]][-filter.fail[[nam]],] else matrix(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]]),ifelse(hasControl,2,1)),
                    avg.counts.matrix[[nam]][-filter.fail[[nam]],],
                    if (!is.null(avg.counts.matrix.norm[[nam]])) avg.counts.matrix.norm[[nam]][-filter.fail[[nam]],] else matrix(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]]),ifelse(hasControl,2,1)),
                    if (!is.null(ratios[[nam]])) ratios[[nam]][-filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])),
                    if (!is.null(ratios.norm[[nam]])) ratios.norm[[nam]][-filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])),
                    all.peaks[[nam]]$significance[-filter.fail[[nam]]],
                    all.peaks[[nam]]$fdr[-filter.fail[[nam]]],
                    #add.matrix[[nam]][-filter.fail[[nam]],], # Avg counts, Ratios, RS
                    if (hasControl && !is.null(iqc.flags[[nam]])) iqc.flags[[nam]] else rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])), # Flags
                    if (ver == 2) all.peaks[[nam]]$pileup[-filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])),
                    #rank.sums[[nam]][-filter.fail[[nam]]],
                    #saturation[[nam]]$max[-filter.fail[[nam]]],
                    saturation[[nam]]$treatment[-filter.fail[[nam]]],
                    if (hasControl) saturation[[nam]]$control[-filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[-filter.fail[[nam]]])),
                    all.peaks[[nam]]$MACS_fe[-filter.fail[[nam]]],
                    stringsAsFactors=FALSE
                )

                if (export.negative)
                {
                    negObj[[nam]] <-
                        data.frame(
                        all.peaks[[nam]]$chromosome[filter.fail[[nam]]],
                        all.peaks[[nam]]$start[filter.fail[[nam]]],
                        all.peaks[[nam]]$end[filter.fail[[nam]]],
                        paste(paste(all.peaks[[nam]]$chromosome[filter.fail[[nam]]],
                                    all.peaks[[nam]]$start[filter.fail[[nam]]],sep=":"),
                                    all.peaks[[nam]]$end[filter.fail[[nam]]],sep="-"),
                        all.peaks[[nam]]$length[filter.fail[[nam]]],
                        if (!is.null(all.peaks[[nam]]$summit))
                            all.peaks[[nam]]$summit[filter.fail[[nam]]]
                        else
                            rep(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]])),
                        count.matrix[[nam]][filter.fail[[nam]],],
                        if (!is.null(count.matrix.norm[[nam]])) count.matrix.norm[[nam]][filter.fail[[nam]],] else matrix(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]]),2),
                        avg.counts.matrix[[nam]][filter.fail[[nam]],],
                        if (!is.null(avg.counts.matrix.norm[[nam]])) avg.counts.matrix.norm[[nam]][filter.fail[[nam]],] else matrix(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]]),2),
                        if (!is.null(ratios[[nam]])) ratios[[nam]][filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]])),
                        if (!is.null(ratios.norm[[nam]])) ratios.norm[[nam]][filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]])),
                        all.peaks[[nam]]$significance[filter.fail[[nam]]],
                        all.peaks[[nam]]$fdr[filter.fail[[nam]]],
                        if (ver == 2) all.peaks[[nam]]$pileup[filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]])),
                        #saturation[[nam]]$max[filter.fail[[nam]]],
                        saturation[[nam]]$treatment[filter.fail[[nam]]],
                        if (hasControl) saturation[[nam]]$control[filter.fail[[nam]]] else rep(NA,length(all.peaks[[nam]]$chromosome[filter.fail[[nam]]])),
                        all.peaks[[nam]]$MACS_fe[filter.fail[[nam]]],
                        stringsAsFactors=FALSE
                    )
                }
            }
            else
            {
                finalObj[[nam]] <- 
                data.frame(
                    #paste("chr",all.peaks[[nam]]$chromosome,sep=""),
                    all.peaks[[nam]]$chromosome,
                    all.peaks[[nam]]$start,
                    all.peaks[[nam]]$end,
                    paste(paste(all.peaks[[nam]]$chromosome,
                                all.peaks[[nam]]$start,sep=":"),
                                all.peaks[[nam]]$end,sep="-"),
                    all.peaks[[nam]]$length,
                    if (!is.null(all.peaks[[nam]]$summit))
                        all.peaks[[nam]]$summit
                    else
                        rep(NA,length(all.peaks[[nam]]$chromosome)),
                    count.matrix[[nam]],
                    if (!is.null(count.matrix.norm[[nam]])) count.matrix.norm[[nam]] else matrix(NA,length(all.peaks[[nam]]$chromosome),ifelse(hasControl,2,1)),
                    avg.counts.matrix[[nam]],
                    if (!is.null(avg.counts.matrix.norm[[nam]])) avg.counts.matrix.norm[[nam]] else matrix(NA,length(all.peaks[[nam]]$chromosome),ifelse(hasControl,2,1)),
                    if (!is.null(ratios[[nam]])) ratios[[nam]] else rep(NA,length(all.peaks[[nam]]$chromosome)),
                    if (!is.null(ratios.norm[[nam]])) ratios.norm[[nam]] else rep(NA,length(all.peaks[[nam]]$chromosome)),
                    all.peaks[[nam]]$significance,
                    all.peaks[[nam]]$fdr,
                    #add.matrix[[nam]][-filter.fail[[nam]],], # Avg counts, Ratios, RS
                    if (hasControl) iqc.flags[[nam]] else rep(NA,length(all.peaks[[nam]]$chromosome)), # Flags
                    if (ver == 2) all.peaks[[nam]]$pileup else rep(NA,length(all.peaks[[nam]]$chromosome)),
                    #saturation[[nam]]$max,
                    saturation[[nam]]$treatment,
                    if (hasControl && !is.null(iqc.flags[[nam]])) saturation[[nam]]$control else rep(NA,length(all.peaks[[nam]]$chromosome)),
                    all.peaks[[nam]]$MACS_fe,
                    stringsAsFactors=FALSE
                )
            }
          
            the.names <- 
                c("chromosome",
                  "start",
                  "end",
                  "id",
                  "length",
                  "summit",
                  "reads_treatment",
                  if (!hasControl) NULL else "reads_control",
                  if (normalize=="none") "na.1" else "reads_treatment_norm",
                  if (normalize=="none" & !hasControl) NULL else if (normalize=="none" & hasControl) "na.2" else "reads_control_norm",
                  "avg_reads_treatment",
                  if (!hasControl) NULL else "avg_reads_control",
                  if (normalize=="none") "na.3" else "avg_reads_treatment_norm",
                  if (normalize=="none" & !hasControl) NULL else if (normalize=="none" & hasControl) "na.4" else "avg_reads_control_norm",
                  if (!hasControl) "na.5" else "fold_enrichment",
                  if (normalize=="none" || !hasControl) "na.6" else "fold_enrichment_norm",
                  "significance",
                  "fdr",
                  if ((!hasControl && is.null(bin.rr)) || !hasControl) "na.7" else "qc_flag",
                  if (ver == 14) "na.8" else "pileup",
                  #"saturation_max",
                  "saturation_treatment",
                  if (hasControl) "saturation_control" else "na.9",
                  "MACS_fe")
            names(finalObj[[nam]]) <- the.names
            if (any(grep("na",names(finalObj[[nam]]))))
                finalObj[[nam]] <- finalObj[[nam]][,-grep("na",names(finalObj[[nam]]))]
            if (export.negative)
                if (any(grep("qc_flag",names(finalObj[[nam]]))))
                    names(negObj[[nam]]) <- names(finalObj[[nam]])[-grep("qc_flag",names(finalObj[[nam]]))]
        }
    }
    else # Just filter for FDR, no bed files for counting..
    {
        if (length(fdr.fail[[nam]]) != 0)
        {
            for (nam in basename(input))
            {
                finalObj[[nam]] <- all.peaks[[nam]][-fdr.fail[[nam]],]
                if (export.negative)
                    negObj[[nam]] <- all.peaks[[nam]][fdr.fail[[nam]],]
            }
        }
        else
        {
            for (nam in basename(input))
            {
                finalObj[[nam]] <- all.peaks[[nam]]
            }
        }
    }
       
    if (write.output)
    {
        new.file <- neg.file <- list()
        if (output=="auto")
        {
            for (nam in names(finalObj))
            {
                f <- unlist(strsplit(nam,".",fixed=TRUE))
                l.f <- length(f)
                n.n <- character(0)
                if (l.f==1)
                    n.n <- f
                else
                {
                  t.n <- f[1:(l.f-1)]
                  n.n <- t.n[1]
                  if (l.f>2)
                    for (j in 2:(l.f-1))
                      n.n <- paste(n.n,t.n[j],sep=".")
                }
                new.file[[nam]] <- paste(file.path(path[[nam]],n.n),".txt",sep="")
                neg.file[[nam]] <- paste(file.path(path[[nam]],n.n),"_negative.txt",sep="")
                cat("Writing file",basename(new.file[[nam]]),"\n")
                flush.console()
                write.table(finalObj[[nam]],file=new.file[[nam]],quote=FALSE,sep="\t",row.names=FALSE)
                if (export.negative)
                {
                    cat("Writing file",basename(neg.file[[nam]]),"\n")
                    flush.console()
                    write.table(negObj[[nam]],file=neg.file[[nam]],quote=FALSE,sep="\t",row.names=FALSE)
                }
            }
        }
        else
        {
            nams <- names(finalObj)
            for (i in 1:length(output))
            {
                cat("Writing file",output[i],"\n")
                flush.console()
                write.table(finalObj[[nams[i]]],file=output[i],quote=FALSE,sep="\t",row.names=FALSE)
            }
            # TODO: Write something about exporting negative files when the output is given
        }
    }
  
    # Plot if required
    if (plot.svr)
    {
        for (nam in names(finalObj))
        {
            plotSatRat(saturation[[nam]]$treatment,ratios[[nam]],sfilt=sat.fail[[nam]],rfilt=fc.fail[[nam]],
                       main.title=paste(nam,"svr",sep="."),out=image.format,path=path[[nam]])
        }
    }
    if (plot.rvr && !is.null(bin.rr))
    {
        for (nam in names(finalObj))
        {
            plotRankRatMulti(finalObj[[nam]]$avg_counts_treatment,finalObj[[nam]]$fold_enrichment,iqc.flags[[nam]],
                             rank.qnt[[nam]],rat.qnt[[nam]],main.title=paste(nam,"rvr",sep="."),out=image.format,path=path[[nam]])
        }
    }
  
    return(finalObj)   
}

# Plot saturation against ratio
plotSatRat <- function(sat,rat,sfilt=NULL,rfilt=NULL,main.title="data",out="x11",path=NA)
{
    if (out=="pdf")
    {
        if (!is.na(path))
            o <- paste(file.path(path,main.title),".pdf",sep="")
        else
            o <- paste(main.title,".pdf",sep="")
        pdf(file=o,width=8,height=6)
        sex <- 0.5
        sex.axis <- 1.8
        sex.main <- 2.5
        sex.lab <- 1.8
        sex.leg <- 1.7
        lsd <- 1.7
    } 
    else
    {
        sex <- 0.7
        sex.axis <- 1.3
        sex.main <- 1.5
        sex.lab <- 1.3
        sex.leg <- 1.3
        lsd <- 2
    }
  
    if (out=="png")
    {
        if (!is.na(path))
            o <- paste(file.path(path,main.title),".png",sep="")
        else
            o <- paste(main.title,".png",sep="")
        png(filename=o,width=640,height=480)
    }
    if (out=="x11") x11()

    leg.text <- c("Remaining")
    leg.col <- c("blue")

    par(cex=sex,cex.axis=sex.axis,cex.main=sex.main,cex.lab=sex.lab,font.lab=2)
    plot(sat,rat,pch=20,col="blue",main=paste("Saturation vs Enrichment plot for",main.title),
        xlab="Tag saturation in peak area",ylab="log2(Enrichment)")
    if (!is.null(rfilt))
    {
        points(sat[rfilt],rat[rfilt],pch=20,col="red")
        leg.text <- c(leg.text,"R filtered")
        leg.col <- c(leg.col,"red")
    }
    if (!is.null(sfilt))
    {
        points(sat[sfilt],rat[sfilt],pch=20,col="orange")
        leg.text <- c(leg.text,"S filtered")
        leg.col <- c(leg.col,"orange")
    }
    grid()
    
    # Legend stuff, put outside plotting region on the right
    legend("topright",legend=leg.text,col=leg.col,text.col=leg.col,cex=sex.leg,pch=20,
          box.lwd=0.5,xjust=0.4,yjust=0.5)

    if (out!="x11") dev.off()
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
    }
    else
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
plotRankRatMulti <- function(ranks,ratios,flags,rank.qnt,rat.qnt,main.title="data",
                             out="x11",path=NA)
{
    if (out=="pdf")
    {
        if (!is.na(path))
            o <- paste(file.path(path,main.title),".pdf",sep="")
        else
            o <- paste(main.title,".pdf",sep="")
        pdf(file=o,width=8,height=6)
        sex <- 0.5
        sex.axis <- 1.8
        sex.main <- 2.5
        sex.lab <- 1.8
        sex.leg <- 1.7
        lsd <- 1.7
    } 
    else
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
  
    if (out=="png")
    {
        if (!is.na(path))
            o <- paste(file.path(path,main.title),".png",sep="")
        else
            o <- paste(main.title,".png",sep="")
        png(filename=o,width=640,height=480)
    }
    if (out=="x11") x11()

    un.flags <- sort(unique(flags),decreasing=TRUE)
    par(cex=sex,cex.axis=sex.axis,cex.main=sex.main,cex.lab=sex.lab,font.lab=2)
    plot.new()
    plot.window(c(range(ranks)[1],range(ranks)[2]+0.025),range(ratios))
    axis(1,at=pretty(ranks))
    axis(2,at=pretty(ratios))
    title(main=paste("RankSum vs Ratio plot for",main.title),xlab="RankSum",ylab="log2(Ratio)")

    for (i in 1:length(un.flags))
        points(ranks[which(flags==un.flags[i])],ratios[which(flags==un.flags[i])],
               pch=20,col=fixed.colors[i])

    for (i in 1:length(ab.rank))
        abline(v=ab.rank[i],lty="dashed",col="green",lwd=lsd)
    for (i in 1:length(ab.rat))
        lines(seq(range(ranks)[1],range(ranks)[2],length.out=100),
              rep(ab.rat[i],100),lty="dashed",col="red",lwd=lsd)

    # Legend stuff, put outside plotting region on the right
    old.par<-par(xpd=TRUE)
    y.pos<-mean(par("usr")[3:4])
    in2x.pos<-(par("usr")[2]-par("usr")[1])/par("pin")[1]
    x.pos<-par("usr")[2]-par("mai")[4]*in2x.pos
    legend(x.pos,y.pos,legend=un.flags,col=fixed.colors,text.col=fixed.colors,
          cex=sex.leg,pch=20,box.lwd=0.5,xjust=0.4,yjust=0.5)
    par(old.par)

    if (out!="x11") dev.off()
}

make.windows <- function(chrlen,winsize,chrname)
{
    names(chrlen) <- chrname
    winObj <- GRangesList()
    rem <- chrlen%%winsize
    len <- floor(chrlen/winsize)
    clen <- len + 1
    snames <- Rle(chrname,clen)
  
    starts <- ends <- vector("list",length(chrlen))
    for (i in 1:length(chrlen))
    { 
        starts[[i]] <- c(seq(1,by=winsize,length=len[i]),len[i]*winsize + 1)
        ends[[i]] <- c(seq(winsize,by=winsize,length=len[i]),len[i]*winsize + rem[i])
    }
    winObj <- GRanges(seqnames=snames,
                    IRanges(start=unlist(starts),end=unlist(ends)),
                    seqlengths=chrlen)
    return(winObj)
}

# x: vector of window counts or a length
make.sliding <- function(x,win,slide,type="list")
{
    len = ifelse(length(x)==1,x,length(x))
    nams = ifelse(is.null(attr(x,"names")),NULL,names(x))
    st <- seq(1,len-win,by=slide)
    en <- seq(win+1,len,by=slide)
    last.slide <- len - en[length(en)]
    st <- c(st,st[length(st)]+last.slide)
    en <- c(en,en[length(en)]+last.slide)

    
    
    if (type=="IRanges")
        return(IRanges(start=st,end=en))
    else
        return(list(start=st,end=en))
}

calc.ksp <- function(x,y,st,en)
{
    if (require(utils))
    {
        pb.available <- TRUE
        pb <- txtProgressBar(min=0,max=length(st),style = 3)
    }
    else pb.available <- FALSE
    for (i in 1:length(st))
    {
        if (pb.available) setTxtProgressBar(pb,i)
        ks[[i]] <- ks.test(x[st[i]:en[i]],y[st[i]:en[i]])
        pvalue[i] <- ks[[i]]$p.value
    }
    if (pb.available) close(pb)
}

getChromInfo <- function(org)
{
    download.file(paste("http://hgdownload.cse.ucsc.edu/goldenPath/",org,"/database/chromInfo.txt.gz",sep=""),
                file.path(tempdir(),"chromInfo.txt.gz"))
    if (.Platform$OS.type == "unix")
        system(paste("gzip -df",file.path(tempdir(),"chromInfo.txt.gz")))
    else
        unzip(file.path(tempdir(),"chromInfo.txt.gz"))
    chrom.info <- read.delim(file.path(tempdir(),"chromInfo.txt"),header=FALSE)
    genome.size <- sum(as.numeric(chrom.info[,2]))
    chrom.info <- chrom.info[-grep("rand|chrM|hap",chrom.info[,1]),]
    attr(chrom.info,"genome.size") <- genome.size
    return(chrom.info)
}

split.data <- function(x,n)
{
    l <- vector("list",n)
    if (!is.null(dim(x))) # x is a matrix/data.frame etc.
    {
        e <- ceiling(nrow(x)/n)
        for (i in 1:n)
            l[[i]] <- rep(i,e)
        f <- unlist(l)
        f <- f[1:nrow(x)]
        y <- split(as.data.frame(x),f)
        for (i in 1:n)
            y[[i]] <- as.matrix(y[[i]])
        return(y)
    }
    else # vector
    {
        e <- ceiling(length(x)/n)
        for (i in 1:n)
            l[[i]] <- rep(i,e)
        f <- unlist(l)
        f <- f[1:length(x)]
        return(split(x,f))
    }
}

find.p.peaks <- function(pvalue,estimator,SNR)
{
    if (!require(MALDIquant))
        stop("Bioconductor package MALDIquant is required!")
    s <- createMassSpectrum(mass=1:length(pvalue),intensity=pvalue,metaData=list(name="KS pvalue"));
    n <- estimateNoiseSmoother(mass(s),intensity(s),estimator=estimator)
    p <- detectPeaks(s,halfWindowSize=1,fun=estimateNoiseSmoother,SNR=SNR,estimator=estimator)
    return(mass(p))
}

estimateNoiseSmoother <- function(x,y,estimator="super",...)
{
    if (estimator=="super")
    {
        s <- supsmu(x=x,y=y,...);
        return(cbind(x,s$y));
    }
    else if (estimator=="mean") # xcms style
    {
        s <- mean(y)
        return(cbind(x,s))
    }
}

#rank.sums[[nam]] <- 
#   1/(avg.counts.matrix[[nam]][,1] + ifelse(avg.counts.matrix[[nam]][,2]==0,1,avg.counts.matrix[[nam]][,2]))

# Chunk for possible use of tagcounter.pl
#if (!is.na(mapfile))
#{
#   # First, import the peak files and store them, as they are small
#   input.bed <- output.bed <- ovObj.treat <- ovObj.control <- ovObj.rawt <- ovObj.rawc <-
#   ovObj.nreads <- ovObj.peakseq <- ovObj.taglen <- list()
#   for (i in 1:n)
#   {
#       cat("Importing peak file ",basename(input[i])," as bed... Please wait...\n",sep="")
#       flush.console()
#       input.bed[[basename(input[[i]])]] <- tempfile()
#       output.bed[[basename(input[[i]])]] <- tempfile()
#       tmp.frame <- cbind(all.peaks[[basename(input[[i]])]][,1:4])
#       write.table(tmp.frame,input.bed[[basename(input[[i]])]],sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#       perl tagcounter.pl --input map$treatment[i] map$control[i] map$rawt[i] map$rawc[i] --region input.bed[[basename(input[[i]])]]
#           --percent 0.5 --ncore 4 --small --output output.bed[[basename(input[[i]])]]

#       nam <- basename(map$filename[i])
#       the.peaks <- read.delim(output.bed[[basename(input[[i]])]])
#       ovObj.treat[[nam]] <- the.peaks[,basename(map$treatment[i])]
#       if (normalize=="balance")
#       {
#           ovObj.rawt[[nam]] <- the.peaks[,basename(map$rawt[i])]
#       }
#       if (hasControl)
#       {
#           ovObj.control[[nam]] <- the.peaks[,basename(map$control[i])]
#           if (normalize=="balance")
#           {
#               ovObj.rawc[[nam]] <- the.peaks[,basename(map$rawc[i])]
#           }
#       }
#   }
#}
