# This R script implements the TIP method presented in the Bioinformatics 2011, 
# 27(23), 3221-3227:
# TIP: A probabilistic method for identifying transcription factor target genes 
# from ChIP-seq binding proﬁles
# Chao Cheng, Renqiang Min and Mark Gerstein
#
# TIP associates each gene with a regulatory score based on the binding profiles 
# of the transcription factors (TFs) of interest. For more details on the 
# algorithm, results and efficiency please see the above publication. 
# In addition, we have extended the algorithm presented in the paper in the 
# following ways:
#  a) Apart from the standard normal distribution used to calculate the 
#  significance of z-scores (as in the paper), we added three additional options,
#  useful for when the promoter region is asymmetric (in the paper they use
#  +/-10kb from TSS while we prefer (-10kb,+1kb)). In the case of asymmetric
#  promoter regions, the normality of the z-scores is violated as the 
#  distribution becomes right-skewed.
#    i)  The significance may be calculated by estimating an F null distribution, 
#    using the non z-score normalized values of the regulatory score to estimate 
#    the parameters of the F distribution (package MASS, function fitdistr).
#    ii) The significance is estimated based on the standard normal distribution 
#    but calculating the z-score on log2 transformed regulatory scores, to 
#    restore normality.
#    iii) The significance is calculated using Monte Carlo (MC) simulations, 
#    where in each MC iteration, a random set of genes is selected from the 
#    organism pool and the binding signal is calucated for these random gene 
#    sets. The final significance is the number of times that the regulatory 
#    score is larger in MC simulations than the original one.
#  b) Apart from the Storey and Tibshirani qvalue multiple testing correaction, 
#  the classical Benajamini-Hochberg FDR correction is also available.
#
# Usage:
# tip <- TIP(
#           peak.files,
#           signal.files,
#           gene.file,
#           signal.type="bed",
#           promoter=c(-10000,1000),
#           frag.size=200,
#           method=c("cmg","normal","f","mc"),
#           correct=TRUE,
#           correct.method=c("qv","bh"),
#           cmg.zcut=1.96,
#           mc.iter=1000,
#           export=TRUE,
#           output=NA
#       )
# or simply TIP(...)
#
# Parameters:
# peak.files    : A set of peak files for each TF, usually the result of a peak 
#                 calling algorithm, in text tab delimited format WITH a header 
#                 line. The tab delimited format should be similar to the UCSC 
#                 BED format, where the first four columns are (chromosome, 
#                 start, end, peak_id), the 5th column contains a peak score 
#                 (e.g. the peak height which is read pileup in nucleotide with 
#                 the highest coverage, usually output of modern peak-calling 
#                 algorithms). The rest of the columns can contain whatever 
#                 additional information. This parameter is optional but AT 
#                 LEAST ONE of peak.files or signal.files MUST be provided. In
#                 case ONLY peak.files is provided, the the score column MUST
#                 contain peak heights and MUST be named "height", otherwise the
#                 script will probably crash or not return anticipated results.
# signal.files  : A set of mapped read files in BED format, where their input 
#                 order corresponds to the input order of peak.files, if 
#                 provided (e.g. if peak.files = c("A.peak","B.peak"), 
#                 signal.files = c("A.bed","B.bed"), where A and B represent the
#                 order and not the naming of the files, e.g. signal.file = 
#                 c("B.bed","A.bed") is WRONG). This parameter is optional but
#                 AT LEAST ONE of peak.files or signal.files MUST be provided.
# gene.file     : One file containing the organism's annotation in BED format 
#                 WITH a header line. The columns should be as follows:
#                 (chromosome, start, end, gene_id, something, strand). All 
#                 additional columns can contain whatever other information
#                 (e.g. expression values, additional annotation) which will be 
#                 appended to the output. The program can benefit from a 7th
#                 column containing HUGO gene symbols. Such files can be easily 
#                 extracted using e.g. Biomart services for Ensembl.
# signal.type:  : The file type of the signal files. It can be "bed" (default)
#                 or "bam".
# promoter      : The promoter region that will be used for the calculation of 
#                 binding signal and gene scores. It defaults to -10kb upstream 
#                 and +1knb downstream (promoter = c(-10000,1000)). See the 
#                 introduction above for variations of this.
# frag.size     : The original estimated fragment size after sonication 
#                 (defaults to 200bp). The reads from the signal files will be
#                 extended to this size.
# method        : The statistical method which is used to estimate the 
#                 significance of the binding signal (see the paper). It can be:
#                 "cmg" for the original method in the paper ("cmg" from 
#                 Cheng-Min-Gerstein). "normal" for using the standard normal 
#                 distribution with z-scores derived from log2 regulation scores 
#                 (see paper). Useful for assymetric around the TSS promoter 
#                 regions (see introduction above). "f" for using the 
#                 F-distribution to calculate statistical significance of raw 
#                 regulation scores (see paper). Useful for assymetric around 
#                 the TSS promoter regions (see introduction above). "mc" to 
#                 use Monte Carlo simulation for the calculation of statistical 
#                 significance. Useful in all cases but VERY time-consuming, 
#                 even with parallel calulations.
# correct       : Correct for multiple testing? Default to TRUE. If FALSE, FDR 
#                 or q-values will not be reported.
# correct.method: Which method is used for multiple testing correction. Default 
#                 to "qv" for the method used in the paper, or "bh"
#                 for the classical Benjamini-Hochberg FDR correction.
# cmg.zcut      : A z-score cutoff to calculate the False Positive Rate of the 
#                 gene regulatory score assignment (see the paper). It
#                 defaults to 1.96 (p-value < 0.05).
# mc.iter       : Number of Monte Carlo iteration for the "mc" method. It 
#                 defaults to 1000.
# export        : Export the results to a text file with gene scores and 
#                 significance attached? If yes, a number of files equal to the
#                 number of TFs will be generated, else, an R list containing 
#                 the results. Defaults to TRUE.
# output        : The output file names in case export = TRUE. If they are not 
#                 given (NA) or they do not meet criteria (e.g. they must be
#                 equal to the number of input TFs), they will be autogenerated.
# filter.output : If TRUE and peaks are present (peak.files is given), then only
#                 the genes that have peaks in the promoter area will be
#                 exported, otherwise the whole gene.file.
# restrict.cores: A number betyween 0 and 1 which denotes the proportion of the
#                 available processor cores to be used for parallel calculations
#                 in a multicore machine and when package parallel is present.
#                 Defaults to 0.6.
#
# In all cases, the program will run in parallel if multiple cores exist in the 
# running machine and the R package "parallel" is present. Otherwise, a single 
# core will be used.
#
# Examples:
# TIP(
#   peak.files=c("ER_MACS_peaks.txt","GR_MACS_peaks.txt"),
#   signal.files=c("ER_signal.gr","GR_signal.gr"),
#   gene.file="ensembl_hg19_genes.txt",
#   promoter=c(-10000,10000),
#   output=c("ER_gene_scores.txt","GR_gene_scores.txt")
#)
#
# TIP(
#   peak.files=c("ER_MACS_peaks.txt","GR_MACS_peaks.txt"),
#   signal.files=c("ER_signal.gr","GR_signal.gr"),
#   gene.file="ensembl_hg19_genes.txt",
#   method="f"
#)
#
# TIP(
#   peak.files=c("ER_MACS_peaks.txt","GR_MACS_peaks.txt"),
#   signal.files=c("ER_signal.gr","GR_signal.gr"),
#   gene.file="ensembl_hg19_genes.txt",
#   method="mc",
#   mc.iter=10000,
#   output=c("ER_gene_scores.txt","GR_gene_scores.txt")
#)
# tip <- TIP(
#   peak.files=c("ER_MACS_peaks.txt","GR_MACS_peaks.txt"),
#   signal.files=c("ER_signal.gr","GR_signal.gr"),
#   gene.file="ensembl_hg19_genes.txt",
#   export=FALSE
#)
# result.ER <- tip$ER_MACS_peaks # or tip[["ER_MACS_peaks"]]

TIP <- function(
    peak.files=NULL,
    signal.files=NULL,
    gene.file,
    signal.type="bed",
    promoter=c(-10000,1000),
    frag.size=200,
    method=c("cmg","normal","f","mc"),
    correct=TRUE,
    correct.method=c("qv","bh"),
    cmg.zcut=1.96,
    mc.iter=1000,
    export=TRUE,
    output=NULL,
    filter.output=TRUE,
    restrict.cores=0.6
)
{
    if (is.null(peak.files) && is.null(signal.files))
        stop("You must provide at least a set of peak files OR a set of signal",
        " files for each transcription factor!")
    if (!is.null(peak.files) & is.null(signal.files)) {
        warning("You did not provide a signal BED file for each peak file!",
            " Peak areas will be used with a possible loss of accuracy...",
            call.=FALSE,immediate. =TRUE)
        peak.present <- TRUE
        bed.present <- FALSE
        name.files <- peak.files
    }
    else if (is.null(peak.files) & !is.null(signal.files)) {
        warning("You did not provide a set of peak files! Only the signal ",
            "will be used with possibility of certain noisy results...",
            call.=FALSE,immediate.=TRUE)
        peak.present <- FALSE
        bed.present <- TRUE
        name.files <- signal.files
    }
    else {
        peak.present <- TRUE
        bed.present <- TRUE
        name.files <- peak.files
    }
    if (bed.present && !signal.type %in% c("bed","bam"))
        stop("signal.type must be one of \"bed\" or \"bam\"!")
    if (is.null(gene.file))
        stop("You must provide a set of genes (in BED format preferably) ",
            "corresponding to the organism of interest!")

    method <- tolower(method[1])
    if (!(method %in% c("cmg","normal","f","mc")))
        stop("method parameter must be one of \"cmg\", \"normal\", \"f\" or ",
            "\"mc\"")
    correct.method <- tolower(correct.method[1])
    if (!(correct.method %in% c("qv","bh")))
        stop("correct parameter must be one of \"qv\" or \"bh\"")
    if (restrict.cores < 0 | restrict.cores > 1)
        stop("restrict.cores parameter must be a value between 0 and 1!")

    if (!require(MASS) && method=="f")
        stop("R package MASS is required with when method is \"f\"!")
    if (correct.method=="qv" && !require(qvalue)) {
        warning("R package qvalue not present! Benjamini-Hochberg FDR will ",
            "not be reported instead which might be more strict...",
            call.=FALSE,immediate.=TRUE)
        qv.present <- FALSE
        correct.method <- "bh"
    }
    else {
        require(qvalue)
        qv.present <- TRUE
    }
    if (!require(GenomicRanges))
        stop("Bioconductor package IRanges is required!")
    if (!require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")
    if (signal.type=="bam" & !require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required!")

    # Check promoter - method combinations
    if (abs(promoter[1])!=abs(promoter[2])) {
        if (method=="cmg")
            warning("Method cmg might not work well with assymetric promoter ",
                "regions! Consider using normal, f or mc...",
                call.=FALSE,immediate.=TRUE)
    }
    else {
        if (method %in% c("normal","f","mc"))
            warning("Method cmg might work better and faster for symetric ",
                "promoter regions! Consider using it or mc but it will be ",
                "slower...",call.=FALSE,immediate.=TRUE)
    }
    # Check method-correction combination
    if (correct && method=="mc") {
        warning("Correction for mc method is highly experimental! Consider ",
            "using another p-value estimation method if you also need ",
            "corrected p-values...",call.=FALSE,immediate.=TRUE)
    }
    # Check requested output files
    if (export & is.null(output)) {
        output <- character(length(name.files))
        for (i in 1:length(name.files))
            output[i] <- file.path(dirname(name.files[i]),
                paste("TIP_out_",sub("^([^.]*).*","\\1",
                basename(name.files[i])),".txt",sep=""))
    }
    if (export & !is.null(output) & length(output)!=length(name.files)) {
        warning("The number of output files is different than the number of ",
            "the input peak files! They will be autogenerated...",
            call.=FALSE,immediate=TRUE)
        output <- character(length(name.files))
        for (i in 1:length(name.files))
            output[i] <- file.path(dirname(name.files[i]),
                paste("TIP_out_",sub("^([^.]*).*","\\1",
                basename(name.files[i])),".txt",sep=""))
    }
    if (!export) {
        out <- vector("list",length(name.files))
        names(out) <- sub("^([^.]*).*","\\1",basename(name.files))
    }
    #format(Sys.time(),format="%Y%m%d%H%M%S")

    # Read gene file
    message("\nReading genes file ",basename(gene.file),
        " to GenomicRanges object...")
    gene.data <- read.delim(gene.file)
    rownames(gene.data) <- as.character(gene.data[,4])
    genes.gr <- makeGRangesFromDataFrame(
        df=gene.data,
        keep.extra.columns=TRUE,
        seqnames.field="chromosome"
    )
    names(genes.gr) <- rownames(gene.data)
    # Get promoters and fix broken intervals
    ranges(genes.gr) <- promoters(ranges(genes.gr),upstream=abs(promoter[1]),
        downstream=abs(promoter[2]))
    ranges(genes.gr) <- restrict(ranges(genes.gr),start=1)

    # Do job with peak files
    for (i in 1:length(peak.files)) {
        peak.presence <- NULL
        if (peak.present) {
            message("Processing peaks file ",basename(peak.files[i]),"...")
            message("   Reading file to GenomicRanges object...")
            peak.data <- read.delim(peak.files[i])
            if (!bed.present && !("height" %in% names(peak.data)))
                stop("There must be a column named \"height\" (case sensitive)",
                    " in the peak file, containing peak summit heights/pileups",
                    " when signal files are not provided!")
            peaks.gr <- makeGRangesFromDataFrame(
                df=peak.data,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome"
            )
            names(peaks.gr) <- as.character(peak.data[,4])
            
            # Get the promoters with peaks to use for coverage calculation
            message("   Getting promoters containing actual peaks...")
            pinp <- findOverlaps(genes.gr,peaks.gr)
            proms.with.peaks <- GRanges(
                seqnames=seqnames(genes.gr)[queryHits(pinp)],
                ranges=IRanges(
                    start=start(genes.gr)[queryHits(pinp)],
                    end=end(genes.gr)[queryHits(pinp)]
                ),
                strand=strand(genes.gr)[queryHits(pinp)],
                peak=names(peaks.gr)[subjectHits(pinp)],
                gene=names(genes.gr)[queryHits(pinp)]
            )
            names(proms.with.peaks) <- as.character(proms.with.peaks$gene)
            peak.presence <- rep("-",length(genes.gr))
            names(peak.presence) <- names(genes.gr)
            peak.presence[unique(names(proms.with.peaks))] <- "+"
        }
        else
            proms.with.peaks <- genes.gr

        if (bed.present) {
            message("   Reading signal file ",basename(signal.files[i]),"...")
            if (signal.type=="bed")
                signal.gr <- import.bed(signal.files[i],trackLine=FALSE)
            else if (signal.type=="bam") {
				sf <- seqinfo(BamFile(signal.files[i]))
				signal.gr <- as(readGAlignments(file=signal.files[i]),"GRanges")
				seqinfo(signal.gr) <- sf
                signal.gr <- trim(signal.gr)
			}
            message("   Extending reads to ",frag.size,"bp...")
            signal.gr <- trim(resize(signal.gr,frag.size,fix="start"))
            
            # Calculate coverages...
            message("   Calculating signal coverages in promoters with peaks...")
            prom.signal <- calc.coverage(signal.gr,proms.with.peaks,
				rc=restrict.cores)
            names(prom.signal) <- names(proms.with.peaks)
        }
        else { # Only peaks are present, it has been checked in the beginning
            the.genes <- unique(elementMetadata(proms.with.peaks)[,"gene"])
            message("      Calculating peaks per gene...")
            peaks.per.gene <- cmclapply(the.genes,ppg,proms.with.peaks,
                rc=restrict.cores)
            names(peaks.per.gene) <- the.genes
            message("      Calculating promoters signals...")
            prom.signal <- cmclapply(peaks.per.gene,pp,rc=restrict.cores)
        }
        
        message("   Calculating signal matrix...")
        signal.matrix <- do.call("rbind",cmclapply(prom.signal,as.numeric,
            rc=restrict.cores))
        message("   Calculating binding profiles...")
        binding.profile <- apply(signal.matrix,2,sum,na.rm=TRUE)/
			sum(signal.matrix,na.rm=TRUE)
        message("   Calculating regulatory score...")
        reg.score <- apply(sweep(signal.matrix,2,binding.profile,"*"),1,sum,
			na.rm=TRUE)

        message("   Calculating peak association score... method: ",method)
        switch(method,
            cmg = {
                z.reg.score <- (reg.score - mean(reg.score))/sd(reg.score)
                p <- 2*pnorm(-abs(z.reg.score))
            },
            normal = {
                z.reg.score <- (log2(reg.score) - mean(log2(reg.score)))/
                    sd(log2(reg.score))
                p <- 2*pnorm(-abs(z.reg.score))
            },
            f = {
                fit <- fitdistr(reg.score,"f",start=list(df1=50,df2=5))
                p <- pf(reg.score,fit$estimate["df1"],fit$estimate["df2"])
            },
            mc = { # Yikes!
                # 1. Select a random set of promoters with size 
                #    length(proms.with.peaks) WITH replacement, as a peak might 
                #    be close to multiple genes
                # 2. For each promoter set, repeat the coverage calculation 
                #    procedure with scores etc.
                message("      Executing MC iterations...")
                mc.set <- vector("list",mc.iter)
                if (multic) {
                    mc.set <- cmclapply(mc.set,function(x,gb,pwp) {
                        tmp <- sample(gb,length(pwp),replace=TRUE)
                        return(tmp)
                    },genes.gr,proms.with.peaks,rc=restrict.cores)
                    mc.scores <- cmclapply(mc.set,mc.scoring,signal.gr,
                        rc=restrict.cores)
                    ptf <- do.call("rbind",cmclapply(mc.scores,function(x,r) { 
                        return(x > r) },rc=restrict.cores))
                }
                else {
                    mc.set <- lapply(mc.set,function(x.gb,pwp) {
                        tmp <- sample(gb,length(pwp),replace=TRUE)
                        return(tmp)
                    },genes.gr,proms.with.peaks)
                    mc.scores <- lapply(mc.set,mc.scoring,signal.gr,multic)
                    ptf <- do.call("rbind",lapply(mc.scores,function(x,r) { 
                        return(x > r) }))
                }
                p <- apply(ptf,1,which)/mc.iter
            }
        )

        switch(correct,
            qv = {
                fdr <- qvalue(p)$qvalues
            },
            bh = {
                fdr <- p.adjust(p,method="BH")
            }
        )

        fpr = length(which(reg.score <= -cmg.zcut))/length(which(reg.score >= 
            -cmg.zcut))
        if (method %in% c("f","mc") | abs(promoter[1])!=abs(promoter[2])) {
            message("\n========== False Positive Rate for ",
                basename(peak.files[i])," : ",fpr," ==========\n")
            warning("WARNING: FPR might not be accurate when method parameter ",
                "is \"f\" or \"mc\" or when promoter regions are assymetric! ",
                "Use with caution...",call.=FALSE,immediate=TRUE)
        }
        else
            message("\n========== False Positive Rate for ",
                basename(peak.files[i])," : ",fpr," ==========\n")

        if (peak.present & !filter.output) {
            reg.score.out <- p.out <- rep(NA,nrow(gene.data))
            names(reg.score.out) <- names(p.out) <- rownames(gene.data)
            reg.score.out[names(reg.score)] <- reg.score
            p.out[names(p)] <- p
            if (correct) {
                fdr.out <- rep(NA,nrow(gene.data))
                names(fdr.out) <- rownames(gene.data)
                fdr.out[names(fdr)] <- fdr
            }
            final <- cbind(gene.data,reg.score.out,p.out)
            if (correct)
                final <- cbind(final,fdr.out)
            if (!is.null(peak.presence))
                final <- cbind(final,peak.presence)
        }
        else {
            final <- cbind(gene.data[names(reg.score),],reg.score,p)
            if (correct)
                final <- cbind(final,fdr)
            if (!is.null(peak.presence)) {
                final <- cbind(final,peak.presence[names(reg.score)])
                names(final)[ncol(final)] <- "peak.presence"
            }
        }
        
        if (export) {
            message("Writing output to ",output[i],"...")
            write.table(final,file=output[i],sep="\t",row.names=FALSE,quote=FALSE)
        }
        else
            out[[sub("^([^.]*).*","\\1",basename(name.files[i]))]] <- final
    }
    
    cat("\n\nDone!\n\n")

    if (!export) return(out)
}

calc.coverage <- function(input,mask,strand=NULL,ignore.strand=TRUE,rc=NULL) {
    if (!is(input,"GRanges"))
        stop("The input argument must be a GenomicRanges object")
    if (!is(mask,"GRanges"))
        stop("The mask argument must be a GenomicRanges object")
    if (!is.null(strand) && !is.list(strand)) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    if (!is.list(input))
        input <- split.by.seqname(input)
    index <- 1:length(mask)
    message("    Calculating coverage...")
    coverage <- cmclapply(index,function(i,mask,input,ignore.strand) {
        x <- mask[i]
        y<-list(
            chromosome=as.character(seqnames(x)),
            start=start(x),
            end=end(x),
            strand=as.character(strand(x)),
            reads=NULL,
            coverage=NULL
        )
        if (!is.null(input[[y$chromosome]])) {
            y$reads <- input[[y$chromosome]][
                subjectHits(findOverlaps(x,input[[y$chromosome]],
                    ignore.strand=ignore.strand))]
        }
        else {
            message(y$chromosome,"not found!")
            y$reads <- NULL
        }
        if (length(y$reads)>0) {
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        }
        else
            return(NULL)
    },mask,input,ignore.strand,rc=rc)
    names(coverage) <- names(mask)
    gc(verbose=FALSE)
    message("    Done!")
    return(coverage) # Rle
}

split.by.seqname <- function(gr) {
    message("    Splitting input regions by seqname...")
    gr.list <- cmclapply(levels(seqnames(gr)),function(x,lib) {
        message("      ",x)
        tmp <- lib[seqnames(lib)==x]
        if (length(tmp)>0) return(tmp) else return(NULL)
    },gr)
    names(gr.list) <- levels(seqnames(gr))
    null <- which(sapply(gr.list,is.null))
    if (length(null)>0)
        gr.list <- gr.list[-null]
    return(gr.list)
}

ppg <- function(x,p) {
    tmp.p <- p[elementMetadata(p)[,"gene"]==x]
    return(list(
        peaks=tmp.p,
        gene=list(
            chr=as.character(seqnames(tmp.p))[1],
            start=start(tmp.p)[1],
            end=end(tmp.p)[1]
        )
    ))
}

ps <- function(x) {
    tmp <- coverage(x$reads)
    tmp <- tmp[[x$gene$chr]]
    len <- length(tmp)
    tmp <- c(tmp[x$gene$start:len],rep(0,x$gene$end-len))
    return(tmp)
}

pp <- function(x) {
    tmp <- numeric(x$gene$end - x$gene$start + 1)
    for (p in x$peaks) {
        tmp[start(p):end(p)] <- p[,"height"]
    }
    return(tmp)
}

mc.scoring <- function(x,s,rc) {
    pr.si <- calc.coverage(s,x)
    names(pr.si) <- names(x)
    sm <- do.call("rbind",cmclapply(pr.si,as.numeric,rc))
    bp <- apply(sm,2,sum)/sum(sm)
    rs <- apply(sweep(sm,2,bp,"*"),1,sum)
    return(rs)
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!require(parallel)) || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.na(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            options(cores=ncores)
        }
    }
    if (m)
        return(mclapply(...,mc.cores=getOption("cores"),mc.set.seed=FALSE))
    else
        return(lapply(...))
}

check.parallel <- function(rc) {
    if (!require(parallel))
            multi <- FALSE
    else {
        multi <- TRUE
        ncores <- parallel:::detectCores()
        if (!is.na(rc))
            ncores <- ceiling(rc*ncores)
        options(cores=ncores)
    }
    return(multi)
}

# Test
#peak.files <- "/media/HD4/Fleming/data/pantelis_nl/chipseq/analysis/TCF4_peaks_filtered.txt"
#signal.files <- "/media/HD4/Fleming/data/pantelis_nl/chipseq/data/TCF4.sample.bed"
#signal.files <- "/media/HD4/Fleming/data/pantelis_nl/chipseq/data/TCF4.bed"
#gene.file <- "/media/HD4/Fleming/data/pantelis_rnaseq/analysis/annotation/hg18_ensembl_genes.txt"
#tip <- TIP(peak.files,signal.files,gene.file)
#TIP <- peak.files,signal.files,gene.file,promoter=c(-10000,1000),frag.size=200,method=c("cmg","normal","f","mc"),correct=TRUE,
#                   correct.method=c("qv","bh"),cmg.zcut=1.96,mc.iter=1000,export=TRUE,output=NA)

### Trash code
## Get the peaks in promoters and expand them a little
#pinp.hits <- GRanges(
#   seqnames=seqnames(peaks.gr)[subjectHits(pinp)],
#   ranges=IRanges(
#       start=start(peaks.gr)[subjectHits(pinp)]-100,
#       end=end(peaks.gr)[subjectHits(pinp)]+100
#   ),
#   strand=strand(peaks.gr)[subjectHits(pinp)],
#   name=elementMetadata(peaks.gr)[subjectHits(pinp),"name"],
#   score=elementMetadata(peaks.gr)[subjectHits(pinp),"score"],
#   gene=elementMetadata(genes.gr)[queryHits(pinp),]
#)
## In case cmg correction could be done for each gene
#if (correct.method=="cmg" && !method %in% c("cmg","normal"))
#{
#   if (qv.present)
#   {
#       warning("Correction method cmg requires cmg or normal method for p-value estimation! Switching to q-value...",
#           call.=FALSE)
#       correct.method <- "qv"
#   }
#   else
#   {
#       warning("Correction method cmg requires cmg or normal method for p-value estimation! Switching to Benjamini-Hochberg...",
#           call.=FALSE)
#       correct.method <- "bh"
#   }
#}
