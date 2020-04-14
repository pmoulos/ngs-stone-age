# 0. Load genes and exons and get introns
# 1. Tile the genome (using GenomeInfoDb as initial genomic ranges)
# 2. Separate reads to strands
# 3. Fix to start position of each read
# 4. For each strand
#   a) Coverage with the 1-length reads
#   b) Views over the tiles
#   c) Find spikes and their positions over the tiles
#   d) If one window has >1 spikes, sum the spikes and maintain as position
#      the position belonging to the highest spike
#   e) Create GRanges with chr, spike start-end, strand, height
#   f) Retain the spikes with height > threshold. These are the I.I.
# 5. Calculate Intronic II in a strand-specific manner
# 6. Calculate Exonic II in a strand-specific manner
# 7. Calculate antisense Intronic II 
# 8. Calculate sense Intronic II
# 9. Recreate the GRangesList objects from the extended exonic and intronic
#    structures created in 5-8.
# 10. Counts are obtained by iteration over II from each dataset for all genes.
#     Insertion counts for introns are obtained in an orientation specific 
#     manner for assessing the trapping function of the splice acceptor of the
#     gene-trap
# 11. Disrupting insertions (D.I.) are calculated as the sum of exonic I.I., 
#     and intronic I.I. with the splice acceptor aligned with transcription
# 12. A bias is calculated as the ratio of intronic insertions in sense over
#     anti-sense orientation for each gene

# frac: fraction of non-zero (II,DI,Bias) triplets to be used as neighbor data
#       for LOF analysis
# doTest: will test if there is difference in reads
# artPct: looks at introns for sequencing artifacts, excludes falling in >0.99
#         of the read counts distribution
hasar <- function(selected,control,name=NULL,gene,exon,winSize=5,threshold=3,
    frac=0.1,normalize=c("none","downsample","sampleto"),doTest=FALSE,
    artPct=0.95,sampleto=1e+6,outList=FALSE,outPath=NULL,output=NULL,
    outFormat=c("txt","xls"),seed=42,rc=NULL) {
    if (!require(parallel))
        stop("R package parallel is required")
    if (!require(Rlof))
        stop("R package Rlof is required")
    if (!require(GenomicRanges))
        stop("Bioconductor package GenomicRanges is required")
    if (!require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required")
    if (!require(GenomeInfoDb))
        stop("Bioconductor package GenomicAlignments is required")
        
    if (!is.character(selected))
        stop("selected must be a vector of BAM files")
    if (!is.character(control))
        stop("control must be a vector of BAM files")
    if (!is.null(name) && !is.character(name))
        stop("name must be a valid analysis name or NULL")
    if (!is.null(outPath) && !is.character(outPath))
        stop("outPath must be a valid analysis name or NULL")
    if (!is.numeric(winSize) || winSize <= 0)
        stop("winSize must be an integer greater than 0")
    if (!is.numeric(threshold) || threshold <= 0)
        stop("threshold must be an integer greater than 0")
    if (!is.numeric(frac) || frac <= 0 || frac > 1)
        stop("frac must be a number greater than 0 and smaller or equal to 1")
    if (!is.null(output) && !is.logical(output) && !is.character(output))
        stop("output must be NULL, TRUE/FALSE or a string naming an out file")
    if (!is.logical(doTest))
        stop("doTest must be TRUE/FALSE")
    normalize <- normalize[1]
    if (!is.character(normalize) || !(normalize %in% c("none",
        "downsample","sampleto")))
        stop("normalize must be one of \"none\", \"downsample\", \"sampleto\"")
    outFormat <- outFormat[1]
    if (!is.character(outFormat) || !(outFormat %in% c("txt","xls")))
        stop("outFormat must be one of \"txt\", \"xls\"")
    if (outFormat == "xls" && !require(openxlsx))
        stop("R package xlsx is required for Excel output")
    if (normalize == "sampleto" && !is.numeric(sampleto))
        stop("sampleto must be a number when normalize is sampleto")
    if (!is.null(rc) && (rc<0 || rc>1)) {
        warning("rc cannot must be between 0 and 1! Ignoring...",
            immediate.=TRUE)
        rc <- NULL
    }
    
    # Check files exist
    if (any(!file.exists(selected)))
        stop("Please check that all selected BAM files exist!")
    if (any(!file.exists(control)))
        stop("Please check that all control BAM files exist!")
    
    # Check outpath
    if (is.null(outPath))
        outPath <- getwd()
    else if (!dir.exists(outPath)) {
        warning("Output path ",outPath," dies not exist! Using current dir...")
        outPath <- getwd()
    }
    
    # Analysis name to create output dir (in the current)
    if (is.null(name))
         name <- paste("hasaR_output_",format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),
            sep="")
    analysisDir <- file.path(outPath,name)
    if (!dir.exists(analysisDir))
        dir.create(analysisDir,recursive=TRUE,mode="0755")
    
    # If more than one file(s) given for any condition, they must be merged 
    # before importing to R for memory issues. Do for every condition separately
    if (length(selected) > 1) {
        dest <- file.path(analysisDir,"selected")
        mergeBam(files=selected,destination=dest)
        selected <- dest
    }
    if (length(control) > 1) {
        dest <- file.path(analysisDir,"control")
        mergeBam(files=control,destination=dest)
        control <- dest
    }
    
    # TODO: Normalization of merged... Could be below...
    
    # Load reference genomic intervals
    message("Loading required reference genomic intervals")
    load(gene)
    gene <- gene
    load(exon)
    exon <- sexon
    
    libSizes <- normFactor <- NULL
    if (normalize != "none") {
        message("Getting library sizes for normalization")
        message("  selected")
        libSizeS <- getLibsize(selected)
        message("  control")
        libSizeC <- getLibsize(control)
        libSize <- c(libSizeS,libSizeC)
        
        if (normalize == "downsample")
            normFactor <- min(libSize)
        else if (normalize == "sampleto") {
            if (any(libSize<sampleto)) {
                m <- min(libSize)
                warning("A library size was found less than sampleto (",m," <",
                    sampleto,"). Raising to ",m,"...",immediate.=TRUE)
                sampleto <- m
            }
            normFactor <- sampleto
        }
    }
    
    # Read BAM files
    message("Reading selected library ",selected)
    reads1 <- trim(as(readGAlignments(file=selected),"GRanges"))
    message("Reading control library ",control)
    reads2 <- trim(as(readGAlignments(file=control),"GRanges"))
    
    # Normalize
    if (!is.null(normFactor)) {
        message("    normalizing")
        set.seed(seed)
        downsampleIndex1 <- sort(sample(length(reads1),normFactor))
        downsampleIndex2 <- sort(sample(length(reads2),normFactor))
        reads1 <- reads1[downsampleIndex1]
        reads2 <- reads2[downsampleIndex1]
    }
    
    # Triplets
    message("Calculating (II,DI,Bias) spaces for selected and control")
    tripletList <- makeTriplets(reads1,reads2,gene,exon,winSize,threshold,
        artPct,normFactor,seed,rc=rc)
    colnames(tripletList$selected) <- c("II_Sum_Selected","DI_Sum_Selected",
        "Bias_Sum_Selected","Reads_Sum_Selected")
    colnames(tripletList$control) <- c("II_Sum_Control","DI_Sum_Control",
        "Bias_Sum_Control","Reads_Sum_Control")
    
    message("Calculating triplet for LOF analysis")
    theTriplet <- suppressWarnings(
        tripletList$selected[,1:3]/tripletList$control[,1:3])
    theTriplet[is.na(theTriplet)] <- 0
    theTriplet[is.infinite(theTriplet)] <- 0
    colnames(theTriplet) <- c("II","DI","Bias")

    # LOF
    ii <- which(apply(theTriplet,1,function(x) return(all(x==0))))
    triF <- theTriplet
    if (length(ii) > 0)
        triF <- triF[-ii,]
    k = ceiling(frac*nrow(triF))
    
    message("Calculating LOF")
    lofFactor <- lof(triF,k)
    names(lofFactor) <- rownames(triF)
    lofAll <- rep(NA,nrow(theTriplet))
    names(lofAll) <- rownames(theTriplet)
    lofAll[rownames(triF)] <- lofFactor
    
    message("Assembling output")
    outputList <- list(
        lof=lofAll,
        triplet=theTriplet,
        sumSelected=tripletList$selected,
        sumControl=tripletList$control#,
        #selected=selectedTrips,
        #control=controlTrips
    )
    
    if (is.null(output) || (is.logical(output) && !output)) {
        if (outList)
            return(outputList)
        invisible(return(FALSE))
    }   
    else if (is.logical(output) && output) { 
        # Output file requested but not provided
        if (outFormat == "txt") {
            output1 <- file.path(analysisDir,paste("hasaR_summary",
                format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),".",outFormat,sep=""))
            #output2 <- file.path(analysisDir,paste("hasaR_detail",
            #    format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),".",outFormat,sep=""))
        }
        else if (outFormat == "xls")
            output <- file.path(analysisDir,paste("hasaR_",
                format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),".xlsx",sep=""))
    }
    else { # Output base provided, must add extensions
        if (outFormat == "txt") {
            output1 <- file.path(analysisDir,paste(output,"_summary","."
                ,outFormat,sep=""))
            #output2 <- file.path(analysisDir,paste(output,"_detail",".",
            #    outFormat,sep=""))
        }
        else if (outFormat == "xls")
            output <- file.path(analysisDir,paste(output,".xlsx",sep=""))
    }
    
    analData <- data.frame(
        gene_id=as.character(gene$gene_id),
        gene_name=as.character(gene$gene_name),
        II_Fold=outputList$triplet[,1],
        DI_Fold=outputList$triplet[,2],
        Bias_Fold=outputList$triplet[,3],
        LOF=outputList$lof
    )
    #export <- cbind(analData,selectedTriplet,controlTriplet,
    #    do.call("cbind",outputList$selected),do.call("cbind",
    #    outputList$control))
    export1 <- cbind(analData,tripletList$selected,tripletList$control)
    #export2 <- cbind(
    #    analData[,1:2],
    #    do.call("cbind",outputList$selected),
    #    do.call("cbind",outputList$control)
    #)
    
    # Order by LOF
    theOrder <- order(export1$LOF,na.last=TRUE,decreasing=TRUE)
    export1 <- export1[theOrder,]
    #export2 <- export2[theOrder,]
    
    if (outFormat == "txt") {
        message("Writing summary output to ",output1)
        write.table(export1,file=output1,sep="\t",quote=FALSE,row.names=FALSE)
        #message("Writing detailed output to ",output2)
        #write.table(export2,file=output2,sep="\t",quote=FALSE,row.names=FALSE)
    }
    else if (outFormat == "xls") {
        #xlsList <- list(summary=export1,detail=export2)
        xlsList <- list(summary=export1)
        message("Writing summary output to ",output)
        write.xlsx(xlsList,file=output,keepNA=TRUE)
    }
    
    if (outList)
        return(outputList)
}

# file1, file2: since original hasappy algorithm ADDS reads among samples, then
# prior to feeding to makeTriple, we merge reads from all samples for the 
# comparison to be performed. So file1, file2 eventually to become reads1, 
# reads2 after prior reading and merging

makeTriplets <- function(reads1,reads2,gene,exon,winSize,threshold,artPct,
    normFactor=NULL,seed=42,rc=NULL) {
    
    # Separate reads to strands
    message("  separate reads per strand")
    message("    selected")
    rp1 <- reads1[strand(reads1)=="+"]
    rm1 <- reads1[strand(reads1)=="-"]
    message("    control")
    rp2 <- reads2[strand(reads2)=="+"]
    rm2 <- reads2[strand(reads2)=="-"]

    # Fix to start position of each read
    message("  getting first base position")
    message("    + strand")
    message("      selected")
    rps1 <- resize(rp1,width=1,fix="start")
    message("      control")
    rps2 <- resize(rp2,width=1,fix="start")
    message("    - strand")
    message("      selected")
    rms1 <- resize(rm1,width=1,fix="start")
    message("      control")
    rms2 <- resize(rm2,width=1,fix="start")

    # For each strand
    # Coverage with the 1-length reads
    message("  calculating coverage of start positions")
    message("    + strand")
    message("      selected")
    covsp1 <- coverageFromRanges(rps1,gene)
    message("      control")
    covsp2 <- coverageFromRanges(rps2,gene)
    message("    - strand")
    message("      selected")
    covsm1 <- coverageFromRanges(rms1,gene)
    message("      control")
    covsm2 <- coverageFromRanges(rms2,gene)
    
    message("  filtering for 0 or 1 coverage")
    message("    + strand")
    message("      selected")
    fp1 <- which(sapply(covsp1,function(x) {
        v <- runValue(x)
        if (sum(v) <= 1)
            return(TRUE)
        return(FALSE)
    }))
    message("      control")
    fp2 <- which(sapply(covsp2,function(x) {
        v <- runValue(x)
        if (sum(v) <= 1)
            return(TRUE)
        return(FALSE)
    }))
    message("    + strand")
    message("      selected")
    fm1 <- which(sapply(covsm1,function(x) {
        v <- runValue(x)
        if (sum(v) <= 1)
            return(TRUE)
        return(FALSE)
    }))
    message("      control")
    fm2 <- which(sapply(covsm2,function(x) {
        v <- runValue(x)
        if (sum(v) <= 1)
            return(TRUE)
        return(FALSE)
    }))
    fu <- Reduce("intersect",list(fp1,fp2,fm1,fm2))
    covsp1 <- covsp1[-fu]
    covsp2 <- covsp2[-fu]
    covsm1 <- covsm1[-fu]
    covsm2 <- covsm2[-fu]
    
    # Tile only covered genes, speeds up process
    message("Getting and tiling genome (window size: ",winSize,")")
    cnams <- names(covsp1)
    subgene <- gene[cnams]
    subexon <- exon[cnams]
    
    # Final tiles
    tiles <- tile(subgene,width=winSize)
    names(tiles) <- names(subgene)
    tileRanges <- unlist(tiles)
    
    message("Formatting genomic regions for rest of the process")
    # Get introns
    geneSplit <- split(subgene,names(subgene))
    geneSplit <- geneSplit[names(subgene)]
    # Voila
    intronSplit <- setdiff(geneSplit,subexon)
    intron <- unlist(intronSplit)
    intron$gene_id <- names(intron)
    names(intron) <- paste(names(intron),".i",1:length(intron),sep="")
    exon <- unlist(subexon)
    
    # Calculate overlaps. hopefully much faster than coverage/Views
    message("  calculating overlaps of start positions")
    message("    + strand")
    message("      selected")
    op1 <- countOverlaps(tileRanges,rps1)
    message("      control")
    op2 <- countOverlaps(tileRanges,rps2)
    message("    - strand")
    message("      selected")
    om1 <- countOverlaps(tileRanges,rms1)
    message("      control")
    om2 <- countOverlaps(tileRanges,rms2)
    
    # If one window has >1 spikes, sum the spikes and maintain as position the
    # position belonging to the highest spike
    message("  filtering overlaps of start positions")
    message("    + strand")
    message("      selected")
    nozerop1 <- which(op1>0)
    spikesp1 <- op1[nozerop1]
    names(spikesp1) <- paste(seqnames(tileRanges[nozerop1]),":",
        start(tileRanges[nozerop1]),"-",end(tileRanges[nozerop1]),sep="")
    message("      control")
    nozerop2 <- which(op2>0)
    spikesp2 <- op1[nozerop2]
    names(spikesp2) <- paste(seqnames(tileRanges[nozerop2]),":",
        start(tileRanges[nozerop2]),"-",end(tileRanges[nozerop2]),sep="")
        
    message("    - strand")
    message("      selected")
    nozerom1 <- which(om1>0)
    spikesm1 <- om1[nozerom1]
    names(spikesm1) <- paste(seqnames(tileRanges[nozerom1]),":",
        start(tileRanges[nozerom1]),"-",end(tileRanges[nozerom1]),sep="")
    message("      control")
    nozerom2 <- which(om2>0)
    spikesm2 <- om2[nozerom2]
    names(spikesm2) <- paste(seqnames(tileRanges[nozerom2]),":",
        start(tileRanges[nozerom2]),"-",end(tileRanges[nozerom2]),sep="")

    # Create GRanges with chr, spike start-end, strand, height
    message("  formatting overlaps to find candidate insertions")
    message("    + strand")
    message("      selected")
    tmpIIp1 <- GRanges(names(spikesp1))
    strand(tmpIIp1) <- Rle("+",length(tmpIIp1))
    tmpIIp1$reads <- spikesp1
    message("      control")
    tmpIIp2 <- GRanges(names(spikesp2))
    strand(tmpIIp2) <- Rle("+",length(tmpIIp2))
    tmpIIp2$reads <- spikesp2
    message("    - strand")
    message("      selected")
    tmpIIm1 <- GRanges(names(spikesm1))
    strand(tmpIIm1) <- Rle("-",length(tmpIIm1))
    tmpIIm1$reads <- spikesm1
    message("      control")
    tmpIIm2 <- GRanges(names(spikesm2))
    strand(tmpIIm2) <- Rle("-",length(tmpIIm2))
    tmpIIm2$reads <- spikesm2
    
    message("  calculating coverage over genomic windows to refine positions")
    message("    + strand")
    message("      selected")
    Cp1 <- coverageFromRanges(rps1,tmpIIp1)
    names(Cp1) <- names(spikesp1)
    message("      control")
    Cp2 <- coverageFromRanges(rps2,tmpIIp2)
    names(Cp2) <- names(spikesp2)
    message("    - strand")
    message("      selected")
    Cm1 <- coverageFromRanges(rms1,tmpIIm1)
    names(Cm1) <- names(spikesm1)
    message("      control")
    Cm2 <- coverageFromRanges(rms2,tmpIIm2)
    names(Cm2) <- names(spikesm2)
       
    message("  refining candidate insertion positions over genomic windows")
    message("    + strand")
    message("      selected")
    preIIp1 <- cmclapply(names(Cp1),procWin,"+",Cp1,rc=rc)
    preIIp1 <- do.call("rbind",preIIp1)
    preIIp1 <- GRanges(preIIp1)
    preIIp1$reads <- tmpIIp1$reads
    message("      control")
    preIIp2 <- cmclapply(names(Cp2),procWin,"+",Cp2,rc=rc)
    preIIp2 <- do.call("rbind",preIIp2)
    preIIp2 <- GRanges(preIIp2)
    preIIp2$reads <- tmpIIp2$reads
    message("    - strand")
    message("      selected")
    preIIm1 <- cmclapply(names(Cm1),procWin,"-",Cm1,rc=rc)
    preIIm1 <- do.call("rbind",preIIm1)
    preIIm1 <- GRanges(preIIm1)
    preIIm1$reads <- tmpIIm1$reads
    message("      control")
    preIIm2 <- cmclapply(names(Cm2),procWin,"-",Cm2,rc=rc)
    preIIm2 <- do.call("rbind",preIIm2)
    preIIm2 <- GRanges(preIIm2)
    preIIm2$reads <- tmpIIm2$reads
    message("  merging")
    preII1 <- c(preIIp1,preIIm1)
    preII2 <- c(preIIp2,preIIm2)
    message("  fltering (threshold: ",threshold,")")
    canII1 <- preII1[preII1$score>threshold] # II candidates
    canII2 <- preII2[preII2$score>threshold]
    
    # Init further structures
    message("  processing candidate insertions")
    intronCounts1 <- intronCounts2 <- intron
    exonCounts1 <- exonCounts2 <- exon
    # Init intronic structure
    message("    initiating...")
    intronCounts1$score <- numeric(length(intronCounts1))
    intronCounts1$reads <- numeric(length(intronCounts1))
    intronCounts1$II <- character(length(intronCounts1))
    intronCounts1$antisense_score <- numeric(length(intronCounts1))
    intronCounts1$antisense_II <- character(length(intronCounts1))
    intronCounts1$sense_score <- numeric(length(intronCounts1))
    intronCounts1$sense_II <- character(length(intronCounts1))
    intronCounts2$score <- numeric(length(intronCounts2))
    intronCounts2$reads <- numeric(length(intronCounts2))
    intronCounts2$II <- character(length(intronCounts2))
    intronCounts2$antisense_score <- numeric(length(intronCounts2))
    intronCounts2$antisense_II <- character(length(intronCounts2))
    intronCounts2$sense_score <- numeric(length(intronCounts2))
    intronCounts2$sense_II <- character(length(intronCounts2))
    # Init exonic structure
    exonCounts1$score <- numeric(length(exonCounts1))
    exonCounts1$reads <- numeric(length(exonCounts1))
    exonCounts1$II <- character(length(exonCounts1))
    exonCounts2$score <- numeric(length(exonCounts2))
    exonCounts2$reads <- numeric(length(exonCounts2))
    exonCounts2$II <- character(length(exonCounts2))

    # Calculate Intronic II in a strand-specific manner
    message("  processing candidate independent insertions")
    message("    intron II")
    message("      selected")
    intronic1 <- findOverlaps(intron,canII1,ignore.strand=FALSE)
    intronCounts1$score[queryHits(intronic1)] <- 
        canII1$score[subjectHits(intronic1)]
    intronCounts1$reads[queryHits(intronic1)] <- 
        canII1$reads[subjectHits(intronic1)]
    intronCounts1$II[queryHits(intronic1)] <- 
        paste(seqnames(canII1)[subjectHits(intronic1)],":",
            start(canII1)[subjectHits(intronic1)],"-",
            end(canII1)[subjectHits(intronic1)],sep="")
    message("      control")
    intronic2 <- findOverlaps(intron,canII2,ignore.strand=FALSE)
    intronCounts2$score[queryHits(intronic2)] <- 
        canII2$score[subjectHits(intronic2)]
    intronCounts2$reads[queryHits(intronic2)] <- 
        canII2$reads[subjectHits(intronic2)]
    intronCounts2$II[queryHits(intronic2)] <- 
        paste(seqnames(canII2)[subjectHits(intronic2)],":",
            start(canII2)[subjectHits(intronic2)],"-",
            end(canII2)[subjectHits(intronic2)],sep="")
            
    # Optionally, filter intronicCounts structure for potential artifacts
    if (!is.na(artPct)) {
        seqArtI1 <- quantile(intronCounts1$reads[intronCounts1$reads>0],artPct)
        seqArtI2 <- quantile(intronCounts2$reads[intronCounts2$reads>0],artPct)
        f1 <- which(intronCounts1$reads<=seqArtI1)
        f2 <- which(intronCounts2$reads<=seqArtI2)
        ff <- union(f1,f2)
        intronCounts1 <- intronCounts1[ff]
        intronCounts2 <- intronCounts2[ff]
    }

    # Calculate Exonic II in a strandless manner
    message("    exon II")
    message("      selected")
    exonic1 <- findOverlaps(exon,canII1,ignore.strand=TRUE)
    exonCounts1$score[queryHits(exonic1)] <- canII1$score[subjectHits(exonic1)]
    exonCounts1$reads[queryHits(exonic1)] <- canII1$reads[subjectHits(exonic1)]
    exonCounts1$II[queryHits(exonic1)] <- 
        paste(seqnames(canII1)[subjectHits(exonic1)],":",
            start(canII1)[subjectHits(exonic1)],"-",
            end(canII1)[subjectHits(exonic1)],sep="")
    message("      control")
    exonic2 <- findOverlaps(exon,canII2,ignore.strand=TRUE)
    exonCounts2$score[queryHits(exonic2)] <- canII2$score[subjectHits(exonic2)]
    exonCounts2$reads[queryHits(exonic2)] <- canII2$reads[subjectHits(exonic2)]
    exonCounts2$II[queryHits(exonic2)] <- 
        paste(seqnames(canII2)[subjectHits(exonic2)],":",
            start(canII2)[subjectHits(exonic2)],"-",
            end(canII2)[subjectHits(exonic2)],sep="")
    
    if (!is.na(artPct)) {
        seqArtE1 <- quantile(intronCounts1$reads[intronCounts1$reads>0],artPct)
        seqArtE2 <- quantile(intronCounts2$reads[intronCounts2$reads>0],artPct)
        f1 <- which(intronCounts1$reads<=seqArtE1)
        f2 <- which(intronCounts2$reads<=seqArtE2)
        ff <- union(f1,f2)
        exonCounts1 <- intronCounts1[ff]
        exonCounts2 <- intronCounts2[ff]
    }
    
    # FIXME: Introns!
    # Calculate antisense Intronic II
#~     message("    antisense intron II")
#~     sin <- which(strand(intron)=="+")
#~     scn <- which(strand(canII)=="-")
#~     if (length(sin) > 0 && length(scn) > 0) {
#~         intronicAsense <- findOverlaps(intron[sin],canII[scn],
#~             ignore.strand=TRUE)
#~         anti_names <- names(intron[sin])[(queryHits(intronicAsense))]
#~         intronCounts[anti_names]$antisense_score <- 
#~             canII[scn]$score[subjectHits(intronicAsense)]
#~         intronCounts[anti_names]$antisense_II <- paste(
#~             seqnames(canII[scn])[subjectHits(intronicAsense)],":",
#~             start(canII[scn])[subjectHits(intronicAsense)],"-",
#~             end(canII[scn])[subjectHits(intronicAsense)],sep="")
#~     }
    message("    antisense intron II")
    message("      selected")
    sin1 <- which(strand(intronCounts1)=="+")
    scn1 <- which(strand(canII1)=="-")
    if (length(sin1) > 0 && length(scn1) > 0) {
        intronicAsense1 <- findOverlaps(intronCounts1[sin1],canII1[scn1],
            ignore.strand=TRUE)
        anti_names1 <- names(intronCounts1[sin1])[(queryHits(intronicAsense1))]
        intronCounts1[anti_names1]$antisense_score <- 
            canII1[scn1]$score[subjectHits(intronicAsense1)]
        intronCounts1[anti_names1]$antisense_II <- paste(
            seqnames(canII1[scn1])[subjectHits(intronicAsense1)],":",
            start(canII1[scn1])[subjectHits(intronicAsense1)],"-",
            end(canII1[scn1])[subjectHits(intronicAsense1)],sep="")
    }
    message("      control")
    sin2 <- which(strand(intronCounts2)=="+")
    scn2 <- which(strand(canII2)=="-")
    if (length(sin2) > 0 && length(scn2) > 0) {
        intronicAsense2 <- findOverlaps(intronCounts2[sin2],canII2[scn2],
            ignore.strand=TRUE)
        anti_names2 <- names(intronCounts2[sin2])[(queryHits(intronicAsense2))]
        intronCounts2[anti_names2]$antisense_score <- 
            canII2[scn2]$score[subjectHits(intronicAsense2)]
        intronCounts2[anti_names2]$antisense_II <- paste(
            seqnames(canII2[scn2])[subjectHits(intronicAsense2)],":",
            start(canII2[scn2])[subjectHits(intronicAsense2)],"-",
            end(canII2[scn2])[subjectHits(intronicAsense2)],sep="")
    }
    
    # Calculate sense Intronic II
#~     message("    sense intron II")
#~     sin <- which(strand(intron)=="-")
#~     scn <- which(strand(canII)=="+")
#~     if (length(sin) > 0 && length(scn) > 0) {
#~         intronicSense <- findOverlaps(intron[sin],canII[scn],
#~             ignore.strand=TRUE)
#~         sense_names <- names(intron[sin])[(queryHits(intronicSense))]
#~         intronCounts[sense_names]$sense_score <- 
#~             canII[scn]$score[subjectHits(intronicSense)]
#~         intronCounts[sense_names]$antisense_II <- paste(
#~             seqnames(canII[scn])[subjectHits(intronicSense)],":",
#~             start(canII[scn])[subjectHits(intronicSense)],"-",
#~             end(canII[scn])[subjectHits(intronicSense)],sep="")
#~     }
    message("    sense intron II")
    message("      selected")
    sin1 <- which(strand(intronCounts1)=="-")
    scn1 <- which(strand(canII1)=="+")
    if (length(sin1) > 0 && length(scn1) > 0) {
        intronicSense1 <- findOverlaps(intronCounts1[sin1],canII1[scn1],
            ignore.strand=TRUE)
        sense_names1 <- names(intronCounts1[sin1])[(queryHits(intronicSense1))]
        intronCounts1[sense_names1]$sense_score <- 
            canII1[scn1]$score[subjectHits(intronicSense1)]
        intronCounts1[sense_names1]$antisense_II <- paste(
            seqnames(canII1[scn1])[subjectHits(intronicSense1)],":",
            start(canII1[scn1])[subjectHits(intronicSense1)],"-",
            end(canII1[scn1])[subjectHits(intronicSense1)],sep="")
    }
    message("      control")
    sin2 <- which(strand(intronCounts2)=="-")
    scn2 <- which(strand(canII2)=="+")
    if (length(sin2) > 0 && length(scn2) > 0) {
        intronicSense2 <- findOverlaps(intronCounts2[sin2],canII2[scn2],
            ignore.strand=TRUE)
        sense_names2 <- names(intronCounts2[sin2])[(queryHits(intronicSense2))]
        intronCounts2[sense_names2]$sense_score <- 
            canII2[scn2]$score[subjectHits(intronicSense2)]
        intronCounts2[sense_names2]$antisense_II <- paste(
            seqnames(canII2[scn2])[subjectHits(intronicSense2)],":",
            start(canII2[scn2])[subjectHits(intronicSense2)],"-",
            end(canII2[scn2])[subjectHits(intronicSense2)],sep="")
    }
    
    # Recreate the GRangesList objects from the extended exonic and intronic
    # structures created in 5-8.
    ## Problem as some transcripts do not have exons
    #intronStru <- split(intronCounts,names(intronCounts))
    ## So we do the following workaround
    message("  reassembling gene structures")
    intronStru1 <- intronStru2 <- intronSplit
    exonStru1 <- exonStru2 <- subexon
    message("    selected")
    # and then empty it
    intronStru1[names(intronStru1)] <- GRangesList(GRanges())
    tmp1 <- split(intronCounts1,intronCounts1$gene_id)
    intronStru1[names(tmp1)] <- tmp1
    intronStru1 <- intronStru1[names(subgene)]
    # Recreate gene names from names of exonCounts (with the same way)
    #exonCounts1 <- unname(exonCounts1)
    #names(exonCounts1) <- exonCounts1$exon_id
    #exonStru1 <- split(exonCounts1,exonCounts1$gene_id)
    #exonStru1 <- exonStru1[names(subgene)]
    exonStru1[names(exonStru1)] <- GRangesList(GRanges())
    tmp1e <- split(exonCounts1,exonCounts1$gene_id)
    exonStru1[names(tmp1e)] <- tmp1e
    exonStru1 <- exonStru1[names(subgene)]
    message("    control")
    # and then empty it
    intronStru2[names(intronStru2)] <- GRangesList(GRanges())
    tmp2 <- split(intronCounts2,intronCounts2$gene_id)
    intronStru2[names(tmp2)] <- tmp2
    intronStru2 <- intronStru2[names(subgene)]
    # Recreate gene names from names of exonCounts
    #exonCounts2 <- unname(exonCounts2)
    #names(exonCounts2) <- exonCounts2$exon_id
    #exonStru2 <- split(exonCounts2,exonCounts2$gene_id)
    #exonStru2 <- exonStru2[names(subgene)]
    exonStru2[names(exonStru2)] <- GRangesList(GRanges())
    tmp2e <- split(exonCounts2,exonCounts2$gene_id)
    exonStru2[names(tmp2e)] <- tmp2e
    exonStru2 <- exonStru2[names(subgene)]

    #  Counts are obtained by iteration over II from each dataset for all genes.
    #  Insertion counts for introns are obtained in an orientation specific 
    #  manner for assessing the trapping function of the splice acceptor of the
    #  gene-trap
    message("  calculating II and DI")
    message("    selected for introns")
    tmp_intronII1 <- as.data.frame(mcols(intronStru1,level="within"))
    tmp_intronII1$strand <- 
        suppressWarnings(as.character(unlist(strand(intronStru1))))
    pre_intronII1 <- vector("list",length(subgene))
    tmp_intronII1 <- split(tmp_intronII1,tmp_intronII1$group_name)
    names(pre_intronII1) <- names(subgene)
    pre_intronII1[names(tmp_intronII1)] <- tmp_intronII1
    intronII1 <- unlist(cmclapply(names(pre_intronII1),function(n,D) {
        x <- D[[n]]
        return(sum(x$score))
    },pre_intronII1,rc=0.5))
    anti_intronII1 <- unlist(cmclapply(names(pre_intronII1),function(n,D) {
        x <- D[[n]]
        if (!is.null(x)) {
            if (x$strand[1] == "-")
                return(sum(x$sense_score))
            return(sum(x$antisense_score))
        }
        return(0)
    },pre_intronII1,rc=0.5))
    names(intronII1) <- names(anti_intronII1) <- names(pre_intronII1)
    message("    control for introns")
    tmp_intronII2 <- as.data.frame(mcols(intronStru2,level="within"))
    tmp_intronII2$strand <- 
        suppressWarnings(as.character(unlist(strand(intronStru2))))
    pre_intronII2 <- vector("list",length(subgene))
    tmp_intronII2 <- split(tmp_intronII2,tmp_intronII2$group_name)
    names(pre_intronII2) <- names(subgene)
    pre_intronII2[names(tmp_intronII2)] <- tmp_intronII2
    intronII2 <- unlist(cmclapply(names(pre_intronII2),function(n,D) {
        x <- D[[n]]
        return(sum(x$score))
    },pre_intronII2,rc=0.5))
    anti_intronII2 <- unlist(cmclapply(names(pre_intronII2),function(n,D) {
        x <- D[[n]]
        if (!is.null(x)) {
            if (x$strand[1] == "-")
                return(sum(x$sense_score))
            return(sum(x$antisense_score))
        }
        return(0)
    },pre_intronII2,rc=0.5))
    names(intronII2) <- names(anti_intronII2) <- names(pre_intronII2)

    # Some exons are now empty because of artifact filtering
    message("    selected for exons")
    pre_exonII1 <- as.data.frame(mcols(exonStru1,level="within"))
    pre_exonII1 <- split(pre_exonII1,pre_exonII1$group_name)
    #pre_exonII1 <- pre_exonII1[names(subgene)]
    exonII1 <- unlist(cmclapply(names(pre_exonII1),function(n,D) {
        x <- D[[n]]
        return(sum(x$score))
    },pre_exonII1,rc=0.5))
    names(exonII1) <- names(pre_exonII1)
    message("    control for exons")
    pre_exonII2 <- as.data.frame(mcols(exonStru2,level="within"))
    pre_exonII2 <- split(pre_exonII2,pre_exonII2$group_name)
    #pre_exonII2 <- pre_exonII2[names(subgene)]
    exonII2 <- unlist(cmclapply(names(pre_exonII2),function(n,D) {
        x <- D[[n]]
        return(sum(x$score))
    },pre_exonII2,rc=0.5))
    names(exonII2) <- names(pre_exonII2)
    
    # Harmonize the resulted vectors of IIs as there are fitlered genes not
    # present in either of intron/exon structures
    theNamesList <- list(names(exonII1),names(exonII2),names(intronII1),
        names(intronII2),names(anti_intronII1),names(anti_intronII2))
    theNames <- Reduce("intersect",theNamesList)
    
    # Harmonize
    exonII1 <- exonII1[theNames]
    exonII2 <- exonII2[theNames]
    intronII1 <- intronII1[theNames]
    intronII2 <- intronII2[theNames]
    anti_intronII1 <- anti_intronII1[theNames]
    anti_intronII2 <- anti_intronII2[theNames]
    
    # And calculate II, DI
    II1 <- exonII1 + intronII1 + anti_intronII1
    II2 <- exonII2 + intronII2 + anti_intronII2

    # Disrupting insertions (D.I.) are calculated as the sum of exonic I.I., and
    # intronic I.I. with the splice acceptor aligned with transcription
    DI1 <- exonII1 + intronII1
    DI2 <- exonII2 + intronII2

    # A bias is calculated as the ratio of intronic insertions in sense over
    # anti-sense orientation for each gene
    message("  calculating Bias")
    message("    selected")
    Bias1 <- unlist(cmclapply(names(pre_intronII1),function(n,D) {
        x <- D[[n]]
        if (!is.null(x)) {
            a <- sum(x$score)
            if (x$strand[1] == "-") {
                b <- sum(x$sense_score)
                if (a >= 0 && b == 0)
                    return(0)
                return(a/b)
            }
            b <- sum(x$antisense_score)
            if (a >= 0 && b == 0)
                return(0)
            return(a/b)
        }
        return(0)
    },pre_intronII1,rc=0.5))
    names(Bias1) <- names(pre_intronII1)
    message("    control")
    Bias2 <- unlist(cmclapply(names(pre_intronII2),function(n,D) {
        x <- D[[n]]
        if (!is.null(x)) {
            a <- sum(x$score)
            if (x$strand[1] == "-") {
                b <- sum(x$sense_score)
                if (a >= 0 && b == 0)
                    return(0)
                return(a/b)
            }
            b <- sum(x$antisense_score)
            if (a >= 0 && b == 0)
                return(0)
            return(a/b)
        }
        return(0)
    },pre_intronII2,rc=0.5))
    names(Bias2) <- names(pre_intronII2)
    
    # Harmonize
    Bias1 <- Bias1[theNames]
    Bias2 <- Bias2[theNames]
    
    message("  calculating overlaping reads to report")
    message("    selected")
    o1 <- findOverlaps(gene,preII1)
    indexByGene1 <- split(subjectHits(o1),queryHits(o1))
    re1 <- unlist(cmclapply(indexByGene1,function(x,R) {
        return(sum(R[x]$reads))
    },preII1,rc=rc))
    names(re1) <- names(gene)[as.numeric(names(re1))]
    message("    control")
    o2 <- findOverlaps(gene,preII2)
    indexByGene2 <- split(subjectHits(o2),queryHits(o2))
    re2 <- unlist(cmclapply(indexByGene2,function(x,R) {
        return(sum(R[x]$reads))
    },preII2,rc=rc))
    names(re2) <- names(gene)[as.numeric(names(re2))]
    
    message("  preparing triplets output")
    out1 <- matrix(0,length(gene),4)
    rownames(out1) <- names(gene)
    out2 <- matrix(0,length(gene),4)
    rownames(out2) <- names(gene)
    
    out1[names(II1),1] <- II1
    out1[names(DI1),2] <- DI1
    out1[names(Bias1),3] <- Bias1
    out1[names(re1),4] <- re1
    out2[names(II2),1] <- II2
    out2[names(DI2),2] <- DI2
    out2[names(Bias2),3] <- Bias2
    out2[names(re2),4] <- re2
    
    colnames(out1) <- colnames(out2) <- c("II","DI","Bias","Reads")

    gc(verbose=FALSE)

    return(list(selected=out1,control=out2))
}

getLibsize <- function(bams) {
    libsize <- numeric(length(bams))
    names(libsize) <- bams
    for (bam in bams) {
        message("    ",basename(bam))
        libsize[[bam]] <- length(readGAlignments(file=bam))
    }
    return(libsize)
}

coverageFromRanges <- function(input,mask,rc=NULL) {
    if (is(mask,"GRanges")) {
        chrs <- as.character(unique(seqnames(input)))
        preCov <- coverage(input)
        preCov <- preCov[chrs]
        maskList <- split(mask,seqnames(mask))
        covs <- cmclapply(names(maskList),function(x,maskList,preCov) {
            return(lazyRangesCoverage(x,maskList,preCov))
        },maskList,preCov,rc=rc)
        covs <- unlist(covs)
        names(covs) <- names(mask)
        return(covs)
    }
}

procWin <- function(n,s,y) {
    # Get member
    x <- y[[n]]
    
    # Get window coordinates from name
    tmp <- strsplit(n,":|-")
    chr <- tmp[[1]][1]
    start <- as.numeric(tmp[[1]][2])
    end <- as.numeric(tmp[[1]][3])
    
    len <- runLength(x)
    if (length(len) > 3) {
        # Values
        val <- runValue(x)
        
        preSummit <- which(x == max(x))
        # Ties? Use first position
        if (length(preSummit) > 1)
            preSummit <- preSummit[1]
            
        # Rest to be added to max
        m <- which(val == max(val))
        len <- len[-m]
        val <- val[-m]
        rest <- which(val > 0)
        add <- sum(len[rest]*val[rest])
        
        # Final
        int <- as.numeric(x[preSummit] + add)
    }
    else {
        preSummit <- which(x == max(x))
        # Ties? Use first position
        if (length(preSummit) > 1)
            preSummit <- preSummit[1]
        int <- as.numeric(x[preSummit])
    }
    pos <- start + preSummit - 1
    return(data.frame(chromosome=chr,start=pos,end=pos,strand=s,score=int))
}

extendGeneRanges <- function(ranges,flank=c(500,500)) {
    w <- width(ranges)
    ranges <- promoters(ranges,upstream=flank[1],downstream=0)
    return(resize(ranges,width=w+flank[1]+flank[2]))
}

lazyRangesCoverage <- function(x,maskList,preCov) {
    message("        processing ",x)
    m <- maskList[[x]]
    pre <- preCov[[x]]
    if (!is.null(m) && !is.null(pre)) { # Sanity...
        V <- Views(pre,ranges(m))
        cot <- unlist(viewApply(V,function(x) x))
        names(cot) <- names(m)
        inv <- which(strand(m)=="-")
        cot[inv] <- lapply(cot[inv],rev)
        return(cot)
    }
    else {
        message("        ",x," not found!")
        return(Rle(NA))
    }
}

getChromInfo <- function(org,
    goldenPath="http://hgdownload.cse.ucsc.edu/goldenPath/") {
    download.file(paste(goldenPath,org,"/database/chromInfo.txt.gz",sep=""),
        file.path(tempdir(),"chromInfo.txt.gz"))
    chromInfo <- read.delim(file.path(tempdir(),"chromInfo.txt.gz"),
        header=FALSE)
    chromInfo <- chromInfo[,1:2]
    chromInfo[,1] <- as.character(chromInfo[,1])
    chromInfo$V3 <- rep(FALSE,nrow(chromInfo))
    m <- grep("M",chromInfo[,1])
    if (length(m) > 0)
        chromInfo$V3[m] <- TRUE
    return(chromInfo)
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
            else 
                m <- FALSE
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}
