source("/media/local/software/hts-tools/R/diagplots.R")

color.unreads <- function(ur,st=TRUE,log.it="FALSE",output="x11",fil=NULL,...)
{
    if (missing(ur))
        stop("Please specify the file containing the unique reads and the respective filenames!")

    r.f <- as.matrix(read.delim(ur,row.names=1,header=FALSE))
    r.f <- as.matrix(r.f[-nrow(r.f),])
    if (log.it)
    r.f <- log2(r.f)

    labels <- basename(rownames(r.f))
    in.labs <- formatC(r.f[,1],digits=2,format="e")
    vals <- 1-r.f[,1]/max(r.f[,1])
    if (st)
    {
        stru <- sort(basename(rownames(r.f)),index.return=TRUE)
        labels <- stru$x
        vals <- vals[stru$ix]
        in.labs <- in.labs[stru$ix]
    } 
  
    openGraphics(output,fil)
  
    flip <- nrow(r.f):1
    in.labs <- in.labs[flip]
    par(mai=c(0.1,3,0.1,1.5),mgp=c(3,1,0.5))
    image(1,1:nrow(r.f),t(vals[flip]),axes=FALSE,xlab="",ylab="",cex=1.5)
    axis(2,at=seq(1,nrow(r.f),by=1),labels=labels[flip],las=2,col="white",cex.axis=0.8)
    axis(4,at=seq(1,nrow(r.f),by=1),labels=r.f[stru$ix[flip],1],las=2,col="white",cex.axis=0.8)
    for (i in 1:nrow(r.f))
        text(1,i,in.labs[i],col="blue",cex=0.6)                            

    closeGraphics(output)
}

cov.bar <- function(btcov,what=c("genome","chromosome"),output="x11",fil=NULL,...)
{
    if (missing(btcov))
        stop("Please specify the file containing coverage data!")
    what <- tolower(what[1])
    if (what!="genome" && what!="chromosome")
        stop("what must be one of \"genome\", \"chromosome\"")

    cov.raw <- read.delim(btcov,header=FALSE)
    names(cov.raw) <- c("chr","depth","nbase","chrlen","frac")

    cov.genome<-cov.raw$frac[grep("genome",cov.raw$chr)]
    names(cov.genome) <- as.character(0:(length(cov.genome)-1))
    barlab <- paste(formatC(100*cov.genome,digits=3,format="f"),"%",sep="")             
    x <- -1/log10(cov.genome)
  
    openGraphics(output,fil)

    if (what=="genome")
    {
        par(mai=c(0.8,0.8,0.5,0.5))
        barx <- barplot(x,space=0,ylim=c(0,max(x)+1),col="red3",border="yellow3",xaxt="n",yaxt="n")
        axis(1,at=seq(0.5,length(x)-0.5,by=1),labels=c(0:(length(x)-1)),cex.axis=0.5,font=2,padj=-2,tcl=-0.3,col="lightgrey")
        axis(2,cex.axis=0.7,padj=1)
        text(x=barx,y=x+0.3,label=barlab,cex=0.6,srt=90,font=2,col="blue")
        title(main="Coverage plot",cex.main=1)
        mtext(side=1,text="coverage depth",line=2,cex=0.9,font=2)
        mtext(side=2,text="-1/log10(coverage)",line=2,cex=0.9,font=2)
        grid()
    }
  
    if (what=="chromosome")
    {
        chromosomes <- unique(as.character(cov.raw$chr[grep("chr",cov.raw$chr)]))
        chromosomes <- suppressWarnings(chromosomes[order(as.numeric(substr(chromosomes,start=4,stop=nchar(chromosomes))))])
        cov.chr <- inv.log.cov <- bar.labs <- vector("list",length(chromosomes))
        names(cov.chr) <- chromosomes
        for (chr in chromosomes)
        {
            cov.chr[[chr]] <- cov.raw$frac[which(cov.raw$chr==chr)]
            names(cov.chr[[chr]]) <- as.character(0:(length(cov.chr[[chr]])-1))
            bar.labs <- paste(formatC(100*cov.chr[[chr]],digits=3,format="f"),"%",sep="")             
            inv.log.cov[[chr]] <- -1/log10(cov.chr[[chr]])
        }

        par(mfrow=c(4,6),mar=c(1,1,1,1),oma=c(1,1,0,0),mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7)
        for (chr in chromosomes)
        {
            x <- inv.log.cov[[chr]]
            barx <- barplot(x,space=0,ylim=c(0,max(x)+1),col="red3",border="red3",xaxt="n",yaxt="n")
            title(main=chr,cex.main=1)
            grid(nx=0,ny=NULL)
        }
    }
  
    closeGraphics(output)
}

htseqtools.qc <- function(input,classes,qctype=c("mds","ssd","gini"),n.cores=1,output="x11",fil=NULL,...)
{
    if (missing(input))
        stop("Please specify input files!")
    if (missing(classes))
        stop("Please specify class vector!")
    # An additional control for qctype using htSeqTools
  
    if (!require(IRanges))
        stop("Bioconductor package IRanges is required!")
    if (!require(htSeqTools))
        stop("Bioconductor package htSeqTools is required!")
    if (!require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")
  
    qctype <- tolower(qctype)
    if (!all(is.element(qctype,c("mds","ssd","gini"))))
        stop("Please provide only supported htSeqTools QC methods (\"mds\", \"ssd\", \"gini\")!")
 
    if (n.cores>1)
    {
        if (!require(multicore))
        {
            warning("R package multicore not present... Switching to single core...",
                    call.=FALSE)
            n.cores <- 1
        }
        else
        {
            require(multicore)
            if (multicore:::detectCores()==1)
            {
                warning("Only one core detected... Switching to single core...",
                        call.=FALSE)
                n.cores <- 1
            }
        }
    }
  
    if (!is.null(fil)) # We will auto-generate some names
        real.fil <- prepareAutoGenName(fil,output)

    nams <- basename(input)
    classes <- as.factor(classes)
    design <- as.numeric(classes)
    colspace <- c("red","blue","yellowgreen","orange","aquamarine2",
                  "pink2","seagreen4","brown","purple","chocolate")
    pchspace <- c(20,17,15,16,8,3,2,0,1,4)
    rgdObj <- RangedDataList()
  
    for (i in 1:length(input))
    {  
        cat("Importing track ",nams[i],"... Please wait...\n",sep="")
        flush.console()
        rgdObj[[i]] <- import.bed(input[i],trackLine=FALSE,asRangedData=TRUE)
    }
    gc(verbose=FALSE)
    names(rgdObj) <- nams
  
    if ("mds" %in% qctype)
    {
        cmdsObj <- cmds(rgdObj,k=2,mc.cores=n.cores)
        real.fil[3] <- "_MDS"
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                           sep=.Platform$file.sep)
        openGraphics(output,fig.fname)
        plot(cmdsObj,col=colspace[1:length(levels(classes))][design],
             pch=pchspace[1:length(levels(classes))][design],
             cex=0.7,cex.axis=0.9,cex.main=0.9)
        grid()
        closeGraphics(output)
    }

    if (is.element(output,c("png","jpg","tiff","bmp","x11")))
        d=c(800,400)
    else if (is.element(output,c("pdf","ps")))
        d=c(12,6)
     
    if ("ssd" %in% qctype)
    {
        ssd <- ssdCoverage(rgdObj)
        lens <- sapply(rgdObj,nrow)
        ssdu <- ssd*lens
        real.fil[3] <- "_SSD"
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                           sep=.Platform$file.sep)
        openGraphics(output,fig.fname,width=d[1],height=d[2])
        par(mai=c(2,0.8,0.8,0.5),lwd=2,mfcol=c(1,2))
        bar.ssd <- barplot(ssd,col="darkred",border="yellow",xaxt="n",yaxt="n")
        axis(1,at=bar.ssd,,labels=nams,cex.axis=0.5,font=2,las=2)
        axis(2,cex.axis=0.7,padj=1)
        title(main="Adjusted coverage deviation",cex.main=1)
        mtext(side=2,text="Coverage AdjSD",line=2,cex=0.9,font=2)
        bar.ssdu <- barplot(ssdu,col="darkblue",border="orange",xaxt="n",yaxt="n")
        axis(1,at=bar.ssdu,,labels=nams,cex.axis=0.5,font=2,las=2)
        axis(2,cex.axis=0.7,padj=1)
        title(main="Coverage deviation",cex.main=1)
        mtext(side=2,text="Coverage SD",line=2,cex=0.9,font=2)
        closeGraphics(output)
    }
  
    if ("gini" %in% qctype)
    {
        gini <- giniCoverage(rgdObj)
        real.fil[3] <- "_GINI"
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                           sep=.Platform$file.sep)
        openGraphics(output,fig.fname,width=d[1],height=d[2])
        par(mai=c(2,0.8,0.8,0.5),lwd=2,mfcol=c(1,2))
        bar.gin <- barplot(gini[,2],col="darkred",border="yellow",xaxt="n",yaxt="n")
        axis(1,at=bar.gin,,labels=nams,cex.axis=0.5,font=2,las=2)
        axis(2,cex.axis=0.7,padj=1)
        title(main="Adjusted coverage Gini index",cex.main=1)
        mtext(side=2,text="Coverage AdjGini",line=2,cex=0.9,font=2)
        bar.ginu <- barplot(gini[,1],col="darkblue",border="orange",xaxt="n",yaxt="n")
        axis(1,at=bar.ginu,,labels=nams,cex.axis=0.5,font=2,las=2)
        axis(2,cex.axis=0.7,padj=1)
        title(main="Coverage Gini index",cex.main=1)
        mtext(side=2,text="Coverage Gini",line=2,cex=0.9,font=2)
        closeGraphics(output)
    }
}

# plot.method: default heatmap, see function plot.cor in diagplots.R for other possibilities
cor.qc <- function(input,classes,cor.type=c("counts","pca","mds","hilbert"),plot.method="heatmap",org="hg18",
                   annot="legend",chr.size.file=NULL,winsize=10000,save.rdata=FALSE,output="x11",fil=NULL,...)
{
    if (missing(input))
        stop("Please specify input files!")
    if (missing(classes))
        stop("Please specify class vector!")
    cor.type <- tolower(cor.type)
  
    if (!require(IRanges))
        stop("Bioconductor package IRanges is required!")
    if (!require(GenomicRanges))
        stop("Bioconductor package GenomicRanges is required!")
    if ("hilbert" %in% cor.type && !require(HilbertVis))
        stop("Bioconductor package HilbertVis is required!")
    if (!require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")
    if (!require(gplots))
        stop("R package gplots is required!")
  
    cor.type <- tolower(cor.type)
    if (!all(is.element(cor.type,c("counts","pca","mds","hilbert"))))
        stop("Please provide one of (\"counts\", \"pca\", \"mds\" or \"hilbert\")!")

    if (is.null(chr.size.file) && (missing(org) || !is.element(org,c("hg18","hg19","mm8","mm9"))))
        stop("Please provide a valid organism or a chromosome size file!")
    if (!is.null(chr.size.file) && is.element(org,c("hg18","hg19","mm8","mm9")))
    {
        warning("Both chromosome size file and organism provided... Will use the file...",call.=FALSE)
        useFile <- TRUE
    } else { useFile <- FALSE }
  
    if (!is.null(fil))
        real.fil <- prepareAutoGenName(fil,output)  
  
    if (useFile)
    {
        chrom.info <- read.delim(chr.size.file,header=FALSE)
        chrom.info <- chrom.info[-grep("rand|chrM|hap",chrom.info[,1]),]
    }
    else
        chrom.info <- getChromInfo(org)
    
    chrname <- as.character(chrom.info[,1])
    chrlen <- as.numeric(chrom.info[,2])
    names(chrlen) <- chrname
  
    w <- make.windows(chrlen=chrlen,winsize=winsize,chrname=chrname)
    wl <- make.windows.list(chrlen=chrlen,winsize=winsize,chrname=chrname)

    nams <- basename(input)
    if (any(c("counts","pca","mds") %in% cor.type))
        window.counts <- matrix(0,length(w),length(input))

    classes <- as.factor(classes)
    design <- as.numeric(classes)
    colspace <- c("red","blue","yellowgreen","orange","aquamarine2",
                  "pink2","seagreen4","brown","purple","chocolate")
    pchspace <- c(20,17,15,16,8,3,2,0,1,4)

    for (i in 1:length(input))
    {
        cat("\nImporting and processing file ",nams[i],"...",sep="")
        flush.console()
        bed <- import.bed(input[i],trackLine=FALSE,asRangedData=FALSE)

        # Find the total window counts
        if ("counts" %in% cor.type)
            window.counts[,i] <- countOverlaps(w,bed)

        # Look for coverage in genomic features
        #if ("coverage" %in% cor.type){ }

        # Try the hilbert curve...
        if ("hilbert" %in% cor.type)
        {
            beds <- split(bed,f=as.factor(seqnames(bed)))
            for (chr in as.character(runValue(seqnames(bed))))
            {
                cat("\n  Creating wiggle vector for ",chr,sep="")
                wigvec <- makeWiggleVector(start=start(beds[[chr]]),end=end(beds[[chr]]),
                            value=rep(1,length(start(beds[[chr]]))),chrlength=chrlen[[chr]])
                cat("\n  Creating Hilbert image for ",chr,sep="")
                hMat <- hilbertImage(wigvec,level=10)     
                real.fil[3] <- paste("_",nams[i],"_HILBERT_",chr,sep="")
                fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                                   sep=.Platform$file.sep)
                openGraphics(output,fig.fname)
                #showHilbertImage(hMat)
                print(showHilbertImage(hMat)) # Damn trellis bug!!!
                closeGraphics(output) 
            }
        }

        gc(verbose=FALSE) 
    }
 
    colnames(window.counts) <- nams

    if (save.rdata)
        save(window.counts,file=file.path(real.fil[1],paste(real.fil[2],".RData",sep="")))

    # Filter completely dead zones to avoid bias
    not.dead <- which(apply(window.counts,1,any))
    window.counts <- window.counts[not.dead,]

    #if (save.rdata)
    #   save(window.counts,file=file.path(real.fil[1],paste(real.fil[2],".RData",sep="")))

    # Calculate correlations and make plots
    if ("counts" %in% cor.type)
    {
        n <- length(input)
        cor.mat <- cor(window.counts,method="spearman")
        colnames(cor.mat) <- colnames(window.counts)
        labs <- matrix(NA,n,n)
        for (i in 1:n)
            for (j in 1:n)
                labs[i,j] <- sprintf("%.2f",cor.mat[i,j])
        real.fil[3] <- paste("_CORHEATMAP",sep="")
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                         sep=.Platform$file.sep)
        plot.cor(cor.mat,type=plot.method,is.cor.mat=TRUE,output=output,fil=fil,...)
    }

    if ("pca" %in% cor.type)
    {
        pcaObj <- prcomp(nat2log(window.counts),scale=TRUE)
        real.fil[3] <- "_PCA"
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                           sep=.Platform$file.sep)
        if (output %in% c("png","jpg","bmp","tiff"))
            openGraphics(output,fil,width=640,height=480)
        else if (output %in% c("ps","pdf"))
            openGraphics(output,fil,width=9,height=7)
        plot(pcaObj$rotation,col=colspace[1:length(levels(classes))][design],
             pch=pchspace[1:length(levels(classes))][design],
             cex=0.7,cex.axis=0.9,cex.main=0.9)
        if (annot=="legend")
            legend(x="topright",cex=0.7,
                   legend=as.character(unique(classes)),
                   col=colspace[1:length(levels(classes))][unique(design)],
                   pch=pchspace[1:length(levels(classes))][unique(design)])
        else if (annot=="text")
            text(pcaObj$rotation,labels=nams,pos=3,cex=0.7)
        grid()
        closeGraphics(output)
    }

    if ("mds" %in% cor.type)
    {
        d <- as.dist(0.5*(1-cor(nat2log(window.counts),method="spearman")))
        mdsObj <- cmdscale(d,eig=TRUE,k=2)
        real.fil[3] <- "_MDS"
        fig.fname <- paste(real.fil[1],paste(real.fil[2],real.fil[3],real.fil[4],sep=""),
                           sep=.Platform$file.sep)
        if (output %in% c("png","jpg","bmp","tiff"))
            openGraphics(output,fil,width=640,height=480)
        else if (output %in% c("ps","pdf"))
            openGraphics(output,fil,width=9,height=7)
        plot(mdsObj$points[,1],mdsObj$points[,2],
             xlab="MDS 1",ylab="MDS 2",
             col=colspace[1:length(levels(classes))][design],
             pch=pchspace[1:length(levels(classes))][design],
             cex=0.7,cex.axis=0.9,cex.main=0.9)
        if (annot=="legend")
            legend(x="topright",cex=0.7,
                   legend=as.character(unique(classes)),
                   col=colspace[1:length(levels(classes))][unique(design)],
                   pch=pchspace[1:length(levels(classes))][unique(design)])
        else if (annot=="text")
            text(mdsObj$points[,1],mdsObj$points[,2],labels=nams,pos=3,cex=0.7)
        grid()
        closeGraphics(output)
    }

    #if (save.rdata)
    #   save(window.counts,file=paste(real.fil[1],paste(real.fil[2],".RData",sep=""),sep=.Platform$file.sep))
}

plot.sat.macs <- function(macs.file,output="x11",fil=NULL,...)
{
    if (missing(macs.file))
        stop("Please provide a MACS diagnostic file!")
  
    main.title <- paste("Peak saturation analysis for",basename(macs.file))

    macs.out <- read.delim(macs.file,check.names=FALSE)
    dims <- dim(macs.out)
    sat.table <- as.matrix(macs.out[,3:dims[2]])
    sat.table <- cbind(rep(100,dims[1]),sat.table)
    xnames <- c("100%","90%",names(macs.out)[4:dims[2]])
    ynames <- c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
    leg.txt <- as.character(macs.out[,1])
    red.dims <- dim(sat.table)
    x <- 1:ncol(sat.table)

    # Generate colors for the points
    color <- rgb(matrix(runif(prod(red.dims),0,1),red.dims[1],3))
    pch.ids <- c(1:4,8,15:25)
    pchs <- runif(red.dims[1],1,length(pch.ids))

    openGraphics(output,fil)

    # Plot series
    par(cex.axis=0.8,cex.main=1.2,cex.lab=1,font.lab=2,col.main="blue",pty="m")
    plot.new()
    plot.window(c(1,length(xnames)+1),c(0,100))
    axis(1,at=x,labels=xnames)
    axis(2,at= seq(0,100,10),labels=ynames)
    title(main=main.title,xlab="Percentage of tags in sample",ylab="Percentage of identified peaks")
    for (i in 1:red.dims[1])
    {
        points(x,sat.table[i,],pch=pchs[i],cex=0.8,col=color[i],bg=color[i])
        lines(x,sat.table[i,],cex=0.8,col=color[i],lty=5)
    }
    grid(col="darkgrey")
  
    # Legend stuff, put outside plotting region on the right
    # Set xpd=TRUE, so plots are clipped to the figure (not the plot) region:
    old.par <- par(xpd=TRUE)
    # Add the legend outside the plot on the right, centered on the plot y axis
    y.pos <- mean(par("usr")[3:4])
    # Calculate the ratio of inches to x user coords
    in2x.pos <- (par("usr")[2]-par("usr")[1])/par("pin")[1]
    # Place x at the right of the margin
    x.pos <- par("usr")[2]-par("mai")[4]*in2x.pos
    # Justify the legend so that is centered in y (yjust=0.5) and on top of x (xjust=0)
    legend(x.pos,y.pos,legend=leg.txt,col=color,text.col=color,lwd=1,lty=1,pch=pchs,
         merge=TRUE,cex=0.7,box.lwd=0.5,xjust=0.2,yjust=0.5)
    par(old.par)

    closeGraphics(output)
}

make.windows.list <- function(chrlen,winsize,chrname)
{
    names(chrlen) <- chrname
    winObj <- GRangesList()

    rem <- chrlen%%winsize
    len <- floor(chrlen/winsize)
    clen <- len + 1
    snames <- Rle(chrname,clen)
  
    for (i in 1:length(chrlen))
    { 
        starts <- seq(1,by=winsize,length=len[i])
        ends <- seq(winsize,by=winsize,length=len[i])
        last.start <- len[i]*winsize + 1
        last.end <- last.start + rem[i] - 1
        winObj[[chrname[i]]] <- GRanges(seqnames=snames[snames==chrname[i]],
                                    IRanges(start=c(starts,last.start),
                                    end=c(ends,last.end)),
                                    seqlengths=chrlen)
    }
  
    return(winObj)
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

getChromInfo <- function(org)
{
    download.file(paste("http://hgdownload.cse.ucsc.edu/goldenPath/",org,"/database/chromInfo.txt.gz",sep=""),
                file.path(tempdir(),"chromInfo.txt.gz"))
    if (.Platform$OS.type == "unix")
        system(paste("gzip -df",file.path(tempdir(),"chromInfo.txt.gz")))
    else
        unzip(file.path(tempdir(),"chromInfo.txt.gz"))
    chrom.info <- read.delim(file.path(tempdir(),"chromInfo.txt"),header=FALSE)
    chrom.info <- chrom.info[-grep("rand|chrM|hap",chrom.info[,1]),]
    return(chrom.info)
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

nat2log <- function(x,base=2)
{
    x[x==0] <- 1
    if (base==2)
        return(log2(x))
    else
        return(log10(x))
}

### This proved to be EXTREMELY slow... We will switch to BEDTools or GenomicTools with
### streaming and real time output parsing
# Find the total window coverages
#cov.bed <- coverage(bed) # This returns a list with chromosomes
#theCounter <- 0
##curr.start <- 0
##curr.stop <- 0
#for (chr in names(wl))
#{
#  cat("\nChromosome ",chr,"...\n",sep="")
#  cov.view <- Views(cov.bed[[chr]],start=start(wl[[chr]]),end=end(wl[[chr]]))
#  pb <- txtProgressBar(max=length(cov.view),char=">",style=3)
#  for (j in 1:length(cov.view))
#  {
#    theCounter <- theCounter + 1
#    window.cov[theCounter,i] <- mean(cov.view[j])
#    setTxtProgressBar(pb,j)
#  }
#  #curr.start <- curr.stop + 1
#  #curr.stop <- curr.start + length(wl[[chr]]) - 2
#  #window.cov[curr.start:curr.stop,i] <- sapply(cov.view[1:(length(cov.view)-1)],mean)
#  close(pb)
#}
