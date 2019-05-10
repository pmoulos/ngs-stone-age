# name : vector of names to be used
# y.lim : how to construct y-axis limits? "default", "auto" or a 2-length vector, defaults
# to "default" (internal R plotting system is taking care of this)
# output : output plotting device (pdf, png etc.), defaults to x11 and if "same", no device
#          will be called so that user can plot many bopxlots using par externally. In this
#          case the second matrix (tam, if given) will be disarded
# ... further arguments to boxplot
boxplot.mat <- function(mat,tam,name=NULL,log.it="auto",y.lim="default",output="x11",fil=NULL,...)
{
    # We have both normalized and un-normalized data?
    two <- TRUE
    if (missing(tam))
        two <- FALSE

    # Need to log?
    if (log.it=="auto")
    {
        if (diff(range(mat,na.rm=TRUE))>1000) # Data would rather be log scaled
        {
            mat <- log2disp(mat)
            # Check the same for tam
            if (two)
                if (diff(range(tam,na.rm=TRUE))>1000)
                    tam <- log2disp(tam)
        }
    }
    else if (log.it=="yes")
    {
        mat <- log2disp(mat)
        if (two)
            tam <- log2disp(tam)
    }

    # Define the axis limits based on user input
    if (!is.numeric(y.lim) && y.lim=="auto")
    {
        if (two)
        {
            min.y <- floor(min(min(mat),min(tam)))
            max.y <- ceiling(max(max(mat),max(tam)))
        }
        else
        {
            min.y <- floor(min(mat))
            max.y <- ceiling(max(mat))
        }
    }
    else if (is.numeric(y.lim))
    {
        min.y <- y.lim[1]
        max.y <- y.lim[2]
    }
    
    if (is.null(name))
        nams <- paste("Sample",1:ncol(mat),sep=" ")
    else if (length(name)==1 && name=="none")
        nams <- rep("",ncol(mat))
    else
        nams <- name
    
    cols <- c("red","green","blue","yellow","aquamarine","orange","burlywood")
    mat.list <- list()
    for (i in 1:ncol(mat))
        mat.list[[i]] <- mat[,i]
    if (two)
    {
        tam.list <- list()
        for (i in 1:ncol(tam))
            tam.list[[i]] <- tam[,i]
    }

    if (output!="same") # So as to be able to use par(mfxxx) from outside
        openGraphics(output,fil)

    if (two)
        par(mfrow=c(1,2))
    if (!is.numeric(y.lim) && y.lim=="default")
    {
        boxplot(mat.list,names=nams,col=cols,las=2,...)
        if (two)
            boxplot(tam.list,names=nams,col=cols,las=2,...)
    }   
    else
    {
        boxplot(mat.list,names=nams,col=cols,ylim=c(min.y,max.y),las=2,...);
        if (two)
            boxplot(tam.list,names=nams,ylim=c(min.y,max.y),col=cols,las=2,...)
    }

    if (output!="same")
        closeGraphics(output)
}

plot.cor <- function(mat,method="pearson",type="circle",log.mat=FALSE,is.cor.mat=FALSE,output="x11",fil=NULL,...)
{
    method <- tolower(method[1])
    if (!is.element(method,c("pearson","spearman")))
        stop("method must be one of \"pearson\" or \"spearman\"!")
    type <- tolower(type[1])
    if (!is.element(type,c("circle","ellipse","heatmap","cgram.shade","cgram.pts")))
        stop("type must be one of \"circle\", \"ellipse\", \"heatmap\", \"cgram.shade\" or \"cgram.pts\"!")
    
    if (!require(corrplot) && (type %in% c("circles","ellipse")))
        stop("R package corrplot is required!")
    if (!require(corrgram) && effect %in% c("cgram.shade","cgram.pts"))
        stop("R package corrgram is required!")

    if (log.mat)
        mat <- log2(mat)
    if (is.cor.mat)
        cor.mat <- mat
    else
        cor.mat <- cor(mat,method=method)
    if (!is.null(colnames(mat)))
        colnames(cor.mat) <- colnames(mat)

    openGraphics(output,fil,...)

    if (type %in% c("circle","ellipse"))
        corrplot(cor.mat,method=type,order="hclust",...)
    else if (type=="cgram.pts")
    {
        corrgram(mat,type="data",order=TRUE,lower.panel=panel.pts,upper.panel=panel.conf,text.panel=panel.txt,diag.panel=panel.minmax,
                 pch=20,col="blue",cex=0.8)
    }
    else if (type=="cgram.shade")
    {
        corrgram(cor.mat,type="cor",order=TRUE,lower.panel=panel.shade,upper.panel=panel.conf,text.panel=panel.txt,diag.panel=panel.minmax)
    }
    else if (type=="heatmap")
    {
        n <- dim(cor.mat)[1]
        labs <- matrix(NA,n,n)
        for (i in 1:n)
            for (j in 1:n)
                labs[i,j] <- sprintf("%.2f",cor.mat[i,j])
        if (n <= 5)
            notecex <- 1.2
        else if (n > 5 & n < 10)
            notecex <- 0.9
        else
            notecex <- 0.7
        #my.heatmap.2(cor.mat,col=colorRampPalette(c("yellow","grey","blue")),revC=TRUE,trace="none",symm=TRUE,Colv=TRUE,#symkey=TRUE,
        #            cellnote=labs,keysize=1,density.info="density",notecex=notecex,cexRow=1,cexCol=1,margins=c(13,13),slimkey=TRUE)
        heatmap.2(cor.mat,col=colorRampPalette(c("yellow","grey","blue")),revC=TRUE,trace="none",symm=TRUE,Colv=TRUE,#symkey=TRUE,
                     cellnote=labs,keysize=1,density.info="density",notecex=notecex,cexRow=0.9,cexCol=0.9,margins=c(13,13))
    }
    closeGraphics(output)
}

plot.pairs <- function(x,type="simple",output="x11",fil=NULL,...)
{
    type <- tolower(type[1])
    if (!is.element(type,c("simple","image")))
        stop("type must be one of \"simple\" or \"image\"!")
    
    if (!require(IDPmisc) && type=="image")
        stop("R package IDPmisc is required!")
    
    n <- ncol(x)
    if (!is.null(colnames(x)))
        nams <- colnames(x)
    else
        nams <- paste("Sample_",1:ncol(x),sep="")

    if (output!="same")
        openGraphics(output,fil)
    
    # Setup the grid
    par(mfrow=c(n,n),mar=c(1,1,1,1),oma=c(1,1,0,0),mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7)

    # Plot
    for (i in 1:n)
    {
        for (j in 1:n)
        {
            #if (i==j && i!=n && j!=n)
            if (i==j)
            {
                plot(0:10,0:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="") # Diagonal
                #text(5,5,textwrap(nams[i],15),cex=1.2)
                #text(c(3,5,3),c(9.5,5,1),c("X-Y plots",textwrap(nams[i],15),"M-D plots"),cex=c(1,1.2,1))
                text(c(3,5,3),c(9.5,5,1),c("X-Y plots",nams[i],"M-D plots"),cex=c(1,1.2,1))
                arrows(6,9.5,9.5,9.5,angle=20,length=0.1,lwd=1.1)
                arrows(0.2,3.2,0.2,0.2,angle=20,length=0.1,lwd=1.1)
            }
            else if (i<j) # XY plot
            {
                if (type=="simple")
                    plot(x[,i],x[,j],pch=20,col="blue",cex=0.4,xlab=nams[i],ylab=nams[j],...)
                else 
                {
                    plot(x[,i],x[,j],type="n")
                    Image(x[,i],x[,j],pixs=0.5)
                }
                lines(lowess(x[,i],x[,j]),col="red")
                grid()
            }
            else if (i>j) # MD plot
            {
                if (type=="simple")
                    plot((x[,i]+x[,j])/2,x[,j]-x[,i],pch=20,col="blue",cex=0.4,...)
                else
                {
                    plot((x[,i]+x[,j])/2,x[,j]-x[,i],type="n")
                    Image((x[,i]+x[,j])/2,x[,j]-x[,i],pixs=0.5)
                }
                lines(lowess((x[,i]+x[,j])/2,x[,j]-x[,i]),col="red")
                grid()
            }
            # TODO: Display somehow the colorbar in when type="image"
            #iplotLegend(IDPcolorRamp,cex.axis=1.5,mar=c(3,6,3,6))
        }
    }

    if (output!="same")
        closeGraphics(output)
}

# slimkey is my own intervention so far to make the color key more compact
my.heatmap.2 <- function (x,Rowv=TRUE,Colv=if (symm) "Rowv" else TRUE,distfun=dist,hclustfun=hclust,
                          dendrogram=c("both","row","column","none"),symm=FALSE,scale=c("none","row","column"),
                          na.rm=TRUE,revC=identical(Colv,"Rowv"),add.expr,breaks,col="heat.colors",
                          symbreaks=min(x < 0,na.rm=TRUE) || scale != "none",colsep,rowsep,sepcolor="white",
                          sepwidth=c(0.05,0.05),cellnote,notecex=1,notecol="cyan",na.color=par("bg"),
                          trace=c("column","row","both","none"),tracecol="cyan",hline=median(breaks),
                          vline=median(breaks),linecol=tracecol,margins=c(5,5),ColSideColors,RowSideColors,
                          cexRow=0.2 + 1/log10(nr),cexCol=0.2 + 1/log10(nc),labRow=NULL,labCol=NULL,key=TRUE,
                          keysize=1.5,density.info=c("histogram","density","none"),denscol=tracecol,
                          symkey=min(x < 0,na.rm=TRUE) || symbreaks,densadj=0.25,main=NULL,xlab=NULL,ylab=NULL,
                          lmat=NULL,lhei=NULL,lwid=NULL,slimkey=FALSE,...) 
{
    if (!require(gplots))
        stop("R library gplots is required!")
    
    scale01 <- function(x,low=min(x),high=max(x))
    {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col,mode="function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.","Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("",ncol=ncol(x),nrow=nrow(x))
    if (!inherits(Rowv,"dendrogram"))
    {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both","row")))
        {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE,while dendrogram is `",dendrogram,"'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv,"dendrogram"))
    {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both","column")))
        {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE,while dendrogram is `",dendrogram,"'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv,"dendrogram"))
    {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv))
    {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr,Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv))
    {
        Rowv <- rowMeans(x,na.rm=na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr,Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else
    {
        rowInd <- nr:1
    }
    if (inherits(Colv,"dendrogram"))
    {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv,"Rowv"))
    {
        if (nr != nc) 
            stop("Colv=\"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr"))
        {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv))
    {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc,Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv))
    {
        Colv <- colMeans(x,na.rm=na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc,Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else
    {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd,colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd,colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row")
    {
        retval$rowMeans <- rm <- rowMeans(x,na.rm=na.rm)
        x <- sweep(x,1,rm)
        retval$rowSDs <- sx <- apply(x,1,sd,na.rm=na.rm)
        x <- sweep(x,1,sx,"/")
    }
    else if (scale == "column")
    {
        retval$colMeans <- rm <- colMeans(x,na.rm=na.rm)
        x <- sweep(x,2,rm)
        retval$colSDs <- sx <- apply(x,2,sd,na.rm=na.rm)
        x <- sweep(x,2,sx,"/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1)
    {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1)
    {
        if (!symbreaks) 
            breaks <- seq(min(x,na.rm=na.rm),max(x,na.rm=na.rm),length=breaks)
        else
        {
        extreme <- max(abs(x),na.rm=TRUE)
        breaks <- seq(-extreme,extreme,length=breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize,4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize,4)
    if (missing(lmat) || is.null(lmat))
    {
        lmat <- rbind(4:3,2:1)
        if (!missing(ColSideColors))
        {
            if (!is.character(ColSideColors) || length(ColSideColors) != nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1,] + 1,c(NA,1),lmat[2,] + 1)
            lhei <- c(lhei[1],0.2,lhei[2])
        }
        if (!missing(RowSideColors))
        {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[,1] + 1,c(rep(NA,nrow(lmat) - 
                1),1),lmat[,2] + 1)
            lwid <- c(lwid[1],0.2,lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length=nrow(lmat)=",nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length=ncol(lmat) =",ncol(lmat))
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    layout(lmat,widths=lwid,heights=lhei,respect=FALSE)
    if (!missing(RowSideColors))
    {
        par(mar=c(margins[1],0,0,0.5))
        image(rbind(1:nr),col=RowSideColors[rowInd],axes=FALSE)
    }
    if (!missing(ColSideColors))
    {
        par(mar=c(0.5,0,0,margins[2]))
        image(cbind(1:nc),col=ColSideColors[colInd],axes=FALSE)
    }
    par(mar=c(margins[1],0,0,margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC)
    {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[,iy]
        cellnote <- cellnote[,iy]
    }
    else iy <- 1:nr
    image(1:nc,1:nr,x,xlim=0.5 + c(0,nc),ylim=0.5 + 
        c(0,nr),axes=FALSE,xlab="",ylab="",col=col,
        breaks=breaks,...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x)))
    {
        mmat <- ifelse(is.na(x),1,NA)
        image(1:nc,1:nr,mmat,axes=FALSE,xlab="",ylab="",
            col=na.color,add=TRUE)
    }
    axis(1,1:nc,labels=labCol,las=2,line=-0.5,tick=0,cex.axis=cexCol)
    if (!is.null(xlab)) 
        mtext(xlab,side=1,line=margins[1] - 1.25)
    axis(4,iy,labels=labRow,las=2,line=-0.5,tick=0,
        cex.axis=cexRow)
    if (!is.null(ylab)) 
        mtext(ylab,side=4,line=margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft=csep + 0.5,ybottom=rep(0,length(csep)),
            xright=csep + 0.5 + sepwidth[1],
            ytop=rep(ncol(x) + 1,csep),lty=1,lwd=1,
            col=sepcolor,border=sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft=0,ybottom=(ncol(x) + 1 - rsep) - 0.5,
            xright=nrow(x) + 1,ytop=(ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],lty=1,lwd=1,
            col=sepcolor,border=sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x),min.scale,max.scale)
    if (trace %in% c("both","column"))
    {
        retval$vline <- vline
        vline.vals <- scale01(vline,min.scale,max.scale)
        for (i in colInd) {
            if (!is.null(vline))
            {
                abline(v=i - 0.5 + vline.vals,col=linecol,lty=2)
            }
            xv <- rep(i,nrow(x.scaled)) + x.scaled[,i] - 0.5
            xv <- c(xv[1],xv)
            yv <- 1:length(xv) - 0.5
            lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
        }
    }
    if (trace %in% c("both","row"))
    {
        retval$hline <- hline
        hline.vals <- scale01(hline,min.scale,max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h=i + hline,col=linecol,lty=2)
            }
            yv <- rep(i,ncol(x.scaled)) + x.scaled[i,] - 0.5
            yv <- rev(c(yv[1],yv))
            xv <- length(yv):1 - 0.5
            lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
        }
    }
    if (!missing(cellnote)) 
        text(x=c(row(cellnote)),y=c(col(cellnote)),labels=c(cellnote),
            col=notecol,cex=notecex)
    par(mar=c(margins[1],0,0,0))
    if (dendrogram %in% c("both","row"))
    {
        plot(ddr,horiz=TRUE,axes=FALSE,yaxs="i",leaflab="none")
    }
    else plot.new()
    par(mar=c(0,0,if (!is.null(main)) 5 else 0,margins[2]))
    if (dendrogram %in% c("both","column"))
    {
        plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main,cex.main=1.5 * op[["cex.main"]])
    if (key)
    {
        par(mar=c(5,4,2,1),cex=0.75)
        tmpbreaks <- breaks
        if (symkey)
        {
            max.raw <- max(abs(c(x,breaks)),na.rm=TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x),na.rm=TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x),na.rm=TRUE)
        }
        else
        {
            min.raw <- min(x,na.rm=TRUE)
            max.raw <- max(x,na.rm=TRUE)
        }

        if (slimkey)
        {
            par(cex.main=1.1)
            cex.m <- 0.9
        } else cex.m <- 1
        
        z <- seq(min.raw,max.raw,length=length(col))
        image(z=matrix(z,ncol=1),col=col,breaks=tmpbreaks,xaxt="n",yaxt="n")
        par(usr=c(0,1,0,1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv),min.raw,max.raw)
        axis(1,at=xv,labels=lv)
        if (scale == "row") 
            mtext(side=1,"Row Z-Score",line=2,cex=cex.m)
        else if (scale == "column") 
            mtext(side=1,"Column Z-Score",line=2,cex=cex.m)
        else mtext(side=1,"Value",line=2,cex=cex.m)
        if (density.info == "density")
        {
            dens <- density(x,adjust=densadj,na.rm=TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x,min.raw,max.raw)
            lines(dens$x,dens$y/max(dens$y) * 0.95,col=denscol,lwd=1)
            axis(2,at=pretty(dens$y)/max(dens$y) * 0.95,pretty(dens$y))
            if (slimkey)
                title("Color Key")
            else
                title("Color Key\nand Density Plot")
            par(cex=0.5)
            mtext(side=2,"Density",line=2,cex=cex.m)
        }
        else if (density.info == "histogram")
        {
            h <- hist(x,plot=FALSE,breaks=breaks)
            hx <- scale01(breaks,min.raw,max.raw)
            hy <- c(h$counts,h$counts[length(h$counts)])
            lines(hx,hy/max(hy) * 0.95,lwd=1,type="s",col=denscol)
            axis(2,at=pretty(hy)/max(hy) * 0.95,pretty(hy))
            if (slimkey)
                title("Color Key")
            else
                title("Color Key\nand Histogram")
            par(cex=0.5)
            mtext(side=2,"Count",line=2,cex=cex.m)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low=retval$breaks[-length(retval$breaks)],
        high=retval$breaks[-1],color=retval$col)
    invisible(retval)
}

findOptGrid <- function(n)
{
    m <- 0
    while (n > m*m)
        m <- m+1
    if (n < m*m)
    {
        k <- m-1
        if (n > m*k)
            k <- k+1
        else
        {
            while (n > m*k)
                k=k-1
        }
    }
    else
        k <- m

    return(c(m,k))
}

openGraphics <- function(o,f,...)
{
    if(!checkGraphicsType(o))
        stop("Invalid graphics output type!")
    if(checkGraphicsFile(o) && is.null(f))
        stop("Please specify an output file name for your plot")
    
    switch(o,
        x11 = { x11() },
        pdf = { pdf(file=f,pointsize=10,...) },
        ps = { postscript(file=f,pointsize=10,...) },
        png = { png(file=f,pointsize=12,...) },
        jpg = { jpeg(file=f,pointsize=12,quality=100,...) },
        bmp = { bmp(file=f,pointsize=12,...) },
        tiff = { tiff(file=f,pointsize=12,...) })
}

closeGraphics <- function(o)
{
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    if (o!="x11")
        dev.off()
}

checkGraphicsType <- function(o)
{
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    else
        return(TRUE)
}

checkGraphicsFile <- function(o)
{
    if (is.element(o,c("png","jpg","tiff","bmp","pdf","ps")))
        return(TRUE)
    else
        return(FALSE)
}

log2disp <- function(mat,base=2)
{
    mat[mat==0] <- 1
    if (base==10)
        return(log10(mat))
    else
        return(log2(mat))
}

textwrap <- function(x,len=15)
{
    n <- nchar(x)
    nlines <- ceiling(n/len)
    splitter <- seq(1,n,by=len)
    if (tail(splitter,1) < n)
        splitter <- c(splitter,n)
    wrapped <- character(nlines)
    for (i in 1:(length(splitter)-2))
        wrapped[i] <- substr(x,splitter[i],splitter[i+1]-1)
    wrapped[length(splitter)-1] <- substr(x,splitter[i+1],splitter[i+2])
    return(paste(wrapped,collapse="\n"))
}

nat2log <- function(x,base=2)
{
    x[x==0] <- 1
    if (base==2)
        return(log2(x))
    else
        return(log10(x))
}
