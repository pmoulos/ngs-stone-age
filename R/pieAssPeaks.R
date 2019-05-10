# Nice pie charts from HyperAssignPeaks...

pieAssPeaks <- function(gff.file,main.title,...)
{
    if (missing(gff.file))
    {
        cat("Select HyperAssignPeaks gff-peak file.\n\n")
        flush.console()
        if (require(tcltk))
        {
            gff.file<-tclvalue(tkgetOpenFile())
            if (!nchar(macs.file))
            {
            	tkmessageBox(message="No file was selected!")
                return()
            }
        } else
            stop("You must provide a file for the script to work!")
    }
    if (missing(main.title))
    {
        v<-strsplit(gff.file,"/")
        main.title<-paste("Distance from peak to TSS for",v[[1]][length(v[[1]])])
    }
    
    # Read gff file from HyperAssignPeaks
    gff.tab<-read.delim(gff.file,header=FALSE)
    dists<-as.numeric(gff.tab[,13])
    edges<-c(-100000,-50000,-20000,-10000,-5000,-2000,-1000,0,1000,2000,5000,10000,20000,50000,100000)
    labs<-c('-100k to -50k','-50k to -20k','-20k to -10k','-10k to -5k','-5k to -2k',
            '-2k to -1k','-1k to TSS','TSS to 1k','1k to 2k','2k to 5k','5k to 10k',
            '10k to 20k','20k to 50k','50k to 100k');
    z<-hist(dists,edges)
    cou<-z$counts
    names(cou)<-labs
    pie(cou,...)
    title(main=main.title)
}