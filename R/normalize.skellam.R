# TODO: Perform some "logical" controls, for example the down.to variable should not be "10" but should be at least some fraction of
#       the original inputs
# TODO: Support for two values of f, one for treatment, one for control
# TODO: Write documentation...
# Apparently, just Poisson is not good enough when only treatment available. To
# make Poisson work, we obviously have to completely exclude enriched regions

normalize.skellam <- function(treatment,control=NA,down.to=NA,bin.size=1000,frag.size=300,chrom.info.file=NULL,f=0.9,pval=1e-6,
                              rr.score="uniform",c=1,write.output=TRUE,org="hg19",output=NA,n.core=6,
                              eff.size=if (is.element(org,c("hg18","hg19","mm8","mm9","mm10","dm2","dm3","ce"))) org else NULL)
{
    # If user does not provide files
    if (missing(treatment))
        stop("You must provide a treatment file for the script to work!")
    
    # Get also destination directories
    path <- dirname(treatment)
    
    if (is.na(output))
        output <- "auto"
    else
    {
        if (!is.na(down.to) && !is.na(control) && output != "auto" && length(output) != 2)
        {
            warning("The output filenames must be a vector of two when treatment and control are given in combination with downscaling! Switching to auto...",
                    .call=FALSE)
            output <- "auto"
        }
    }
    if (!require(GenomicRanges))
        stop("Bioconductor package IRanges is required!")
    if (!require(rtracklayer))
        stop("Bioconductor package rtracklayer is required!")
    if (!require(skellam))
        stop("R package skellam is required!")
    if (f<=0 || f>1)
        stop("f must be between 0 and 1!")
    if (pval<=0 || pval>1)
        stop("pval must be between 0 and 1!")
    if (2*frag.size > bin.size)
        stop("bin.size must be at least twice the fragment size")
    if (!is.na(down.to) && (down.to <= 0 || !is.numeric(down.to)))
        stop("down.to must be a positive integer!")
    rr.score <- tolower(rr.score)
    if (!(rr.score %in% c("uniform","linear","exponential")))
        stop("rr.score must be one of \"uniform\", \"linear\" or \"exponential\"")
    if (is.na(control) && is.na(down.to))
        stop("Either a control file or a number of reads to downscale to must be provided!")
    
    if (n.core > 1)
    {
        if (!require(multicore))
        {
            warning("R package multicore not present... Switching to single core...",call.=FALSE)
            n.core <- 1
        }
        else
        {
            require(multicore)
            cores <- multicore:::detectCores()
            if (cores==1)
            {
                warning("Only one core detected... Switching to single core...",call.=FALSE)
                n.core <- 1
            }
        }
    }
    
    cat("\n")
    cat("Gathering chromosome information...\n")
    if (!is.null(chrom.info.file))
    {
        chrom.info <- read.delim(chrom.info.file,header=FALSE)
        if (length(grep("rand|chrM|hap",chrom.info[,1])) > 0)
            chrom.info <- chrom.info[-grep("rand|chrM|hap",chrom.info[,1]),]
        genome.size <- sum(as.numeric(chrom.info[,2]))
        attr(chrom.info,"genome.size") <- genome.size
        rownames(chrom.info) <- chrom.info[,1]
    }
    else
        chrom.info <- get.chrom.info(org)
    
    cat("Creating genomic bins of ",bin.size,"bp...\n",sep="")
    gen.win <- make.windows.list(chrlen=as.numeric(chrom.info[,2]),winsize=bin.size,chrname=as.character(chrom.info[,1]))
    #gen.win <- make.sliding(x=as.numeric(chrom.info[,2]),win=bin.size-1,slide=2*frag.size,"chr18","GRanges") # Test
    
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
    
    # Start doing the job
    cat("Importing track ",basename(treatment),"... Please wait...\n",sep="")
    flush.console()
    treat.bed <- import.bed(treatment,trackLine=FALSE)
    treat.depth <- length(treat.bed)

    if (!is.na(control))
    {
        cat("Importing track ",basename(control),"... Please wait...\n",sep="")
        flush.console()
        control.bed <- import.bed(control,trackLine=FALSE)
        control.depth <- length(control.bed)
    }
    
    if (!is.na(down.to)) # We have to define a number of reads to remove from each chromosome
    {
        cat("Defining number of reads to remove from each chromosome in treatment...\n")
        d.t <- treat.depth - down.to
        if (d.t < 0)
            stop(paste("The number of reads (",down.to,") to downscale the treatment to, must be smaller than the number of reads in the treatment file (",treat.depth,")!",sep=""))
        else
        {
            treat.pop <- numeric(length(seqlevels(treat.bed)))
            names(treat.pop) <- seqlevels(treat.bed)
            for (n in seqlevels(treat.bed))
                treat.pop[n] <- length(treat.bed[seqnames(treat.bed)==n])
            n.remove.t <- round(d.t*treat.pop/treat.depth)
            #n.remove.t <- round(d.t*chrom.info[,2]/attr(chrom.info,"genome.size"))
            names(n.remove.t) <- seqlevels(treat.bed)
        }
        if (!is.na(control))
        {
            cat("Defining number of reads to remove from each chromosome in control...\n")
            d.c <- control.depth - down.to
            if (d.c < 0)
                stop(paste("The number of reads (",down.to,") to downscale the control to, must be smaller than the number of reads in the control file (",control.depth,")!",sep=""))
            else
            {
                control.pop <- numeric(length(seqlevels(control.bed)))
                names(control.pop) <- seqlevels(control.bed)
                for (n in seqlevels(control.bed))
                    control.pop[n] <- length(control.bed[seqnames(control.bed)==n])
                n.remove.c <- round(d.c*control.pop/control.depth)
                #n.remove.c <- round(d.c*chrom.info[,2]/attr(chrom.info,"genome.size"))
                names(n.remove.c) <- seqlevels(control.bed)
            }
        }
        
        cat("\nDistribution of existing treatment reads per chromosome:\n")
        print(treat.pop)
        cat("Distribution of treatment reads to remove per chromosome:\n")
        print(n.remove.t)
        if (!is.na(control))
        {
            cat("\nDistribution of existing control reads per chromosome:\n")
            print(control.pop)
            cat("Distribution of control reads to remove per chromosome:\n")
            print(n.remove.c)
        }
        cat("\n")
    }
    
    bin.counts <- bin.index <- alive <- xe <- ye <- xne <- yne <- remove.index <- list()
    gen.win.e <- gen.win.ne <- GRangesList()
    
    cat("Counting tags in the treatment track ",basename(treatment)," over genomic windows of ",bin.size,"bp...\n",sep="")
    if (n.core > 1)
        bin.counts$treatment <- mclapply(as.list(gen.win),countOverlaps,treat.bed,mc.cores=n.core)
    else
    {
        for (chr in rownames(chrom.info))
        {
            cat("  ",chr,"...\n",sep="")
            bin.counts$treatment[[chr]] <- countOverlaps(gen.win[[chr]],treat.bed[seqnames(treat.bed)==chr])
        }
    }
    
    if (!is.na(control))
    {
        cat("Counting tags in the control track ",basename(control)," over genomic windows of ",bin.size,"bp...\n",sep="")
        if (n.core > 1)
            bin.counts$control <- mclapply(as.list(gen.win),countOverlaps,control.bed,mc.cores=n.core)
        else
        {
            for (chr in rownames(chrom.info))
            {
                cat("  ",chr,"...\n",sep="")
                bin.counts$control[[chr]] <- countOverlaps(gen.win[[chr]],control.bed[seqnames(control.bed)==chr])
            }
        }
    }
    
    gc(verbose=FALSE)
    
    # 1. Remove dead zones for treatment (and control)
    # 2. Calculate local lambdas MACS-like for treatment (and control)
    # 3. Calculate Skellam null distribution and expected reads in each window
    # 4. Separate the windows in enriched and non-enriched ones by
    #    a. Skellam p-value
    #    b. For the ones that pass the Skellam p-value, log2(T/C)<1 to avoid garbage like repeats
    # 5. In order to avoid peak bleeding
    #    a. Expand the enriched regions one fragment length upstream and downstream
    #    b. Shrink the non-enriched regions one fragment length upstream and downstream
    # 6. Merge (IRanges:::reduce) the resulting windows
    # 7. Using the new regions, find the reads that overlap with the new enriched and non-enriched windows
    # 8. Remove reads from non-enriched windows with probability p and from enriched with probability 1-p (notes for more complex)
    
    if (is.na(down.to)) # Launch first main function
    {
        
        total.remove <- 0
        n.remove <- n.remove.e <- n.remove.ne <- n.fall.e <- n.fall.ne <- l.alive <- numeric(nrow(chrom.info))
        names(n.remove) <- names(n.remove.e) <- names(n.remove.ne) <- names(n.fall.e) <- names(n.fall.ne) <- names(l.alive) <- rownames(chrom.info)
        
        cat("Constructing read index...\n")
        if (treat.depth > control.depth)
            read.index <- make.read.index(treat.bed)
        else
            read.index <- make.read.index(control.bed)
    
    #     if (n.core > 1)
    #     {
    #         cat("Parallely processing chromosomes...\n")
    #         cat("  Calculating enriched regions...\n")
    #         if (!is.na(control))
    #         {
    #             alive.c <- mclapply(bin.counts$control,function(x) which(x != 0))
    #             alive.t <- mclapply(bin.counts$treatment,function(x) which(x != 0))
    #             for (chr in rownames(chrom.info))
    #             {
    #                 alive[[chr]] <- union(alive.c[[chr]],alive.t[[chr]])
    #             }
    #         }
    #     }
    #     else
    #     {
            for (chr in rownames(chrom.info))
            {
                cat("Processing ",chr,"...\n",sep="")
                cat("  Calculating enriched regions...\n")
                if (!is.na(control))
                {
                    alive <- which(bin.counts$control[[chr]] != 0 | bin.counts$treatment[[chr]] != 0)
                    l.alive[chr] <- length(alive)
                    x <- bin.counts$control[[chr]][alive]
                    y <- bin.counts$treatment[[chr]][alive]
                    #lamda1 <- (length(treat.bed[seqnames(treat.bed)==chr])*(y/bin.size))/(eff.gen.size*bin.size/attr(chrom.info,"genome.size"))
                    #lamda2 <- (length(control.bed[seqnames(control.bed)==chr])*(x/bin.size))/(eff.gen.size*bin.size/attr(chrom.info,"genome.size"))
                    lamda1 <- (length(treat.bed[seqnames(treat.bed)==chr])*(y/bin.size))/(eff.gen.size*bin.size/chrom.info[chr,2])
                    lamda2 <- (length(control.bed[seqnames(control.bed)==chr])*(x/bin.size))/(eff.gen.size*bin.size/chrom.info[chr,2])
                    if (treat.depth > control.depth)
                    {
                        n.remove[chr] <- length(treat.bed[seqnames(treat.bed)==chr]) - length(control.bed[seqnames(control.bed)==chr])
                        if (n.remove[chr] < 0) n.remove[chr] <- 0
                        exp.tag <- qskellam(1-pval,lambda1=lamda1,lambda2=lamda2)
                        ie <- which(y - x > exp.tag & log2(ifelse(y==0,1,y)/ifelse(x==0,1,x)) > 1)
                    }
                    else # See notes why I am doing this... Enrichment fold does not apply here...
                    {
                        n.remove[chr] <- length(control.bed[seqnames(control.bed)==chr]) - length(treat.bed[seqnames(treat.bed)==chr])
                        if (n.remove[chr] < 0) n.remove[chr] <- 0
                        exp.tag <- qskellam(1-pval,lambda1=lamda2,lambda2=lamda1)
                        ie <- which(x - y > exp.tag)
                    }
                    xe[[chr]] <- x[ie]
                    ye[[chr]] <- y[ie]
                }
                
                gen.win.e[[chr]] <- gen.win[[chr]][alive[ie]]
                
                # Expanding regions is very sensitive as we might cross chromosome limits
                # First we check the start
                cat("  Adjusting bounds for enrichment bleeding...\n")
                if ( start(gen.win.e[[chr]][1]) + (end(gen.win.e[[chr]][1]) - start(gen.win.e[[chr]][1]))/2 - (bin.size+2*frag.size)/2 < 1 )
                {
                    # The first position
                    ranges(gen.win.e[[chr]][1]) <- resize(ranges(gen.win.e[[chr]][1]),bin.size+frag.size,fix="start")
                    # The second position
                    if ( start(gen.win.e[[chr]][2]) + (end(gen.win.e[[chr]][2]) - start(gen.win.e[[chr]][2]))/2 - (bin.size+2*frag.size)/2 < 1 )
                        ranges(gen.win.e[[chr]][2]) <- resize(ranges(gen.win.e[[chr]][2]),bin.size+frag.size,fix="start")
                    else
                        ranges(gen.win.e[[chr]][2]) <- resize(ranges(gen.win.e[[chr]][2]),bin.size+2*frag.size,fix="center")
                }
                else
                    ranges(gen.win.e[[chr]][1:2]) <- resize(ranges(gen.win.e[[chr]][1:2]),bin.size+2*frag.size,fix="center")
                
                # Then the rest
                e <- length(gen.win.e[[chr]])
                if ( start(gen.win.e[[chr]][e]) + (end(gen.win.e[[chr]][e]) - start(gen.win.e[[chr]][e]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
                {
                    # Up to e-2 should be safe
                    ranges(gen.win.e[[chr]][3:(e-2)]) <- resize(ranges(gen.win.e[[chr]][3:(e-2)]),bin.size + 2*frag.size,fix="center")
                    # Check up to e-1
                    if ( start(gen.win.e[[chr]][e-1]) + (end(gen.win.e[[chr]][e-1]) - start(gen.win.e[[chr]][e-1]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
                        ranges(gen.win.e[[chr]][e-1]) <- resize(ranges(gen.win.e[[chr]][e-1]),bin.size+frag.size,fix="end")
                    else
                        ranges(gen.win.e[[chr]][e-1]) <- resize(ranges(gen.win.e[[chr]][e-1]),bin.size+2*frag.size,fix="center")
                }
                else
                    ranges(gen.win.e[[chr]][3:e]) <- resize(ranges(gen.win.e[[chr]][3:e]),bin.size+2*frag.size,fix="center")
                
                # The complement of the final enriched regions should comprise the non-enriched regions
                gen.win.ne[[chr]] <- gaps(gen.win.e[[chr]])
                gen.win.ne[[chr]] <- gen.win.ne[[chr]][-c(1:2)]
                                  
                # Merge the expanded intervals
                cat("  Merging adjusted regions...\n")
                gen.win.e[[chr]] <- reduce(gen.win.e[[chr]])
                
                # ...and recalculate the overlaps, required for the probability step
                if (rr.score %in% c("linear","exponential"))
                {
                    cat("  Recounting tags...\n")
                    ye[[chr]] <- countOverlaps(gen.win.e[[chr]],treat.bed)
                    if (!is.na(control))
                        xe[[chr]] <- countOverlaps(gen.win.e[[chr]],control.bed)
                }
                
                cat("  Indexing treatment reads over genomic windows...\n")
                bin.index$treatment$e[[chr]] <- as.list(findOverlaps(gen.win.e[[chr]],treat.bed[seqnames(treat.bed)==chr]))
                bin.index$treatment$ne[[chr]] <- as.list(findOverlaps(gen.win.ne[[chr]],treat.bed[seqnames(treat.bed)==chr]))
                if (!is.na(control))
                {
                    cat("  Indexing control reads over genomic windows...\n")
                    bin.index$control$e[[chr]] <- as.list(findOverlaps(gen.win.e[[chr]],control.bed[seqnames(control.bed)==chr]))
                    bin.index$control$ne[[chr]] <- as.list(findOverlaps(gen.win.ne[[chr]],control.bed[seqnames(control.bed)==chr]))
                }
                
                cat("  Sampling reads to remove... Selected method: ",rr.score,"\n",sep="")
                
                # The next 2 blocks control what may happen if a region is assigned with more reads to remove than its actual 
                # read content
                if (treat.depth > control.depth)
                {
                    n.fall.e[chr] <- length(unlist(bin.index$treatment$e[[chr]]))
                    n.fall.ne[chr] <- length(unlist(bin.index$treatment$ne[[chr]]))
                }
                else
                {
                    n.fall.e[chr] <- length(unlist(bin.index$control$e[[chr]]))
                    n.fall.ne[chr] <- length(unlist(bin.index$control$ne[[chr]]))
                }
                n.remove.e[chr] <- floor((1-f)*n.remove[chr])
                n.remove.ne[chr] <- ceiling(f*n.remove[chr])
                
                ne.deplete <- FALSE
                if (n.remove.ne[chr] > n.fall.ne[chr]) # Must resize appropriately, removing more from enriched regions...
                {
                    d <- n.remove.ne[chr] - n.fall.ne[chr]
                    n.remove.ne[chr] <- n.fall.ne[chr]
                    n.remove.e[chr] <- n.remove.e[chr] + d
                    ne.deplete <- TRUE
                }
                if (n.remove.e[chr] > n.fall.e[chr]) # This should never happen! Here, only to avoid crash
                {
                    d <- n.remove.e[chr] - n.fall.e[chr]
                    n.remove.e[chr] <- n.fall.e[chr]
                    n.remove.ne[chr] <- n.remove.ne[chr] + d
                }
                
                #cat("    ",n.fall.e[chr]," reads fall into enriched regions\n",sep="")
                #cat("    ",n.fall.ne[chr]," reads fall into non-enriched regions\n",sep="")
                #cat("    ",n.remove.e[chr]," reads will be removed from enriched regions\n",sep="")
                #cat("    ",n.remove.ne[chr]," reads will be removed from non-enriched regions\n",sep="")
                
                if (treat.depth > control.depth)
                {
                    if (rr.score == "exponential")
                    {
                        r <- log2(ifelse(ye[[chr]]==0,1,ye[[chr]])/ifelse(xe[[chr]]==0,1,xe[[chr]]))
                        p <- (1/exp(1-c*r))/((1/(1-f))*max(1/exp(1-c*r)))
                        wei <- rep(p,sapply(bin.index$treatment$e[[chr]],function(x) length(x)))
                        i.remove=c(
                            sample(unlist(bin.index$treatment$e[[chr]]),n.remove.e[chr],prob=wei),
                            sample(unlist(bin.index$control$ne[[chr]]),n.remove.ne[chr])
                        )
                    }
                    else if(rr.score == "linear")
                    {
                        r <- log2(ifelse(ye[[chr]]==0,1,ye[[chr]])/ifelse(xe[[chr]]==0,1,xe[[chr]]))
                        p <- r/((1/(1-f))*max(r))
                        p[p<0] <- 0
                        wei <- rep(p,sapply(bin.index$treatment$e[[chr]],function(x) length(x)))
                        i.remove=c(
                            sample(unlist(bin.index$treatment$e[[chr]]),n.remove.e[chr],prob=wei),
                            sample(unlist(bin.index$control$ne[[chr]]),n.remove.ne[chr])
                        )
                    }
                    else if (rr.score == "uniform")
                        i.remove <- c(
                            sample(unlist(bin.index$treatment$ne[[chr]]),n.remove.ne[chr]),
                            sample(unlist(bin.index$treatment$e[[chr]]),n.remove.e[chr])
                        )
                }
                else # uniform removal should be OK in this case
                {
                    i.remove <- c(
                        sample(unlist(bin.index$control$ne[[chr]]),n.remove.ne[chr]),
                        sample(unlist(bin.index$control$e[[chr]]),n.remove.e[chr])
                    )
                }
                
                ul <- length(unique(i.remove))
                cc <- 1
                
                if (ne.deplete) # Sample a bit more from enriched (does not happen often)
                {
                    while(ul != n.remove[chr]) 
                    {
                        cat("    Iteration ",cc,"...\n",sep="")
                        d <- n.remove[chr] - ul
                        if (treat.depth > control.depth)
                        {
                            sample.space <- setdiff(unlist(bin.index$treatment$e[[chr]]),i.remove)
                            cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                            i.more <- sample(sample.space,ceiling(f*d))
                        }
                        else
                        {
                            sample.space <- setdiff(unlist(bin.index$control$e[[chr]]),i.remove)
                            cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                            i.more <- sample(sample.space,ceiling(f*d))
                        }
                        i.remove <- c(i.remove,i.more)
                        ul <- length(unique(i.remove))
                        cc <- cc + 1
                    }
                }
                else # Sample a bit more from non-enriched
                {
                    while(ul != n.remove[chr]) 
                    {
                        cat("    Iteration ",cc,"...\n",sep="")
                        d <- n.remove[chr] - ul
                        if (treat.depth > control.depth)
                        {
                            sample.space <- setdiff(unlist(bin.index$treatment$ne[[chr]]),i.remove)
                            cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                            i.more <- sample(sample.space,ceiling(f*d))
                        }
                        else
                        {
                            sample.space <- setdiff(unlist(bin.index$control$ne[[chr]]),i.remove)
                            cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                            i.more <- sample(sample.space,ceiling(f*d))
                        }
                        i.remove <- c(i.remove,i.more)
                        ul <- length(unique(i.remove))
                        cc <- cc + 1
                    }
                }
                
                remove.index[[chr]] <- read.index[[chr]][i.remove]
                total.remove <- total.remove + n.remove[chr]
            } # for chr ends here
        #}
        
        cat("Removing actual reads...\n")
        if (treat.depth > control.depth)
            treat.bed <- treat.bed[-unlist(remove.index)]
        else
            control.bed <- control.bed[-unlist(remove.index)]
        
        if (write.output)
        {
            if (output=="auto")
            {
                if (!is.na(control))
                {                
                    if (treat.depth > control.depth)
                    {
                        name <- sub("(.+)[.][^.]+$", "\\1", basename(treatment))
                        o <- file.path(dirname(treatment),paste(name,"norm.bed",sep="."))
                        cat("Writing output to ",o,"...\n",sep="")
                        export(treat.bed,o)
                    }
                    else
                    {
                        name <- sub("(.+)[.][^.]+$", "\\1", basename(control))
                        o <- file.path(dirname(treatment),paste(name,"norm.bed",sep="."))
                        cat("Writing output to ",o,"...\n",sep="")
                        export(control.bed,o)
                    }
                }
                else
                {
                    name <- sub("(.+)[.][^.]+$", "\\1", basename(treatment))
                    o <- file.path(dirname(treatment),paste(name,"norm.bed",sep="."))
                    cat("Writing output to ",o,"...\n",sep="")
                    export(treat.bed,o)
                }
            }
            else
            {
                if (!is.na(control))
                {
                    if (treat.depth > control.depth)
                    {
                        cat("Writing output to ",output,"...\n",sep="")
                        export(treat.bed,output)
                    }
                    else
                    {
                        cat("Writing output to ",output,"...\n",sep="")
                        export(control.bed,output)
                    }
                }
                else
                {
                    cat("Writing output to ",output,"...\n",sep="")
                    export(treat.bed,output)
                }
            }
        }
        
        # Create a nice report!
        cat("\n----- REPORT -----\n")
        cat("Treatment file ",basename(treatment)," has ",treat.depth," reads\n",sep="")
        if (!is.na(control))
            cat("Control file ",basename(control)," has ",control.depth," reads\n",sep="")
        if (treat.depth > control.depth)
        {
            cat("Treatment file ",basename(treatment)," was shrinked by ",total.remove," reads\n",sep="")
            cat("Of these:\n")
            for (chr in rownames(chrom.info))
            {
                cat("  --- ",chr,"\n",sep="")
                cat("  ",length(gen.win[[chr]]) - l.alive[chr]," dead regions were found treatment and control\n",sep="")
                cat("  ",length(gen.win.e[[chr]])," enriched regions were found in treatment over control\n",sep="")
                cat("  ",l.alive[chr] - length(gen.win.e[[chr]])," non-enriched regions were found in treatment over control\n",sep="")
                cat("  In these:\n")
                cat("    ",n.fall.e[chr]," reads fall into enriched regions\n",sep="")
                cat("    ",n.fall.ne[chr]," reads fall into non-enriched regions\n",sep="")
                cat("    Of these:\n")
                cat("      ",n.remove.e[chr]," reads were removed from enriched regions\n",sep="")
                cat("      ",n.remove.ne[chr]," reads were removed from non-enriched regions\n",sep="")
                cat("  Finally ",n.remove[chr]," reads were removed from treatment\n",sep="")
            }
        }
        else
        {
            cat("Control file ",basename(control)," was shrinked by ",total.remove," reads\n",sep="")
            cat("Of these:\n")
            for (chr in rownames(chrom.info))
            {
                cat("  --- ",chr,"\n",sep="")
                cat("  ",length(gen.win[[chr]]) - l.alive[chr]," dead regions were found treatment and control\n",sep="")
                cat("  ",length(gen.win.e[[chr]])," enriched regions were found in control over treatment\n",sep="")
                cat("  ",l.alive[chr] - length(gen.win.ne[[chr]])," non-enriched regions were found in control over treatment\n",sep="")
                cat("  In these:\n")
                cat("    ",n.fall.e[chr]," reads fall into enriched regions\n",sep="")
                cat("    ",n.fall.ne[chr]," reads fall into non-enriched regions\n",sep="")
                cat("    Of these:\n")
                cat("      ",n.remove.e[chr]," reads were removed from enriched regions\n",sep="")
                cat("      ",n.remove.ne[chr]," reads were removed from non-enriched regions\n",sep="")
                cat("  Finally ",n.remove[chr]," reads were removed from control\n",sep="")
            }
        }
        cat("\n")
        
    }
    
    else # Second main function
    
    {
        
        total.remove.t <- total.remove.c <- 0
        n.remove.et <- n.remove.net <- n.remove.ec <- n.remove.nec <-
            n.fall.et <- n.fall.net <- n.fall.ec <- n.fall.nec <- l.alive <- numeric(nrow(chrom.info))
        names(n.remove.et) <- names(n.remove.net) <- names(n.remove.ec) <- names(n.remove.nec) <- 
            names(n.fall.et) <- names(n.fall.net) <- names(n.fall.ec) <- names(n.fall.nec) <- names(l.alive) <- rownames(chrom.info)
        gen.win.et <- gen.win.net <- gen.win.ec <- gen.win.nec <- GRangesList()
        remove.index.t <- remove.index.c <- list()
        
        cat("Constructing read index for treatment...\n")
        read.index.t <- make.read.index(treat.bed)
        if (!is.na(control))
        {
            cat("Constructing read index for control...\n")
            read.index.c <- make.read.index(control.bed)
        }
        
        for (chr in rownames(chrom.info))
        {
            cat("Processing ",chr,"...\n",sep="")
            cat("  Calculating enriched regions...\n")
            if (!is.na(control))
            {
                alive <- which(bin.counts$control[[chr]] != 0 | bin.counts$treatment[[chr]] != 0)
                l.alive[chr] <- length(alive)
                x <- bin.counts$control[[chr]][alive]
                y <- bin.counts$treatment[[chr]][alive]
                if (n.remove.t[chr] < 0) n.remove.t[chr] <- 0
                if (n.remove.c[chr] < 0) n.remove.c[chr] <- 0
                lamda1 <- (length(treat.bed[seqnames(treat.bed)==chr])*(y/bin.size))/(eff.gen.size*bin.size/chrom.info[chr,2])
                lamda2 <- (length(control.bed[seqnames(control.bed)==chr])*(x/bin.size))/(eff.gen.size*bin.size/chrom.info[chr,2])
                exp.tag.t <- qskellam(1-pval,lambda1=lamda1,lambda2=lamda2)
                iet <- which(y - x > exp.tag.t & log2(ifelse(y==0,1,y)/ifelse(x==0,1,x)) > 1)
                exp.tag.c <- qskellam(1-pval,lambda1=lamda2,lambda2=lamda1)
                iec <- which(x - y > exp.tag.c)
                xe[[chr]] <- x[iec]
                ye[[chr]] <- y[iet]
                gen.win.et[[chr]] <- gen.win[[chr]][alive[iet]]
                gen.win.ec[[chr]] <- gen.win[[chr]][alive[iec]]
            }
            else # HERE!!!
            {
                alive <- which(bin.counts$treatment[[chr]] != 0)
                l.alive[chr] <- length(alive)
                y <- bin.counts$treatment[[chr]][alive]
                lamda <- (length(treat.bed[seqnames(treat.bed)==chr])*(y/bin.size))/(eff.gen.size*bin.size/attr(chrom.info,"genome.size"))
                ###***
                exp.tag <- qpois(1-pval,lambda=lamda)
                ###***
                ie <- which(y > exp.tag)
                ye[[chr]] <- y[ie]
                gen.win.et[[chr]] <- gen.win[[chr]][alive[ie]]
            }
            
            # First we check the start of treatment
            cat("  Adjusting treatment bounds for enrichment bleeding...\n")
            if ( start(gen.win.et[[chr]][1]) + (end(gen.win.et[[chr]][1]) - start(gen.win.et[[chr]][1]))/2 - (bin.size+2*frag.size)/2 < 1 )
            {
                cat("    Done first\n")
                # The first position
                ranges(gen.win.et[[chr]][1]) <- resize(ranges(gen.win.et[[chr]][1]),bin.size+frag.size,fix="start")
                cat("    Done second\n")
                # The second position
                if ( start(gen.win.et[[chr]][2]) + (end(gen.win.et[[chr]][2]) - start(gen.win.et[[chr]][2]))/2 - (bin.size+2*frag.size)/2 < 1 )
                    ranges(gen.win.et[[chr]][2]) <- resize(ranges(gen.win.et[[chr]][2]),bin.size+frag.size,fix="start")
                else
                    ranges(gen.win.et[[chr]][2]) <- resize(ranges(gen.win.et[[chr]][2]),bin.size+2*frag.size,fix="center")
                cat("    Done third\n")
            }
            else
                ranges(gen.win.et[[chr]][1:2]) <- resize(ranges(gen.win.et[[chr]][1:2]),bin.size+2*frag.size,fix="center")
            
            # Then the rest of treatment
            e <- length(gen.win.et[[chr]])
            if ( start(gen.win.et[[chr]][e]) + (end(gen.win.et[[chr]][e]) - start(gen.win.et[[chr]][e]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
            {
                # Up to e-2 should be safe
                ranges(gen.win.et[[chr]][3:(e-2)]) <- resize(ranges(gen.win.et[[chr]][3:(e-2)]),bin.size + 2*frag.size,fix="center")
                # Check up to e-1
                if ( start(gen.win.et[[chr]][e-1]) + (end(gen.win.et[[chr]][e-1]) - start(gen.win.et[[chr]][e-1]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
                    ranges(gen.win.et[[chr]][e-1]) <- resize(ranges(gen.win.et[[chr]][e-1]),bin.size+frag.size,fix="end")
                else
                    ranges(gen.win.et[[chr]][e-1]) <- resize(ranges(gen.win.et[[chr]][e-1]),bin.size+2*frag.size,fix="center")
            }
            else
                ranges(gen.win.et[[chr]][3:e]) <- resize(ranges(gen.win.et[[chr]][3:e]),bin.size+2*frag.size,fix="center")
            
            gen.win.net[[chr]] <- gaps(gen.win.et[[chr]])
            gen.win.net[[chr]] <- gen.win.net[[chr]][-c(1:2)]
                        
            # Then we check the start of control
            if (!is.na(control))
            {    
                cat("  Adjusting control bounds for enrichment bleeding...\n")
                if ( start(gen.win.ec[[chr]][1]) + (end(gen.win.ec[[chr]][1]) - start(gen.win.ec[[chr]][1]))/2 - (bin.size+2*frag.size)/2 < 1 )
                {
                    # The first position
                    ranges(gen.win.ec[[chr]][1]) <- resize(ranges(gen.win.ec[[chr]][1]),bin.size+frag.size,fix="start")
                    # The second position
                    if ( start(gen.win.ec[[chr]][2]) + (end(gen.win.ec[[chr]][2]) - start(gen.win.ec[[chr]][2]))/2 - (bin.size+2*frag.size)/2 < 1 )
                        ranges(gen.win.ec[[chr]][2]) <- resize(ranges(gen.win.ec[[chr]][2]),bin.size+frag.size,fix="start")
                    else
                        ranges(gen.win.ec[[chr]][2]) <- resize(ranges(gen.win.ec[[chr]][2]),bin.size+2*frag.size,fix="center")
                }
                else
                    ranges(gen.win.ec[[chr]][1:2]) <- resize(ranges(gen.win.ec[[chr]][1:2]),bin.size+2*frag.size,fix="center")
                
                # Then the rest of control
                e <- length(gen.win.ec[[chr]])
                if ( start(gen.win.ec[[chr]][e]) + (end(gen.win.ec[[chr]][e]) - start(gen.win.ec[[chr]][e]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
                {
                    # Up to e-2 should be safe
                    ranges(gen.win.ec[[chr]][3:(e-2)]) <- resize(ranges(gen.win.ec[[chr]][3:(e-2)]),bin.size + 2*frag.size,fix="center")
                    # Check up to e-1
                    if ( start(gen.win.ec[[chr]][e-1]) + (end(gen.win.ec[[chr]][e-1]) - start(gen.win.ec[[chr]][e-1]))/2 + (bin.size+2*frag.size)/2 > chrom.info[chr,2] )
                        ranges(gen.win.ec[[chr]][e-1]) <- resize(ranges(gen.win.ec[[chr]][e-1]),bin.size+frag.size,fix="end")
                    else
                        ranges(gen.win.ec[[chr]][e-1]) <- resize(ranges(gen.win.ec[[chr]][e-1]),bin.size+2*frag.size,fix="center")
                }
                else
                    ranges(gen.win.ec[[chr]][3:e]) <- resize(ranges(gen.win.ec[[chr]][3:e]),bin.size+2*frag.size,fix="center")
                
                gen.win.nec[[chr]] <- gaps(gen.win.ec[[chr]])
                gen.win.nec[[chr]] <- gen.win.nec[[chr]][-c(1:2)]
            }
            
            # Merge the expanded intervals
            cat("  Merging adjusted regions...\n")
            gen.win.et[[chr]] <- reduce(gen.win.et[[chr]])
            #gen.win.net[[chr]] <- reduce(gen.win.net[[chr]]) # Unlikely but never know...
            if (!is.na(control))
            {
                gen.win.ec[[chr]] <- reduce(gen.win.ec[[chr]])
                #gen.win.nec[[chr]] <- reduce(gen.win.nec[[chr]]) # Unlikely but never know...
            }
            
            cat("  Indexing treatment reads over genomic windows...\n")
            bin.index$treatment$e[[chr]] <- as.list(findOverlaps(gen.win.et[[chr]],treat.bed[seqnames(treat.bed)==chr]))
            bin.index$treatment$ne[[chr]] <- as.list(findOverlaps(gen.win.net[[chr]],treat.bed[seqnames(treat.bed)==chr]))
            if (!is.na(control))
            {
                cat("  Indexing control reads over genomic windows...\n")
                bin.index$control$e[[chr]] <- as.list(findOverlaps(gen.win.ec[[chr]],control.bed[seqnames(control.bed)==chr]))
                bin.index$control$ne[[chr]] <- as.list(findOverlaps(gen.win.nec[[chr]],control.bed[seqnames(control.bed)==chr]))
            }
            
            n.fall.et[chr] <- length(unlist(bin.index$treatment$e[[chr]]))
            n.fall.net[chr] <- length(unlist(bin.index$treatment$ne[[chr]]))
            n.remove.et[chr] <- floor((1-f)*n.remove.t[chr])
            n.remove.net[chr] <- ceiling(f*n.remove.t[chr])
            if (!is.na(control))
            {
                n.fall.ec[chr] <- length(unlist(bin.index$control$e[[chr]]))
                n.fall.nec[chr] <- length(unlist(bin.index$control$ne[[chr]]))
                n.remove.ec[chr] <- floor((1-f)*n.remove.c[chr])
                n.remove.nec[chr] <- ceiling(f*n.remove.c[chr])
            }
            
            cat("  Sampling reads to remove from treatment... Selected method: ",rr.score,"\n",sep="")
            
            #cat("\n\n")
            #cat("    ",n.fall.et[chr]," reads fall into enriched regions\n",sep="")
            #cat("    ",n.fall.net[chr]," reads fall into non-enriched regions\n",sep="")
            #cat("    ",n.remove.et[chr]," reads will be removed from enriched regions\n",sep="")
            #cat("    ",n.remove.net[chr]," reads will be removed from non-enriched regions\n",sep="")
            #cat("\n\n")
            
            ne.deplete <- FALSE
            if (n.remove.net[chr] > n.fall.net[chr]) # Must resize appropriately, removing more from enriched regions...
            {
                d <- n.remove.net[chr] - n.fall.net[chr]
                n.remove.net[chr] <- n.fall.net[chr]
                n.remove.et[chr] <- n.remove.et[chr] + d
                ne.deplete <- TRUE
            }
            if (n.remove.et[chr] > n.fall.et[chr]) # This should never happen! Here, only to avoid crash
            {
                d <- n.remove.et[chr] - n.fall.et[chr]
                n.remove.et[chr] <- n.fall.et[chr]
                n.remove.net[chr] <- n.remove.net[chr] + d
            }
            
            cat("    ",n.fall.et[chr]," reads fall into enriched regions\n",sep="")
            cat("    ",n.fall.net[chr]," reads fall into non-enriched regions\n",sep="")
            cat("    ",n.remove.et[chr]," reads will be removed from enriched egions\n",sep="")
            cat("    ",n.remove.net[chr]," reads will be removed from non-enriched regions\n",sep="")
            cat("\n\n")
            
            # Exponential and linear score are not available because of different lengths of enriched and non-enriched regions in
            # treatment and control
            i.remove <- c(
                sample(unlist(bin.index$treatment$ne[[chr]]),n.remove.net[chr]),
                sample(unlist(bin.index$treatment$e[[chr]]),n.remove.et[chr])
            )
            
            ul <- length(unique(i.remove))
            cc <- 1
            
            if (ne.deplete) # Sample a bit more from enriched (does not happen often)
            {
                while(ul != n.remove.t[chr]) 
                {
                    cat("    Iteration ",cc,"...\n",sep="")
                    d <- n.remove.t[chr] - ul
                    sample.space <- setdiff(unlist(bin.index$treatment$e[[chr]]),i.remove)
                    if (length(sample.space)<ceiling(f*d)) { # Nothing to sample, already reached
                        cat("      Empty sample space! Goal reached.\n")
                        break
                    }
                    cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                    i.more <- sample(sample.space,ceiling(f*d))
                    i.remove <- c(i.remove,i.more)
                    ul <- length(unique(i.remove))
                    cc <- cc + 1
                }
            }
            else # Sample a bit more from non-enriched
            {
                while(ul != n.remove.t[chr])
                {
                    cat("    Iteration ",cc,"...\n",sep="")
                    d <- n.remove.t[chr] - ul
                    sample.space <- setdiff(unlist(bin.index$treatment$ne[[chr]]),i.remove)
                    if (length(sample.space)<ceiling(f*d)) { # Nothing to sample, already reached
                        cat("      Empty sample space! Goal reached.\n")
                        break
                    }
                    cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                    i.more <- sample(sample.space,ceiling(f*d))
                    i.remove <- c(i.remove,i.more)
                    ul <- length(unique(i.remove))
                    cc <- cc + 1
                }
            }
            
            remove.index.t[[chr]] <- read.index.t[[chr]][i.remove]
            total.remove.t <- total.remove.t + n.remove.t[chr]
            
            if (!is.na(control))
            {
                cat("  Sampling reads to remove from control... Selected method: ",rr.score,"\n",sep="")
                
                #cat("\n\n")
                #cat("    ",n.fall.ec[chr]," reads fall into enriched regions\n",sep="")
                #cat("    ",n.fall.nec[chr]," reads fall into non-enriched regions\n",sep="")
                #cat("    ",n.remove.ec[chr]," reads will be removed from enriched regions\n",sep="")
                #cat("    ",n.remove.nec[chr]," reads will be removed from non-enriched regions\n",sep="")
                #cat("\n\n")
                
                ne.deplete <- FALSE
                if (n.remove.nec[chr] > n.fall.nec[chr]) # Must resize appropriately, removing more from enriched regions...
                {
                    d <- n.remove.nec[chr] - n.fall.nec[chr]
                    n.remove.nec[chr] <- n.fall.nec[chr]
                    n.remove.ec[chr] <- n.remove.ec[chr] + d
                    ne.deplete <- TRUE
                }
                if (n.remove.ec[chr] > n.fall.ec[chr]) # This should never happen! Here, only to avoid crash
                {
                    d <- n.remove.ec[chr] - n.fall.ec[chr]
                    n.remove.ec[chr] <- n.fall.ec[chr]
                    n.remove.nec[chr] <- n.remove.nec[chr] + d
                }
                
                #cat("    ",n.fall.ec[chr]," reads fall into enriched regions\n",sep="")
                #cat("    ",n.fall.nec[chr]," reads fall into non-enriched regions\n",sep="")
                #cat("    ",n.remove.ec[chr]," reads will be removed from enriched regions\n",sep="")
                #cat("    ",n.remove.nec[chr]," reads will be removed from non-enriched regions\n",sep="")
                #cat("\n\n")
    
                i.remove <- c(
                    sample(unlist(bin.index$control$ne[[chr]]),n.remove.nec[chr]),
                    sample(unlist(bin.index$control$e[[chr]]),n.remove.ec[chr])
                )
                
                ul <- length(unique(i.remove))
                cc <- 1
                
                if (ne.deplete) # Sample a bit more from enriched (does not happen often)
                {
                    while(ul != n.remove.c[chr])
                    {
                        cat("    Iteration ",cc,"...\n",sep="")
                        d <- n.remove.c[chr] - ul
                        sample.space <- setdiff(unlist(bin.index$control$e[[chr]]),i.remove)
                        if (length(sample.space)<ceiling(f*d)) { # Nothing to sample, already reached
                            cat("      Empty sample space! Goal reached.\n")
                            break
                        }
                        cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                        i.more <- sample(sample.space,ceiling(f*d))
                        i.remove <- c(i.remove,i.more)
                        ul <- length(unique(i.remove))
                        cc <- cc + 1
                    }
                }
                else # Sample a bit more from non-enriched
                {
                    while(ul != n.remove.c[chr])
                    {
                        cat("    Iteration ",cc,"...\n",sep="")
                        d <- n.remove.c[chr] - ul
                        sample.space <- setdiff(unlist(bin.index$control$ne[[chr]]),i.remove)
                            if (length(sample.space)<ceiling(f*d)) { # Nothing to sample, already reached
                            cat("      Empty sample space! Goal reached.\n")
                            break
                        }
                        cat("      Sample space: ",length(sample.space),"\tSample size: ",ceiling(f*d),"\n",sep="")
                        i.more <- sample(sample.space,ceiling(f*d))
                        i.remove <- c(i.remove,i.more)
                        ul <- length(unique(i.remove))
                        cc <- cc + 1
                    }
                }
                
                remove.index.c[[chr]] <- read.index.c[[chr]][i.remove]
                total.remove.c <- total.remove.c + n.remove.c[chr]
            }
        
        } # for chr ends here
        
        cat("Removing actual reads...\n")
        treat.bed <- treat.bed[-unlist(remove.index.t)]
        if (!is.na(control))
            control.bed <- control.bed[-unlist(remove.index.c)]
        
        if (write.output)
        {
            if (output=="auto")
            {
                name <- sub("(.+)[.][^.]+$", "\\1", basename(treatment))
                o <- file.path(dirname(treatment),paste(name,"norm.bed",sep="."))
                cat("Writing output to ",o,"...\n",sep="")
                export(treat.bed,o)
                if (!is.na(control))
                {                
                    name <- sub("(.+)[.][^.]+$", "\\1", basename(control))
                    o <- file.path(dirname(treatment),paste(name,"norm.bed",sep="."))
                    cat("Writing output to ",o,"...\n",sep="")
                    export(control.bed,o)
                }
            }
            else
            {
                cat("Writing treatment output to ",output[1],"...\n",sep="")
                export(treat.bed,output[1])
                if (!is.na(control))
                {
                    cat("Writing output to ",output[2],"...\n",sep="")
                    export(control.bed,output[2])
                }
            }
        }
        
        # Create a nice report!
        cat("\n----- REPORT -----\n")
        cat("Treatment file ",basename(treatment)," has ",treat.depth," reads\n",sep="")
        if (!is.na(control))
            cat("Control file ",basename(control)," has ",control.depth," reads\n",sep="")
        cat("Requested final depth: ",down.to,"\n",sep="")
        cat("Treatment file ",basename(treatment)," was shrinked by ",total.remove.t," reads\n",sep="")
        cat("Of these:\n")
        for (chr in rownames(chrom.info))
        {
            cat("  --- ",chr,"\n",sep="")
            cat("  ",length(gen.win[[chr]]) - l.alive[chr]," dead regions were found in treatment and control\n",sep="")
            cat("  ",length(gen.win.et[[chr]])," enriched regions were found in treatment over control\n",sep="")
            cat("  ",l.alive[chr] - length(gen.win.et[[chr]])," non-enriched regions were found in treatment over control\n",sep="")
            cat("  In these:\n")
            cat("    ",n.fall.et[chr]," reads fall into enriched regions\n",sep="")
            cat("    ",n.fall.net[chr]," reads fall into non-enriched regions\n",sep="")
            cat("    Of these:\n")
            cat("      ",n.remove.et[chr]," reads were removed from enriched regions\n",sep="")
            cat("      ",n.remove.net[chr]," reads were removed from non-enriched regions\n",sep="")
            cat("  Finally ",n.remove.t[chr]," reads were removed from treatment\n",sep="")
        }
        if (!is.na(control))
        {            
            cat("Control file ",basename(control)," was shrinked by ",total.remove.c," reads\n",sep="")
            cat("Of these:\n")
            for (chr in rownames(chrom.info))
            {
                cat("  --- ",chr,"\n",sep="")
                cat("  ",length(gen.win[[chr]]) - l.alive[chr]," dead regions were found in treatment and control\n",sep="")
                cat("  ",length(gen.win.ec[[chr]])," enriched regions were found in control over treatment\n",sep="")
                cat("  ",l.alive[chr] - length(gen.win.ec[[chr]])," non-enriched regions were found in control over treatment\n",sep="")
                cat("  In these:\n")
                cat("    ",n.fall.ec[chr]," reads fall into enriched regions\n",sep="")
                cat("    ",n.fall.nec[chr]," reads fall into non-enriched regions\n",sep="")
                cat("    Of these:\n")
                cat("      ",n.remove.ec[chr]," reads were removed from enriched regions\n",sep="")
                cat("      ",n.remove.nec[chr]," reads were removed from non-enriched regions\n",sep="")
                cat("  Finally ",n.remove.c[chr]," reads were removed from control\n",sep="")
            }
        }
        cat("\n")
        
    }

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

# x: vector of window counts or a length
make.sliding <- function(x,win,slide,chr=NA,type="list")
{
    len = ifelse(length(x)==1,x,length(x))
    if(is.null(attr(x,"names")))
       nams <- NULL
    else 
        nams <- names(x)
    if (is.na(chr) && type=="GRanges")
    {
        warning("chr must be provided for a GRanges object! Will return IRanges instead...")
        type <- "IRanges"
    }
    st <- seq(1,len-win,by=slide)
    en <- seq(win+1,len,by=slide)
    last.slide <- len - en[length(en)]
    st <- c(st,st[length(st)]+last.slide)
    en <- c(en,en[length(en)]+last.slide)
    
    if (type=="IRanges")
        return(IRanges(start=st,end=en))
    else if (type=="GRanges")
        return(
            GRanges(seqnames=Rle(chr,length(st)),
                IRanges(start=st,end=en))
        )
    else
        return(list(start=st,end=en))
}

make.read.index <- function(grange)
{
    read.index <- vector("list",length(seqlevels(grange)))
    for (chr in seqlevels(grange))
    {
        read.index[[chr]] <- which(seqnames(grange)==chr)
    }
    return(read.index)
}

to.GRangesList <- function(grange)
{
    glist <- GRangesList()
    for (n in seqlevels(grange))
    {
        glist[[n]] <- grange[seqnames(grange)==n]
    }
    return(glist)
}

get.chrom.info <- function(org)
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
    rownames(chrom.info) <- chrom.info[,1]
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
