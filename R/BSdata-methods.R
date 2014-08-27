# show method
setMethod('show', 'BSdata', function(object) {
    message("S4 Object of class BSdata")
    message()
    message("TABIX indexed file for this BSdata:")
    message(object@file)
    message()
    message('The first lines of the uncovered regions of this BSdata:')
    print(head(object@uncov))
    message("Associated organism genome:")
    message(organism(object@org))
    message()

                                        # requires Rsamtools TabixFile function...
    tabixRef <- TabixFile(object@file, yieldSize=1000L)
    message('Chromosomes available:')
    print(seqnamesTabix(tabixRef))
    tbdata <- scanTabix(tabixRef)
    tbdata <- unlist(tbdata)
                                        # tabixdata2GR in Allfunctions.R
    tbdata <- tabixdata2GR(tbdata)
    message('The first lines of the data:')
    print(head(tbdata))
})




setGeneric('mCsmoothing', function(Object, refgr, Scorefun='sum', Nbins=20,
                                   Context='CG', plot=TRUE)
           standardGeneric('mCsmoothing'))
setMethod('mCsmoothing','BSdata', function(Object, refgr, Scorefun='sum',
                                           Nbins=20, Context='CG', plot=TRUE) {
    chrs <- seqnames(Object@org)
    chr <- unique(as.character(seqnames(refgr)))
    if(any(!chr %in% chrs))
        stop('This is not a valid chromosome .. ')
    if(!is(refgr, "GRanges"))
        stop('refgr has to be of class GRanges ..')
    if(Scorefun != 'sum' && Scorefun != 'mean')
        stop('Scorefun has to be either sum or mean .. ')
    if(length(which(!(Context %in% c('all','CG','CHG','CHH')))) > 0)
        stop('Context has to be either all or a combination of CG, CHG, and CHH ..')
    if(!is.logical(plot))
        stop('plot has to be of class logical ..')
    if(!is.numeric(Nbins))
        stop('Nbins has to be of class numeric')
    if(!(Nbins >5))
      stop('Nbins has to be of greater than 5')

                                        # a generic function to determing data binning
                                        # through the binning C function included in the package
    chrdensity <- function(GenoRanges, pos, score, scorefun, nbins) {
        bins <- NULL
        binsize <- width(GenoRanges)/nbins
        for (bin in 1:(nbins+1))
          bins[bin] <- round(start(GenoRanges) + (bin-1)*binsize)
        if(scorefun == 'sum') res <- .C(.binning,score=as.double(score),
               pos=as.double(pos), as.integer(length(pos)),
               bins=bins, as.integer(length(bins)),
               binout=as.double(rep(0,length(bins))), doavg=as.integer(0),
               PACKAGE='methylPipe')
        if(scorefun == 'mean') res <- .C(.binning,score=as.double(score),
               pos=as.double(pos), as.integer(length(pos)),
               bins=bins, as.integer(length(bins)),
               binout=as.double(rep(0,length(bins))), doavg=as.integer(1),
               PACKAGE='methylPipe')
        return(res$binout)
    }
    mcdata <- mapBSdata2GRanges(GenoRanges=refgr, Sample=Object,
                                context=Context, mC=0, depth=0, pValue=1)
    ind <- which(is.na(mcdata)==TRUE)
    if(length(ind)!=0)
    {
      mcdata <- mcdata[-c(ind)]
      refgr <- refgr[-c(ind)]
    }

                                        # using the methylation level (proportion of C/(C+T) reads for a given position)
    Score <- sapply(mcdata, function(x) mcols(x)$C/(mcols(x)$C+mcols(x)$T))
    if(!is(Score, "list"))
        Score <- list(Score)

    mCdChr <- NULL
    Pos <- sapply(mcdata, function(x) start(x))
    if(!is(Pos, "list"))
        Pos <- list(Pos)

    for (i in 1:length(mcdata))
        {
        temp <- chrdensity(GenoRanges=refgr[i], pos=c(Pos[[i]]),
                           score=c(Score[[i]]), scorefun=Scorefun,
                           nbins=Nbins)
        mCdChr[[i]] <- temp[-c(length(temp))]
    }

    mCdChr <- matrix(unlist(mCdChr), ncol = Nbins, byrow = TRUE)
    mCdChr <- apply(mCdChr, 2, mean)

                                        # data are smoothed
    Cs <- smooth.spline(x=1:length(mCdChr),
                        y=mCdChr, spar=0.5, keep.data=FALSE)
    totalwidth <- width(refgr)
    bin <- round(Nbins*.1)
    binsize <- totalwidth/bin
    bins <- seq(start(refgr),totalwidth,binsize)
    if(plot)
    {
      plot(x=Cs$x, y=100*Cs$y/max(Cs$y), axes=FALSE,xlab=NA,
           ylab=NA, type='l', lwd=2, ylim=c(0,100), col="darkgreen",
           main="Smoothed methylation levels")
      box()
      axis(side = 2)
      axis(side = 1, at=c(seq(1,Nbins,round(Nbins/bin)),Nbins),
           labels=c(bins,end(refgr)))
      mtext(side = 1, "Genomic co-ordinates", line = 3)
      mtext(side = 2, "Methylation Level (mC/total reads)", line = 3)
  }

    return(list(pos=c(1:Nbins),
                score=mCdChr, smoothed=Cs))
})


setGeneric('findPMDs', function(Object, Nproc=10, Chrs=NULL)
           standardGeneric('findPMDs'))
setMethod('findPMDs', 'BSdata', function(Object, Nproc=10, Chrs=NULL){
    if(!is.null(Chrs) && !is.character(Chrs))
        stop("Chrs has to be either NULL or class character")
    ind <- Chrs %in% seqnamesTabix(Object@file)
    if(any(ind==FALSE))
        stop("Chrs not found in the sample...")
    if(is.null(Chrs))
        Chrs <- seqnamesTabix(Object@file)

    PMDchr <- function(Ind, object, Blocks, org)
        {
            refgr <- GRanges(Blocks$chr[Ind],IRanges(Blocks$start[Ind],Blocks$end[Ind]))
            meth.gr <- mapBSdata2GRanges(GenoRanges= refgr, Sample= object, context= "CG")[[1]]
            mcols(meth.gr) <- mcols(meth.gr)[3:2]
            mcols(meth.gr)$T <- mcols(meth.gr)$T+mcols(meth.gr)$C
            chr <- Blocks[Ind,1]
            chrL <- seqlengths(org)[chr]
            options(menu.graphics=FALSE)
            PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel=chr, seqLengths=chrL, num.cores=1)
        }

    blocks <- splitChrs(chrs=Chrs, org=Object@org)
    cl <- makeCluster(Nproc, 'PSOCK')
    clusterEvalQ(cl, library(MethylSeekR))
    clRes <- clusterApplyLB(cl, 1:nrow(blocks), PMDchr, object=Object,
                            Blocks=blocks, org=Object@org)
    PMDscomb <- GRangesList(clRes)
    return(PMDscomb)
})

