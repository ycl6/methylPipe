# show method
setMethod('show', 'BSdataSet', function(object) {
    message("S4 Object of class BSdataSet")
    message()
    message("BSdata objects contained:")
    print(names(object))
    message()
    message("Groups of objects:")
    print(object@group)
    message()
    message("Associated organism genome:")
    message(organism(object@org))
    message()
})

                                        # subset
setMethod('[[', 'BSdataSet', function(x, i, j = "missing") {
    bsdata <- new('BSdata', file=x@Objlist[[i]]@file,
                  uncov=x@Objlist[[i]]@uncov, org=x@org)
    bsdata
})

setMethod('[[<-', 'BSdataSet', function(x, i, j = "missing", value) {
    x@Objlist[[i]]@file <- value@file
    x@Objlist[[i]]@uncov <- value@uncov
    x
})

setMethod('[','BSdataSet', function(x, inds, i= 'missing', j= 'missing',
                                    drop = "missing") {
    x@Objlist = x@Objlist[inds]
    x
})

                                        # length
setMethod('length', 'BSdataSet', function(x) {
    length(x@Objlist)
})


setMethod("$", "BSdataSet", function(x, name){
    x@Objlist[[name]]
})


                                        # get DMR GRanges
setGeneric('findDMR',
           function(object, Nproc=10, ROI=NULL,
                    pmdGRanges=NULL, MCClass='mCG',
                    dmrSize=10, dmrBp=1000, binsize=0,
                    eprop=0, coverage=1, Pvalue=NULL, SNPs=NULL)
           standardGeneric('findDMR'))
setMethod('findDMR', 'BSdataSet', function(object, Nproc=10, ROI=NULL,
                                           pmdGRanges=NULL, MCClass='mCG',
                                           dmrSize=10, dmrBp=1000,
                                           binsize=0, eprop=0, coverage=1, Pvalue=NULL, SNPs=NULL) {
    if(!is.numeric(Nproc))
        stop('Nproc has to be of class numeric ..')
    if(!is.null(pmdGRanges) && !is(pmdGRanges, "GRanges"))
        stop('pmdGRanges has to be either NULL or an object of class GRanges ..')
    if(!is.null(ROI) && !is(ROI, "GRanges"))
        stop('ROI has to be either NULL or an object of class GRanges ..')
    if(!(MCClass %in% c('mCG','mCHG','mCHH', 'all')))
        stop('MCClass has to be one of mCG, mCHG, mCHH or all')
    if(!is.numeric(dmrSize))
        stop('dmrSize has to be an object of class numeric ..')
    if(dmrSize < 5)
        stop('dmrSize has to be atleast 5 ..')
    if(!is.numeric(dmrBp))
        stop('dmrBp has to be an object of class numeric ..')
    if(!is.numeric(binsize))
        stop('binsize has to be an object of class numeric ..')
    if(!is.numeric(eprop))
        stop('eprop has to be an object of class numeric ..')
    if((eprop < 0) || (eprop > 1))
        stop('eprop has to be from 0 to 1 ..')
    if(!is.null(SNPs) && !is(SNPs,"GRanges"))
        stop('SNPs has to be either NULL or an object of class GRanges ..')
    if(!is.null(Pvalue) && !is.numeric(Pvalue))
        stop('Pvalue has to be either NULL or an object of class numeric ..')
    if(!is.numeric(coverage))
        stop('coverage has to be an object of class numeric ..')

                                        # parallelize tasks based on each chromosome

    DMRchr <- function(Ind, Blocks=blocks, PMDs=pmdGRanges,
                       samples=object, MCClass=MCClass,
                       DmrSize=dmrSize, DmrBp=dmrBp,
                       Binsize=binsize, Eprop=eprop, Coverage=coverage, pValue=Pvalue, SNP=SNPs) {
                                        # Ind is the requested region of this call of the DMRchr function
        Chr <- as.character(Blocks[Ind,1])
        start <- round(Blocks[Ind,2])
        refgr <- GRanges(seqnames = Chr, ranges = IRanges(Blocks[Ind,2], end = Blocks[Ind,3]))
        minC <- DmrSize

                                        # mC calls are retrieved for each sample
        resList <- list()
        Cposind <- NULL

        for(i in 1:length(samples)) {
            message('retreiving methylation data from sample.. ')
            resList[[i]] <- mapBSdata2GRanges(GenoRanges= refgr, Sample= samples[[i]], depth=0,
                                              context= sub('m', '', MCClass))[[1]]
            if(!is(resList[[i]],"GRanges")) return(NULL)
            Cposind <- c(Cposind, start(resList[[i]]))
        }


        if(Binsize > 0){
            message('getting mc data .. ')
            Cposind <- NULL
            Cposind <- unlist(getCpos(refgr, seqContext = sub('m', '', MCClass), nbins = 1,
                                      org = samples@org))
            if(length(Cposind)==0) return(NULL)
            Cposind <- round(unique(Cposind))
            Cposind <- sort(Cposind)
            Cposind <- (Cposind+start-1)
        }
        else
            {
                Cposind <- round(unique(Cposind))
                Cposind <- sort(Cposind)
    }


        if(Coverage > 1)
            {
                covrm <- NULL
                for(i in 1:length(samples)) {
                    message('removing minimum coverage cutoff cytosines.. ')
                    inds <- which((mcols(resList[[i]])$C+mcols(resList[[i]])$T) < Coverage)
                    covrm <- c(covrm, start(resList[[i]])[inds])
                }
                torm <- which((Cposind %in% covrm)==TRUE)
                Cposind <- Cposind[-c(torm)]
            }

        message('building ratioGR .. ')
        ratioMat <- matrix(0, length(Cposind), length(samples))
        ratioGR <- GRanges(Rle(Chr), IRanges(Cposind, Cposind), elementmetaData=ratioMat)
        rm(ratioMat)
        names(mcols(ratioGR)) <- samples@group

                                        # removing SNPs if provided
        if(!is.null(SNP))
            {
                message("removing SNPs")
                ov <- findOverlaps(ratioGR, SNP)
                if(length(ov)>0) ratioGR <- ratioGR[-unique(queryHits(ov))]
            }

                                        # fixing unmethylated C, fill ratioGR

        for(i in 1:length(samples)) {
            message('fixing unmethylated C .. ')
                                        # removing cytosines positions not covered by sequencing
            uncov <- samples[[i]]@uncov
            uncov <- uncov[seqnames(uncov) == Chr]
            ov <- findOverlaps(ratioGR, uncov)
            if(length(ov)>0) ratioGR <- ratioGR[-unique(queryHits(ov))]
            if(length(ratioGR)== 0) return(NULL)

            message('filling ratioGR .. ')
            ind <- findOverlaps(resList[[i]], ratioGR)
            mc <- mcols(resList[[i]])$C/(mcols(resList[[i]])$C+ mcols(resList[[i]])$T)
            mcols(ratioGR)[,i][subjectHits(ind)] <- round(mc[queryHits(ind)],3)

            if(!is.null(pValue))
                {
                    message('keeping only those mC less than Pvalue threshold..')
                    pV <- -10*log10(pValue)
                    ind <- which(mcols(resList[[i]])$Significance < pV)
                                        # only if some mC crosses pvalue threshold those are removed from ratioGR
                    if(length(ind) > 0)
                        {
                            ov <- findOverlaps(ratioGR, resList[[i]][ind])
                            if(length(ov)>0) ratioGR <- ratioGR[-unique(queryHits(ov))]
                        }
                }

        }
        rm(resList)

                                        # if PMD are passed, DMR are not searched within them

        if(!is.null(PMDs))
            {
                message('filtering PMDs out ..')
                genoranges <- PMDs[seqnames(PMDs) == Chr]
                if(length(genoranges)!=0)
                    {
                        ov <- findOverlaps(ratioGR, genoranges)
                        if(length(ov)>0) ratioGR <- ratioGR[-unique(queryHits(ov))]
                    }
            }

        Cposind <- as.numeric(start(ratioGR))
                                        # smoothing ratioMat (useful for nonCG in human for example ..)
        if(Binsize > 0){
            message('applying smoothing for nonCG .. ')
            smbins <- seq(min(Cposind)-1, max(Cposind), Binsize)
            if(smbins[length(smbins)] != max(Cposind))
                smbins <- c(smbins, max(Cposind))
            ncl <- length(mcols(ratioGR))
            smMat <- matrix(NA, length(smbins)-1, ncl)
                                        # determinig the average met level within each bin
            for(icol in 1:ncl) {
                ratios <- unlist(as.data.frame(mcols(ratioGR[,icol])))
                                        # -100 identify NA in C code ..
                ratios[is.na(ratios)] <- -100
                res <- .C(.binning, score= as.double(ratios),
                          pos= as.double(Cposind),
                          poslength=as.integer(length(Cposind)),
                          binvec= smbins, binlength=as.integer(length(smbins)),
                                        # binout contains the bin avg
                          binout= as.double(rep(-100, length(smbins)-1)),
                          doavg= as.integer(1), PACKAGE='methylPipe')$binout
                res[res == (-100)] <- NA
                smMat[,icol] <- res
            }
                                        # ratioGR is replaced by the smoothed matrix
            Cposind <- smbins[-length(smbins)]
            allNAinds <- which(apply(smMat, 1, function(x) all(is.na(x))))
                                        # all NA rows are discarded
            if(length(allNAinds) > 0) {
                smMat <- smMat[-allNAinds,]
                Cposind <- Cposind[-allNAinds]
            }
            ratioGR <- GRanges(Rle(Chr), IRanges(Cposind,Cposind), elementmetaData=smMat)
            rm(smMat)
            names(mcols(ratioGR)) <- samples@group
            minC <- 6
        }
                                        #ratioMat <- as.data.frame(ratioMat)
        maxR <- length(ratioGR)
                                        # Eprop can be used to check only those mC showing
                                        # the highest empirical differential methylation between
                                        # samples to save time
        if(Eprop > 0) {
            mat <- as.data.frame(mcols(ratioGR))
            maxs <- apply(mat, 1, max, na.rm=T)
            mins <- apply(mat, 1, min, na.rm=T)
                                        # the difference between the min and the max met level for a
                                        # given mC over all the samples is used to determine the
                                        # empirical differential methylation
            diffs <- maxs-mins
            diffInds <- which(diffs >= Eprop)
            if(length(diffInds) == 0) return(NA)
            rm(mat)
        }
                                        # else all the mC are considered
        else diffInds <- 1:maxR

        NAvec=rep(NA, length(diffInds))
        dmrDf <- data.frame(chr=rep(Chr, length(diffInds)),
                            start=NAvec, end=NAvec, pValue=NAvec,
                            MethDiff_Perc=NAvec, log2Enrichment=NAvec, stringsAsFactors=FALSE)
        counter <- 0
        oldPos <- min(Cposind)-DmrBp
        message('looking for DMRs .. ')
        for(i in 1:length(diffInds)) {
            counter <- counter+1
            if(counter/100 == trunc(counter/100)) message(paste(counter, '', sep=''))
            ind <- diffInds[i]
            if(as.logical((Cposind[ind]- oldPos) < DmrBp/2)) next
                                        # looking from the current position ind till ind + dmrSize
                                        # (the number of mC in each evaluated block)
            indNext <- ind+DmrSize-1
                                        #if(as.logical((Cposind[indNext]- oldPos) < (DmrBp/2))) next
                                        # checking if over the considered genomic region
            if(indNext > maxR) break
                                        # checking if the distance between the end and the
                                        # begining of the block is within DmrBp
            if((Cposind[indNext]-Cposind[ind]) < DmrBp) {
                df <- as.data.frame(mcols(ratioGR)[ind:indNext,])
                                        # filtering the df based on having some data > 0 there ..
                df <- df[rowSums(df,na.rm=TRUE) > 0,]
                minN <- max(5, round(0.8*minC)) # minimum number of available data required
                if(nrow(df) < minN) next
                cs <- apply(df, 2, function(x) length(which(x > 0)))
                                        # at least one sample with at least minN positions with data > 0
                if(max(cs, na.rm=T) < minN) next
                else {
                                        # non parametric kruskall-wallis test for multi-sample comparison
                    if(length(samples) > 2) {
                        colinds <- which(apply(df, 2, function(x) !all(is.na(x))))
                                        # checking that at least 2 samples remain
                        if(length(colinds) == 1) Pv=NA
                                        # determining Pvalue
                        else Pv <- kruskal.test(df[,colinds])$p.value
                        groups <- unique(names(mcols(ratioGR)))
                        if(any( !groups %in% c("C","E")))
                            {
                                dmrDf$MethDiff_Perc[i] <- NA
                                dmrDf$log2Enrichment[i] <- NA
                            }
                        else
                            {
                                means <- colMeans(df, na.rm=T)
                                means <- t(as.matrix(means))
                                colnames(means) <- names(mcols(ratioGR))
                                ind_C <- which(colnames(means)=="C")
                                ind_E <- which(colnames(means)=="E")
                                        # the methylation difference is determined
                                Methdiff <- (mean(means[ind_E])-mean(means[ind_C]))*100
                                dmrDf$MethDiff_Perc[i] <- round(Methdiff, 3)
                                log2Enrichment <- log2(mean(means[ind_E]) / mean(means[ind_C]))
                                dmrDf$log2Enrichment[i] <- round(log2Enrichment,3)
                            }
                    }
                    else {
                                        # non parametric wilcoxon.test for 2 sample comparison
                        if(any(apply(df, 2,function(x) all(is.na(x))))) Pv <- NA
                        else {
                            noNArows <- which(apply(df,1,
                                                    function(x) length(which(is.na(x))))==0)
                            if(length(noNArows) < minN)  Pv <- NA
                            else
                                {
                                    Pv <- wilcox.test(df$C, df$E,alternative='two',
                                                      paired=TRUE, exact=FALSE)$p.value
                                }
                            Methdiff <- (mean(df$E, na.rm=T)-mean(df$C, na.rm=T))*100
                            dmrDf$MethDiff_Perc[i] <- round(Methdiff, 3)
                            log2Enrichment <- log2(mean(df$E, na.rm=TRUE) / mean(df$C, na.rm=TRUE))
                            dmrDf$log2Enrichment[i] <- round(log2Enrichment,3)
                        }
                    }
                    dmrDf$start[i] <- Cposind[ind]
                    dmrDf$end[i] <- Cposind[indNext]+Binsize
                    dmrDf$pValue[i] <- round(Pv,3)
                }
                oldPos <- Cposind[ind]
            }
        }
        NAinds <- which(is.na(dmrDf$start))
        if(length(NAinds)>0) dmrDf <- dmrDf[-NAinds,]
        if(nrow(dmrDf)==0) dmrDf <- NA
        return(dmrDf)
    }

    Chrs <- list()
    for (i in 1:length(object))
        {
            tabixRef <- TabixFile(object[[i]]@file)
            Chrs[[i]] <- as.character(seqnamesTabix(tabixRef))
        }
    Chrs <- unique(unlist(Chrs))

    if(!is.null(ROI))
        {
            blocks <- as.data.frame(ROI)
            blocks <- blocks[,-c(4:5)]
        }
    else
        blocks <- splitChrs(chrs=Chrs, org=object@org)

    cl <- makeCluster(Nproc, type='PSOCK')
    clRes <- clusterApplyLB(cl, 1:nrow(blocks),
                            DMRchr, Blocks=blocks, PMDs=pmdGRanges,
                            samples=object, MCClass=MCClass,
                            DmrSize=dmrSize, DmrBp=dmrBp,
                            Binsize=binsize, Eprop=eprop, Coverage=coverage, pValue=Pvalue, SNP=SNPs)
    DmrDf <- NULL
    for(clind in 1:length(clRes)) {
        if(is.null(clRes[[clind]])) next
        if(is.null(nrow(clRes[[clind]]))) next
        DmrDf <- rbind(DmrDf, clRes[[clind]])
    }
    if(is.null(DmrDf)) return(NULL)
    NAinds <- which(is.na(DmrDf$pValue))
    if(length(NAinds) > 0) DmrDf <- DmrDf[-NAinds,]
    DmrGR <- GRanges(DmrDf[,1], IRanges(DmrDf[,2], DmrDf[,3]), pValue= DmrDf[,4],
                     MethDiff_Perc= DmrDf[,5], log2Enrichment= DmrDf[,6])
    DmrGR <- sort(DmrGR)
})


setGeneric('methstats',
           function(object, chrom, mcClass='mCG', minC=1, coverage=1, pval=1, Nproc=1)
             standardGeneric('methstats'))
setMethod('methstats', 'BSdataSet', function(object, chrom, mcClass='mCG', minC=1, 
                                             coverage=1, pval=1, Nproc=1) {
  if(!is.null(chrom) && !is.character(chrom))
    stop('chrom has to be either NULL or an object of class character ..')
  if(!(mcClass %in% c('mCG','mCHG','mCHH', 'all')))
    stop('mcClass has to be one of mCG, mCHG, mCHH or all..')
  if(!is.numeric(minC) && minC < 0 )
    stop('minC has to be of class numeric and positive..')
  if(!is.numeric(coverage))
    stop('coverage has to be of class numeric ..')
  if(!is.numeric(pval))
    stop('pval has to be of class numeric ..')
  if(pval > 1)
    stop('pval has to be lower than 1 ..')
  if(!is.numeric(Nproc))
    stop('Nproc has to be of class numeric ..')
  
  
  # parallelize tasks based on each chromosome
  Methchr <- function(Ind, Blocks, samples, mcContext, mCs, cov, Pval)
  {
    resList <- list()
    Cposind <- NULL
    Chr <- as.character(Blocks[Ind,1])
    gr <- GRanges(seqnames = Chr, ranges = IRanges(Blocks[Ind,2], end = Blocks[Ind,3]))
    
    for(j in 1:length(samples)) {
      message('retreiving methylation data from sample.. ')
      resList[[j]] <- mapBSdata2GRanges(GenoRanges= gr, Sample= samples[[j]], depth=cov,
                                        context= sub('m', '', mcContext), mC=mCs, pValue=Pval)[[1]]
      if(!is(resList[[j]],"GRanges")) return(NULL)
      Cposind <- c(Cposind, start(resList[[j]]))    
    }
    
    Cposind <- round(unique(Cposind))
    Cposind <- sort(Cposind)
    
    message('building ratioGR .. ')
    ratioMat <- matrix(0, length(Cposind), length(samples))
    ratioGR <- GRanges(Rle(Chr), IRanges(Cposind, Cposind), elementmetaData=ratioMat)
    rm(ratioMat)
    names(mcols(ratioGR)) <- names(samples)
    
    for(i in 1:length(samples)) {
      message('fixing unmethylated C .. ')
      # removing cytosines positions not covered by sequencing
      uncov <- samples[[i]]@uncov
      uncov <- uncov[seqnames(uncov) == Chr]
      ov <- findOverlaps(ratioGR, uncov)
      if(length(ov)>0) ratioGR <- ratioGR[-unique(queryHits(ov))]
      if(length(ratioGR)== 0) return(NULL)
      
      message('filling ratioGR .. ')
      ind <- findOverlaps(resList[[i]], ratioGR)
      mc <- mcols(resList[[i]])$C/(mcols(resList[[i]])$C+ mcols(resList[[i]])$T)
      mcols(ratioGR)[,i][subjectHits(ind)] <- round(mc[queryHits(ind)],3)    
    }
    ratioGR
  }
  
  Chrs <- list()
  for (i in 1:length(object))
  {
    tabixRef <- TabixFile(object[[i]]@file)
    Chrs[[i]] <- as.character(seqnamesTabix(tabixRef))
  }
  Chrs <- unique(unlist(Chrs))
  
  if(!is.null(chrom))
  {
    if(!chrom %in% Chrs) 
    {
      stop('chromosome not available in samples .. ')
    }
    blocks <- splitChrs(chrs=Chrs, org=object@org)
    ind <- which(blocks[,1]==chrom)
    blocks <- blocks[ind,]
  }
  else
    blocks <- splitChrs(chrs=Chrs, org=object@org)
  
  cl <- makeCluster(Nproc, type='PSOCK')
  #clusterEvalQ(cl, library(GenomicRanges))
  clRes <- clusterApplyLB(cl, 1:nrow(blocks),
                          Methchr, Blocks=blocks, samples=object, 
                          mcContext=mcClass, mCs=minC, cov=coverage, Pval=pval)
  ind <- which(is.na(clRes)==TRUE)
  if(length(ind)!=0)
    clRes <- clRes[-c(ind)]
  if(length(clRes)==1)
    clRes <- clRes[[1]]
  else
  {
    clRes_comb <- GRangesList(clRes)
    clRes <- unlist(clRes_comb)
    rm(clRes_comb)
  }
  mc_table <- as.data.frame(mcols(clRes))
  rm(clRes)
  
  
  ########## statistics and plots ##########
  
  ######### descriptive stats ##########
  statistics_results <- list()
  statistics_results[[1]] <- summary(mc_table)
  statistics_results[[2]] <- cor(mc_table)
  names(statistics_results) <- c("descriptive_stats","correlation_mat")
  
  ##### correlation plots #########
  
  panel.cor <- function(x, y, digits=2, cex.cor)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    test <- cor.test(x,y)
    Signif <- ifelse(round(test$p.value,3)<0.01,"p<0.01",paste("p=",round(test$p.value,3)))  
    text(0.5, 0.25, paste("r=",txt))
    text(.5, .75, Signif)
  }
  
  panel.smooth<-function (x, y, col = "orange", bg = NA, pch = 46, 
                          cex = 1.5, col.smooth = "red", span = 2/3, iter = 3, ...) 
  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
  }
  
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    opar=par(ps=10)
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="brown", ...)
  }
  
  pairs(mc_table,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist, 
          cex.labels = 1,font.labels = 2,main="Correlation Matrix")
  ######### clustering plots ####
  
  clust_mat <- t(mc_table)
  dist_mat <- dist(clust_mat, method = "euclidean")
  x <- hclust(dist_mat)
  dev.new()
  plot(x, main="Methylation Clustering", xlab="Samples")
  return(statistics_results)
}
)