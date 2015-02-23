                              # function to determine hmC level from mlml output file

process.hmc <- function(file, output_folder, Coverage){
    if(!is.character(file))
        stop('file has to be of class character ..')
    if(!is.character(output_folder))
        stop('output_folder has to be of class character ..')
    if(!is(Coverage,"GRanges"))
        stop('Coverage has to be of class GRanges with column name coverage..')

    output_folder <- normalizePath(output_folder)  
    sample_name <- unlist(strsplit(file, split=".txt"))
    output.files <- list()
    temp_data <- fread(file, sep="\t", header=FALSE)
    temp_data_gr <- GRanges(temp_data[,1],IRanges(temp_data[,2], temp_data[,2]))
    ov <- findOverlaps(temp_data_gr, Coverage)
    mcols(temp_data_gr)$coverage <- 0
    ind1 <- queryHits(ov)
    ind2 <- subjectHits(ov)
    mcols(temp_data_gr)$coverage[ind1] <- mcols(Coverage)$coverage[ind2]
    temp_data[,4] <- round((temp_data[,4]*mcols(temp_data_gr)$coverage))
    temp_data[,5] <- round((temp_data[,5]*mcols(temp_data_gr)$coverage))
    temp_data[,6] <- round((temp_data[,6]*mcols(temp_data_gr)$coverage))
    ind <- which(temp_data[,7]==0)
    temp_data <- temp_data[ind,]
    temp_data_CpG <- temp_data[,c(1,2,4,6)]
    temp_data_hmC <- temp_data[,c(1,2,5,6)]
    filename <- file.path(output_folder,
                          paste(sample_name, "_", "CpG", ".txt", sep=""))
    write.table(temp_data_CpG, file=filename)
    filename <- file.path(output_folder,
                          paste(sample_name, "_", "hmC", ".txt", sep=""))
    write.table(temp_data_hmC, file=filename)
}

                                        # function to determine the methylation

meth.call <- function(files_location, output_folder, no_overlap, read.context, Nproc){
    if(!is.character(files_location))
        stop('files_location has to be of class character ..')
    if(!is.character(output_folder))
        stop('output_folder has to be of class character ..')
    if(!is.logical(no_overlap))
        stop('no_overlap has to be of class logical ..')
    if(!is.character(read.context))
        stop('read.context has to be of class character ..')
    if(!is.numeric(Nproc))
        stop('Nproc has to be of class numeric ..')
    
    files_location <- normalizePath(files_location)
    output_folder <- normalizePath(output_folder)

                                        # read.type
    if( !( read.context %in% c("CpG","All")))
        stop("wrong 'read.context' argument supplied, only 'CpG' or 'All' accepted")

                                        # list of input sam files
    all_files <- list.files(path = files_location, pattern = ".sam")
    sample_name <- unlist(strsplit(all_files, split=".sam"))

                                        # list of output files
    output.files <- list()
    output.files[[read.context]] <- file.path(output_folder,
                                              paste(sample_name, "_", read.context, ".txt", sep=""))
                                        # creation of system command
    read.files <- file.path(files_location, all_files)
    perl.arguments <- paste("--file", read.files, "--sam_type paired_sam")
    if(read.context == "CpG" )
        {perl.arguments <- paste(perl.arguments, "--CpG", output.files[["CpG"]] )}
    if(read.context == "All" )
        {perl.arguments <- paste(perl.arguments, "--All", output.files[["All"]] )}
    if(no_overlap){perl.arguments <- paste(perl.arguments, "--no_overlap" )}
                                        # location of the perl script
    perl.loc <- (system.file("exec", "methinfo.pl", package="methylPipe"))
    cmd <- paste("perl", perl.loc, perl.arguments)
                                        # then call perl to process the file
    message(    )
    message("Extracting methylation information from input SAM files and creating output files for each sample...")
    message()
    runperl <- function(x)
        {
            status <- try(system(x))
            if(status != 0)
                {stop("\nError in methylation calling...\n
                Check if it is a sorted Bismark SAM file\n")}
    }

    parallel::mclapply(cmd, runperl, mc.cores=Nproc, mc.preschedule=TRUE)
    
    message("Creating temporary BAM files from input SAM files...")
    message()
    setwd(output_folder)
    system('mkdir tempBAM')
    temp_folder <- paste0(output_folder,"/tempBAM")

    for (i in 1:length(all_files))
        asBam(file.path(files_location, all_files[i]), destination=file.path(temp_folder, sample_name[i]), overwrite=TRUE)

    message("Creating uncovered regions objects for each sample from BAM files...")
    message()
    bam.files <- file.path(path = temp_folder, paste(sample_name, ".bam", sep=""))
    for (i in 1:length(all_files))
        {
            aln_file <- bam.files[i]
            aln <- readGAlignments(aln_file)
            cov <- coverage(aln)
            ind <- sapply(cov, function(x)(length(which(runValue(x)== 0))))
            cov <- cov[which(ind > 1)]
            uncov_GR <- as(cov, "GRanges")
            uncov_GR <- uncov_GR[mcols(uncov_GR)$score== 0]
            filename <- file.path(output_folder, paste0(sample_name[[i]],"_uncov", ".Rdata"))
            Objectout <- paste0(sample_name[[i]],"_uncov")
            assign(Objectout, uncov_GR)
            save(list=Objectout, file=filename)
        }
    message('Removing all temporary BAM files...')
    message()
    str <- paste("rm -r","tempBAM")
    system(str)
    message('Methylation info and Uncovered regions output files created for each sample in output_folder...')
    message()
    message("Processing done Successfully...")
    message()
}


tabixdata2GR <- function(x) {
    if(length(x) == 0) return(NA)
    temp <- strsplit(x, split="\t")
    mat  <- matrix(unlist(temp), ncol=7, byrow=TRUE)
    df   <- as.data.frame(mat, stringsAsFactors=FALSE)
    df_gr <- GRanges(df[,1], IRanges(as.numeric(df[,2]), as.numeric(df[,2])),
                     strand=df[,3], Context=df[,4], C=as.numeric(df[,5]), T=as.numeric(df[,6]), Significance=as.numeric(df[,7]))
    rm(df)
    rm(mat)
    return(df_gr)
}

BSdata <- function(file, uncov, org) {
    bsdata <-  new('BSdata', file=file, uncov=uncov, org=org)
    return(bsdata)
}


BSdataSet <-  function(org, group, ...) {
    Objname <- names(list(...))
    if(!is.null(Objname))
        new("BSdataSet", org=org, group=group, Objlist=list(...), names=names(list(...)))
    else stop('define the names for the Objects in the arguments..')
}


GElist <- function(...)
{
    Objname <- names(list(...))
    if(!is.null(Objname))
        new("GElist", Objlist=list(...), names=names(list(...)))
    else stop('define the names for the Objects in the arguments..')
}


BSprepare <-  function(files_location, output_folder, tabixPath, bc=1.5/100) {
    if(!is.character(files_location))
        stop('files_location has to be of class character ..')
    if(!is.character(output_folder))
        stop('output_folder has to be of class character ..')
    if(!is.character(tabixPath))
        stop('tabixPath has to be of class character ..')
    if(!file.exists(paste(tabixPath, '/tabix', sep='')))
        stop('tabix not found at tabixPath ..')
    if(!file.exists(paste(tabixPath, '/bgzip', sep='')))
        stop('bgzip not found at tabixPath ..')
    
    files_location <- normalizePath(files_location)
    output_folder <- normalizePath(output_folder)
    tabixPath <- normalizePath(tabixPath)
    binomTestMulti <- function(mat, maxN=50, p=bc, bh=TRUE) {
        bmat <- matrix(NA, maxN, maxN)
        for(i in 1:maxN) {
            for(j in i:maxN) bmat[i,j] <- binom.test(x=i, n=j, p=p,
                                                     alternative='greater')$p.value
        }
        pVfun <-  function(x) {
            if(x[2] < (maxN+1)) return(bmat[x[1], x[2]])
            else return(binom.test(x=x[1], n=x[2], p=p,
                                   alternative='greater')$p.value)
        }
        pValues <-  apply(mat, 1, pVfun)
        if(bh) pValues <-  p.adjust(pValues, method='BH')
        pValues <-  -round(10*log10(pValues))
        return(pValues)
    }
    
    
    ####### checking for the chromosome column ##########
    all_files <- list.files(path = files_location, pattern = ".txt")
    for (i in 1:length(all_files))
    {
      path <- paste0(files_location,"/",all_files[[i]])
      temp_data <- fread(path,nrows=10)
      if(length(grep("chr",temp_data$V1))==0)
      {
        cmd <- paste("sed -i 's/\"//g'",path)
        system(cmd)
        cmd <- paste("sed -i 's/^/chr/'",path) 
        system(cmd)
      }
    }
    
    ############### sorting chromosome files ############
    
    sample_name <- unlist(strsplit(all_files, split=".txt"))
    cmd <- paste("sort -k1,1 -k2,4n", paste0(files_location, "/", all_files) , ">", 
                 paste0(files_location, "/", sample_name,"_sort.txt"))
    sapply(cmd,system)
    cmd <- paste("mv", paste0(files_location, "/", sample_name,"_sort.txt"), 
                 paste0(files_location, "/", sample_name,".txt"))
    sapply(cmd,system)
    
    ######################################################
    
    for(i in 1:length(all_files)) {
        filecg <- all_files[i]
        message(filecg)
        filemC <- paste0(files_location, "/", filecg, '.mC')
        str <- paste('cut -f5,6', paste0(files_location, "/", filecg), '>', filemC)
        system(str)
                                        # computing binomial pvalues
        mCdat <- fread(filemC)
        mCdat <- as.matrix(mCdat)
        mCdat[,2] <- mCdat[,1]+ mCdat[,2]
        pV <- binomTestMulti(mCdat, p=bc)
                                        # attaching pvalues

        filepV <- paste0(filemC, '.pvalues')
        fileTabix <- paste0(files_location, "/", sample_name[i], '_tabix.txt')
        write(pV, file=filepV, ncolumns=1)
        str <- paste('paste', paste0(files_location, "/", filecg), filepV, '>', fileTabix)
        system(str)
        system(paste('rm', filepV, filemC))
    }
                                        # concatenating, compressing and building tabix index
    message('postprocessing ..')
    Tabix_files <- list.files(path = files_location, pattern = "_tabix.txt")
                                        
                                        # generating the TABIX compressed file
    for(filetb in Tabix_files) {
        filetb_name <- unlist(strsplit(filetb,split=".txt"))
        str <- paste('cat', paste0(files_location, "/", filetb), '>', 
                     paste0(output_folder, "/", filetb_name,"_out.txt"))
        system(str)
        str <- paste0(tabixPath, '/bgzip ', output_folder, "/", filetb_name,"_out.txt")
        system(str)
                                        # generating the TABIX index file
        fileoutgz <- paste0(output_folder, "/", filetb_name,"_out.txt", '.gz')
        str <- paste(tabixPath, '/tabix -s 1 -b 2 -e 2 -f ', fileoutgz, sep='')
        system(str)
    }
}


                                        # splitting chromosomes in N Mb regions ..
splitChrs <- function(chrs, org) {
    if(!is(org,"BSgenome") && !is(org,"list"))
        stop('org has to be either a list or an object of class BSgenome ..')
    chrs_org <- seqnames(org)
    if(any(!chrs %in% chrs_org))
        stop('This is not a valid chromosome .. ')
    chrL <- seqlengths(org)[chrs]
    splitDf <- data.frame(chrs, 1, chrL, row.names = NULL, stringsAsFactors=FALSE)
    colnames(splitDf) <- c('chr','start','end')
    return(splitDf)
}

                                        # fisher's method to determine a Pvalue based on a set of Pvalues
chiCombP <- function(Pvalues) {
    chistat <- -2*(sum(log(Pvalues)))
    chiP <- 1-pchisq(chistat, df=2*length(Pvalues))
    return(chiP)
}


                                        # consolidating DMRs ..

consolidateDMRs <- function(DmrGR, pvThr=0.05, MethDiff_Thr=NULL, log2Er_Thr=NULL,
                            GAP=0, type=NULL, correct=FALSE){
    if(!is(DmrGR, "GRanges"))
        stop('DmrGR has to be of class GRanges ..')
    if(!is.numeric(pvThr))
        stop('pvThr has to be of class numeric ..')
    if(!is.null(MethDiff_Thr) && !is.numeric(MethDiff_Thr))
        stop('MethDiff_Thr has to be either NULL or of class numeric ..')
    if(!is.null(log2Er_Thr) && !is.numeric(log2Er_Thr))
      stop('log2Er_Thr has to be either NULL or of class numeric ..')
    if(!is.numeric(GAP))
        stop('GAP has to be of class numeric ..')
    if(!is.null(MethDiff_Thr) && any(!type %in% c("hyper","hypo")))
        stop('type has to be either NULL or one of hyper and hypo ..')
    if(!is.logical(correct))
        stop('correct has to be of class logical ..')

    if(correct) mcols(DmrGR)$pValue= p.adjust(mcols(DmrGR)$pValue, method='BH')
    sInds <- which(mcols(DmrGR)$pValue <= pvThr)
    if(!is.null(MethDiff_Thr)) {
                                        # selecting based on methylation diff
        mInds <- which(abs(mcols(DmrGR)$MethDiff_Perc) >= MethDiff_Thr)
        sInds <- IRanges::intersect(sInds, mInds)
    }
    if(!is.null(log2Er_Thr)) {
                                        # selecting based on log2enrichment
        eInds <- which(abs(mcols(DmrGR)$log2Enrichment) >= log2Er_Thr)
        sInds <- IRanges::intersect(sInds, eInds)
    }
    if(length(sInds) == 0) return(NULL)
    DmrGR <- DmrGR[sInds,]
    if(!is.null(type))
    {
      if(type=="hypo")
        DmrGR <- DmrGR[mcols(DmrGR)$MethDiff_Perc < 0]
      else
        DmrGR <- DmrGR[mcols(DmrGR)$MethDiff_Perc > 0]
    }

                                        # building a GRanges for selected DMRs

    joinDMR <- function(query, Gap=0, fun=chiCombP) {
        if(!is.numeric(Gap))
            stop('Gap has to be of class numeric .. ')
        if(class(fun) != 'function')
            stop('fun has to be of class function .. ')
        commonChrs <- unique(as.character(seqnames(query)))
        resIR_gr <- GRanges()
        for(chrom in commonChrs) {
            queryF <- query[seqnames(query) == chrom]
            Pv <- mcols(queryF)$pValue
            MethDiff <-  mcols(queryF)$MethDiff_Perc
            Enrichment <-  mcols(queryF)$log2Enrichment
                                        # extracting IRanges from the dataframe
            IR1 <- IRanges(start(queryF), end(queryF))
            IR2 <- IR1
            if(Gap > 0) {
                                        # Gap are added on each side of the genomic regions
                                        # they will be taken back after the intersection
                start(IR1) <- start(IR1) - Gap
                end(IR1) <- end(IR1) + Gap
                start(IR2) <- start(IR2) - Gap
                end(IR2) <- end(IR2) + Gap
            }
            iIR <- IRanges::intersect(IR1, IR2)
            if(length(iIR) == 0) next
            if(Gap > 0) {
                start(iIR) <- start(iIR) + Gap
                end(iIR) <- end(iIR) - Gap
            }
            combPv <- rep(NA, length(iIR))
            combMethDiff <- rep(NA, length(iIR))
            combEnrichment <- rep(NA, length(iIR))
                                        # determining which IR elements were combined in iIR
            result1 <- findOverlaps(iIR, IR1)
            res1 <- cbind(queryHits(result1), subjectHits(result1))
            result2 <- findOverlaps(iIR, IR2)
            res2 <- cbind(queryHits(result1), subjectHits(result1))
            for(uind in unique(res2[,1])) {
                inds <- which(res2[,1] == uind)
                IR2ind <- res2[inds,2]
                Pvs <- Pv[IR2ind]
                Meth <- MethDiff[IR2ind]
                Enr <- Enrichment[IR2ind]
                combPv[uind] <- do.call('fun', list(Pvs))
                combMethDiff[uind] <- do.call(mean, list(Meth))
                combEnrichment[uind] <- do.call(mean, list(Enr))
            }
            suppressWarnings(resIR_gr <- append(resIR_gr, GRanges(Rle(chrom), IRanges(start(iIR), end(iIR)),
                                                                  pValue=round(combPv,3), MethDiff_Perc=round(combMethDiff,3),
                                                                  log2Enrichment=round(combEnrichment,3))))
        }
                                        # appending the new iIR to those generated for the other chromosomes
        zeroWinds <- which(width(resIR_gr) == 1)
        if(length(zeroWinds) > 0) resIR_gr <- resIR_gr[-zeroWinds,]
        resIR_gr
    }

    dmrGRanges <- joinDMR(DmrGR, Gap=GAP)
    return(dmrGRanges)
}


                                        # retrieve mC calls for genomic regions given a BSdata object for a sample

mapBSdata2GRanges <- function(GenoRanges, Sample, context='all', mC=1, depth=0,
                              pValue=1) {
    if(!is(GenoRanges, "GRanges"))
        stop('GenoRanges has to be of class GRanges ..')
    if(!is(Sample, "BSdata"))
        stop('Sample has to be of class BSdata ..')
    if(length(which(!(context %in% c('all','CG','CHG','CHH')))) > 0)
        stop('context has to be either all or a combination of CG, CHG, and CHH ..')
    if(!is.numeric(mC) && mC < 0 )
        stop('mC has to be of class numeric and positive..')
    if(!is.numeric(depth))
        stop('depth has to be of class numeric ..')
    if(!is.numeric(pValue))
        stop('pValue has to be of class numeric ..')
    if(pValue > 1)
        stop('pValue has to be lower than 1 ..')


                                        # querying TABIX indexed files based on the tabix seqnames,
                                        # assigning NA to GRanges regions in chrs not represented
    tabixChrs <- seqnamesTabix(Sample@file)
    representedChr <- which(as.character(seqnames(GenoRanges)) %in% tabixChrs)
    if(length(representedChr) == 0)
      {
      res <- as.list(rep(NA, length(GenoRanges))) 
      return(res)
      }
    otherChr <- which(!(as.character(seqnames(GenoRanges)) %in% tabixChrs))
    if(length(otherChr) > 0) res <- as.list(rep(NA, length(GenoRanges)))

    regions <- GenoRanges[representedChr,]
    resScan <- scanTabix(Sample@file, param=regions)
    resScan <- lapply(resScan, tabixdata2GR)
    if(length(otherChr) > 0) res[representedChr] <- resScan
    else res <- resScan
                                        # filtering results
    filterFun <- function(gr, context, mC, depth, pValue) {
        indsList <- list()
        if(context[1] != 'all') indsList$cl <- which(mcols(gr)$Context %in% context)
        if(mC > 0) indsList$mC <- which(mcols(gr)$C > mC)
        if(depth > 0) indsList$d <- which((mcols(gr)$C+mcols(gr)$T) > depth)
        pV <- -10*log10(pValue)
        if(pValue < 1) indsList$p <- which(mcols(gr)$significance > pV)
        tableInds <- table(unlist(indsList))
        inds <- names(tableInds)[tableInds == length(indsList)]
        if(is.null(inds) || length(inds) == 0) return(NA)
        return(gr[as.numeric(inds),])
    }


    if(context != 'all' || mC > 1 || depth > 0 || pValue < 1) {
        NAinds <- as.numeric(which(is.na(res)))
        if(length(NAinds) > 0) {
            res[-NAinds] <- lapply(res[-NAinds], filterFun, context, mC, depth, pValue)

        }
        else res <- lapply(res, filterFun, context, mC, depth, pValue)
    }
    return(res)
}


                                        # for  GenomicRanges with N bins, extract the information
                                        # of the Genomic ranges for a given bin

extractBinGRanges <- function(GenoRanges, bin, nbins) {
    if(!is(GenoRanges, "GRanges"))
        stop('GenoRanges has to be of class GRanges ..')
    if(!is.numeric(bin) || length(bin) != 1 || is.na(bin))
        stop(' bin has to be a single non-NA integer ')
    if(!is.numeric(nbins) || length(nbins) != 1 || is.na(nbins))
        stop(' nbins has to be a single non-NA integer ')
    if( bin < 1 || bin > nbins)
        stop(' bin has to be between 1 and nbins ')
    binsize <- round(width(GenoRanges)/nbins)
    starts <- start(GenoRanges) + (bin-1)*binsize
    gnew <- GRanges(as.character(seqnames(GenoRanges)),
                    IRanges(starts, width=binsize),
                    strand=as.character(strand(GenoRanges)))
    return(gnew)
}


                                        # As in mapBSdata2GRanges, but returns mC calls
                                        # within each bin of each Genomic region

mapBSdata2GRangesBin <- function(GenoRanges, Sample,
                                 context='all', mC=1, depth=0, pValue=1, nbins){
    if(!is(GenoRanges, "GRanges"))
        stop('GenoRanges has to be of class GRanges ..')
    if(!is(Sample,'BSdata'))
        stop('Sample has to be of class BSdata ..')
    if(length(which(!(context %in% c('all','CG','CHG','CHH')))) > 0)
        stop('context has to be either all or a combination of CG, CHG, and CHH ..')
    if(!is.numeric(mC) && mC < 0 )
        stop('mC has to be of class numeric and positive..')
    if(!is.numeric(depth))
        stop('depth has to be of class numeric ..')
    if(!is.numeric(pValue))
        stop('pValue has to be of class numeric ..')
    if(pValue > 1)
        stop('pValue has to be lower than 1 ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')
    resList <- list()
    if(nbins == 1)
        resList[[1]] <- mapBSdata2GRanges(GenoRanges, Sample, context)
    else
        {
            for(bin in 1:nbins) {
                gel <- extractBinGRanges(GenoRanges=GenoRanges, bin=bin, nbins=nbins)
                resList[[bin]] <- mapBSdata2GRanges(GenoRanges=gel, Sample=Sample,
                                                    context, mC, depth, pValue)
          }
        }
    resrev <- list()
    for(GEind in 1:length(GenoRanges)) resrev[[GEind]] <-
        lapply(resList, function(x) x[[GEind]])
    return(resrev)
}


                                        # get genomic Cxx positons for a series of genomic regions
getCpos <- function(GenoRanges, seqContext='all', nbins, org) {
    if(!is(GenoRanges, "GRanges"))
         stop('GenoRanges has to be of class GRanges ..')
    if(!(seqContext %in% c('all','CG','CHG','CHH')))
        stop('seqContext has to be one of all, CG, CHG or CHH ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')
    if(!is(org,"BSgenome"))
        stop('org has to be of class BSgenome ..')
    chrs <- as.character(seqnames(GenoRanges))
    if(length(which(!(unique(chrs) %in% names(org)))) > 0)
        stop('Object covers chromosomes not available in org ..')

    resList <- as.list(rep(NA, length(GenoRanges)))
    Nbins <- nbins
    for(Chr in unique(chrs)) {
        chrinds <- which(as.character(seqnames(GenoRanges)) == Chr)
        GenoRangesChr <- GenoRanges[seqnames(GenoRanges) == Chr]
        chrseq <- unmasked(org[[Chr]])
        if(seqContext != 'all') resList[chrinds] <-
            getCposChr(GenoRanges= GenoRangesChr, seqContext=seqContext,
                       chrseq=chrseq, Nbins)
        else {
            resCG <- resList[chrinds] <- getCposChr(GenoRanges= GenoRangesChr,
                                                    seqContext='CG', chrseq=chrseq, Nbins)
            resCHG <- resList[chrinds] <- getCposChr(GenoRanges= GenoRangesChr,
                                                     seqContext='CHG', chrseq=chrseq, Nbins)
            resCHH <- resList[chrinds] <- getCposChr(GenoRanges= GenoRangesChr,
                                                     seqContext='CHH', chrseq=chrseq, Nbins)
            resALL <- resCG
            for(i in 1:length(resALL)) {
                if(Nbins > 1) { for(bin in 1:Nbins)
                                    resALL[[i]][[bin]] <- sort(c(unlist(resCG[[i]][[bin]]),
                                                                 unlist(resCHG[[i]][[bin]]), unlist(resCHH[[i]][[bin]]))) }
                else resALL[[i]][[1]] <- sort(c(unlist(resCG[[i]]),
                                                unlist(resCHG[[i]]), unlist(resCHH[[i]])))
            }
            resList <- resALL
        }
    }
    return(resList)
}


                                        # get genomic Cxx positons for a series of genomic regions on a given chr seq

getCposChr <- function(GenoRanges, seqContext, chrseq, nbins) {
    if(length(unique(as.character(seqnames(GenoRanges)))) > 1)
        stop('GenoRanges maps to several chromosomes ..')
    if(!(seqContext %in% c('CG','CHG','CHH')))
        stop('seqContext has to be one of CG, CHG or CHH ..')
    if(!is(chrseq,"DNAString"))
        stop('chrseq has to be of class DNAString (unmasked) ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')
    typeD <- DNAString(seqContext)
    typeDr <- reverse(DNAString(seqContext))
    Nbins <- nbins

    chrv <- Views(chrseq, start=start(GenoRanges), end=end(GenoRanges))
    vls <- as(chrv,'DNAStringSet')
    vlsC <- Biostrings::complement(vls)
    Cpos <- startIndex(vmatchPattern(typeD, vls, fixed='subject'))
    CposC <- endIndex(vmatchPattern(typeDr, vlsC, fixed='subject'))
    for(i in 1:length(Cpos)) {
        res <- c(Cpos[[i]], CposC[[i]])
        if(!is.null(res)) Cpos[[i]] <- sort(res)
    }

                                        # restore bin structure (Nbins entry for each genomic region)
    CposList <- list()
    for(i in 1:length(Cpos))
        {
            binList <- list()
            startPos <- 1
            binsize <- round(width(GenoRanges)[i])/Nbins
            endPos <- binsize
            for(bin in 1:Nbins)
                {
                    ind <- which(Cpos[[i]] > startPos & Cpos[[i]] < endPos)
                    binList[[bin]] <- Cpos[[i]][ind]
                    startPos <- endPos+1
                    endPos <- endPos+binsize
                }
            if(length(binList) == 0)
                CposList[[i]] <- NA
            else
                CposList[[i]] <- binList
        }
    return(CposList)
}


getCposDensity<- function(GenoRanges, Cpos, nbins) {
    allwidths <- width(GenoRanges)
    Nbins <- nbins
    for(i in 1:length(Cpos)) {
        widths <- rep(round(allwidths[i]/Nbins), Nbins)
        if(is.na(Cpos[[i]]))
          Cpos[[i]] <- rep(0,nbins)
        else
          Cpos[[i]] <- sapply(Cpos[[i]], length)/widths
    }
    return(Cpos)
}

                                        # assign binmC, binC, binrC for GEcollection using Tabix based BSdata ...

profileDNAmetBin <- function(GenoRanges, Sample,
                             mcCLASS='mCG',
                             mC=1, depthThr=0, mCpv=1,
                             minCoverage=0.75, nbins=2) {
    if(!is(GenoRanges, "GRanges"))
        stop('GenoRanges has to be of class GRanges ..')
    if(!is(Sample,"BSdata"))
        stop('Sample has to be of class BSdata ..')
    if(length(which(!(mcCLASS %in% c('mCG','mCHG','mCHH')))) > 0)
        stop('mcCLASS has to be one of mCG, mCHG, and mCHH ..')
    if(!is.numeric(mC) && mC < 0 )
        stop('mC has to be of class numeric and positive..')
    if(!is.numeric(depthThr))
        stop('depthThr has to be of class numeric ..')
    if(!is.numeric(mCpv))
        stop('mCpv has to be of class numeric ..')
    if(!is.null(minCoverage) && !is.numeric(minCoverage))
        stop('minCoverage has to be either NULL or of class numeric ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')

    Nbins <- nbins
                                        # extracting the mC methylation level 2sec
    bsdat <- mapBSdata2GRangesBin(GenoRanges=GenoRanges, Sample= Sample,
                                  context=sub('m','',mcCLASS), mC=mC,
                                  depth=depthThr, pValue=mCpv, nbins=Nbins)
                                        # binmC
    sumFun <- function(x) {
        suppressWarnings(if(is.na(x)[[1]]) return(NA))
        sum(mcols(x)$C/(mcols(x)$C+mcols(x)$T))
    }
    mlsum <- lapply(bsdat, sapply, sumFun)
    binmC <- matrix(unlist(mlsum), length(mlsum), Nbins, byrow=T)
    unmethList <- list()
    uncov <- Sample@uncov
    for(bin in 1:Nbins) {
      gel <- extractBinGRanges(GenoRanges=GenoRanges, bin=bin, nbins=Nbins)
      ov <- findOverlaps(gel,uncov)
      if(length(ov) > 0)
      {
        val <- rep(0,length(GenoRanges))
        ind <- unique(queryHits(ov))
        val[ind] <- NA
        unmethList[[bin]] <- val 
      }
      else 
        unmethList[[bin]] <- 0
    }
    
    unmethList <- matrix(unlist(unmethList), Nbins, length(GenoRanges), byrow=T)
    unmethList <- t(unmethList)
    ind <- which(is.na(binmC))
    binmC[ind] <- unmethList[ind]
    
    wbp <- round(width(GenoRanges)/Nbins)
    binmC <- apply(binmC, 2, function(x) x/wbp)

                                        # binC
    Cpos <- getCpos(GenoRanges= GenoRanges,
                    seqContext= sub('m','', mcCLASS),
                    org= Sample@org, nbins=Nbins)
    suppressWarnings(CposD <- getCposDensity(GenoRanges, Cpos=Cpos, nbins=Nbins))
    binC <- matrix(unlist(CposD), length(CposD), Nbins, byrow=T)

                                        # reversing for minus strand if strand is defined
    minusStrand <- which(as.character(strand(GenoRanges)) == '-')
    if(length(minusStrand) > 0) {
        binmC[minusStrand,] <- binmC[minusStrand, Nbins:1]
        binC[minusStrand,] <- binC[minusStrand, Nbins:1]
    }

                                        # adding relative methylation (C/mC) expressed on [0,100] range

    if(is.null(ncol(binmC)))
        binmC <- matrix(binmC, ncol=length(binmC))

    binrC <- binmC
    for(i in 1:ncol(binrC)) binrC[,i] <- 100*binmC[,i]/binC[,i]
    binrC[is.infinite(binrC)] <- NA

                                        # adding uncovered region information
    maskUncovered <- Sample@uncov
    seqlengths(maskUncovered) <- NA
    seqlengths(GenoRanges) <- NA
    suppressWarnings(nonCoveredInds <- findOverlaps(GenoRanges, maskUncovered))
    orgInds <- queryHits(nonCoveredInds)
    binrC[orgInds,] <- NA

                                        # cutting some decimals and specifying row/column names for the slots
    row_names <- paste(as.character(seqnames(GenoRanges)),
                       ":", start(GenoRanges),"-", end(GenoRanges), sep="")
    col_names <- c(1:Nbins)
    binmC <- signif(binmC, 3)
    colnames(binmC) <- col_names
    binC <- signif(binC, 3)
    colnames(binC) <- col_names
    binrC <- signif(binrC, 3)
    colnames(binrC) <- col_names
    binscore <- matrix(NA, nrow(binmC),
                       ncol(binmC), dimnames=list(row_names, col_names))
                                        # saving chr results
    Object <- new("GEcollection",
                  SummarizedExperiment(assays=list(binmC=binmC,
                                           binC=binC, binrC=binrC,
                                           binscore=binscore),
                                       rowData=GenoRanges))
    rownames(Object) <- row_names
    Object
}


profileDNAmetBinParallel <- function(GenoRanges, Sample, mcCLASS='mCG',
                                     mC=1, depthThr=0, mCpv=1,
                                     minCoverage=0.75, Nproc=1, nbins=2){

    if(!is(GenoRanges, "GRanges"))
        stop('GenoRanges has to be of class GRanges ..')
    if(!is(Sample,"BSdata"))
        stop('Sample has to be of class BSdata ..')
    if(length(which(!(mcCLASS %in% c('mCG','mCHG','mCHH')))) > 0)
        stop('mcCLASS has to be one of mCG, mCHG, and mCHH ..')
    if(!is.numeric(mC) && mC < 0 )
        stop('mC has to be of class numeric and positive..')
    if(!is.numeric(depthThr))
        stop('depthThr has to be of class numeric ..')
    if(!is.numeric(mCpv))
        stop('mCpv has to be of class numeric ..')
    if(!is.null(minCoverage) && !is.numeric(minCoverage))
        stop('minCoverage has to be either NULL or of class numeric ..')
    if(!is.numeric(Nproc) && Nproc < 1)
        stop('Nproc has to be of class numeric ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')

    Chrs <- unique(as.character(seqnames(GenoRanges)))
    tabixChrs <- seqnamesTabix(Sample@file)
    otherChr <- which(!(as.character(seqnames(GenoRanges)) %in% tabixChrs))
    otherChr <- which(!(as.character(seqnames(GenoRanges)) %in% tabixChrs))
    if(length(otherChr) > 0) GenoRanges <- GenoRanges[-c(otherChr)]
    Chrs <- unique(as.character(seqnames(GenoRanges)))
    Nproc <- min(Nproc, length(Chrs))
    cl <- makeCluster(Nproc, 'PSOCK')
                                        # a load balanced parallel execution of the profileDNAmetBin method is run
    profileDNAmetBinChr <- function(CHR, GenoRanges,
                                    Sample, mcCLASS, mC,
                                    depthThr, mCpv,
                                    minCoverage, nbins) {
        GenoRangesChr <- GenoRanges[seqnames(GenoRanges) == CHR]
        res <- profileDNAmetBin(GenoRanges=GenoRangesChr, Sample=Sample,
                                mcCLASS=mcCLASS, mC=mC, depthThr=depthThr,
                                mCpv=mCpv, minCoverage=minCoverage, nbins=nbins)
        res
    }

    clRes <- clusterApplyLB(cl, Chrs, profileDNAmetBinChr,
                            GenoRanges=GenoRanges, Sample=Sample,
                            mcCLASS=mcCLASS, mC=mC, depthThr=depthThr,
                            mCpv=mCpv, minCoverage=minCoverage, nbins=nbins)
    chrcomb <- NULL
    chrcomb <- do.call("rbind", clRes)
    Object <- new("GEcollection", SummarizedExperiment(assays=assays(chrcomb),
                                                       rowData(chrcomb)))
    Object
}



                                  # function to plot gene locus visualization with methylation/omics data

plotMeth <- function(grl, colors=NULL, datatype, yLim, brmeth=NULL, mcContext="CG", annodata=NULL, Datatrackname,
                     transcriptDB, chr, start=NULL, end=NULL, org){
    if(!is(grl,'list') && !is(grl,'GElist'))
        stop('grl has to be of class list or GElist...')
    if(is(grl,'list'))
        {
            if(any(!(sapply(grl, class) %in% c("GRanges", 'GEcollection'))))
                stop('grl has to be a list of either GRanges or GEcollection objects...')
        }
    if(!is.null(colors) && !is.character(colors))
      stop('colors has to be either NULL or of class character ...')
    if(is.character(colors))
    {
      if(length(colors) != length(grl))
        stop('length of colors has to be same as length of grl ...')
    }
    for(i in 1:length(datatype)) {
      dti <- datatype[i]
      if(!(dti %in% c('density', 'C','mC', 'rC','cols', 'gr')))
        stop('datatype has to be an array containing
           one of: density, C, mC, rC, cols or gr')
    }
    if(!is.numeric(yLim))
      stop('yLim has to be of class numeric ..')
    if(!is.null(brmeth) && !is(brmeth,'BSdataSet'))
      stop('brmeth has to be either NULL or of class BSdataSet ...')
    if(length(which(!(mcContext %in% c('CG','CHG','CHH')))) > 0)
      stop('mcContext has to be one of CG, CHG, and CHH ..')
    if(!is.null(annodata) && !is(annodata,'GRangesList'))
      stop('annodata has to be either NULL or of class GRangesList ...')
    if(!is.character(Datatrackname))
      stop('Datatrackname has to be of class character ...')
    if(!is.null(annodata) && !is.null(brmeth))
    {
      if(length(Datatrackname) != length(grl)+length(annodata)+length(brmeth))
        stop('Datatrackname has to be equal to length of grl, brmeth and annodata...')
      grl_trackname <- Datatrackname[1:length(grl)]
      brmeth_trackname <- Datatrackname[length(grl)+1:length(brmeth)]
      anno_trackname <- Datatrackname[length(grl)+length(brmeth)+1:length(annodata)]
    }
    else
    {
      if(!is.null(brmeth))
      {
        if(length(Datatrackname) != length(grl)+length(brmeth))
          stop('Datatrackname has to be equal to length of grl and brmeth')
        grl_trackname <- Datatrackname[1:length(grl)]
        brmeth_trackname <- Datatrackname[length(grl)+1:length(brmeth)]
      }
      if(!is.null(annodata))
      {
        if(length(Datatrackname) != length(grl)+length(annodata))
          stop('Datatrackname has to be equal to length of grl and annodata...')
        grl_trackname <- Datatrackname[1:length(grl)]
        anno_trackname <- Datatrackname[length(grl)+1:length(annodata)]
      }
      if(is.null(annodata) && is.null(brmeth))
      {
        if(length(Datatrackname) != length(grl))
          stop('Datatrackname has to be equal to length of grl...')
        grl_trackname <- Datatrackname[1:length(grl)]
      }
    }

    if(length(datatype) != length(grl) && length(grl) != length(yLim))
      stop('grl, datatype and yLim has to be of same length ..')
    if(!is(transcriptDB,"TxDb") && !is.null(transcriptDB))
        stop('transcriptDB has to be either NULL or an object of class TxDb ..')
    if(!is.character(chr))
        stop('chr has to be of class character ..')
    if(!is.null(start) && !is.numeric(start))
      stop('start has to be of either NULL or of class numeric ..')
    if(!is.null(end) && !is.numeric(end))
        stop('end has to be of either NULL or of class numeric ..')
    if(!is(org, "BSgenome"))
      stop('org has to be of class BSgenome ..')

    gen <- org@provider_version
    itrack <- IdeogramTrack(genome = gen, chromosome = chr)
    axisTrack <- GenomeAxisTrack()
    matlist <- list()
    dTrack <- list()
    for (i in 1:length(grl))
        {
            if(datatype[i] == 'mC') {
                refgr <- rowData(grl[[i]])
                matlist[[i]] <- t(binmC(grl[[i]]))
            }
            if(datatype[i] == 'C') {
                refgr <- rowData(grl[[i]])
                matlist[[i]] <- t(binC(grl[[i]]))
            }
            if(datatype[i] == 'rC') {
                refgr <- rowData(grl[[i]])
                matlist[[i]] <- t(binrC(grl[[i]]))
            }
            if(datatype[i] == 'density') {
                refgr <- rowData(grl[[i]])
                matlist[[i]] <- t(binscore(grl[[i]]))
            }
            if(datatype[i] == 'cols') {
                refgr <- rowData(grl[[i]])
                matlist[[i]] <- t(as.matrix(mcols(grl[[i]])[1]))
            }
            if(datatype[i] == 'gr') {
              refgr <- grl[[i]]
              matlist[[i]] <- t(as.matrix(mcols(grl[[i]])[1]))
            }
            if(is.null(colors))
              dTrack[[i]] <- DataTrack(range=refgr, data=matlist[[i]] , name=grl_trackname[i], chromosome=chr,
                                     fill="darkgreen", ylim=c(0, yLim[i]))
            else
              dTrack[[i]] <- DataTrack(range=refgr, data=matlist[[i]] , name=grl_trackname[i], chromosome=chr,
                                       fill=colors[i], ylim=c(0, yLim[i]))
        }

    if(!is.null(start) && !is.null(end))
    {
      brmeth_gr <- GRanges(chr,IRanges(start,end))
      if(abs(start-end) < 50)
      {
        strack <- SequenceTrack(org)
        anntrack <- list(itrack, axisTrack, strack)
      }
      else
        anntrack <- list(itrack, axisTrack)
    }
    else
    {
      start <- min(start(refgr[seqnames(refgr)==chr,]))
      end <- max(end(refgr[seqnames(refgr)==chr,]))
      anntrack <- list(itrack, axisTrack)
      brmeth_gr <- GRanges(chr,IRanges(start,end))
    }


    brTrack <- list()
    bsdat <- NULL
    brlist <- list()
    plotype_length <- length(grl)
    if(!is.null(brmeth))
    {
      for (i in 1:length(brmeth))
      {
        bsdat[i] <- unlist(mapBSdata2GRanges(GenoRanges=brmeth_gr, Sample= brmeth[[i]],
                                      context=mcContext))
        mc <- mcols(bsdat[[i]])$C/(mcols(bsdat[[i]])$C+ mcols(bsdat[[i]])$T)
        mc <- as.matrix(mc)
        temp_gr <- bsdat[[i]]
        strand(temp_gr) <- "*"
        brTrack[[i]] <- DataTrack(range=temp_gr, data=t(mc) , name=brmeth_trackname[i], chromosome=chr,
                                 col.histogram="black", ylim=c(0, 1))
      }
      plotype_length <- plotype_length+length(brmeth)
    }

    plottype <- rep("histogram", plotype_length)
#### adding annotation track

    AnnoTrack <- list()
    if(!is.null(annodata))
    {
      col <- rainbow(length(annodata))
      for (i in 1:length(annodata))
      {
        AnnoTrack[[i]] <- AnnotationTrack(annodata[[i]], name=anno_trackname[i], fill=col[i])
      }
    }

    if(!is.null(transcriptDB))
     {
    	txdb <- transcriptDB
    	txTrack <- GeneRegionTrack(txdb, chromosome=chr, name="Transcripts", showId=TRUE,
    	                           geneSymbol=TRUE, background.panel = "#FFFEDB", background.title = "brown")
    	combtrack <- c(anntrack, dTrack, brTrack, AnnoTrack, txTrack)
     }
    else combtrack <- c(anntrack, dTrack, brTrack, AnnoTrack)

    plotTracks(combtrack, type=plottype, chromosome=chr, from=start, to=end,
               reverseStrand = FALSE, background.panel = "#FFFEDB", background.title = "brown", col="black")
}
