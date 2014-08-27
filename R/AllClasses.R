setClass(Class = 'BSdata',
         slots = c(file = 'character',
             uncov = 'GRanges', org = 'BSgenome'))
setValidity("BSdata", function(object) {
    tabixRef <- TabixFile(object@file, yieldSize=1000L)
    tbdata <- scanTabix(tabixRef)
    tbdata <- unlist(tbdata)
    tbdata <- tabixdata2GR(tbdata)
    chrs <- seqnames(object@org)
    if(length(which( !(as.character(seqnames(tbdata)) %in% chrs))) > 0)
        return("the 1st column of the BSdata has to contain
chr assignments in the form chr1, chr2, ...
included in the org reference sequence")
    if(length(which( !(mcols(tbdata)$Context %in% c('CG','CHG', 'CHH')))) > 0)
        return("the 4th column of the BSdata
         has to contain one of CG, CHG or CHH ..")
    if(!is.numeric(mcols(tbdata)$C) || length(which(mcols(tbdata)$C < 1)) > 0)
        return("the 5th column of the BSdata
         has to contain only positive integers >=1 ..")
    if(!is.numeric(mcols(tbdata)$T) || length(which(mcols(tbdata)$T < 0)) > 0)
        return("the 6th column of the BSdata
         has to contain only positive integers ..")
    if(!is.numeric(mcols(tbdata)$Significance) ||
       length(which(mcols(tbdata)$Significance < 0)) > 0)
        return("the 7th column of the BSdata
         has to contain only positive integers ..")
    NULL
})



setClass(Class = 'BSdataSet', slots = c(Objlist='list', names='character',
                                  group='character', org='BSgenome'))
setValidity("BSdataSet", function(object) {
    if(length(object@Objlist) != length(object@group))
        return("length of group has to be equal to number of objects")
    if(length(object@Objlist) == 1)
        return("more then one BSdata object has to be provided")
    if(!is(object@names,"character"))
        return("names has to be of class character")
    if(!is(object@group,"character"))
        return("group has to be of class character")
    if(length(object@Objlist)==2)
    {
      if( any( !object@group %in% c("C","E")))
        return("if only two groups they can only be 'C' or 'E'")
    }
    for(i in 1:length(object@Objlist))
        if(!is(object@Objlist[[i]],"BSdata"))
            return("the object has to be of class BSdata")
    for(i in 1:length(object@Objlist))
        {
            object@Objlist[[i]]@file <- (object@Objlist[[i]])@file
            object@Objlist[[i]]@uncov <- (object@Objlist[[i]])@uncov
        }
    NULL
})


setClass(Class = 'GEcollection', contains = "SummarizedExperiment",
         prototype = prototype(SummarizedExperiment(assays = SimpleList(binmC = matrix(0, 0, 0),
                                                        binC = matrix(0, 0, 0),
                                                        binrC = matrix(0, 0, 0),
                                                        binscore = matrix(0, 0, 0)))))
setValidity("GEcollection", function(object) {
    requiredAssays <- c("binmC", "binC", "binrC", "binscore")
    if ( !identical(requiredAssays, c(names(assays(object)))))
        return("assays slots have to in the order: binmC, binC, binrC, binscore")
    if ( !is.na(assays(object)[["binmC"]]) &&
        min(assays(object)[["binmC"]], na.rm = TRUE) < 0)
        return("The binmC values have to either NA or positive")
    if ( !is.na(assays(object)[["binC"]]) &&
        min(assays(object)[["binC"]], na.rm = TRUE) < 0)
        return("The binC values have to either NA or positive")
    if (!is.na(assays(object)[["binrC"]]) &&
        min(assays(object)[["binrC"]], na.rm = TRUE) < 0 &&
        max(assays(object)[["binrC"]], na.rm = TRUE) > 1)
        return("The binrC values has to be between 0 and 1")
    if (!is.na(assays(object)[["binscore"]]) &&
        min(assays(object)[["binscore"]], na.rm = TRUE) < 0)
        return("The binscore values have to either NA or positive")
    NULL
})


setClass(Class = 'GElist', slots = c(Objlist="list",
                               names="character"))
setValidity("GElist", function(object) {
    for(i in 1:length(object@Objlist))
        if(!is(object@Objlist[[i]],"GEcollection"))
            return("object has to be of class GEcollection")
    NULL
})


