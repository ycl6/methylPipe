                                        # chr
setGeneric('chr', function(object) standardGeneric('chr'))
setMethod('chr','GEcollection', function(object)
          as.character(seqnames(object)))
                                        # binmC
setGeneric('binmC', function(object)
           standardGeneric('binmC'))
setMethod('binmC', 'GEcollection', function(object)
          assays(object)[["binmC"]])
                                        # binC
setGeneric('binC', function(object)
           standardGeneric('binC'))
setMethod('binC', 'GEcollection', function(object)
          assays(object)[["binC"]])

                                        # binrC
setGeneric('binrC', function(object)
           standardGeneric('binrC'))
setMethod('binrC', 'GEcollection', function(object)
          assays(object)[["binrC"]])

                                        # binscore
setGeneric('binscore', function(object)
           standardGeneric('binscore'))
setMethod('binscore', 'GEcollection', function(object)
          assays(object)[["binscore"]])

                                        # set binscore
setGeneric('binscore<-', function(object, value)
           standardGeneric('binscore<-'))
setMethod('binscore<-', 'GEcollection', function(object, value) {
    if(!is.numeric(value))
        stop("The binscore values has to be numeric")
    assays(object)[["binscore"]] <- value
    object
})

                                        # nbins
setGeneric('nbins', function(object)
           standardGeneric('nbins'))
setMethod('nbins', 'GEcollection', function(object)
          ncol(object))

                                        # length
setMethod('length', 'GEcollection', function(x)
          length(rowData(x)))

                                        # show
setMethod('show', 'GEcollection', function(object) {
    cat("S4 Object of class GEcollection; ")
    message()
    print(rowData(object))
    cat("\nbinmC : ")
    if(!is.null(assays(object)[["binmC"]])) cat('ok')
    else cat('NA')
    cat("\nbinC : ")
    if(!is.null(assays(object)[["binC"]])) cat('ok')
    else cat('NA')
    cat("\nbinrC : ")
    if(!is.null(assays(object)[["binrC"]])) cat('ok')
    else cat('NA')
    cat("\nbinscore : ")
    if(length(which(is.na(assays(object)[["binscore"]])))
       != length(assays(object)[["binscore"]])) cat('ok')
    else cat('NA')
    cat('\n')
})

