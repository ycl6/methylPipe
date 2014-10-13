                                        # subset
setMethod('[[', 'GElist', function(x, i, j = "missing"){
    GEl = x@Objlist[[i]]
    GEl
})

setMethod('[[<-', 'GElist', function(x, i, j = "missing", value) {
    x@Objlist[[i]] = value
    x
})

setMethod('[', 'GElist', function(x, inds, i='missing', j='missing',
                                  drop ="missing") {
    x@Objlist = x@Objlist[inds]
    x
})

                                        # length
setMethod('length', 'GElist', function(x) {
    length(x@Objlist)
})

setMethod("$", "GElist", function(x, name){
    x@Objlist[[name]]
})
