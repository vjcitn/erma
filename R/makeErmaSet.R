setClass("ErmaSet", contains="GenomicFiles")

setGeneric("cellTypes", function(x, ...) standardGeneric("cellTypes"))
setMethod("cellTypes", "ErmaSet", function(x, ...)
   colData(x)[,15])

makeErmaSet = function() {
# fonly = dir(system.file("bed_tabix", package="erma")) # no path
 ermaset = 
  GenomicFiles(  
    files=dir(system.file("bed_tabix", package="erma"), 
                full.names=TRUE, pattern="gz$")
  )
 mm = mapmeta()
 rownames(mm) = mm[,2]
 ftags = gsub("_.*", "", basename(files(ermaset)))
# cd = as(mm[ftags,], "DataFrame")
# rownames(cd) = ftags
 ermaset@colData = mm[ftags,]
 new("ErmaSet", ermaset)
}

setMethod("show", "ErmaSet", function(object) {
    callNextMethod()
    cat("cellTypes() for type names; data(short_celltype) for abbr.\n")
})


