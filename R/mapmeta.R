
mapmeta = function() {
  tab = read.csv(
   dir(system.file("metadata_csv",package="erma"), full=TRUE), 
         stringsAsFactors=FALSE)
 sexcol = grep("SEX", colnames(tab))
 colnames(tab)[sexcol] = "SEX"
 narrDF(fixNonASCII(DataFrame(tab[-c(1,2),])))
}

narrDFcols = function() {
   basic = c(2,4,6,17,18,21)
   cand = options()$narrDFcols
   if (is.null(cand)) return(basic)
   cand
}
 
setClass("narrDF", contains="DataFrame")
setMethod("show", "narrDF", function(object) {
cat("(showing narrow slice of ", nrow(object), "x", ncol(object), " DataFrame)\   ")
if (ncol(object)==95) callNextMethod(object[,narrDFcols()])
else callNextMethod(object)
cat("use colnames() for full set of metadata attributes.\n")
})

narrDF = function(x) new("narrDF", x)
