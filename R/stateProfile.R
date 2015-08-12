
#liberalImport = function(con, which, pad=1e5, genome, ...) {
 # intention is to import request plus padding so
 # that result can be trimmed to fill out the 
 # where region with any relevant metadata
# stopifnot(is(which, "GRanges"))
# tmp = import(con, which=which+pad, genome=genome, ...)
 #restrict(tmp, start=start(which), end=end(which),
 #   keep.all.ranges=TRUE)
# subsetByOverlaps(tmp, which)
#}

csProfile = function(ermaset, symbol, flanksize=10000, useShiny=FALSE,
    ctsize=10, shortCellType=TRUE, orient5p=TRUE) {
#
# chromatin state profile: 5p end of gene txstart to flanksize nt upstream
#
     if (!useShiny) stateProfile( ermaset=ermaset, 
             symbol=symbol, ctsize=ctsize, width=flanksize, 
             shortCellType=shortCellType, orient5p=orient5p )
     else stateProf(ermaset=ermaset, shortCellType=shortCellType, ctsize=ctsize)
}

stateProfile = function(ermaset, symbol="IL33", width=50000,
    ctsize=10, shortCellType=TRUE, orient5p = TRUE) {
   mod = try(genemodel(symbol))
   if (inherits(mod, "try-error")) stop("can't resolve symbol")
   uil = flank(resize(range(mod), 1), start=orient5p, width=width)
   
   ## ----bind----------------------------------------------------------------
   ermaset@rowRanges = uil
   ## ----getcss--------------------------------------------------------------
   csstates = lapply(reduceByFile(ermaset, MAP=function(range, file) {
     imp = liberalImport(file, which=range, genome=genome(range))
     seqlevels(imp) = seqlevels(range)
     imp$rgb = rgbByState(imp$name)
     imp
   }), "[[", 1) 
   tys = cellTypes(ermaset)  # need to label with cell types
   ## ---- annotate
   csstates = lapply(1:length(csstates), function(x) {
      csstates[[x]]$celltype = tys[x]
      csstates[[x]]
      })
   
   ## ----doviz, fig=TRUE-----------------------------------------------------
   cssgr = unlist(GRangesList(csstates))
   if (!exists("short_celltype")) data(short_celltype)
#   short_celltype = get(load(dir(system.file("data",package="erma"),full.names=TRUE, pattern="short")))
#   states_25 = get(load(dir(system.file("data",package="erma"),full.names=TRUE, pattern="states_25")))
    if (!exists("states_25")) data(states_25)
#   data(states_25)
   mycol = states_25$rgb
   names(mycol) = paste0(1:25, "_", states_25$MNEMONIC)
   cssd = as.data.frame(cssgr)
   chrn = as.character(seqnames(cssgr))[1]
   if (shortCellType) cssd$celltype = short_celltype[ cssd$celltype ]
   n2 = ggplot(cssd)
   n2 + geom_rect(aes(ymin=0, ymax=1, xmin=start, xmax=end, fill=name)) + 
       facet_grid(celltype~.) +
       theme(strip.text.y = element_text(size = ctsize, angle = 0)) +
        ylim(0,1) + scale_y_continuous(breaks=NULL, limits=c(0,1)) +
        scale_fill_manual( name="state", values=mycol ) + xlab(chrn) +
        ggtitle(paste0(symbol, " [", as.character(strand(mod)[1]), "] ", "(", floor(width/1000), "kb ", ifelse(orient5p, "upstream", "downstream"), ")"))
}
   

liberalImport = function(con, which, pad=1e5, genome, ...) {
 # intention is to import request plus padding so
 # that result can be trimmed to fill out the 
 # where region with any relevant metadata
 stopifnot(is(which, "GRanges"))
 tmp = import(con, which=which+pad, genome=genome, ...)
 tmp = restrict(tmp, start=start(which), end=end(which),
    keep.all.ranges=TRUE)
 tmp[ which(end(tmp) >= start(which) & start(tmp) <= end(which)) ]
# subsetByOverlaps(tmp, which)
}
