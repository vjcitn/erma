
csProfile = function(ermaset, symbol, upstream=2000, downstream=200,
    useShiny=FALSE,
    ctsize=10, shortCellType=TRUE, tsswidth=3) {
#
# chromatin state profile: use IRanges 'promoters' concept to pick
#   states and colors from ChromImpute resources
#
     if (!useShiny) stateProfile( ermaset=ermaset, 
             symbol=symbol, ctsize=ctsize, upstream=2000, downstream=200,
             shortCellType=shortCellType)
     else stateProf(ermaset=ermaset, shortCellType=shortCellType, ctsize=ctsize)
}

stateProfile = function(ermaset, symbol="IL33", upstream=2000,
     downstream=200, ctsize=10, shortCellType=TRUE, tsswidth=3) {
   mod = try(genemodel(symbol))
   tss = start(resize(range(mod),1))
   if (inherits(mod, "try-error")) stop("can't resolve symbol")
   uil = promoters(range(mod), upstream=upstream, downstream=downstream) 
   ## ----bind----------------------------------------------------------------
   ermaset@rowRanges = uil  # binding this in allows reduceByFile to work
   ## ----getcss--------------------------------------------------------------
#   imps = reduceByFile(ermaset, MAP=function(range, file) {
#     imp = liberalImport(file, which=range, genome=genome(range))
#     seqlevels(imp) = seqlevels(range)
#     imp$rgb = rgbByState(imp$name)
#     imp
#   })
   range = rowRanges(ermaset)
   fe = files(ermaset)
   csstates = foreach(i = 1:length(fe)) %dopar% {
#     path = files(ermaset)[i]
#     con = file(path)
#     open(con, type="r")
#     on.exit(close(con))
     imp = liberalImport(fe[i], which=range, genome=genome(range))
     seqlevels(imp) = seqlevels(range)
     imp$rgb = rgbByState(imp$name)
     imp
   }
#   csstates = lapply(imps, "[[", 1) # only one range
   tys = cellTypes(ermaset)  # need to label with cell types
   ## ---- annotate
   csstates = lapply(1:length(csstates), function(x) {
      csstates[[x]]$celltype = tys[x]
      len = length(csstates[[x]])
      csstates[[x]] = c(csstates[[x]], csstates[[x]][len]) # dummy
      ranges(csstates[[x]][len+1]) = IRanges(tss,width=tsswidth)
      csstates[[x]]$name[len+1] = "TSS"
      csstates[[x]]$rgb[len+1] = "#000000"
      #print(range(mod))
      #print(csstates[[x]])
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
   mycol = c(mycol, "TSS"="#000000")
   cssd = as.data.frame(cssgr)
   chrn = as.character(seqnames(cssgr))[1]
   if (shortCellType) cssd$celltype = short_celltype[ cssd$celltype ]
   n2 = ggplot(cssd)
   n2 + geom_rect(aes(ymin=0, ymax=1, xmin=start, xmax=end, fill=name)) + 
       facet_grid(celltype~.) +
       theme(strip.text.y = element_text(size = ctsize, angle = 0)) +
        ylim(0,1) + scale_y_continuous(breaks=NULL, limits=c(0,1)) +
        scale_fill_manual( name="state", values=mycol ) + xlab(chrn) +
        ggtitle(paste0(symbol, " [", as.character(strand(mod)[1]), "]"))
}
   

liberalImport = function(con, which, pad=1e5, genome, ...) {
 # intention is to import request plus padding so
 # that result can be trimmed to fill out the 
 # 'which' region with any relevant metadata
 stopifnot(is(which, "GRanges"))
 tmp = import(con, which=which+pad, genome=genome, ...)
 tmp = restrict(tmp, start=start(which), end=end(which),
    keep.all.ranges=TRUE)
 tmp[ which(end(tmp) >= start(which) & start(tmp) <= end(which)) ]
# subsetByOverlaps(tmp, which)
}

subsetByRanges = function( ermaset, range ) {
   rowRanges(ermaset) = range
   reduceByFile(ermaset, MAP=function(range, file) {
   imp = import(file, which=range, genome=genome(range)[1])
   seqlevels(imp) = seqlevels(range)
   imp$rgb = rgbByState(imp$name)
   imp })
}

statesByRange = function(ermaset, rangeToUse,
     ctsize=10, shortCellType=TRUE, tsswidth=3) {
   if (inherits(mod, "try-error")) stop("can't resolve symbol")
   uil = promoters(range(mod), upstream=upstream, downstream=downstream) 
   ## ----bind----------------------------------------------------------------
   ermaset@rowRanges = uil  # binding this in allows reduceByFile to work
   ## ----getcss--------------------------------------------------------------
#   imps = reduceByFile(ermaset, MAP=function(range, file) {
#     imp = liberalImport(file, which=range, genome=genome(range))
#     seqlevels(imp) = seqlevels(range)
#     imp$rgb = rgbByState(imp$name)
#     imp
#   })
   range = rangeToUse
   fe = files(ermaset)
   csstates = foreach(i = 1:length(fe)) %dopar% {
#     path = files(ermaset)[i]
#     con = file(path)
#     open(con, type="r")
#     on.exit(close(con))
     imp = liberalImport(fe[i], which=range, genome=genome(range))
     seqlevels(imp) = seqlevels(range)
     imp$rgb = rgbByState(imp$name)
     imp
   }
   tys = cellTypes(ermaset)  # need to label with cell types
   ## ---- annotate
   csstates = lapply(1:length(csstates), function(x) {
      csstates[[x]]$celltype = tys[x]
      len = length(csstates[[x]])
      csstates[[x]] = c(csstates[[x]], csstates[[x]][len]) # dummy
      ranges(csstates[[x]][len+1]) = IRanges(tss,width=tsswidth)
      csstates[[x]]$name[len+1] = "TSS"
      csstates[[x]]$rgb[len+1] = "#000000"
      #print(range(mod))
      #print(csstates[[x]])
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
   mycol = c(mycol, "TSS"="#000000")
   cssgr
}
