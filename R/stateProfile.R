stateProfile = function(ermaset, symbol="IL33", width=50000,
    ctsize=7, shortCellType=TRUE) {
   mod = try(genemodel(symbol))
   if (inherits(mod, "try-error")) stop("can't resolve symbol")
   uil = flank(resize(range(mod), 1), width=width)
   
   ## ----bind----------------------------------------------------------------
   rowRanges(ermaset) = uil
   ## ----getcss--------------------------------------------------------------
   csstates = lapply(reduceByFile(ermaset, MAP=function(range, file) {
     imp = import(file, which=range)
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
   data(short_celltype)
   data(states_25)
   mycol = states_25$rgb
   names(mycol) = paste0(1:25, "_", states_25$MNEMONIC)
   cssd = as.data.frame(cssgr)
   if (shortCellType) cssd$celltype = short_celltype[ cssd$celltype ]
   n2 = ggplot(cssd)
   n2 + geom_rect(aes(ymin=0, ymax=1, xmin=start, xmax=end, fill=name)) + 
       facet_grid(celltype~.) +
       theme(strip.text.y = element_text(size = ctsize, angle = 0)) +
        ylim(0,1) + scale_y_continuous(breaks=NULL, limits=c(0,1)) +
        scale_fill_manual( name="state", values=mycol )
}
   
