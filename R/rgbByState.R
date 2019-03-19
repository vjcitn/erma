
rgbByState = function(sts) {
  rgbmap = structure(c("#FE0000", "#FE4500", "#FE4500", "#FE4500", "#008000", 
"#008000", "#008000", "#009500", "#C1E005", "#C1E005", "#C1E005", 
"#C1E005", "#FEC24D", "#FEC24D", "#FEC24D", "#FEFE00", "#FEFE00", 
"#FEFE00", "#FEFE66", "#66CCA9", "#8990CF", "#E5B7B6", "#70309F", 
"#808080", "#FEFEFE"), .Names = c("TssA", "PromU", "PromD1", 
"PromD2", "Tx5'", "Tx", "Tx3'", "TxWk", "TxReg", "TxEnh5'", "TxEnh3'", 
"TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", "EnhAc", 
"DNase", "ZNF/Rpts", "Het", "PromP", "PromBiv", "ReprPC", "Quies"))
 names(rgbmap) = paste0(1:25, "_", names(rgbmap))
 stopifnot(all(sts %in% names(rgbmap)))
 rgbmap[sts]
}

