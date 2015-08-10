
geneTxRange = function (sym, annoResource = Homo.sapiens)
{
    stopifnot(is.atomic(sym) && (length(sym) == 1))
    if (!exists(dsa <- deparse(substitute(annoResource)))) 
        require(dsa, character.only = TRUE)
    num = select(annoResource, keys = sym, keytype = "SYMBOL", 
        columns = c("TXCHROM", "TXSTART", "TXEND"))
    si = seqinfo(annoResource)
    ans = GRanges(as.character(num$TXCHROM[1]), IRanges(min(num$TXSTART), max(num$TXEND)))
    mcols(ans)[["SYMBOL"]] = sym
    seqinfo(ans) = si[as.character(num$TXCHROM[1]),]
    ans
}

