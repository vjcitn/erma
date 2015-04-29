
genemodel = function (sym, genome = "hg19", annoResource = Homo.sapiens, 
    getter = exonsBy, byattr = "gene") 
{
    stopifnot(is.atomic(sym) && (length(sym) == 1))
    if (!exists(dsa <- deparse(substitute(annoResource)))) 
        require(dsa, character.only = TRUE)
    num = select(annoResource, keys = sym, keytype = "SYMBOL", 
        columns = c("ENTREZID", "SYMBOL"))$ENTREZID
    getter(annoResource, by = byattr)[[num]]
}

