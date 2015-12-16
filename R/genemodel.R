
genemodelOLD = function (sym, genome = "hg19", annoResource = Homo.sapiens, 
    getter = exonsBy, byattr = "gene") 
{
    stopifnot(is.atomic(sym) && (length(sym) == 1))
    if (!exists(dsa <- deparse(substitute(annoResource)))) 
        require(dsa, character.only = TRUE)
    num = AnnotationDbi::select(annoResource, keys = sym, keytype = "SYMBOL", 
        columns = c("ENTREZID", "SYMBOL"))$ENTREZID
    getter(annoResource, by = byattr)[[num]]
}

genemodel = function(key, keytype="SYMBOL", annoResource=Homo.sapiens) {
#
# purpose is to get exon addresses given a symbol
# propagate seqinfo from annotation resource
#
  oblig=c("EXONCHROM", "EXONSTART", "EXONEND", "EXONSTRAND", "EXONID")
  addrs = AnnotationDbi::select(annoResource, keys=key, keytype=keytype, columns=oblig)
  ans = GRanges(addrs$EXONCHROM, IRanges(addrs$EXONSTART, addrs$EXONEND),
      strand=addrs$EXONSTRAND, exon_id=addrs$EXONID)
  mcols(ans)[[keytype]] = key
  useq = unique(as.character(seqnames(ans)))
  si = seqinfo(annoResource)
  seqinfo(ans) = si[useq,]
  ans
}

map2range = function(maptag="17q12", annoResource=Homo.sapiens) {
  t1 = AnnotationDbi::select(annoResource, keys=maptag, keytype="MAP", columns=
      c("TXSTART", "TXEND", "TXCHROM"))
  ans = GRanges(t1$TXCHROM[1], IRanges(min(t1$TXSTART), max(t1$TXEND)))
  si = seqinfo(annoResource)
  seqinfo(ans) = si[t1$TXCHROM[1],]
  ans
}
