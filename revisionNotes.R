# in place of transcripList2GElist ...
	library(GenomicFeatures)
	library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	library(org.Mm.eg.db)
	txdb= TxDb.Mmusculus.UCSC.mm9.knownGene
	egdb= org.Mm.eg.db

	isActiveSeq(txdb)
	cols(txdb)
	keytypes(txdb)

	keys(txdb, keytype='TXNAME')[1:5]
	select(txdb, keys='uc007afg.1', keytype='TXNAME',cols=cols(txdb))
	select(txdb, keys='18777', cols=cols(txdb), keytype='GENEID')
	select(txdb, keys='18777', cols='TXCHROM', keytype='GENEID')
	select(egdb, keys='18777', cols='GENENAME', keytype='ENTREZID')

	# to select gene bodies
	transcripts(txdb, list(tx_name='uc007aew.1'))
	# to select promoters
	promoters(txdb, list(tx_name='uc007aew.1'), upstream=1000, downstream=500)
	# to select exons
	exons(txdb, list(tx_name='uc007aew.1'))
	# to select cds
	cds(txdb, list(tx_name='uc007aew.1'))
	# to select TSS positions
	tss = resize(transcripts(txdb), width=1)


# within a GEList object (to be associated with binscore, bin*C etc ..):
	gr= GRanges(seqnames= Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),		ranges= IRanges(1:10, end = 7:16, names = head(letters, 10)),		strand= Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)))
	# or GRangesList ???


# to determine uncovered regions based on a BAM files of alignments
readsS= sample(1:1e8, 1e7)
system.time(grbig<- GRanges(seqnames=Rle('chr1'), ranges=IRanges(readsS, readsS+100))) # this is a GRanges collecting 10 million randomly generated reads in a 1e8 bases chromosome
  user  system elapsed 
 0.436   0.277   0.713 
system.time(grcov<- coverage(grbig))
  user  system elapsed 
23.113   1.967  25.082 
grcov
SimpleRleList of length 1
$chr1
integer-Rle of length 100000098 with 18000018 runs
 Lengths:  1 41  3  6 22 28  1 14 13  4 10 ...  2 11  2 17  8  4 45  2  4  1
 Values :  0  1  2  3  4  5  6  5  6  7  8 ... 10  9  8  7  6  5  4  3  2  1

values is the coverage; this means that your coverage vector can be read as one 0, followed by forty-one 1, three 2 etc.. :
as.numeric(grcov$chr1)[1:100]
 [1] 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
[38] 1 1 1 1 1 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5
[75] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
system.time(inds0<- which(grcov$chr1==0)) # uncovered genomic pos
  user  system elapsed 
 4.785   0.991   5.778
system.time(gr0<- findRange(x=inds0, vec=grcov$chr1))
  user  system elapsed 
 0.192   0.054   0.246
gr0
IRanges of length 2631
         start      end width
[1]           1        1     1
[2]      288575   288575     1
[3]      399486   399494     9
[4]      399486   399494     9
[5]      399486   399494     9
[6]      399486   399494     9
[7]      399486   399494     9
[8]      399486   399494     9
[9]      399486   399494     9
...         ...      ...   ...
[2623] 99889592 99889607    16
[2624] 99889592 99889607    16
[2625] 99889592 99889607    16
[2626] 99889592 99889607    16
[2627] 99889592 99889607    16
[2628] 99889592 99889607    16
[2629] 99889592 99889607    16
[2630] 99889592 99889607    16
[2631] 99889592 99889607    16




