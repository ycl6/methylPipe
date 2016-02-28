# methylPipe for custom genomes
This is a fork of the methylPipe Bioconductor package (https://github.com/Bioconductor-mirror/methylPipe)
Some of the functions in the original package is not compatible with non-UCSC genomes, as a result it is not most suitable for custome genomes. I have make some changes to fix this.

### Changes
##### R/BSdataSet-methods.R
1. methstats()
* line 548 - disable dev.new() when using pdf()

##### R/Allfunctions.R
1. BSprepare()
* line 187 - add "addchr" switch
* line 226 - check "addchr"

2. plotMeth()
* line 825 - add "ucscOrg" switch and band dataframe required if ucscOrg=FALSE
* line 886-894 - use custom band dataframe to create Ideogram track with custom genome

