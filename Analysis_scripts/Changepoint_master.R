



library(cn.mops)
library(changepoint)
library(Rsamtools)
library(parallel)
#ncores <- detectCores()
#print(ncores)
 

print(sessionInfo())

# source file of necessary functions
source(file="../Analysis_scripts/Changepoint_functions.R")


# Use defaults for required parameters if not specified in initiation file

if ( ! exists("minScaffReads") ) { 
	warning("Minimum scaffold reads per sample not specified.  Using default value: 100")
	minScaffReads  <- 100 
}
if ( ! exists("window.size") ) { 
	warning("Read count window size not specified.  Using default value: 500 bp")
	window.size <- 500 
}
if ( ! exists("minReads") ) {
	warning("Minimum window read count not speficied.  Using default value: 10 reads")
	minReads <- 10
}
if ( ! exists("base.name") ) {
	warning("No file basename specified. Using \"chgpt-tmp\" as default ")
	base.name <- "chgpt-tmp"
}
if (! exists("minSegment") ) {
	warning("No minimum segment size was specified. Using default of 2000 bp")
	minSegment <- 2000
}



############################
# Read data from bam files #
############################

# If reading from bam files on server, get read counts and scaffold details. Otherwise they are loaded from image
if (startWithBams) {
	# reading in data from bam files
	message("Reading data from bam files")
	genomeReadsTotal <- getTotalReads(bamFiles)
	
	scaffLengths <- getScaffoldLengths(bamFiles[1])
	
	# Create Granges object containing read counts from each sample for each window.
	bamDataRanges <- countReads(bamFiles, names(scaffLengths) )
	
	message("Completed reading data from bam files; saving workspace image")
	# Save environment for local analysis off of cluster
	save.image(file = paste0(outpath, ".Rdata"))
}


message("Generating total read counts per scaffold")
# Get total reads for each scaffold in each sample
scaffoldTotalReads <- lapply(bamDataRanges, function(x) { apply( X = mcols(x), MARGIN = 2, FUN = sum ) }  )   
# Partition scaffolds for sufficient reads in each sample
informativeScaffs <-  sapply(scaffoldTotalReads, function(x) {all(x >= minScaffReads) } )
# Get M:F ratios based on summed readcounts 
ratiosTotalReads <- sapply(scaffoldTotalReads, normRatios)


# Initiate dataframe to track various levels of filtering applied to scaffolds.
scaffoldFilters <- data.frame("scaffold" = names(scaffLengths), "lengths" = scaffLengths, "informative" = informativeScaffs, "ratios.total"= ratiosTotalReads, stringsAsFactors = F )



# Exclude scaffolds with insufficient read counts
# bamDataRanges <- bamDataRanges[informativeScaffs]
# scaffolds <- names(scaffLengths)[informativeScaffs]  # create simple vector of scaffold names


message("Filtering windows with minimum read count")
# for each scaffold, create indexing vector for windows with sufficient reads.
countsFilterList <- lapply(bamDataRanges, minReadsFilter, minimReads = minReads )  

# filter each GRanges object to remove low-read windows
bamDataFiltered <- lapply(1:length(bamDataRanges), function(i) {bamDataRanges[[i]][ countsFilterList[[i]] ] }) 
names(bamDataFiltered) <- scaffoldFilters$scaffold
	

# extract read count data for each scaffold, then calculate median readcount per window
counts <- lapply(bamDataFiltered, function(x){elementMetadata(x)})
countsMedian <- sapply(X = counts, FUN = function(x) { apply(x, 2, median)}) 

# Set minimum segment size for detection of changepoints
minSeg <- minSegment/window.size  # Min segment in basepairs divided by window size

# Determine scaffolds without sufficient windows for using changepoint analysis. 
# scaffolds must have twice as many windows as the minimum segment length, or cgpt.mean function chokes.
Nwindows <- sapply(counts, nrow)
scaffoldFilters$windows <-  Nwindows > minSeg * 2  # logical vector of scaffolds with sufficient windows for changepoint analysis

# These steps not needed if using scaffoldFilters dataframe.
#enoughWindows <- which(Nwindows > minSeg * 2)
#tooFewWindows <- sum(! Nwindows > minSeg * 2) 
#counts <- counts[enoughWindows]


########################
# changepoint analysis #
######################## 
message("Preparing for changepoint analysis")
# only analyze scaffolds with sufficient windows
ratiosList <- lapply(counts[scaffoldFilters$windows], normRatios)  # get normalized M:F ratio for each window in each scaffold
meanMFratio <- sapply(ratiosList, mean)  # get the mean M:F ratio across windows in each scaffold
ratiosFrame <- data.frame( "scaffold" = scaffoldFilters$scaffold[scaffoldFilters$windows], "MFratio" = meanMFratio, "scaffL" = scaffoldFilters$length[scaffoldFilters$windows], stringsAsFactors = F)

message("Running changepoint function on scaffolds")
# apply changepoint detection to all scaffolds
changeList <- lapply(ratiosList, function(x) {cpt.mean(x, method = "PELT", Q = 5, penalty= "SIC", minseglen = minSeg) })
#names(changeList) <- scaffolds


message("Parsing results of changepoint")
# some scaffolds don't have changepoints detected, so need to separate constant from chpt scaffolds
changeMeansList <- lapply(changeList, function(x) {return(x@param.est$mean)})  # extract changepoint mean vectors
changeCounts <- sapply(changeMeansList, length)  # count number of changes
scaffsWithChanges <- which(changeCounts > 1) # gives us indexing vector for which scaffolds have changes


# For scaffolds where there are detected changes, parse data in Granges objects, and then a dataframe.
ZAranges  <- lapply(names(scaffsWithChanges),  function(i) {parseChanges( changes = changeList[[i]], grange = bamDataFiltered[[i]])})
names(ZAranges) <- names(scaffsWithChanges)
ZAframe.list <- lapply(ZAranges, as.data.frame)

# check whether changes are within chromosome-type or are actually chimeric.
chimeraValidxn <- lapply(ZAranges, changeNotChimera)


#############################
# Classify scaffold linkage #
#############################

message("Assigning scaffolds to linkage")
# categorize scaffolds as Z, W, or chimeric
# apply classification to scaffolds analyzed via changepoint
ratiosFrame$linkage <- cut(x = ratiosFrame$MFratio, breaks = c( -Inf, cutoffW, cutoffZ, Inf), labels= c("W", "Autosome", "Z")) 
levels(ratiosFrame$linkage) <- c(levels(ratiosFrame$linkage), "Chimera") # adding "chimera" so that it is included as a factor
ratiosFrame[names(chimeraValidxn), "linkage"] <- unlist(chimeraValidxn) # update linkage assignment based on changepoint analysis.

# apply classifications to ALL scaffolds, based on Total read counts
scaffoldFilters$linkage <- cut(x = scaffoldFilters$ratios.total, breaks = c( -Inf, cutoffW, cutoffZ, Inf), labels= c("W", "Autosome", "Z")) 
levels(scaffoldFilters$linkage) <- c(levels(scaffoldFilters$linkage), "Chimera")
scaffoldFilters[names(chimeraValidxn), "linkage"] <- unlist(chimeraValidxn)


# for "true" chimeric scaffolds, get portions of A,Z,W
chimericSegments <- sapply(ZAranges[ rownames(ratiosFrame[ ratiosFrame$linkage == "Chimera", ]) ] , evalChimeras)
chimericSegments <- t(chimericSegments)
chimericProporxns <- chimericSegments/rowSums(chimericSegments)



#############################
#     Write out results     #
#############################
message("Writing out summary and tabular results")

cat("Assembly version: ", assembly, "\n\nRUN PARAMETERS\n" , 
	minScaffReads, "\tMinimum reads required across scaffold (in each sample)\n",
	window.size, "\tWindow size for counting reads\n",
	minReads, "\tMinimum reads required per window for inclusion in Changepoint analysis\n",
	minSegment, "\tMinimum segment size required for detecting changepoints\n\n\n",
	file = paste0(outpath, "_summary.txt"), sep = ""
	)

cat("RESULTS OVERVIEW\n",
	length(scaffLengths), "\tTotal scaffolds in assembly\n", 
	sum(scaffoldFilters$windows), "\tScaffolds analyzed via Changepoint\n", 
	sum(!scaffoldFilters$informative), "\tScaffolds with samples having below threshold read counts\n\n\n", 
	file = paste0(outpath, "_summary.txt"), sep = "", append = T)

scaffTable <- table(scaffoldFilters$linkage)

cat("CHROMOSOMAL ASSIGNMENT\n",
	scaffTable["Autosome"], "\tAutosomal scaffolds\n",
	scaffTable["Z"], "\tZ-linked scaffolds\n",
	scaffTable["W"], "\tW-linked scaffolds\n",
	scaffTable["Chimera"], "\tChimeric scaffolds\n",
 	append = T, file=paste0(outpath, "_summary.txt"), sep=""
 	)


# Write table of changepoint results
ZAframe <- do.call(what = rbind, args = ZAframe.list)  # only report change coords for scaffolds with detected changes.
write.csv(x=ZAframe, file = paste0(outpath, "_changepoint_segments.csv"), quote = F, row.names = T)

# Write table of analyzed scaffolds with linkage assignments.
# write.csv(x=scaffoldFilters, file = paste0(outpath, "_scaffold_assignments.csv"), quote = F, row.names = F)
# write.csv(x=ratiosFrame, file = paste0(outpath, "_scaffold_assignments.csv"), quote = F, row.names = F)
# Write table of chimeric scaffold regions
write.csv(x= chimericSegments, file = paste0(outpath, "_chimeric_portions-bp.csv"), quote = F, row.names = T)
write.csv(x= round(chimericProporxns, 8), file = paste0(outpath, "_chimeric_portions-ppxn.csv"), quote = F, row.names = T)


comboRatio <- merge(x = scaffoldFilters, y = ratiosFrame, by="scaffold", all.x=T)
names(comboRatio)[names(comboRatio) == "linkage.x"] <- "linkage.total"
names(comboRatio)[names(comboRatio) == "linkage.y"] <- "linkage.windows"
names(comboRatio)[names(comboRatio) == "MFratio"] <- "ratios.windows"

# Sort dataframe, and exclude redundant "scaffL" column
comboRatio <- comboRatio[rev(order(comboRatio$lengths)), grep(pattern = "scaffL", x = names(comboRatio), invert=T)]

write.csv(x=comboRatio, file = paste0(outpath, "_scaffold_assignments.csv"), quote = F, row.names = F)

comboRatio[which(comboRatio$linkage.total != comboRatio$linkage.windows), ]

#######################
#    PLOT DATA        #
#######################
message("Starting data plotting")
# Give distributions of median read counts values.
pdf(file = paste0(outpath, "_MedianReadCounts.pdf"), width = 8, height = 6 )
boxplot(data.frame(t(countsMedian)), col = "sky blue", xlab = "Sample", ylab = "Scaffold median reads per window", sub = paste("Window size:", window.size, "bp"), outline = F , main = base.name )
dev.off()


message("Plotting Lengh-Ratio")
pdf(file = paste0(outpath, "_ScaffoldMeanRatios-Windows.pdf"), width = 8, height = 6 )
plotLengthRatio(lengths = ratiosFrame$scaffL, ratios = ratiosFrame$MFratio, linkage = ratiosFrame$linkage, main = assembly, sub="Ratio Method: Mean Windows", ylim = c(-5,2))
dev.off()

pdf(file = paste0(outpath, "_ScaffoldMeanRatios-Total.pdf"), width = 8, height = 6 )
plotLengthRatio(lengths = scaffoldFilters$lengths, ratios = scaffoldFilters$ratios.total, linkage = scaffoldFilters$linkage, main = assembly, sub="Ratio Method: summed Reads", ylim = c(-5,2))
dev.off()


message("Plotting correlation of M:F ratios between total and window read counts")
#
pdf(file = paste0(outpath, "_ratios-correlation.pdf"), width = 6, height = 6 )
plot(comboRatio$ratios.total ~ comboRatio$ratios.window, xlim = c(-5, 2), ylim = c(-5,2), pch=19, cex = .3,
xlab = "M:F ratios (Windows Based)", ylab="M:F ratios (Total counts)", main = assembly)
abline(a = 0, b = 1, col = "red")
dev.off()


### get unfiltered data to support plotting
counts.unfilt <- lapply(bamDataRanges, function(x){elementMetadata(x)})
ratiosList.unfilt <- lapply(counts.unfilt, normRatios)



message("Plotting individual scaffolds, by groupings")
# Plot scaffolds with reported changepoints
pdf(file = paste0(outpath, "_ChangepointScaffs.pdf"), width = 8, height = 6 )
for (i in names(scaffsWithChanges) ) {  # only plot scaffolds with detected changes
	plotBoth(changes = ratiosList.unfilt[[i]], unfilt = ratiosList.unfilt[[i]], zaframe = ZAframe.list[[i]], ratiosFrame = ratiosFrame[i,], main = i, ylim = c(-6.5,6), ySeg = -5.8 , sub=assembly )
}
dev.off()

# Plot scaffolds that are chimeric (will be a subset of "changepoint" scaffolds)
pdf(file = paste0(outpath, "_ChimericScaffs.pdf"), width = 8, height = 6 )
for (i in rownames(ratiosFrame[ratiosFrame$linkage == "Chimera",]) ) {  # only plot scaffolds with detected changes
	plotBoth(changes = ratiosList.unfilt[[i]], unfilt = ratiosList.unfilt[[i]], zaframe = ZAframe.list[[i]], ratiosFrame = ratiosFrame[i,], main = i, ylim = c(-6.5,6), ySeg = -5.8 , sub=assembly )
}
dev.off()

# Plot Z scaffolds
pdf(file = paste0(outpath, "_Z-Scaffs.pdf"), width = 8, height = 6 )
for (i in scaffoldFilters$scaffold[scaffoldFilters $linkage == "Z"]) {
		plotBoth(changes = ratiosList.unfilt[[i]], unfilt = ratiosList.unfilt[[i]], zaframe = ZAframe.list[[i]], ratiosFrame = ratiosFrame[i,], main = i, ylim = c(-6.5,6), ySeg = -5.8 , sub=assembly   )
}
dev.off()

# Plot W scaffolds
pdf(file = paste0(outpath, "_W-Scaffs.pdf"), width = 8, height = 6 )
for (i in scaffoldFilters $scaffold[scaffoldFilters $linkage == "W"]) {
		plotBoth(changes = ratiosList.unfilt[[i]], unfilt = ratiosList.unfilt[[i]], zaframe = ZAframe.list[[i]], ratiosFrame = ratiosFrame[i,], main = i , ylim = c(-6.5,6), ySeg = -5.8 , sub=assembly  )
}
dev.off()

save.image(file = paste0(outpath, ".Rdata"))

