


# scan bam files to get total counts of aligned reads for normalization.
getTotalReads <- function (bamFiles) {  # bamFiles is vector of bam files, with full filepath
	samidxList <- lapply(bamFiles, idxstatsBam)  # read in alignment stats from bam index to get all reads
	names(samidxList) <- sampleIDs
	genomeReads <- sapply(samidxList, function(x) (return(x$mapped)))  # extract read counts into matrix
	rownames(genomeReads) <- samidxList[[1]]$seqname  
	genomeReadsTotal <- colSums(genomeReads)
	return(genomeReadsTotal)
}


# Extract scaffold lengths from Bam, and return length-sorted vector
getScaffoldLengths <- function(bam.1 = bamFiles[1]) {
	header <- scanBamHeader(files=bam.1)
	scaffLengths <- header[[1]][[1]]
	return (sort(scaffLengths, decreasing = T) )
	
	# BELOW: legacy code that filtered and returned ONLY selected scaffold. This will be done in main code now, so accounting of length-filtered scaffolds can be made. 
	# scaffolds <- names(scaffLengths[scaffLengths > minScaffL]) # Use this to line to analyze all scaffolds greather than L bases...
	# target.scaffLengths <- scaffLengths[scaffolds] 
	# scaffolds <- names(sort(target.scaffLengths, decreasing = T))  # order selected scaffold names by scaffold length
	# return(scaffolds)
}


countReads <- function(bamFiles, scaffolds) {
	bamDataRanges <- lapply(scaffolds, function(scaff) { suppressMessages(getReadCountsFromBAM(bamFiles, parallel = ncores, sampleNames = sampleIDs, refSeqName = scaff, WL= window.size) ) } ) #
	names(bamDataRanges) <- scaffolds # assign names so they trickle down through analysis
	return(bamDataRanges)
}

# create and apply a filter function to remove windows with too few reads.
minReadsFilter <- function(x, minimReads ) { 
	countsFrame <- elementMetadata(x)
	selectWindows <- which(apply(countsFrame, 1, function(x) {all(x >  minimReads)} ) )
	return(selectWindows)
}

# function to normalize read counts per library and generate log2 transformed M:F ratios
normRatios <- function(cD) {
	libSizes <- genomeReadsTotal/1e8

	if ( is.null(nrow(cD)) ) {   # operating on a vector
		counts.norm <- cD/libSizes
		dat.m <- mean(counts.norm[males])
		dat.f <- mean(counts.norm[females])
	} else { # operating on matrix or dataframe
		counts.norm <- t(apply(cD, 1, function(x) {x/libSizes} ))
		if (nrow(counts.norm) == 1 ) { # for some bewildering reason, subsetting a single-row matrix returns a vector, where "apply()" won't work
			dat.m <- mean(counts.norm[, males])
			dat.f <- mean(counts.norm[, females])
		} else {
		#	dat.m <- apply(X = counts.norm[, males], MARGIN = 1, FUN = mean)  
		#	dat.f <- apply(X = counts.norm[, females], MARGIN = 1, FUN = mean)
			dat.m <- counts.norm[, males]  
			dat.f <- counts.norm[, females]
		}
	}
	ratios <- log2(dat.m/dat.f)
	return(ratios)
}




# process results from changepoint analysis into Granges object
parseChanges <- function(changes, grange, za.cutoff = cutoffZ, wa.cutoff = cutoffW) {
	myScaff <- seqnames(grange)[1]  # get current scaffold being analyzed
	chPoints <- grange[cpts(changes)] # get GRanges corresponding to changepoints
	starts <- c(1, start(ranges(chPoints)))  # start points for any Z|A segments
	lastBase <- end(ranges(tail(grange,1)))  # get last base of scaffold
	ends <- c(starts[-1] - 1 ,  lastBase)  # end points for any segments Z|A segments
	linkage <- cut(x = changes@param.est$mean, breaks = c( -1000, wa.cutoff, za.cutoff, 1000), labels= c("W", "Autosome", "Z"))  # assign linkage for region based on coverage. 

	# build new Granges object
	zaRanges <- GRanges(seqnames = myScaff, 
						ranges = IRanges(start = starts, end = ends),
						ZAlinkage = linkage,
						chgpt.mean = changes@param.est$mean
						)
	return(zaRanges)
}



# Evaluate changepoint-positive scaffolds for consistency of classification, e.g. are all segments of the same type...
changeNotChimera <- function( grNow ) {
	segmentLinkage <- elementMetadata(grNow)$ZAlinkage  # extract vector of chromosomal assignment for each changepoint segment detected.
	segmentUniq <- unique(as.character(segmentLinkage)) # linkage initially a factor, so transform to character
	ifelse(length(segmentUniq) == 1, segmentUniq, "Chimera")
}


evalChimeras <- function( grNow ) {
	df <- as.data.frame(grNow)
	chr <- c("Autosome", "Z", "W")
	segments <- c(0,0,0)
	names(segments) <- chr
	for (i in chr) {
		segments[i] <- sum(df$width[df$ZAlinkage == i])
	}
	return(segments)
}


# Set default plotting colors
zpoint.col <- "dodger blue"  # Z linked scaffolds
zline.col <- "dodger blue" # Theoretical line for Z linkage (Log2 = 1)
apoint.col <- "black"  # color for plotting autosomal scaffolds
aline.col <- "grey"  # Theoretical line for Autosomes (log2 = 0)
wpoint.col <- "green" # point color for W-linked scaffolds
cpoint.col <- "deep pink" # point color for Chimeric scaffolds



# Plot M:F ratio as function of Log10(scaffold length), color coded as linkage assignment.
plotLengthRatio <- function(lengths, ratios, linkage, ... ) 	
	{
		plot(x = log10(lengths), y = ratios, pch = 19, cex = 0.5, 
			ylab = "Log2(Mean M:F read counts)", xlab = "Log10(Scaffold Length)",
			col = c(wpoint.col, apoint.col, zpoint.col, cpoint.col)[linkage], ...
		)
		abline(h = 0, col = aline.col)
		abline(h = 1, col = zline.col, lty = 2)
		legend("bottomright", pch = 19, legend = c("W", "Autosomal", "Z", "Chimera"), 
			col = c(wpoint.col, apoint.col, zpoint.col, cpoint.col) )
	}





# Function to plot scaffolds with changepoint results, with segments differentiated.
# Inputs:   "changes" : unfiltered M:Fratios for the target scaffold
#			"zaframe" : dataframe version of Granges object defining changepoint results
#			"Yseg"	  : y-value for plotting segment labels

# Test code for plotCPza function
# myScaff <- "chr1"
# plotCPza(changes = ratiosList.unfilt[[myScaff]], zaframe = ZAframe.list[[myScaff]], main = myScaff, ylim = c(-7,7), ySeg = -6)

plotCPza <- function(changes , zaframe, ySeg = -2, ...) {
	plot(x = ( 1:length(changes) * window.size )/1000 , y = changes, 
		pch = 20 , col = "grey", type = "p",  
		xlab = "Scaffold position (Kbp)", ylab = "Log2(Male:Female Read Counts)", ...
		)
	yval <- ySeg
	ySegLab <- ySeg + 0.1 * ySeg
	offset <- tail(zaframe$end/1000, 1) * 0.01  # create slight space between adjacent lines.
	edgeL <- 0.05
	
	# delineate Autosomal sections 
	if (any(zaframe$ZAlinkage == "Autosome")) {
		segments(x0 = zaframe$start[zaframe$ZAlinkage == "Autosome"]/1000, 
			 x1 = zaframe$end[zaframe$ZAlinkage == "Autosome"]/1000, 
			 y0 = zaframe$chgpt.mean[zaframe$ZAlinkage == "Autosome"], 
			 lend = "round", lwd = 1
		)
		arrows(	x0 = zaframe$start[zaframe$ZAlinkage == "Autosome"]/1000 + offset, 
				x1 = zaframe$end[zaframe$ZAlinkage == "Autosome"]/1000, 
				y0= yval, lwd = 3, angle = 90, code = 3, length = edgeL
		)
		Alab <- rowMeans(zaframe[zaframe$ZAlinkage == "Autosome", c("start", "end")])
		text(x = Alab/1000, y = ySegLab, labels= "A")
	}
	
	# delineate Z sections
	if (any(zaframe$ZAlinkage == "Z")) {
		segments(x0 = zaframe$start[zaframe$ZAlinkage == "Z"]/1000 , 
				 x1 = zaframe$end[zaframe$ZAlinkage == "Z"]/1000, 
				 y0 = zaframe$chgpt.mean[zaframe$ZAlinkage == "Z"], col = zline.col 
			)
		arrows( x0 = zaframe$start[zaframe$ZAlinkage == "Z"]/1000 + offset, 
				x1 = zaframe$end[zaframe$ZAlinkage == "Z"]/1000, 
				y0= yval, lwd = 3, angle = 90, code = 3, col = "dodger blue", lty = 1, length = edgeL 
			)
		Zlab <- rowMeans(zaframe[zaframe$ZAlinkage == "Z", c("start", "end")])
		text(x = Zlab/1000, y = ySegLab, labels= "Z")
	}	
	
	# delineate W sections
	if (any(zaframe$ZAlinkage == "W")) {
		segments(x0 = zaframe$start[zaframe$ZAlinkage == "W"]/1000 , 
				 x1 = zaframe$end[zaframe$ZAlinkage == "W"]/1000, 
				 y0= zaframe$chgpt.mean[zaframe$ZAlinkage == "W"], col = wpoint.col, lty = 1 
			)
		arrows( x0 = zaframe$start[zaframe$ZAlinkage == "W"]/1000 + offset, 
				x1 = zaframe$end[zaframe$ZAlinkage == "W"]/1000, 
				y0= yval, lwd = 3, angle = 90, code = 3, col = wpoint.col, lty = 1, length = edgeL 
			)
	Zlab <- rowMeans(zaframe[zaframe$ZAlinkage == "W", c("start", "end")])
	text(x = Zlab/1000, y = ySegLab, labels= "W")
	}	
}



# Plot windows and average M:F ratio for chromosomes without changepoints.

# test code for plotScaff
# plotScaff(unfilt = ratiosList.unfilt[[30]], ratiosFrame[30,])
plotScaff <- function(unfilt, ratiosFrame, ...) {
	plot(	x = ( 1:length(unfilt) * window.size )/1000, 
			y = unfilt, 
			pch = 20 , col = "grey", type = "p",  
			xlab = "Scaffold position (Kbp)", ylab = "log2(Male:Female) Read Counts", ... )
		abline(h = ratiosFrame[,"MFratio"], col =  c(wpoint.col, apoint.col, zpoint.col, cpoint.col)[ratiosFrame[,"linkage"]])
}


# Function to plot either contstant scaffolds, or those with changepoints detected.
plotBoth <- function (zaframe, changes, unfilt, ratiosFrame, ySeg, ... ) {
	if ( is.null(zaframe) ) {
		plotScaff(unfilt = unfilt, ratiosFrame = ratiosFrame, ...)
	} else {
		plotCPza(changes = changes, zaframe = zaframe, ySeg, ...)
	}
}
