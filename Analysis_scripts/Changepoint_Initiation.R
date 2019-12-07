 
# The starting point for the changepoint analysis will be 
# 1.	List of bam files
# 2.	Path to list of bamfiles
# 3.	Sample identifiers -- sampleIDs
# 4.	Vectors identifying sex of samples. Must be indexing vectors --> males, females
# 5.  	Minimum numer of reads required over entire scaffold, for each sample --> minScaffReads
# 6.  	Window size (in bp) for changepoint analysis --> window.size
# 7.  	Minimum reads required in window for inclusion in CP --> minReads
# 8.  	File Versions (shared portion of bam file name) --> file_version
# 9.  	base name for analysis --> base.name
# 10. 	Minimum segment size for detecting changepoints --> minSegment
# 11. 	log2(M:F) ratio cutoffs for Z and W detection: cutoffZ   cutoffW
# 12. 	Assembly version being analyzed --> assembly


ncores <- 6 # how many cores to run in parallel when reading bam file
base.name <- "FALCON_HiFi_Q30"
outpath <- paste0("./", base.name)

startWithBams <- FALSE  # Set to TRUE for analyzing from scratch, starting with reading bams.

# If running on server and starting from counting reads in bam, then run this segment. Otherwise, load saved workspace.
if (startWithBams) {
	
	datapath <- "/scratch/l338g110/Danaus_Trio/Bowtie_FALCON_HiFi_Q30" # data directory where relevant bams are temporarily stored.
	file_version <- "1.bam"
	bams <- list.files(pattern = paste0(file_version, "$"), path=datapath )  # bam files to read 
	bamFiles <- paste0(datapath,bams)
	
	# parse file names to get sample ID and sex.
	sampleIDs <- sub(pattern = paste0("(\\w{4})", file_version), replacement ="\\1", perl = T, x = bams)
	males <- grep(pattern = c("M"),  x = sampleIDs)
	females <-  grep(pattern = c("F"),  x = sampleIDs)
} else {
	load(paste0(outpath, ".Rdata"))
    startWithBams <- FALSE   # Do not change this value
}

# Specify parameters for analysis
assembly <- "FALCON_HiFi_Q30"
window.size <- 1000		# How many basepairs in each window for counting reads
minScaffReads <- 50    # Each sample must have at least this many reads for each sample to be analyzed.
minReads <- 6            # Minimum reads required in window for inclusion in CP
minSegment <- 5000        # Minimum segment size for detecting changepoints; valid windows must provide twice this for changepoint analysis
cutoffZ <- 0.585            # Threshold of Log2(M:F) above which scaffolds/regions are considered Z-linked; default 0.585 correspondong to M:F ratio of 1.5
cutoffW <- -0.585        # Threshold of Log2(M:F) below which scaffolds/regions are considered W-linked

# Calling master script also calls all functions and runs entire analysis.
source(file="../Analysis_scripts/Changepoint_master.R")
