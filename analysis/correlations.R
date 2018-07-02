## Rscript basedir insubdir infile outsubdir

args <- commandArgs(trailingOnly=TRUE)

# Setup variable parameters.
WORKPATH <- getwd()
DATAPATH <- 'combined_counts'
DATAFILE <- args[[1]]
OUTPATH <- "results"
Headers <- TRUE
index <- "Geneid"

# Import data.
DATA <- read.table(file.path(WORKPATH, DATAPATH, DATAFILE), header=Headers)
#View(DATA)

# Shorten the column labels (file designations)
names(DATA) <- gsub('__Aligned.out.bam.exon_..1', '', names(DATA), fixed = TRUE)

# Do all the pairwise correlations between the columns, ignoring a non-numerical column by name.
CORS <- cor(DATA[names(DATA) != index])

# Save the correlation matrix.
write.csv(CORS, file=paste0(file.path(WORKPATH, OUTPATH, DATAFILE), "_cors.csv"), dim(CORS)[1])

#EOF
