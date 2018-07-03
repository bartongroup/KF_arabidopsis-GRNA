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

# Do all the pairwise correlations between the columns, ignoring a non-numerical column by name.
CORS <- cor(DATA[names(DATA) != index])

# Save the correlation matrix.
write.table(CORS, paste0(file.path(WORKPATH, OUTPATH, DATAFILE), '_cors.tsv'), col.names=Headers, append=FALSE, sep="\t")

#EOF
