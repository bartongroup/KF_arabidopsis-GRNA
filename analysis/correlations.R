# Setup variable parameters.
WORKPATH <- "~/Projects/AtWTdist/"
DATAPATH <- "combined_counts/"
OUTPATH <- "results/"
DATAFILE <- "AtWT-unlev.tsv"
Headers <- TRUE
index <- "Geneid"

# Import data.
DATA <- read.table( paste(WORKPATH, DATAPATH, DATAFILE, sep=""), header=Headers)
#View(DATA)

# Do all the pairwise correlations between the columns, ignoring a non-numerical column by name.
CORS <- cor(DATA[names(DATA) != index])
#View(CORS)

# Remove the symmetric bottom-left half of the matrix
# Start with an empty matrix of same size as CORS
nc = length(CORS[1,])
nr = length(CORS[,1])
DIAGON <- matrix(ncol=nc, nrow=nr)
rownames(DIAGON) <- rownames(CORS)
colnames(DIAGON) <- colnames(CORS)
# Loop through rows (of increasing column offset)
for(r in seq(1,nr)){
    DIAGON[r,r:nc] <- CORS[r,r:nc] # include the diagonal of 1s
}
# OR
#for(r in seq(1,nr-1)){
#}
#  DIAGON[r,(r+1):nc] <- CORS[r,(r+1):nc] # leave out the diagonal of 1s
#View(DIAGON)

# Save the correlation matrix and non-redundant matrix to files.
write.csv(CORS,
file=paste(WORKPATH, OUTPATH, DATAFILE, "_cors.csv", sep=""),
dim(CORS)[1])
write.csv(DIAGON,
file=paste(WORKPATH, OUTPATH, DATAFILE, "_corsnr.csv", sep=""),
dim(DIAGON)[1])


# Descriptive stats
library(pastecs)
metrics <- stat.desc(as.vector(DIAGON))
#metrics
write.table(metrics, file=paste(WORKPATH, OUTPATH, DATAFILE, "_cors.csv_desc.txt", sep=""))


#EOF
