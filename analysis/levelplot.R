# Setup variable parameters.
WORKPATH <- "~/Projects/AtWTdist/results/"
VALFILE <- "AtWT-unlev.tsv_cors.csv"
Headers <- TRUE
RowLabels <- TRUE
Title <- "\nInter-Replicate Correlation"
XLabel <- YLabel <- "Replicate"

# Import data.
VALS <- read.csv( paste(WORKPATH, VALFILE, sep=""), header=Headers)
# Set up row labels
if (RowLabels == TRUE) {
    rn <- VALS[,1]
    # Get rid of the useless first column.
    VALS <- VALS[,2:length(VALS[1,])]
    rownames(VALS) <- rn
}
#View(VALS)

# OPTIONAL log tranform of values
#VALS <- log(VALS)
#View(VALS)

library(RColorBrewer)
library(lattice)

# Set up colours and range.
# EITHER monotone
myPalette <- colorRampPalette(c("darkred", "white"))(40)
# OR two-tone
#myPalette <- c(colorRampPalette(c("darkblue","white"))(20),
#               colorRampPalette(c("white", "darkred"))(20))


# Draw the image.
pdf(paste(WORKPATH, VALFILE, "_heat.pdf", sep=""))
levelplot(as.matrix(VALS), col.regions=myPalette, scales=list(x=list(rot=90)),
          xlab=XLabel, ylab=YLabel, main=Title, cuts= length(myPalette)-1)
dev.off()

# Draw without outliers.
outliers = c(11)
outliernames = c("AtWT_11")
pdf(paste(WORKPATH, VALFILE, "_heat_fltr.pdf", sep=""))
levelplot(as.matrix(VALS[-outliers, ! names(VALS) %in% outliernames]), col.regions=myPalette, scales=list(x=list(rot=90)),
          xlab=XLabel, ylab=YLabel, main=paste(Title, "- filtered", seq=""), cuts= length(myPalette)-1)
dev.off()


#EOF
