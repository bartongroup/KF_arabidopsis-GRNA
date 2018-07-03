args <- commandArgs(trailingOnly=TRUE)

# Setup variable parameters.
WORKPATH <- getwd()
VALFILE <- args[[1]]
Headers <- TRUE
RowLabels <- TRUE
Title <- "\nInter-Replicate Correlation"
XLabel <- YLabel <- "Replicate"
outliers = c(11)
outliernames = c("AtWT_11")

# Import data.
VALS <- read.table(file.path(WORKPATH, VALFILE), header=Headers, sep="\t")

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
pdf(file=paste0(file.path(WORKPATH, VALFILE), '_heat.pdf'), onefile=TRUE)
levelplot(as.matrix(VALS), col.regions=myPalette, scales=list(x=list(rot=90)),
          xlab=XLabel, ylab=YLabel, main=Title, cuts= length(myPalette)-1)
dev.off()

# Draw without outliers.
pdf(file=paste0(file.path(WORKPATH, VALFILE), "_heat_fltr.pdf", sep=""))
levelplot(as.matrix(VALS[-outliers, ! names(VALS) %in% outliernames]), col.regions=myPalette, scales=list(x=list(rot=90)),
          xlab=XLabel, ylab=YLabel, main=paste(Title, "- filtered", seq=""), cuts= length(myPalette)-1)
dev.off()


#EOF
