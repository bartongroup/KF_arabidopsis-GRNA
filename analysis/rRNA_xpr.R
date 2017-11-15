library(data.table)
library(ggplot2)

ribids <- c("ATMG01390", "AT2G01010", "AT2G01020", "AT3G41768", "AT3G41979", "ATCG00920",
            "ATCG00950", "ATCG00960", "ATCG00970", "ATCG01160", "ATCG01170", "ATCG01180",
            "ATCG01210", "ATMG00020", "ATMG01380")


#xpr <- as.data.table( read.csv("/Volumes/kfroussios/PROJECTS/AtWTdist/combined_counts/AtWT-lev.tsv", header = TRUE, sep = "\t") )
#xpru <- as.data.table( read.csv("/Volumes/kfroussios/PROJECTS/AtWTdist/combined_counts/AtWT-unlev.tsv", header = TRUE, sep = "\t") )
xpr <- as.data.table( read.csv("./combined_counts/AtWT-lev.tsv", header = TRUE, sep = "\t") )
xpru <- as.data.table( read.csv("./combined_counts/AtWT-unlev.tsv", header = TRUE, sep = "\t") )

xpr[, rRNA := Geneid %in% ribids]
xpr[, dataset := rep("re-sampled", length(rRNA))]
xpru[, rRNA := Geneid %in% ribids]
xpru[, dataset := rep("non-normalized", length(rRNA))]
setkey(xpr, Geneid)
setkey(xpru, Geneid)

mydata <- rbind(xpr[(rRNA), ], xpru[(rRNA), ])
setkey(mydata, Geneid)

df <- data.frame( "exp"= c( mydata[, AtWT_1], mydata[, AtWT_2], mydata[, AtWT_3], mydata[, AtWT_4], 
                            mydata[, AtWT_5], mydata[, AtWT_6], mydata[, AtWT_7], mydata[, AtWT_8], 
                            mydata[, AtWT_9], mydata[, AtWT_10], mydata[, AtWT_11], mydata[, AtWT_12],
                            mydata[, AtWT_13], mydata[, AtWT_14], mydata[, AtWT_15], mydata[, AtWT_16], 
                            mydata[, AtWT_17] ),
                  "rRNA"= rep(mydata[, Geneid], times= 17),
                  "dataset" = rep(mydata[, dataset], times= 17),
                  "sample"= rep(names(mydata)[2:18], each= 30) )

g <- ggplot(df, aes(x=sample, y=exp, fill=rRNA) ) +
  facet_grid(rRNA ~ dataset, scales="free") +
  geom_bar(stat="identity") +
  guides(fill="none") +
  labs(title= "rRNA abundance across replicates", y="Reads", x= "Replicate") +
  theme(title= element_text(size= rel(1.5)),
        axis.text.x= element_text(angle= 90, size= rel(1.5)),
        axis.text.y= element_text(size= rel(1.5)),
        strip.text.y= element_text(size= rel(1.5)),
        strip.text.x= element_text(size= rel(1.8)),
        strip.background= element_rect(fill= "grey85"),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey90"))

#pdf("~/rRNA_xpr.pdf", height=18, width=12, title="Supplementary Figure 1")
#pdf("/Volumes/kfroussios/PROJECTS/AtWTdist/results/rRNA_xpr.pdf", height=18, width=12, , title="rRNA abundance across replicates")
pdf("./results/rRNA_xpr.pdf")
print(g)
dev.off()

#mydata[, rRNA := NULL]
#setkey(mydata, dataset, Geneid)
#write.csv(mydata, "~/ribosomal_xpr.csv")
