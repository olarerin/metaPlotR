### Read in file
dist <- read.delim ("dist.measures.txt", header = T)

## Select the longest isoforms
trx_len <- dist$utr5_size + dist$cds_size + dist$utr3_size
dist <- dist[order(dist$gene_name, trx_len),] # sort by gene name, then transcript length
dist <- dist[duplicated(dist$gene_name),] # select the longest isoform

## A simple histogram
hist(dist$rel_location, breaks = 200, col = "black")
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## Rescale regions
## Determine scale factor
utr5.SF <- median(dist$utr5_size, na.rm = T)/median(dist$cds_size, na.rm = T)
utr3.SF <- median(dist$utr3_size, na.rm = T)/median(dist$cds_size, na.rm = T)

# assign the regions to new dataframes
utr5.dist <- dist[dist$rel_location < 1, ]
cds.dist <- dist [dist$rel_location < 2 & dist$rel_location >= 1, ]
utr3.dist <- dist[dist$rel_location >= 2, ]

# rescale 5'UTR and 3'UTR
library("scales")
utr5.dist$rel_location <- rescale(utr5.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.dist$rel_location <- rescale(utr3.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))

# Combine the regions in a new dataframe and plot
all.regions <- c(utr5.dist$rel_location, cds.dist$rel_location, utr3.dist$rel_location)
hist.data <- hist(all.regions, breaks = 200, col = "black", xlim = c(0,3)) # plot and save to variable
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## Alternate representations of the metagene
## A line plot
plot(hist.data$breaks[1:length(hist.data$breaks)-1], hist.data$counts, type = 'l')
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## A smooth density plot
plot(density(all.regions))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")

## An absolute distance plot
hist(dist$utr3_st, xlim = c(-500,500), breaks = 2000, col = "black")
abline(v=0, col = "red")