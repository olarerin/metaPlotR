# metaPlotR
A Perl/R pipeline for plotting metagenes

---
title: "Visualizing metagenes with MetaPlotR"
output: html_document
---

```test this block of code
  and see what happens```
  
### Understanding the data
This is the final script in the MetaPlotR pipeline.  The input for this file is 
the metagene coordinates file outputted from rel_and_abs_dist_calc.pl

Read in file
```{r}
setwd (dir = '~/Projects/metagene_Bioinformatics/')
dist <- read.delim ("dist.measures.txt", header = T)
```

View the number of rows and columns in the dataset
```{r}
dim(dist)
```
View the first few lines
```{r}
head(dist)
```
The input file contains `r dim(dist)[1]` rows and `r dim(dist)[2]` columns.  Each row represents a single site (in this example an m6A site).  The column headers for the first four columns are self explanatory.  The fifth column "rel_location" (for relative location) contains the calculated metagene coordinates.  In its simplest form (i.e. non-normalized), the metagene coordinates from 0 to 1 represent the 5'UTR with 0 being closer to the beginnig of the 5'UTR and 1 closer to the end.  Similarly, 1 to 2 represents the CDS and 2 to 3 the 3'UTR.  A histogram/density plot of the "rel_location" value gives the standard metagene.

In addition to the standard metagene which is based on the relative location of sites in transcripts, this next six columns (utr5_st, utr5_end, cds_st, cds_end, utr3_st, utr3_end) contain information for plotting the absolute distance of sites from several points of interest. For example, in this dataset the third row has a value of +63 under column header "utr3_st".  That means the site is 63 nucleotides upstream of the 3'UTR start site.

The last three columns contain the lengths of the 5'UTRs, coding sequences and 3'UTRs.

### Selecting gene isoforms for metagene analysis
The dataset is redundant -- a given site is represented by multiple transcript isoforms.  The choice of which isoforms to choose should be informed by the underlying biology.  For example, if a gene expression dataset is available, one option may be to pick the highest expressed isoform.  Another option is to pick the longest isoform, which is likely to capture more sites.  Below is sample code for picking the largest isoforms
```{r}

trx_len <- dist$utr5_size + dist$cds_size + dist$utr3_size
dist <- dist[order(dist$gene_name, trx_len),] # sort by gene name, then transcript length
dist <- dist[duplicated(dist$gene_name),] # select the longest transcript
dim(dist)
```

### Visualizing the metagene
#### A simple histogram
```{r}
hist(dist$rel_location, breaks = 200, col = "black")
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
```

In this plot, the range 0 to 1 represents the 5'UTR, 1 to 2 the CDS, and 2 to 3 the 3'UTR (as delineated by the red vertical lines).  From this figure, one may conclude that the events (in this case m6A sites) occur througout the gene body with a peak around the stop codon and a precipitous transition from the 5'UTR to the CDS.  However, one caveat is that the three regions of interest are drawn with equal widths.  On average, this is not the case.  We can view the average lenghts in this dataset:

```{r}
summary(data.frame(dist$utr5_size, dist$cds_size, dist$utr3_size))
```

The median lengths are `r median(dist$utr5_size, na.rm = T)`, `r median(dist$cds_size, na.rm = T)`, and `r median(dist$utr3_size, na.rm = T)` for the 5'UTR, CDS and 3'UTR, respectively.

To account for these length differences in the metagene, we can re-scale the widths of the 5'UTR and 3'UTR relative to the CDS (which is set constant to a width of 1 unit).  So first we calculate a simple scale factor (SF):

```{r}
utr5.SF <- median(dist$utr5_size, na.rm = T)/median(dist$cds_size, na.rm = T)
utr3.SF <- median(dist$utr3_size, na.rm = T)/median(dist$cds_size, na.rm = T)
```

The SF for the 5'UTR is `r round(utr5.SF,2)` and for the 3'UTR is `r round(utr3.SF,2)`.  The followign code rescales these regions accordingly:

```{r}
# assign the regions to new dataframes
utr5.dist <- dist[dist$rel_location <= 1, ]
cds.dist <- dist [dist$rel_location <= 2 & dist$rel_location > 1, ]
utr3.dist <- dist[dist$rel_location > 2, ]

# rescale 5'UTR and 3'UTR
library("scales")
utr5.dist$rel_location <- rescale(utr5.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.dist$rel_location <- rescale(utr3.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))

```

Finally, plot the metagene with the rescaled UTRs
```{r}
# Combine and plot
## Histogram
all.regions <- c(utr5.dist$rel_location, cds.dist$rel_location, utr3.dist$rel_location)
hist.data <- hist(all.regions, breaks = 200, col = "black", xlim = c(0,3)) # plot and save to variable
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
```

### Alternate representations of the metagene
#### A line plot
```{r}
plot(hist.data$breaks[1:length(hist.data$breaks)-1], hist.data$counts, type = 'l')
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
```

A smooth density plot
```{r}
plot(density(all.regions))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
```

### Mapping the absolute distance of sites from fixed features
An alternative to the metagene plot, which shows the relative location of sites along a virtual transcipt, is a feature distance plot which shows the absolute distance of sites from a given transcriptomic feature (e.g. stop codon, transcription start site, etc).  As discussed earlier columns 6-11 of the dataset contains the absolute distance data.

```{r}
head(dist[,6:11])
```

For example, we can view the distribution of sites within 100 nucleotides of the stop codon

