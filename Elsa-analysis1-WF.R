## Load required libraries
library(ggnewscale) # this package allows you to have multiple fill and color scales when using ggplot2, a common plotting tool in R. 
library(ggtree) # this package creates visual representations of phylogenetic trees
library(R.utils) # this is a utility package, it helps R run smoothly
library(tidyverse) # this is a utility package, it helps R run smoothly
library(ape) # APE stands for 'Analyses of Phylogenetics and Evolution,' thus this package is used for phylogenetic analysis
library(devtools) # this is a utility package, it helps R run smoothly
library(ggplot2) # ggplot2 is a helpful package for making graphs in R. 
library(hillR) # this package is used to calculate Hill numbers, a metric for 'combined' diversity 
library(spaa) # this package holds tools used for the analysis of species associations, and niche overlap
library(vegan) # Vegan is a common package used for the analysis of ecology data. Here, we apply it for its tools on diversity analysis. 
library(hilldiv) # this package is used for the calculation of Hill numbers (see also hillR)

### library(phyloseq) ERROR "installation of package ‘phyloseq’ had non-zero exit status"

library(phytools) #phytools is a package that focuses on comparative phylogenetic analysis 
### library(microbiome) ERROR "installation of package ‘microbiome’ had non-zero exit status," source package is phyloseq

library(matrixStats) # this package is used for opperations performed on matrices

# library(microbiomeutilities) #ERROR "installation of package ‘microbiomeutilities’ had non-zero exit status," source package is phyloseq

library(lme4) # this package is used for fitting and analyzing linear and nonlinear models/ 
library(MuMIn) # MuMIn = Multi-model inference, and contains function for information-theoretic model selection, and model averaging 
library(nlme) # this package contains packages related to nonlinear mixed-effects models
library(knitr) # Knitr is a package for turning R script into reports, used here for the generation of an R markdown document 
library(kableExtra) # Kable builds on Knitr, and adds extra functions
library(pairwiseAdonis) # this is a function of (vegan) and is used for multi-level pairwise comparison
library(sjPlot) # this package is used for visualizing statistical analysis 

# library(distillR) #ERROR package ‘DistillR’ is not available for this version of R

library(RColorBrewer) # this is a series of preset color palettes
library(reshape2) # this package is used for reformatting and "reshaping" data to fit ideal formats 
library(ggpubr) # this package is used for formatting ggplot2 figures to be publication ready. 
library(ggdendro) # this package allows you to apply ggplot2 to make Dendrograms and tree diagrams. 
library(grid) # this package adds gridlines to existing plots
library(gplots)
library(dendextend) # this package allows you to extend dendrograms made in R, for the purposes of visualization and comparison. 
library(stringr) # this is a utility package
library(Rtsne) # this is a utility package 
library(glue) # this is a utility package

## Download the tables from, change the working directory and files names if needed

## in this section you define the names and locations of the files you will be using, and the working directory you will be using for the rest of the script.

## Declare directories and files
workingdir="/Users/elsa/Desktop/THESIS/R-analysis1"
counts_file="/Users/elsa/Desktop/THESIS/R-analysis1/DMB0024_counts.tsv"
tree_file="/Users/elsa/Desktop/THESIS/R-analysis1/DMB0024.tree"
taxonomy_file="/Users/elsa/Desktop/THESIS/R-analysis1/DMB0024_mag_info.tsv"
metadata_file="/Users/elsa/Desktop/THESIS/R-analysis1/DMB0024_metadata.tsv"
coverage_file="/Users/elsa/Desktop/THESIS/R-analysis1/DMB0024_coverage.tsv"

# ELSA NOTE: note, my files are just .tsv, not .tsv.gz

## Load data 
setwd(workingdir) # this line sets the working directory of the script to the location established in the block above, where i defined 'workingdir'

# the lines below take the tables imported above, and turn them into readable table files we can use for the rest of our analysis.  
counts_table <- read.table(counts_file, sep="\t",row.names=1,header=T)
coverage_table <- read.table(coverage_file, sep="\t",row.names=1,header=T)
metadata <- read.table(metadata_file,sep="\t",header=T)%>%
  rename(sample=EHI_plaintext)
mags_table <- read.table(taxonomy_file,sep="\t",header=T)
tree <- read.tree(tree_file)

# ELSA NOTE:  here I changed the lines a bit because my files did not need to be unzipped. 

# Load EHI taxonomy colours. Download it first if you have not them in your computer
#colours_URL="https://raw.githubusercontent.com/earthhologenome/EHI_taxonomy_colour/main/ehi_phylum_colors.tsv"
#download.file(colours_URL, "ehi_phylum_colors.tsv")

# in this section we upload the EHI's standard color palatte. This is so that all EHI analyses follow a similar format, for ease of readability and comprehension. 

ehi_phylum_colors <- read.table("~/Desktop/Thesis/R-analysis1/ehi_phylum_colors.tsv",sep="\t",header=T,comment.char = "")

#Delete *__ from the taxonomy names
ehi_phylum_colors1 <- ehi_phylum_colors %>%
  mutate_at(vars(phylum), ~ str_replace(., "[dpcofgs]__", ""))
taxonomyclean <- mags_table%>%
  mutate_at(vars(domain,phylum,class,order,family,genus,species), ~ str_replace(., "[dpcofgs]__", ""))

nsamples <- ncol(counts_table) # this line defines our number of samples, by reading from the counts table
metagenomic_bases <- sum(metadata$metagenomic_bases) #this line defines the number of metagenomic bases, by summing each sample from the metadata file. 
host_bases <- sum(metadata$host_bases) #this line defines the number of host bases, from the sum from each sample from the metadata file
discarded_bases <- sum(round(((metadata$metagenomic_bases+metadata$host_bases)/(1-metadata$bases_lost_fastp_percent))-(metadata$metagenomic_bases+metadata$host_bases))) # this line defines the number of discarded bases 
total_bases <- discarded_bases + host_bases + metagenomic_bases # this line defines the number of total bases
singlem_bases <- sum(metadata$metagenomic_bases * metadata$singlem_fraction) 
nmags <- nrow(counts_table) # this line defines the number of MAGs
new_species <- mags_table %>%
  filter(species == "s__") %>%
  nrow()

sequencing_depth <- colSums(counts_table) # this line defines a variable for sequencing depth
sequencing_depth_sum <- sum(sequencing_depth) # this line defines a variable for total sequencing depth sum 
sequencing_depth_mean <- mean(sequencing_depth) # this line defines the variable for mean sequencing depth
sequencing_depth_sd <- sd(sequencing_depth) # this line defines a variable for the standard deviation for the mean sequencing depth. 

cat(nsamples) # this tells you the number of samples included in analysis = 58

# cat() is used to concatenate strings. Thus, cat(nsamples) will tell us the number of samples in a readable integer form

cat(nmags) # this tells you the number of metagenome-assembled genomes (MAGs) included in analysis = 472

totalgb <- round(total_bases / 1000000000,2) #this line defines a variable, 'totalgb,' to be the total number of bases, written in the unit Gigabases. The ',2' denotes that the output is round to 2 decimal points. 

cat(totalgb) # total DNA sequenced = 333.71 gigabases

discardgb <- round(discarded_bases / 1000000000,2)
cat(discardgb) # total amount of discarded data (low quality, lack of info) = 10.37 Gigabases

discarddata <- round(discarded_bases / total_bases * 100,2) #similar to the above calculations that convert to Gigabases, this line converts to a percentage value. 
cat(discarddata) # amount of discarded data = 3.11% of total data. 3.11<5%, therefore good! and within expected range. 

hostGB <- round(host_bases / 1000000000,2)
cat(hostGB) # amount of host data = 6.49 gigabases

hostdata <- round(host_bases / (total_bases-discarded_bases) * 100,2)
cat(hostdata) # host data represents 2.01% of quality filtered data

prokaGB <- round(singlem_bases / 1000000000,2)
cat(prokaGB) # an estimated amount of prokaryotic data = 293.14 gigabases of data

prokadata <- round(singlem_bases / (metagenomic_bases) * 100,2)
cat(prokadata) # prokaryotic data is estimated to be 92.51% of metagenomic data. 

metaGB <- round(metagenomic_bases / 1000000000,2)
cat(metaGB) # =316.85 gigabases of metagenomic data

metaperce <- round(metagenomic_bases / (total_bases-discarded_bases) * 100,2)
cat(metaperce) # = 97.99%  metagenomic data

totalreads <- round(sequencing_depth_sum / 1000000,2)
cat(totalreads) # = 1833.71 million reads

mappedGB <- round(sequencing_depth_sum / 1000000000 * 143,2)
cat(mappedGB) # total mapped sequencing depth = 262.22 gigabases

meanreads <- round(sequencing_depth_mean / 1000000,2)
cat(meanreads) # average mapped sequencing depth = 31.62 million reads

meanGB <- round(sequencing_depth_mean / 1000000000 * 143,2)
cat(meanGB) # average mapped sequencing depth = 4.52 gigabases 

phyla <- ehi_phylum_colors1 %>%
  right_join(taxonomyclean, by=join_by(phylum == phylum)) %>% 
  arrange(match(genome, tree$tip.label)) %>% 
  select(phylum, colors) %>%
  unique()

heatmap <- ehi_phylum_colors1 %>%
  right_join(taxonomyclean, by=join_by(phylum == phylum)) %>%
  arrange(match(genome, tree$tip.label)) %>%
  select(genome,phylum) %>%
  mutate(phylum = factor(phylum, levels = unique(phylum))) %>%
  column_to_rownames(var = "genome")

colors_alphabetic <- ehi_phylum_colors1 %>%
  right_join(taxonomyclean, by=join_by(phylum == phylum)) %>%
  arrange(match(genome, tree$tip.label)) %>%
  select(phylum, colors) %>%
  unique() %>%
  arrange(phylum) %>%
  select(colors) %>%
  pull()

circular_tree <- force.ultrametric(tree,method="extend") %>%
  ggtree(., layout = 'circular', size = 0.3)

circular_tree <- gheatmap(circular_tree, heatmap, offset=0.65, width=0.1, colnames=FALSE) +
  scale_fill_manual(values=colors_alphabetic) +
  geom_tiplab2(size=1, hjust=-0.1) +
  theme(plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0))
