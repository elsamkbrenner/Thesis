
## Upset visualization: MAGs in different locations and shared among locations

locationcolors=c('#c4d7d1','#408892','#2d3749','#c04062','#6b3a59','#e08683')

MAGrel <- genome_counts_rel
MAGrel_pa=1*(MAGrel_red>0)
#MAGrel_pa[1:6,1:6]
table_upset_analysis_cont=t(aggregate(t(MAGrel_pa),by=list(cat_metadata_red$Location),FUN=sum)[,-1])
colnames(table_upset_analysis_cont)=levels(as.factor(cat_metadata_red$Location))
table_upset_analysis=(table_upset_analysis_cont>0)*1
table_upset_analysis=data.frame(table_upset_analysis)
table_upset_analysis=apply(table_upset_analysis,2,as.integer)
rownames(table_upset_analysis) <- rownames(MAGrel_pa)

#pdf("figures/MAG_intersection.pdf",width=8,height=6, onefile=F)
upset(as.data.frame(table_upset_analysis),
      keep.order = T,
      sets = rev(c("Aruba","Brazil","CaboVerde","Denmark","Malaysia","Spain")),
      sets.bar.color= rev(locationcolors),
      mb.ratio = c(0.55, 0.45), order.by = "freq")
#dev.off()




## RDA

count_table <- data.frame(physeq_all@otu_table)
count_table_t <- data.frame(t(count_table))
hel1 <- decostand(count_table_t, "hellinger")

physeq_clr <- microbiome::transform(physeq_all, 'clr')

count_crl <- data.frame(physeq_clr@otu_table)
count_crl_t <- data.frame(t(count_crl))
metadata_crl <- data.frame(physeq_clr@sam_data)
metadata_crl <- rownames_to_column(metadata_crl, "sample")
metadata_crl <- metadata_crl[match(rownames(count_crl_t),metadata_crl$sample),]
#mean(rownames(count_crl_t)==metadata_crl$sample)
#metadata_crl$Day=factor(metadata_crl$Day)
design <- metadata_crl[, c("sample", "Origin","Location")]

rda.model <- rda(count_crl_t ~ Origin, data = design)# Run the model
rda.model.hel <- rda(hel1 ~ Origin, data = design)# Run the model
#summary(rda.model)
# Anova
anova(rda.model)
anova(rda.model.hel)
# anova(hel_rda,by="axis")
# summary(hel_rda)
RsquareAdj(rda.model)
RsquareAdj(rda.model.hel)
anova.cca(rda.model, step = 1000)
anova.cca(rda.model, step = 1000, by = "term")
anova.cca(rda.model, step = 1000, by = "axis")
perc <- round(100*(summary(rda.model)$cont$importance[2, 1:2]), 2)
## extract scores - these are coordinates in the RDA space
sc_si <- vegan::scores(rda.model, display="sites", choices=c(1,2), scaling=1)
sc_sp <- vegan::scores(rda.model, display="species", choices=c(1,2), scaling=1)
sc_bp <- vegan::scores(rda.model, display="bp", choices=c(1, 2), scaling=1)

# Prepare data for plotting
rda_sc_wa <- vegan::scores(rda.model,
                           display = "wa",
                           scaling = 1
)
rda_sc_cn <- data.frame(vegan::scores(rda.model,
                                      display = "cn",
                                      scaling = 1
))
rda_sc_sp <- data.frame(vegan::scores(rda.model,
                                      display = "sp",
                                      scaling = 2
))

# Select MAGs from quantile
ASV_Scores_RDA1quantiles <- quantile(rda_sc_sp$RDA1,
                                     probs = c(0.01, 0.99)
)
HighFit_ASVs <- rda_sc_sp[
  rda_sc_sp$RDA1 < ASV_Scores_RDA1quantiles[1] |
    rda_sc_sp$RDA1 > ASV_Scores_RDA1quantiles[2],
]

HighFit_ASVNames <- rownames(HighFit_ASVs)
species_index <- data.frame(
  Species = HighFit_ASVNames,
  Index = 1:length(HighFit_ASVNames)
)

# combine dataframes
all_dataframe <- data.frame(rda_sc_wa, design)

# Plot
ggplot(all_dataframe, aes(y = PC1, x = RDA1)) +
  geom_point(aes(y = PC1, x = RDA1, colour = Origin)) +
  # text
  xlab(paste0("RDA1 (", perc[1], "%)")) +
  ylab(paste0("PC1 (", perc[2], "%)")) +
  #  ggtitle("Caecum microbiome development (AdjR2 = 0.20)") +
  # segment
  geom_segment(data = HighFit_ASVs,
               aes(x = 0, xend = RDA1, y = 0, yend = PC1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "grey") +
  geom_point(data = HighFit_ASVs,
             aes(y = PC1, x = RDA1),
             show.legend = FALSE,
             color = "grey") +
  geom_label_repel(data = HighFit_ASVs,
                   aes(y = PC1, x = RDA1,label = species_index$Species),
                   size = 3,
                   show.legend = FALSE,
                   colour = "grey",
                   max.overlaps = 20) +
  #  scale_colour_manual(values = c("red", "blue")) +
  scale_shape_discrete(solid = TRUE) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    legend.position = "none"
  )

