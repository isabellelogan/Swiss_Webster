# Load packages.

library (ape)
library (ecodist)
library (ggplot2)
library (vegan)
library (dplyr)
library (scales)
library (grid) 
library (reshape2)
library (knitr)
library (BiocStyle)
library (gridExtra)
library (phyloseq)
library (caret)
library (igraph)
library (ergm)
library (randomForest)
library (ggnetwork)
library (phyloseqGraphTest)
library (intergraph)
library (gridExtra)
library (flextable)
library (ggraph)
library (igraph)
library (tidyverse)
library (ggpubr)
library (rstatix)
library (readxl)
library (writexl)
library (RColorBrewer)

# Read in the phyloseq data sets.
phyloObjR_Second <- readRDS("/Users/isabelle/16S/Bbduk/phyloObjR_Second.rds")

###############################################################################################
### Kendall Tau.
###############################################################################################

###
# For the start of the second feeding study
###

# Now we need to keep only ASVs that are present in at least 100% of samples

# agglomerate

phylo.glom <- tax_glom (phyloObjR_Start_Sec, taxrank = "Genus")

# Use genus as taxa names

taxa_names(phylo.glom) <- tax_table(phylo.glom)[,6]

# Export ASV table

phylo.otu <- as.data.frame(otu_table(phyloObjR_Start_Sec))

# Make binary

phylo.otu[phylo.otu >= 1] <- 1

# Remove any not present in 100% of samples

phylo.hundred.start <- phylo.otu[,colSums(phylo.otu != 0) >= 29] # 29 is 100% of 29

# Extract relevant ASV names

keep.otu <- colnames(phylo.hundred.start)

# Generate new phyloseq object

phyloObjR.Start_100 <- pop_taxa(phyloObjR_Start_Sec, keep.otu)

phyloObjR.Start_100

my.Taxa = names(sort(taxa_sums(phyloObjR.Start_100), decreasing = TRUE))
tax1 = prune_taxa(my.Taxa, phyloObjR.Start_100)
replace <- c(rep(3, 9), rep (1, 10), rep (2, 10))

###
# How about ASVs in 100% of samples for the second feeding study by relative abundance
###

my.Taxa = names(sort(taxa_sums(phyloObjR.hundred.rel), decreasing = TRUE))
tax1 = prune_taxa(my.Taxa, phyloObjR.hundred.rel)
replace <- c(rep(3, 9), rep (1, 10), rep (2, 10))

###
# How about ASVs in 100% of samples for the second feeding study, rel. abundance, each family by itself
###

phylo.hundred.fam = tax_glom(phyloObjR_Second_Rel, "Family")
my.Taxa = names(sort(taxa_sums(phylo.hundred.fam), decreasing = TRUE))
tax1 = prune_taxa(my.Taxa, phylo.hundred.fam)
replace <- c(rep(3, 9), rep (1, 10), rep (2, 10))

tax1.r <- tax1

# Change taxa names to something more palatable

taxa_names(tax1.r) <- paste0("ASV_", seq(ntaxa(tax1.r)))

# Change Diet factor to numbers

replace <- as.factor(replace)
sample_data(tax1.r)$Diet = factor(replace)
sample_data(tax1.r)$Diet = as.numeric (sample_data(tax1.r)$Diet)

# Check correlations

OTU1_ent = as(otu_table(tax1.r), "matrix")
# if(taxa_are_rows(tax1.r)){SV1_ent <- t(OTU1_ent)}
OTU1_ent_df = as.data.frame(OTU1_ent)
# View(OTU1_ent_df)
OTU1_ent_df <- t(OTU1_ent_df)

# Now perform the test

kruskal.results <- list(0)

for (i in 1:nrow(OTU1_ent_df)) {
  kruskal.results[[i]] <- kruskal.test (OTU1_ent_df[i,], sample_data(tax1.r)$Diet)
}

# Now we need to pull out all the p-values

kruskal.w.p <- c(rep(0, length(kruskal.results)))

for (i in 1:length(kruskal.results)) {
  kruskal.w.p[i] <- kruskal.results[[i]]$p.value
}

# False discovery rate adjustment
kruskal.w.p <- p.adjust(kruskal.w.p, method = "fdr")

# Which taxa have a p-value smaller than 0.05?

kruskal.select <- NULL
kruskal.select <- which (kruskal.w.p < 0.05)
kruskal.select

# Prepare data to put into the plotting function

tax1.r.psdf <- data.table::data.table(psmelt(tax1.r))

tax1.r.psdf <- as_tibble (tax1.r.psdf)

# Select OTUs, this only works for all samples

otu.select <- NULL
for (i in 1:length(kruskal.select)) {
  otu.select[i] <- c(paste0("ASV_", kruskal.select[i]))
}
otu.select

# Extract the actual p-values for the OTUs

kruskal.p <- kruskal.w.p[kruskal.select]
kruskal.p

# Perform Mann-Whitney U tests
# From here: https://rcompanion.org/rcompanion/d_06.html
# You have to do each one at a time for each pairwise comparison

p.h <- c(1:3)
p.h.hfd <- as.list(p.h)

data.test.lfd <- OTU1_ent_df[, 1:9]
data.test.hfd <- OTU1_ent_df[, 10:19]

for (i in 1:length(kruskal.select)) {
  p.h.hfd[[i]] <- wilcox.test (data.test.lfd[kruskal.select[i],], 
                               data.test.hfd[kruskal.select[i],], 
                               # comparisons = list(c("HFD", "HXN")), #c("LFD", "HFD")), 
                               p.adjust.method = "fdr", paired = FALSE, exact = FALSE)
}
p.h.values <- c(0,0,0)
for (i in 1:length(kruskal.select)) {
  p.h.values[i] <- p.h.hfd[[i]]$p.value
}
p.h.values
matrix.out <- rbind (kruskal.select, p.h.values)
matrix.out.hfd <- t(matrix.out)
matrix.out.hfd

data.test.hxn <- OTU1_ent_df[,20:29]

for (i in 1:length(kruskal.select)) {
  p.h.hxn[[i]] <- wilcox.test (data.test.hfd[kruskal.select[i],], 
                               data.test.hxn[kruskal.select[i],],
                               # comparisons = list(c("HFD", "HXN")), #c("LFD", "HFD")), 
                               p.adjust.method = "fdr", paired = FALSE, exact = FALSE)
}
p.h.values <- c(0,0,0)
for (i in 1:length(kruskal.select)) {
  p.h.values[i] <- p.h.hxn[[i]]$p.value
}
p.h.values

matrix.out <- rbind (kruskal.select, p.h.values)
matrix.out.hxn <- t(matrix.out)
matrix.out.hxn

# Build a table to add the p-values to the graph for the second study

group1 <- c(rep(c("HFD", "HFD", "HXN"),length(p.h.hfd)))
group2 <- c(rep(c("HXN", "LFD", "LFD"),length(p.h.hfd)))
p.val.hfd <- list()
for (i in 1:length(p.h.hfd)) {
  # p.val [[i]] <- c(p.h[[i]]$res$P.adj) # old way 
  p.val.hfd [[i]] <- c(p.h.hfd[[i]]$p.value) # mann whitney way
}
p.val.hfd <- unlist (p.val.hfd, use.names = FALSE)

p.val.hxn <- list()
for (i in 1:length(p.h.hfd)) {
  # p.val [[i]] <- c(p.h[[i]]$res$P.adj) # old way 
  p.val.hxn [[i]] <- c(p.h.hxn[[i]]$p.value) # mann whitney way
}
p.val.hxn <- unlist (p.val.hxn, use.names = FALSE)

# Combine vectors
p.combined <- c(rep(1, (length(kruskal.select)*3)))
# HXN first
p.combined [seq(1, length(p.combined), 3)] <- c(p.val.hxn)
# HFD next
p.combined [seq(2, length(p.combined), 3)] <- c(p.val.hfd)
# No need to add dummies!

# This is our table to plot p values onto the graph!

dun.table <- cbind (group1, group2, p.combined)
colnames(dun.table) <- c("group1", "group2", "p.adj")
dun.table <- as_tibble(dun.table)
dun.table$p.adj <- as.numeric(dun.table$p.adj)
dun.table$p.adj <- round(dun.table$p.adj, 3)

# Adjust p values (prints better this way)
for (x in 1:(length(kruskal.select)*3)) {
  if (dun.table$p.adj[x] < 0.001) {
    dun.table$p.adj[x] = 0.001
  }
}

# Split the table into chunks of 3, change as needed.
dun.list <- split(dun.table, rep (1:(nrow(dun.table)/3), each = 3))

# Build a function that plots the graphs

plot_kruskal = function (taxon.table, OTU.in, replace = c(rep("LFD", 9), rep ("HFD", 10), rep ("HXN", 10)),
                         levels = c("HFD", "HXN", "LFD"), p.value = 0.05, fill = c("grey53", "green3", "white"),
                         line = c("black"), dunn.p) {
  keep.otu <- taxon.table %>% filter (OTU == OTU.in)
  keep.otu <- as.tibble (keep.otu)
  keep.otu <- keep.otu %>% arrange(Sample)
  # Set the factors for the diets back
  replace <- as.factor(replace)
  keep.otu$Diet = factor(replace)
  # This line sets the order of the diets
  keep.otu <- keep.otu %>% dplyr::mutate (Diet = factor(Diet, levels = levels))
  # Set y-axis limit
  y.limit = (1.6*round(max(keep.otu$Abundance), 2))
  # Set the p-value text
  mylabel = bquote(Kruskal~Wallis~italic(p) == .(round(p.value, digits = 3)))
  # mytitle = bquote(.(keep.otu$Phylum[1])~italic(keep.otu$Family[1])~italic(keep.otu$Genus[1]))
  # Add the location of the labels to the dunn test
  dunn.p <- as.data.frame (dunn.p)
  # Change this to the number of comparisons made (6 for all combined)
  y.location <- c((0.925*y.limit), (0.825*y.limit), (0.725*y.limit) )#, (0.80*y.limit), (0.75*y.limit), (0.70*y.limit))
  dunn.p <- dunn.p %>% add_column (y.location)
  dunn.p$p.adj <- round (dunn.p$p.adj, 3)
  colnames(dunn.p) <- c("group1", "group2", "p.adj", "y.position")
  # Change the 6 to 3 for the individual studies
  for (i in 1:3) {
    if (dunn.p$p.adj[i] >= 0.05) {
      dunn.p$p.adj[i] <- NA
    }
  }
  dunn.p <- dunn.p %>% drop_na()
  if (nrow(dunn.p) == 0) {
    p <- ggplot(keep.otu, aes(x = Diet, y = Abundance)) + 
      geom_boxplot (fill = fill, colour = line, alpha = 0.7, outlier.shape = NA, shape = 21) +
      geom_jitter (position = position_jitter(0.2), shape = 21) +
      theme_classic() +
      scale_x_discrete (name = "Diet") +
      scale_y_continuous (name= "Relative Abundance (%)\n", breaks = seq(0, y.limit, y.limit/10), 
                          labels = waiver(), limits = c(0, y.limit)) +
      guides (fill = guide_legend(keywidth = 1, keyheight = 1)) +
      ggtitle (paste (keep.otu$Phylum[1], keep.otu$Family[1], keep.otu$Genus[1], sep = " ")) +
      annotate ("text", x = "HFD", y = (y.limit), label = mylabel) 
  }
  else {
    p <- ggplot(keep.otu, aes(x = Diet, y = Abundance)) + 
      geom_boxplot (fill = fill, colour = line, alpha = 0.7, outlier.shape = NA, shape = 21) +
      geom_jitter (position = position_jitter(0.2), shape = 21) +
      theme_classic() +
      scale_x_discrete (name = "Diet") +
      scale_y_continuous (name= "Relative Abundance (%)\n", breaks = seq(0, y.limit, y.limit/10), 
                          labels = waiver(), limits = c(0, y.limit)) +
      guides (fill = guide_legend(keywidth = 1, keyheight = 1)) +
      ggtitle (paste (keep.otu$Phylum[1], keep.otu$Family[1], keep.otu$Genus[1], sep = " ")) +
      annotate ("text", x = "HFD", y = (y.limit), label = mylabel) +
      stat_pvalue_manual(data = dunn.p, label = "p.adj", xmin = "group1", xmax = "group2", 
                         y.position = "y.position")#, tip.length = 0)
  }
  p <- p + theme (axis.text = element_text (size = rel (1.75))) +
    theme (axis.title = element_text (size = rel (1.75))) +
    theme (legend.text = element_text (size = rel (1.50))) +
    theme (legend.title = element_text (size = rel (1.50))) +
    scale_y_continuous(name= "Relative Abundance (%)\n",
                       labels = scales::number_format(accuracy = 0.1), breaks = scales::breaks_extended(n = 10)) +
    # scale_y_continuous(breaks = integer_breaks()) + # not a great solution
    # scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + # this makes y-axis pretty
    rremove ("xlab")
  plot(p)
}

# And plot (at genus level)
for (i in 1:length(kruskal.p)) {
  plot_kruskal(taxon.table = tax1.r.psdf, OTU.in = otu.select[i], p.value = kruskal.p[i], dunn.p = dun.list[[i]])
  ggsave(filename = paste0("Sec_Kruskal_Eighty_Gen_New",i,".png"))
}

# And plot (at family level)
for (i in 1:length(kruskal.p)) {
  plot_kruskal(taxon.table = tax1.r.psdf, OTU.in = otu.select[i], p.value = kruskal.p[i], dunn.p = dun.list[[i]])
  ggsave(filename = paste0("Sec_Kruskal_Eighty_Fam_New",i,".png"))
}




###############################################################################################
### Alpha Diversity
###############################################################################################


# Now we plot alpha diversity metrics (all the extra colors are not needed)
# Second study only
# With dots:
plot_richness(phyloObjR_Second, x="Diet", measures=c("Shannon"), 
              color="Diet") + theme_classic() + 
  scale_color_manual (values = c("grey53", "green3", "black"))

# Boxplot:
# Came from https://tips.cgrb.oregonstate.edu/posts/phyloseq-bug-meeting-presentation-fall-2019/
richness.rare <- cbind(estimate_richness(phyloObjR_Second, 
                                         measures = c('Shannon', 'Simpson', 'Chao1')),
                       sample_data(phyloObjR_Second)$Diet)
colnames(richness.rare) <- c('Shannon', 'Simpson', 'Chao1', 'Diet')
richness.rare$Labels <- rownames(richness.rare)

ggplot(data = richness.rare, aes(x = Shannon, y = Simpson)) + 
  geom_point()

# I don't really see an outlier, so I'm not going to remove any samples

ad.test.df <- richness.rare[,c('Shannon', 'Simpson', 'Chao1')]
ad.test.df <- cbind(ad.test.df,
                    sample_data(phyloObjR_Second)[,'Diet'])
colnames(ad.test.df) <- c('Shannon', 'Simpson', 'Chao1', 'Diet')
alpha_test <- kruskal.test(Shannon ~ Diet, data=ad.test.df)

# Kruskal-Wallis rank sum test

# data:  Shannon by Diet
# Kruskal-Wallis chi-squared = 11.029, df = 2, p-value = 0.004028

kruskal.test(Simpson ~ Diet, data = ad.test.df)

# Kruskal-Wallis rank sum test

# data:  Simpson by Diet
# Kruskal-Wallis chi-squared = 12.303, df = 2, p-value = 0.002131

# p value label
mylabel = bquote(Kruskal~Wallis~italic(p) == .(round(alpha_test$p.value, digits = 3)))
y.limit = 1.15 * max (ad.test.df$Shannon)

# Now we need to perform pairwise comparisons, first for shannon

krusk.table <- ad.test.df %>%
  dunn_test (Shannon ~ Diet, p.adjust.method = "none")
krusk.table <- krusk.table %>% add_xy_position(x = "Diet")

# First LFD-HFD
ad.test.df <- as_tibble(ad.test.df)
HFD.test <- ad.test.df %>%
  filter (Diet == "LFD" | Diet == "HFD")
HFD.test$Diet <- as.numeric(HFD.test$Diet)
lfd.krusk.results <- wilcox.test (HFD.test$Shannon, 
                                  HFD.test$Diet, 
                                  p.adjust.method = "none", paired = FALSE, exact = FALSE)
#summary(lfd.krusk.results)
lfd.krusk.results$p.value
# Then HFD-HXN
HXN.test <- ad.test.df %>%
  filter (Diet == "HFD" | Diet == "HXN")
HXN.test$Diet <- as.numeric(HXN.test$Diet)
hxn.krusk.results <- wilcox.test (HXN.test$Shannon, 
                                  HXN.test$Diet, 
                                  p.adjust.method = "none", paired = FALSE, exact = FALSE)
#summary(hxn.krusk.results)
hxn.krusk.results$p.value

# Make the table for the p values for the plot

krusk.table[1,8] <- lfd.krusk.results$p.value
krusk.table[2,8] <- hxn.krusk.results$p.value
krusk.table[3,8] <- 1
krusk.table[2,9] <- "***"
krusk.table$p.adj <- signif (krusk.table$p.adj, digits = 3)
krusk.table <- krusk.table %>%
  mutate (to.plot = p.adj)
krusk.table$to.plot <- round (krusk.table$to.plot, digits = 3)

# Fix p values
for (x in 1:3) {
  if (krusk.table$to.plot[x] < 0.001) {
    krusk.table$to.plot[x] = "<0.001"
  }
}

shannon.plot <- ggboxplot(
  ad.test.df, x = "Diet", y = "Shannon",
  color = "black",
  fill = "Diet", palette = c("grey53", "grey72", "white"), ylab = ("Shannon Index (a.u.)"),
  add = "jitter", outlier.shape = NA, shape = 21) + #, shape = "Diet") # add this in for different shapes
  theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  rremove ("xlab")
shannon.plot
# add p values
shannon.plot <- shannon.plot + stat_pvalue_manual(krusk.table, label = "to.plot", hide.ns = TRUE)
# add overall value
shannon.plot <- shannon.plot + theme(legend.position = 'top') +
  annotate ("text", x = "HFD", y = y.limit, label = mylabel)

shannon.plot
ggsave("Alpha_diversity_Shannon.tiff")
ggsave("Alpha_diversity_Shannon.png")
ggsave("Alpha_diversity_Shannon_color.png")

# Now we need to perform pairwise comparisons, now for chao1
alpha_test <- kruskal.test(Chao1 ~ Diet, data=ad.test.df)

# p value label
mylabel = bquote(Kruskal~Wallis~italic(p) == .(round(alpha_test$p.value, digits = 3)))
y.limit = 1.15 * max (ad.test.df$Chao1)

krusk.table <- ad.test.df %>%
  dunn_test (Chao1 ~ Diet, p.adjust.method = "none")
krusk.table <- krusk.table %>% add_xy_position(x = "Diet")

# First LFD-HFD
ad.test.df <- as_tibble(ad.test.df)
HFD.test <- ad.test.df %>%
  filter (Diet == "LFD" | Diet == "HFD")
HFD.test$Diet <- as.numeric(HFD.test$Diet)
lfd.krusk.results <- wilcox.test (HFD.test$Chao1, 
                                  HFD.test$Diet, 
                                  p.adjust.method = "none", paired = FALSE, exact = FALSE)
#summary(lfd.krusk.results)
lfd.krusk.results$p.value
# Then HFD-HXN
HXN.test <- ad.test.df %>%
  filter (Diet == "HFD" | Diet == "HXN")
HXN.test$Diet <- as.numeric(HXN.test$Diet)
hxn.krusk.results <- wilcox.test (HXN.test$Chao1, 
                                  HXN.test$Diet, 
                                  p.adjust.method = "none", paired = FALSE, exact = FALSE)
#summary(hxn.krusk.results)
hxn.krusk.results$p.value

# Make the table for the p values for the plot

#krusk.table <- read_excel("/Users/isabelle/16S/Raw/Kruskal.xlsx")
#krusk.table <- as_tibble (krusk.table)
#krusk.table
krusk.table[1,8] <- lfd.krusk.results$p.value
krusk.table[2,8] <- hxn.krusk.results$p.value
krusk.table[3,8] <- 1
krusk.table[2,9] <- "***"
krusk.table$p.adj <- signif (krusk.table$p.adj, digits = 3)
krusk.table <- krusk.table %>%
  mutate (to.plot = p.adj)
krusk.table$to.plot <- round (krusk.table$to.plot, digits = 3)

# Fix p values
for (x in 1:3) {
  if (krusk.table$to.plot[x] < 0.001) {
    krusk.table$to.plot[x] = "<0.001"
  }
}

chao1.plot <- ggboxplot(
  ad.test.df, x = "Diet", y = "Chao1",
  color = "black",
  fill = "Diet", palette = c("grey53", "grey72", "white"), ylab = ("Chao1 Index (a.u.)"),
  add = "jitter", outlier.shape = NA, shape = 21) + # shape = "Diet") + # add this in for different shapes
  theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  rremove ("xlab")
chao1.plot
# add p values
chao1.plot <- chao1.plot + stat_pvalue_manual(krusk.table, label = "to.plot", hide.ns = TRUE)
# add overall value
chao1.plot <- chao1.plot + theme(legend.position = 'top') +
  annotate ("text", x = "HFD", y = 4, label = mylabel)
chao1.plot
ggsave("Alpha_diversity_Chao1.png")
ggsave("Alpha_diversity_Chao1_new.png")
ggsave("Alpha_diversity_Chao1_color.png")


###############################################################################################
### Beta Diversity (un- and weighted unifrac)
###############################################################################################


# Plot unweighted unifrac for the second study
unifrac_dist_second = phyloseq::distance(phyloObjR_Second, method="unifrac", weighted=F)
u.ordination_second = ordinate(phyloObjR_Second, method="PCoA", distance=unifrac_dist_second)
q1 <- plot_ordination(phyloObjR_Second, u.ordination_second, color="Diet", shape = "Diet", # added shape = "Diet"
                      title="Unweighted Unifrac PCoA of Rarefied Data")+#,
  #shape="Study_number") +
  theme_classic() + 
  geom_point (size = 3) + # removed shape = 16, 
  scale_color_manual(values=c("grey53", "green3", "black")) + # original color was grey72
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 10)) +
  stat_ellipse(geom = "polygon", type="norm", level = 0.95, alpha = 0.001)
q1$data$Diet <- factor(q1$data$Diet, levels = desired.order)
q1 <- q1 + theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  theme (plot.title = element_text (size = rel (1.75)))
plot(q1)
ggsave('Unwei_PCoA_second_new.png')
ggsave('Unwei_PCoA_second_new_color.png')

# Let's do some statistics!
# Adonis
u.second <- adonis (unifrac_dist_second ~ sample_data(phyloObjR_Second)$Diet)
u.second

# Anosim on all the data
u.sim.second <- anosim (unifrac_dist_second, sample_data(phyloObjR_Second)$Diet, permutations = 999)
u.sim.second
# R = 0.7668
# p = 0.001

# Anosim on LFD - HFD only
phylo.lfd <- prune_samples (sample_names(phyloObjR_Second) < "121", phyloObjR_Second)
unifrac_dist_lfd = phyloseq::distance(phylo.lfd, method="unifrac", weighted=F)
u.ordination_lfd = ordinate(phylo.lfd, method="PCoA", distance=unifrac_dist_lfd)
u.sim.lfd <- anosim (unifrac_dist_lfd, sample_data(phylo.lfd)$Diet, permutations = 999) 
u.sim.lfd
# R = 0.4483
# p = 0.001

# Anosim on HFD - HXN only
phylo.hfd <- prune_samples (sample_names(phyloObjR_Second) > "110", phyloObjR_Second)
unifrac_dist_hfd = phyloseq::distance(phylo.hfd, method="unifrac", weighted=F)
u.ordination_hfd = ordinate(phylo.hfd, method="PCoA", distance=unifrac_dist_lfd)
u.sim.hfd <- anosim (unifrac_dist_hfd, sample_data(phylo.hfd)$Diet, permutations = 999) 
u.sim.hfd
# R = 0.9147
# p = 0.001

# Plot weighted unifrac for the second study
desired.order <- c("HFD", "HXN", "LFD")
unifrac_dist_second = phyloseq::distance(phyloObjR_Second, method="unifrac", weighted=T)
u.ordination_second = ordinate(phyloObjR_Second, method="PCoA", distance=unifrac_dist_second)
q1 <- plot_ordination(phyloObjR_Second, u.ordination_second, color="Diet", shape = "Diet", # added shape = "Diet"
                      title="Weighted Unifrac PCoA of Rarefied Data") +
  theme_classic() + 
  geom_point (size = 3) + # removed shape = 16
  scale_color_manual(values=c("grey53", "green3", "black")) + # old was grey72
  theme(legend.position="bottom") +
  stat_ellipse(geom = "polygon", type="norm", level = 0.95, alpha = 0.001)
q1$data$Diet <- factor(q1$data$Diet, levels = desired.order)
q1 <- q1 + theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  theme (plot.title = element_text (size = rel (1.75)))
plot(q1)
ggsave('Wei_PCoA_second_new.png')
ggsave('Wei_PCoA_second_new_color.png')

# Let's do some statistics!
# Adonis
u.second <- adonis (unifrac_dist_second ~ sample_data(phyloObjR_Second)$Diet)
u.second

# Anosim on all the data
u.sim.second <- anosim (unifrac_dist_second, sample_data(phyloObjR_Second)$Diet, permutations = 999)
u.sim.second
# R = 0.4977
# p = 0.001

# Anosim on LFD - HFD only
phylo.lfd <- prune_samples (sample_names(phyloObjR_Second) < "121", phyloObjR_Second)
unifrac_dist_lfd = phyloseq::distance(phylo.lfd, method="unifrac", weighted=T)
u.ordination_lfd = ordinate(phylo.lfd, method="PCoA", distance=unifrac_dist_lfd)
u.sim.lfd <- anosim (unifrac_dist_lfd, sample_data(phylo.lfd)$Diet, permutations = 999) 
u.sim.lfd
# R = 0.3506
# p = 0.001

# Anosim on HFD - HXN only
phylo.hfd <- prune_samples (sample_names(phyloObjR_Second) > "110", phyloObjR_Second)
unifrac_dist_hfd = phyloseq::distance(phylo.hfd, method="unifrac", weighted=T)
u.ordination_hfd = ordinate(phylo.hfd, method="PCoA", distance=unifrac_dist_lfd)
u.sim.hfd <- anosim (unifrac_dist_hfd, sample_data(phylo.hfd)$Diet, permutations = 999) 
u.sim.hfd
# R = 0.4560
# p = 0.001


###############################################################################################
### Conditional Correspondence Analysis
###############################################################################################


phylo.cca.diet <- ordinate (phyloObjR_Second ~ Diet, "CCA")

a <- plot_ordination(phyloObjR_Second, phylo.cca.diet, color="Diet", shape = "Diet", 
                     title="Conditional Correspondence Analysis, Constrained by Diet") +
  theme_classic() + 
  geom_point(size = 3) + # removed shape = 16, 
  stat_ellipse(geom = "polygon", type="norm", level = 0.95, alpha = 0.001) +
  scale_color_manual(values=c("grey53", "green3", "black")) + #original was green3 now grey72
  theme(legend.position="bottom")
a$data$Diet <- factor(a$data$Diet, levels = desired.order)
a <- a + theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  theme (plot.title = element_text (size = rel (1.50)))
plot(a)
ggsave('CCA_Second_Study.png')
ggsave('CCA_Second_Study_Color.png')

# Get a p-value for the beta-dispersion
metadata <- as (sample_data(phyloObjR_Second), "data.frame")
distances_data <- vegdist(phyloseq::distance(phyloObjR_Second, method="bray"))
anova(betadisper(distances_data, metadata$Diet))

# Perform capscale function
cap <- capscale (phyloseq::distance(phyloObjR_Second, method="bray") ~ Diet,
                 data = metadata)
plot(cap)
summary(cap)
anova(cap)

# Perform adonis function
adonis (phyloseq::distance(phyloObjR_Second, method="bray") ~ Diet,
        data = metadata)
# p = 0.001
# R squared = 0.3411


###############################################################################################
### CCA statistics
###############################################################################################


##############################
## Convert phyloseq to vegan
##############################

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Do this for the second study
pss <- pssd2veg(phyloObjR_Second)
pss
psotu <- psotu2veg(phyloObjR_Second)
psotu

## Now use vegan!!!!
# Good info https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html

mod0 <- cca(psotu ~ 1, pss)
mod1 <- cca(psotu ~ ., pss)
mod1 <- cca(psotu ~ Diet, pss)
mod1 <- cca(psotu ~ Diet + Study_number + Diet:Study_number, pss)
m <- ordistep(mod0, scope = formula(mod1), perm.max = 200, direction ="forward")
m

#check for variance of inflation in model

vif.cca(m)

# Test of CCA result

anova.cca(m, step = 1000)

# Test of all canonical axes

anova.cca(m, by='axis', step=1000)

# let's make some plots

plot (m, scaling = 1, display = c("lc", "cn"), main = "Biplot CCA-scaling 1")
plot (mod1, scaling = 1, display = c("lc", "cn"), main = "Biplot CCA-scaling 1")

screeplot(m)

plot(m, scaling=1, display=c('sp', 'lc', 'cn'), main='Triplot CCA matrix ~ env -scaling 1')

# Scaling 2 (default): site scores scaled to relative eigenvalues, species are weighted averages of the sites  
plot(m, display=c('sp', 'lc', 'cn'), main="Triplot CCA matrix ~ env -scaling 2")

# set the seed for consistent results

set.seed (998)

# Test the significance of the global model

anova (mod1, step = 1000)

# Test the significance of each axis

anova (mod1, by = 'axis', step = 1000)
# The first two axes are statistically significant

# Test the significance of each explanatory variable

anova (mod1, by = 'terms')

# Diet is significant

# Let's do a forward selection using the function ordistep()
# Comparing variables is based on AIC criteria and p-values from the Monte Carlo
# permutation test

ordi <- ordistep(cca(psotu ~ 1, data = pss), scope = formula (mod1), direction = "forward",
                 pstep = 1000)

# Let's plot the eigenvalues from the ordistep results

evplot <- function (ev) { # broken stick model
  n <- length (ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n+1-i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow = c(2,1))
  barplot (ev, main = "Eigenvalues", col = "bisque", las =2)
  abline (h=mean(ev), col="red")
  legend ("topright", "Average eigenvalue", lwd =1, col=2, bty="n")
  barplot (t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside = TRUE,
           main ="% variation", col = c("bisque", 2), las = 2)
  legend("topright", c("% eigenvalue", "Broken stick model"),
         pch = 15, col = c("bisque", 2), bty="n")
  par(op)
}

# We can do this either for the CA or the CCA eigenvalues
ev <- ordi$CA$eig
evplot (ev)


###############################################################################################
### NMDS on the PICRUSt2 analysis
###############################################################################################


# Import pathway data

pathways <- read_excel("/Users/isabelle/16S/Picrust/picrust2_out_pipeline/pathways_out/Pathways.xlsx")
pathways <- as.data.frame (pathways)
pathways2 <- pathways[,-1]
rownames(pathways2) <- pathways[,1]
pathways2 <- as.matrix (pathways2)
pathways <- pathways2
View(pathways)

# Now make a PCoA plot for the PICRUSt2 analysis
# Calculate distance matrix
# distance.matrix <- dist (scale (t (pathways), center = TRUE, scale = TRUE), method = "euclidian")
# distance.matrix <- distance (t (pathways), method = "bray")
distance.matrix <- vegdist (t (pathways), method = "bray")

# calculate PCoA from the distance matrix
sw.pcoa <- pcoa(distance.matrix, correction = "cailliez")

# Perform multidimensional scaling
mds.stuff <- cmdscale (distance.matrix, eig = TRUE, x.ret = TRUE)

# Round
mds.var.per <- round (mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

# Format data for ggplot
pcoa.data <- data.frame (Sample = rownames (sw.pcoa$vectors),
                         X = sw.pcoa$vectors[,1],
                         Y = sw.pcoa$vectors[,2])
mds.values <- mds.stuff$points
mds.data <- data.frame (Sample = rownames (mds.values),
                        X = mds.values[,1],
                        Y = mds.values[,2])

# Add diet column
Diet <- c(rep("LFD", 9), rep("HFD", 10), rep("HXN", 10))
pcoa.data <- cbind (pcoa.data, Diet)
mds.data <- cbind (mds.data, Diet)
mds.data$Diet <- as.factor(mds.data$Diet)
pcoa.data$Diet <- as.factor(pcoa.data$Diet)

# And plot
p <- ggplot (data = pcoa.data, aes (x=X, y=Y, color=Diet, shape = Diet)) +  # added shape = Diet
  geom_point (size = 3) +
  scale_color_manual(values=c("grey53", "green3", "black")) + #original was grey72
  theme_classic() + 
  xlab (paste ("Axis.1   [62%]")) + # number comes from mds.var.per
  ylab (paste ("Axis.2   [16.4%]")) +
  #  xlab (paste ("PCoA1 - ", mds.var.per[1], "%", sep = "")) +
  #  ylab (paste ("PCoA2 - ", mds.var.per[2], "%", sep = "")) + add this for percentages
  ggtitle ("PCoA of Metabolic Pathways (Bray-Curtis distance)") +
  theme(legend.position="bottom") +
  stat_ellipse(geom = "polygon", type="norm", level = 0.95, alpha = 0.001)
p$data$Diet <- levels (c("HFD", "HXN", "LFD"))
p <- p + theme (axis.text = element_text (size = rel (1.75))) +
  theme (axis.title = element_text (size = rel (1.75))) +
  theme (legend.text = element_text (size = rel (1.50))) +
  theme (legend.title = element_text (size = rel (1.50))) +
  theme (plot.title = element_text (size = rel (1.75)))
plot (p)
ggsave("Metabolic_PCoA_new.png")
ggsave("Metabolic_PCoA_new_color.png")

# Let's do some statistics!
# Adonis
pcoa.ad <- adonis (distance.matrix ~ sample_data(phyloObjR_Second)$Diet)
pcoa.ad

# Anosim on LFD - HFD only
phylo.lfd <- prune_samples (sample_names(phyloObjR_Second) < "121", phyloObjR_Second)
pathways.lfd <- pathways [, 1:19]
dist_lfd = vegdist (t (pathways.lfd), method = "bray")
pcoa.a.lfd <- anosim (dist_lfd, sample_data(phylo.lfd)$Diet, permutations = 999) 
pcoa.a.lfd
# R = 0.3150
# p = 0.006

# Anosim on HFD - HXN only
phylo.hfd <- prune_samples (sample_names(phyloObjR_Second) > "110", phyloObjR_Second)
pathways.hfd <- pathways [, 10:29]
dist_hfd = vegdist (t (pathways.hfd), method = "bray")
pcoa.a.hfd <- anosim (dist_hfd, sample_data(phylo.hfd)$Diet, permutations = 999) 
pcoa.a.hfd
# R = 0.2844
# p = 0.006


###############################################################################################
### Heatmaps
###############################################################################################


############################################
### Second feeding study only, by counts ###
############################################

ps.sec <- phyloObjR_Second

ps.abun <- merge_samples(ps.sec, "Diet", fun=median)

# Plot at phylum level
ps.abun.phylum <- tax_glom(ps.abun, taxrank="Phylum")
plot_heatmap(ps.abun.phylum, sample.label="Diet", 
             taxa.label="Phylum", sample.order=c("HFD", "HXN", "LFD"),
             low="#66CCFF", high="#FF0000", na.value="white")+
  labs (fill = "Abundance (counts)") +
  scale_x_discrete(labels=c("HFD", "HXN", "LFD")) +
  scale_x_discrete(name ="Diet")
ggsave('Heatmap_phylum_sec_new.tiff', antialias="none")



# Plot at family level
ps.abun.fam <- tax_glom(ps.abun, taxrank="Family")
plot_heatmap(ps.abun.fam, sample.label="Diet", 
             taxa.label="Family", sample.order=c("HFD", "HXN", "LFD"), taxa.order = "Phylum",
             low="#66CCFF", high="#FF0000", na.value="white")+
  labs (fill = "Abundance (counts)") +
  scale_x_discrete(labels=c("HFD", "HXN", "LFD")) +
  scale_x_discrete(name ="Diet")
ggsave('Heatmap_family_sec_new.tiff', antialias="none")


# Plot at genus level
ps.abun.gen <- tax_glom(ps.abun, taxrank="Genus")
taxaOrder = unique(taxa_names(ps.abun.gen))
plot_heatmap(ps.abun.gen, 
             sample.label = "Diet", 
             taxa.label = "Genus", 
             sample.order = c("HFD", "HXN", "LFD"), 
             taxa.order = taxaOrder, # original was "Family"
             low = "#66CCFF", high = "#FF0000", na.value = "white") +
  labs (fill = "Abundance (counts)") +
  scale_x_discrete(labels=c("HFD", "HXN", "LFD")) +
  scale_x_discrete(name ="Diet")
ggsave('Heatmap_genus_sec_new.tiff', antialias="none")


###############################################################################################
### Random Forest
###############################################################################################


# First we define the data set, taking our data from the respective 
# phyloseq objects

ps <- phyloObjR_Second

# Change the ASVs to plain numbers

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))

# Now add the AUC column to the phyloseq object

# For the second feeding study
add.auc <- data.frame (AUC = sw2.auc$AUC)

sample_data(ps)$AUC <- add.auc$AUC

pslog <- transform_sample_counts(ps, function(x) log(1 + x))

# order of plotting
desired.order <- c("HFD", "HXN", "LFD")

set.seed (998)

dataMatrix <- data.frame (diet = sample_data(pslog)$Diet,
                          auc = sample_data(pslog)$AUC,
                          mouse = sample_data(pslog)$Mouse,
                          study = sample_data(pslog)$Study_number,
                          otu_table(pslog))

# While keeping AUC as a continuous variable works just fine, it doesn't provide a confusion matrix.
# Let's try setting that to a categorical variable.

dataMatrix$auc <- as.numeric(dataMatrix$auc)
dataMatrix$auc
dataMatrix$auc[1:9] <- mean(dataMatrix$auc[1:9])
dataMatrix$auc[10:19] <- mean(dataMatrix$auc[10:19])
dataMatrix$auc[20:29] <- mean(dataMatrix$auc[20:29])
dataMatrix$auc <- as.factor(dataMatrix$auc)
levels (dataMatrix$auc) <- c("low", "med", "high")

# Define a 75/25% train/test split of the dataset
inTraining <- createDataPartition(sample_data(pslog)$Mouse, p = .75, list = FALSE)
training <- dataMatrix[inTraining,]
testing <- dataMatrix[-inTraining,]

# Next we can predict class labels on the test set using the predict function and 
# compare to the truth.

# The following came from the tutorial
# https://www.youtube.com/watch?v=HeTT73WxKIc

rf <- randomForest(formula = auc ~ ., data = training, mtry = 4, ntree = 2001, importance = TRUE)
# OOB estimate of  error rate: 4.17%

rf
plot(rf)

result <- data.frame(testing$auc, predict(rf, testing, type = "response"))
result

plot(result)

# This is what I actually used.

rfFit <- train(auc ~ diet, data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table.all <- table(rfClasses, testing$diet)
table.all
# Result:
# rfClasses HFD HXN LFD
#      low    0   0   2
#      med    0   1   0
#      high   2   0   0

# To interpret these PLS and random forest results, it is standard to produce 
# biplots and proximity plots, respectively. The code below extracts coordinates 
# and supplies annotation for points to include on the PLS biplot.

## Draw biplot

## Draw Random Forest biplot
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTraining, ])
a <- ggplot(rf_prox) +
  theme_classic() +
  geom_point(aes(x = X1, y = X2, col = Diet),
             size = 2, alpha = 0.6) +
  scale_color_manual(values = c("black", "grey", "blue", "green")) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(col = "Diet", x = "Axis1", y = "Axis2")
a$data$Diet <- factor(a$data$Diet, levels = desired.order)
plot(a)

ggsave(filename = "Rando_AUC_Second.tiff")

# To further understand the fitted random forest model, we identify the 
# microbe with the most influence in the random forest prediction. 
as.vector(tax_table(pslog)[which.max(importance(rfFit$finalModel)), 
                           c("Family", "Genus")])







################################################################################
### End of script.
################################################################################
