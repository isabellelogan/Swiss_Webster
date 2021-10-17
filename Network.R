setwd ("/Users/isabelle/Downloads/flannery_stagaman_analysis-master/")

require (cowplot)
require (cplm)
require (factoextra)
require (Hmisc)
library (doBy)
include=FALSE
source("packages_and_sources.R")

# Set options
options(knitr.table.format = "html")
.pardefault <- par(no.readonly = T)
theme_set(theme_cowplot())
figure_theme <- theme_update(
  plot.background = element_blank(),
  plot.title = element_text(size = 12),
  plot.subtitle = element_text(face = "italic", size = 8),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 10),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 5),
  legend.box.background = element_blank(),
  legend.key = element_blank()
)
opts_chunk$set(
  echo = FALSE,
  tidy = TRUE,
  warning = FALSE,
  message = FALSE,
  # out.width = 650,
  # out.height = 650 * 0.6,
  dpi = 300
)


###### SET VARIABLES ######
main.seed <- 1337
set.seed(main.seed)
pnas.2col <- 7
blue_to_red <- c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")
#kegg.map.files <- c(
#  "Static_data/map_kos_mod.rds",
#  "Static_data/my_kos_mods_clpsd.rds"
#)

# Read metadata
meta.data.file <- "Static_data/Sharpton_Compendium_R_CV_only.csv"

data.date <- {
  strsplit(meta.data.file, split = ".", fixed = T)[[1]][1] %>% strsplit(split = "_")
}[[1]] %>%
  tail(1) # grab date from metadata file name

# name and create directories for saving objects and plots
saveDir <- paste("Saved_objs", data.date, "data", sep = "_")
if (!dir.exists(saveDir)) {
  dir.create(saveDir)
  cat(paste(saveDir, "created"))
}

#########################
###### IMPORT DATA ######
#########################
merged.abund.stats <- read.csv("Static_data/stats.csv", header = TRUE) # abundance stats file

# Here you need to pick whichever file you want
# Genes
merged.rel.abunds <- readRDS("Static_data/picrust2_rel_abunds.rds") # KO relative abundances data
# Pathways
merged.rel.abunds <- read.csv("Static_data/picrust2_rel_pathways.csv", header = TRUE) # pathways relative abund.
#
metaphlan.tax.tbl <- readRDS("Static_data/otu_tax_tbl.rds") # taxon abundance table

# Remove weird characters
metaphlan.tax.tbl$c <- gsub('_', '-', metaphlan.tax.tbl$c)

# Delete columns that sum to zero for the taxon abundance table.

metaphlan.tax.tbl <- metaphlan.tax.tbl[,-(which(colSums(metaphlan.tax.tbl) < 10))]
metaphlan.tax.tbl <- metaphlan.tax.tbl[, order (colSums(-metaphlan.tax.tbl))]

View(metaphlan.tax.tbl)
# metaphlan.tax.tbl <- metaphlan.tax.tbl[,1:40]

# Remove columns that have a zero in it from the taxon and pathway tables
metaphlan.tax.tbl[metaphlan.tax.tbl == 0] <- NA
metaphlan.tax.tbl <- t(metaphlan.tax.tbl)
metaphlan.tax.tbl <- metaphlan.tax.tbl [complete.cases(metaphlan.tax.tbl),]
metaphlan.tax.tbl <- t(metaphlan.tax.tbl)

# Round the pathways data frame to integers
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

merged.rel.abunds <- round_df(merged.rel.abunds, 0)

# Remove pathways that have zeroes in them

merged.rel.abunds[merged.rel.abunds == 0] <- NA
merged.rel.abunds <- t(merged.rel.abunds)
merged.rel.abunds <- merged.rel.abunds [complete.cases(merged.rel.abunds),]
merged.rel.abunds <- t(merged.rel.abunds)



#################################################3
# Remove the low fat diet samples
#################################################3
merged.rel.abunds <- merged.rel.abunds[-(1:9),]
metaphlan.tax.tbl <- metaphlan.tax.tbl[-(1:9),]
merged.abund.stats <- merged.abund.stats[-(1:9),]
#################################################3
#################################################3


#################################################3
# Remove the high fat diet samples
#################################################3
merged.rel.abunds <- merged.rel.abunds[-(1:10),]
metaphlan.tax.tbl <- metaphlan.tax.tbl[-(1:10),]
merged.abund.stats <- merged.abund.stats[-(1:10),]
#################################################3
#################################################3


# Check filetype of metadata file
if (grepl(".txt", meta.data.file)) {
  data.file.sep <- "\t"
} else if (grepl(".csv", meta.data.file)) {
  data.file.sep <- ","
} else {
  stop("Meta-data file type not recognized")
}


vars.legend <- read.table(
  "Static_data/variables_legend.txt",
  sep = "\t",
  header = TRUE,
  na.strings = "N/A",
  stringsAsFactors = FALSE
) %>% as.data.table()
setkey(vars.legend, var.name)
# a data.frame with covariate categories
varIDs <- data.table(
  "group" = c(
    "intr",           "intr",       "intr", "intr",        "intr",  "intr",        "intr",  "intr"
  ),
  "varID" = c(
    "Weight_",   "Diet_",      "GTT_", "plasma_",        "LV_", "SM_",        "IL_", "CN_"
  ),
  "name" =  c(
    "Weight", "Diet Consumed", "Glucose Tolerance",  "Plasma", "Liver",  "Muscle", "Ileum", "Colon"
  ),
  stringsAsFactors = FALSE
)



# Import metadata
stool.metadata <- read.table(
  meta.data.file,
  header = TRUE,
  row.names = 1,
  sep = data.file.sep,
  na.strings = "NA",
  stringsAsFactors = FALSE
)

stool.tax.tbl <- metaphlan.tax.tbl

test <- as.data.frame(test)
rownames (test) <- c(111:130)
stool.tax.tbl <- test
# stool.tax.tbl <- readRDS ("Static_data/otu_tax_tbl.rds")
stool.tax.tbl <- stool.tax.tbl[!duplicated(stool.tax.tbl)]
# create a phyloseq object of all reads and their taxon assignments
stool.tax.phy <- phyloseq(
  sample_data(stool.metadata),
  otu_table(stool.tax.tbl, taxa_are_rows = FALSE)
)

# grab only stool samples
stool.kos.rel.abunds <- merged.rel.abunds
# make sample numbers rownames
stool.kos.rel.abunds2 <- stool.kos.rel.abunds[,-1]
rownames(stool.kos.rel.abunds2) <- stool.kos.rel.abunds[,1]
stool.kos.rel.abunds <- stool.kos.rel.abunds2

#(stool.kos.phy <- phyloseq(
#  sample_data(stool.metadata),
#  otu_table(as.matrix(stool.kos.rel.abunds), taxa_are_rows = FALSE)
#))

#stool.tax.tbl <- readRDS ("Static_data/otu_tax_tbl.rds")
#stool.tax.tbl <- stool.tax.tbl[!duplicated(stool.tax.tbl)]






#################################################3
#################################################3
# Remove low fat diet
#################################################3
stool.tax.tbl <- stool.tax.tbl[-(1:9),]
#################################################3
#################################################3





#################################################3
#################################################3
# Remove high fat diet
#################################################3
stool.tax.tbl <- stool.tax.tbl[-(1:10),]
#################################################3
#################################################3



stool.tax.phy <- phyloseq(
  sample_data(stool.metadata),
  otu_table(stool.tax.tbl, taxa_are_rows = FALSE)
)

taxa.names <- taxa_names(stool.tax.phy)
tax.tbl <- do.call(
  "rbind",
  lapply(
    taxa.names,
    function(taxon) {
      tax.split <- strsplit(taxon, "|", fixed = T)[[1]]
      if (length(tax.split) < 8) {
        tax.split <- c(tax.split, rep(NA, (8 - length(tax.split))))
      }
      tax.row <- sub("\\w__", "", tax.split, perl = T)
    }
  )
)
colnames(tax.tbl) <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"
)
rownames(tax.tbl) <- taxa.names
View(tax.tbl)

# Remove weird characters
tax.tbl[,6] <- gsub('/Shigella', '', tax.tbl[,6])

tax_table(stool.tax.phy) <- tax_table(tax.tbl)



saveRDS(stool.tax.phy, file = file.path(saveDir, "phyloseq_stool_metagenome_taxa.rds"))

######################################################
## Compound Poisson Generalized Linear Models
######################################################

# These chunks are best run on a machine with many cores and a good amount of storage.
# I usually run these overnight

{set.seed(main.seed)
  terms.file <- file.path(saveDir, "list_analysis_terms_envfit.rds")

  # Change the following as needed, once for ASV data, once for PICRUSt2 data

  #type <- "tax" # ASVs
   type <- "kos" # PICRUSt2

  dfsDir <- paste("CondsInModel", type, "DFs", sep = "_")
  if (!dir.exists(file.path(saveDir, dfsDir))) {
    dir.create(file.path(saveDir, dfsDir))
  }

  phy.name <- paste("stool", type, "phy", sep = ".")

  phy.obj <- eval(parse(text = phy.name))
  # phy.obj <- metaphlan.tax.tbl
  vars <- c("GTT_AUC", "GTT_HOMA_IR", "plasma_insulin", "Weight_gain", "Weight_Gain_Daily",
            "plasma_leptin", #"plasma_HDL", "Diet_Total_eaten",
            "Weight_Food_Efficiency", "plasma_LDL", "plasma_TAG", "LV_TAG", "LV_CHO")
  analysis.terms <- vars
  #  vars <- c("GTT_AUC") # Used for testing
  conds <- c("Diet")
  pruned.phy <- prune_samples(
    row.names(sample_data(phy.obj)[complete.cases(sample_data(phy.obj)[, vars]), ]),
    phy.obj
  )
  pruned.phy <- prune_taxa(taxa_sums(pruned.phy) > 0, pruned.phy)
  # pruned.phy <- phy.obj # for all asvs

  # metadata <- as(sample_data(pruned.phy)[, c(vars, conds)], "data.frame")
  # Make the metadata the compendium data!
  metadata <- stool.metadata



  ##################################################################
  # For HFD - HXN only
   metadata <- metadata[-(1:9),]
  ##################################################################


  ##################################################################
  # For HXN only, run after removing LFD samples
  metadata <- metadata[-(1:10),]
  ##################################################################



  #metadata$Diet <- c(rep("LFD", 9), rep("HFD", 10), rep("HXN", 10))
  comm.tbl <- as(otu_table(pruned.phy), "matrix")

  # Remove columns not occurring in all samples
  comm.tbl[comm.tbl == 0] <- NA
  comm.tbl <- t(comm.tbl)
  comm.tbl <- comm.tbl [complete.cases(comm.tbl),]
  comm.tbl <- t(comm.tbl)

  # Get Akkermansia only (second column) (not used)
  #test <- comm.tbl[,1:2]
  #comm.tbl <- comm.tbl[,1:49]#test

  if (type == "tax") {
    colnames(comm.tbl) <- gsub("|", ".", colnames(comm.tbl), fixed = T)
  }

  all.units <- colnames(comm.tbl)

  all.pvals <- list (length(vars))
  section.pvals <- c(1:length(all.units))
  for (j in 1:length(vars)) {

    for (i in 1:length(all.units)) {
      # for (i in 1:50) { # for testing

      section.pvals [i] <- sapply(
        vars [j],
        function(var) {
          # important!!
          # change var to vars when testing
          ##################################
          model.data <- cbind(metadata[, c(var, conds),
                                       drop = F], comm.tbl[, all.units[i], drop = F])

          # Use this line with diet
          # unit.frm <- paste(all.units[i], " ~ ", paste(var), " + ", paste (conds))

          # This line doesn't work
          # unit.frm <- paste(all.units[i], "~", paste(c(vars[j], conds)), collapse = " + ")

          # Using the dot works
          # unit.frm <- paste(all.units[i], "~ .")

          # Use this line without diet
          unit.frm <- paste(all.units[i], "~", paste(c(var)))

          unit.model <- try( #cpglm was original, can also use glm.nb
            cpglm(

              as.formula(unit.frm),
              data = model.data,
              link = "log", # can also use identity, log works
              #optimizer = "L-BFGS-B" # can also use L-BFGS-B
              # bobyqa works here
              optimizer = "bobyqa"
            ),
            silent = F
          )

          if (class(unit.model) == "try-error") {
            unit.pval <- NA
            intercept <- NA
            slope <- NA
          } else {
            sink(file = "tmp.txt")
            model.sum <- summary(unit.model)
            sink()
            unit.pval <- model.sum$coefficients[var, "Pr(>|t|)"] # change to "Pr(>|z|)" for glm.nb,
            # original is "Pr(>|t|)"
            intercept <- unname(unit.model$coefficients[1])
            slope <- unname(unit.model$coefficients[var])
          }
          ### ^ There is no r2 value for this method (AIC-based), so going to use slope instead
          unit.df <- data.frame(
            "unit" = unit.pval, # this was unit originally
            "b" = intercept,
            "m" = slope
          )
          saveRDS(
            unit.df,
            file.path(saveDir, dfsDir, paste("df_", var, "-", all.units[i], ".rds", sep = ""))
          )
          return(unit.pval)
        }
      )
    }
    all.pvals[[j]] <- section.pvals
  }

}
# For ASVs save it using the following code
saveRDS(
  all.pvals,
  file = file.path(
    saveDir,
    paste("list_", type, "_childSesVars_allUnits_cpglms_condsInModel_unadjPvals2.rds",
          sep = ""
    )
  )
)
# For Pathways or Genes save it using the following code
saveRDS(
  all.pvals,
  file = file.path(
    saveDir,
    paste("list_", type, "_cpglms_pathways.rds",
          sep = ""
    )
  )
)

# Added this to get the names of the bacteria
taxon <- colnames (comm.tbl)

# <!-- ### ASV-Covariate Relationships -->


dfsDir <- "CondsInModel_tax_DFs"
dfsTar <- paste(dfsDir, "tgz", sep = ".")
varPlotDir <- file.path(saveDir, "cpglm_condsInModel_tax_plots")
if (!dir.exists(varPlotDir)) {
  dir.create(varPlotDir)
}


tax.unadj.pvals <- readRDS(
  file = file.path(saveDir, "list_tax_childSesVars_allUnits_cpglms_condsInModel_unadjPvals2.rds")
)

tax.unadj.pvals.df <- data.frame(tax.unadj.pvals, stringsAsFactors = F)
# The column names need to match those from line 1532 (where we run the regression model)
colnames (tax.unadj.pvals.df) <- c("GTT_AUC", "GTT_HOMA_IR", "plasma_insulin", "Weight_gain", "Weight_Gain_Daily",
                                   "plasma_leptin", #"plasma_HDL", "Diet_Total_eaten",
                                   "Weight_Food_Efficiency", "plasma_LDL", "plasma_TAG", "LV_TAG", "LV_CHO")
tax.going.in <- comm.tbl # Need to get this from section 5
tax.unadj.pvals.df$unit <- colnames (tax.going.in)
View(tax.unadj.pvals.df)

for (i in 1:(ncol(tax.unadj.pvals.df)-1)) {
  tax.unadj.pvals.df[,i] <- p.adjust(tax.unadj.pvals.df[,i], method = "fdr")
}
# For all diets combined
#tax.adj.pvals <- data.frame(
#  "var" = c(rep ("GTT_AUC", 14), rep ("GTT_HOMA_IR", 14), rep("plasma_insulin", 14), rep("Weight_gain", 14),
#            rep ("Weight_Gain_Daily", 14),rep ("plasma_leptin", 14),
#            rep ("Weight_Food_Efficiency", 14), rep ("plasma_LDL", 14),rep ("plasma_TAG", 14),
#            rep ("LV_TAG", 14), rep ("LV_CHO", 14)),
tax.adj.pvals <- data.frame(
  "var" = c(rep ("GTT_AUC", nrow(tax.unadj.pvals.df)), rep ("GTT_HOMA_IR", nrow(tax.unadj.pvals.df)), rep("plasma_insulin", nrow(tax.unadj.pvals.df)), rep("Weight_gain", nrow(tax.unadj.pvals.df)),
            rep ("Weight_Gain_Daily", nrow(tax.unadj.pvals.df)),rep ("plasma_leptin", nrow(tax.unadj.pvals.df)),
            #rep ("plasma_HDL", nrow(tax.unadj.pvals.df)), rep ("Diet_Total_eaten" ,nrow(tax.unadj.pvals.df)),
            rep ("Weight_Food_Efficiency", nrow(tax.unadj.pvals.df)), rep ("plasma_LDL", nrow(tax.unadj.pvals.df)),rep ("plasma_TAG", nrow(tax.unadj.pvals.df)),
            rep ("LV_TAG", nrow(tax.unadj.pvals.df)), rep ("LV_CHO", nrow(tax.unadj.pvals.df))),
  "unit" = c(rep (tax.unadj.pvals.df$unit, (ncol (tax.unadj.pvals.df) - 1 ))),
  "adj.pval" = c(tax.unadj.pvals.df$GTT_AUC, tax.unadj.pvals.df$GTT_HOMA_IR, tax.unadj.pvals.df$plasma_insulin,
                 tax.unadj.pvals.df$Weight_gain, tax.unadj.pvals.df$Weight_Gain_Daily,
                 tax.unadj.pvals.df$plasma_leptin, #tax.unadj.pvals.df$plasma_HDL, tax.unadj.pvals.df$Diet_Total_eaten,
                 tax.unadj.pvals.df$Weight_Food_Efficiency, tax.unadj.pvals.df$plasma_LDL, tax.unadj.pvals.df$plasma_TAG,
                 tax.unadj.pvals.df$LV_TAG, tax.unadj.pvals.df$LV_CHO),
  #"adj.pval" = c(tax.unadj.pvals.df$GTT_AUC, tax.unadj.pvals.df$GTT_HOMA_IR, tax.unadj.pvals.df$plasma_insulin),
  stringsAsFactors = F
)





###################################################
# For HFD HXN only # might not need this after fixing the code above
tax.adj.pvals <- data.frame(
  "var" = c(rep ("GTT_AUC", 17), rep ("GTT_HOMA_IR", 17), rep("plasma_insulin", 17), rep("Weight_gain", 17),
            rep ("Weight_Gain_Daily", 17),rep ("plasma_leptin", 17),
            rep ("Weight_Food_Efficiency", 17), rep ("plasma_LDL", 17),rep ("plasma_TAG", 17),
            rep ("LV_TAG", 17), rep ("LV_CHO", 17)),
  "unit" = c(rep (tax.unadj.pvals.df$unit, 11)),
  "adj.pval" = c(tax.unadj.pvals.df$GTT_AUC, tax.unadj.pvals.df$GTT_HOMA_IR, tax.unadj.pvals.df$plasma_insulin,
                 tax.unadj.pvals.df$Weight_gain, tax.unadj.pvals.df$Weight_Gain_Daily, tax.unadj.pvals.df$plasma_leptin,
                 tax.unadj.pvals.df$plasma_HDL, tax.unadj.pvals.df$Diet_Total_eaten,
                 tax.unadj.pvals.df$Weight_Food_Efficiency, tax.unadj.pvals.df$plasma_LDL, tax.unadj.pvals.df$plasma_TAG,
                 tax.unadj.pvals.df$LV_TAG, tax.unadj.pvals.df$LV_CHO),
  #"adj.pval" = c(tax.unadj.pvals.df$GTT_AUC, tax.unadj.pvals.df$GTT_HOMA_IR, tax.unadj.pvals.df$plasma_insulin),
  stringsAsFactors = F
)
###################################################





tax.sig.pvals <- droplevels(subset(tax.adj.pvals, adj.pval < 0.05))

if (nrow(tax.sig.pvals) == 0) {
  print("No Taxa have a significant relationship with any co-variates")
} else {
  # metadata <- readRDS(file.path(saveDir, "df_allVars_condsResids.rds"))
  metadata <- as(sample_data(stool.tax.phy), "data.frame")
  comm.tbl <- as(otu_table(stool.tax.phy), "matrix")
  comm.tbl <- comm.tbl[row.names(metadata), ]
  colnames(comm.tbl) <- gsub("|", ".", colnames(comm.tbl), fixed = TRUE)
  # Remove columns not occurring in all samples
  comm.tbl[comm.tbl == 0] <- NA
  comm.tbl <- t(comm.tbl)
  comm.tbl <- comm.tbl [complete.cases(comm.tbl),]
  comm.tbl <- t(comm.tbl)
#  comm.tbl <- comm.tbl [,1:49]
  tax.sig.pvals$spearman.r <- apply(
    tax.sig.pvals[, 1:2], 1,
    function(row) {
      print(row)
      cor(as.numeric(metadata[, as.character(row["var"])]),
          comm.tbl[, as.character(row["unit"])],
          method = "spearman",
          use = "complete.obs"
      )
    }
  )
  tax.plot.list <- list()
  unit.label.list <- list()
  # v <- "child_behav_CBCL_extern_t1"
  for (v in sort(unique(tax.sig.pvals$var))) {
    print(v)
    var.sig.pvals <- subset(tax.sig.pvals, var == v) %>% droplevels()
    if (nrow(var.sig.pvals) == 0) {
      NULL
    } else {
      var.data <- cbind(
        metadata[, v, drop = F],
        comm.tbl[, as.character(var.sig.pvals$unit), drop = F]
      ) %>%
        reshape2::melt(id.var = v, variable.name = "unit", value.name = "rel.abund")
      var.data$unit.short <- strsplit(as.character(var.data$unit), ".", fixed = T) %>%
        sapply(tail, 2) %>%
        apply(2, paste, collapse = "\n")
      unit.med.relAbund <- doBy::summaryBy(rel.abund ~ unit, data = var.data, FUN = median)
      good.units <- subset(unit.med.relAbund, rel.abund.median > 0)
      if (nrow(good.units) == 0) {
        NULL
      } else {
        var.data <- subset(var.data, unit %in% good.units$unit) %>% droplevels()
        cor.data <- subset(var.sig.pvals, var == v & unit %in% unique(var.data$unit)) %>%
          droplevels()
        unit.label.list[[v]] <- data.frame(
          "var" = v,
          "relationship" = ifelse(cor.data$spearman.r < 0, "-",
                                  ifelse(cor.data$spearman.r > 0, "+", "flat")
          ),
          "taxon" = unique(var.data$unit)
        )
        reg.lines <- data.frame(NULL)
        # u <- unique(as.character(var.data$unit))[1]
        reg.slopes <- data.frame(
          matrix(ncol = 2, nrow = length(unique(as.character(var.data$unit))))
        ) %>% setNames(c("unit", "reg.slope"))
        row.names(reg.slopes) <- unique(as.character(var.data$unit))
        for (u in unique(as.character(var.data$unit))) {
          print(u)
          model.data <- subset(var.data, unit == u)
          file.vu <- paste(dfsDir, "/df_", v, "-", u, ".rds", sep = "")
          setwd("/Users/isabelle/Downloads/flannery_stagaman_analysis-master/Saved_objs_only_data")
          #untar(file.path(saveDir, dfsTar),
          #      files = file.vu,
          #      compressed = "gzip"
          #)
          unit.df <- readRDS(file.vu)
          #file.remove(file.vu)
          #file.remove(dfsDir)
          unit.df$unit.short <- strsplit(as.character(u), ".", fixed = T) %>%
            sapply(tail, 2) %>%
            apply(2, paste, collapse = "\n")
          reg.slopes[u, c("unit", "reg.slope")] <- unit.df[1, c("unit.short", "m")]
          reg.lines <- rbind(reg.lines, unit.df)
        }
        unit.label.list[[v]]$reg.slope <- reg.slopes$reg.slope
        var.data$log.rel.abund <- log(var.data$rel.abund)
        var.plot <- ggplot(data = var.data, aes_string(x = v, y = "log.rel.abund")) +
          geom_abline(
            data = reg.lines, aes(intercept = b, slope = m),
            color = "blue", size = 0.5
          ) +
          geom_point(aes(shape = ifelse(rel.abund == 0, "zero", "not zero")),
                     alpha = 0.6
          ) +
          scale_shape_manual(name = "Relative abundance", values = c(16, 1)) +
          labs(
            y = "log(Unit Relative Abundance)",
            x = mgsub(vars.legend$var.name, vars.legend$description, v),
            title = "Compound Poisson GLMs on variable residuals"
          ) +
          facet_wrap(~unit.short, scales = "free_y") +
          theme(
            axis.text.y = element_blank(),
            legend.position = "top",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 8)
          )
        tax.plot.list[[v]] <- var.plot
        nCol <- ceiling(sqrt(length(unique(var.data$unit))))
        nRow <- ceiling(length(unique(var.data$unit)) / nCol)
        setwd("/Users/isabelle/Downloads/flannery_stagaman_analysis-master/")
        ggsave(
          var.plot,
          file = file.path(varPlotDir, paste(gsub(".", "_", v, fixed = T), "png", sep = ".")),
          height = nRow * 2,
          width = nCol * 2.5,
          dpi = 300
        )
      }
    }
  }

  # unit.label.df <- do.call("rbind", e$unit.label.list)
  unit.label.df <- do.call("rbind", unit.label.list)
  row.names(unit.label.df) <- NULL
  unit.label.df$var <- mgsub(
    c(vars.legend$var.name, ".resid"),
    c(vars.legend$description, " Residuals"),
    unit.label.df$var
  )
  taxa.split <- strsplit(as.character(unit.label.df$taxon), ".", fixed = TRUE)
  taxa.tbl <- sapply(taxa.split, function(i) {
    if (length(i) < 8) {
      i <- c(i, "t_NA")
    }
    i
  }) %>%
    t() %>%
    as.data.frame()
  taxa.tbl <- taxa.tbl[, -8]
  names(taxa.tbl) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  names(unit.label.df) <- c("Variable", "Relationship", "lab", "Slope")
  unit.label.df <- cbind(unit.label.df[, c(1, 2, 4)], taxa.tbl)
  unit.label.df[, -c(1:3)] <- apply(
    unit.label.df[, -c(1:3)], 2,
    function(col) sapply(strsplit(col, "__"), `[`, 2)
  )
  # write.table(
  #   unit.label.df,
  #   file = file.path(saveDir, "df_sigVarTaxRels_condsInModel.txt"),
  #   sep = "\t", quote = F, row.names = F
  # )
  saveRDS(unit.label.df, file = file.path(saveDir, "df_sigVarTaxRels_condsInModel.rds"))
  unit.label.dt <- as.data.table(unit.label.df)
  saveRDS(unit.label.dt, file = file.path(saveDir, "dt_sigVarTAXrels_condsInModel.rds"))
  tax.plot.list <- tax.plot.list[!sapply(tax.plot.list, is.null)]
  nRows <- sapply(tax.plot.list, function(plot) {
    panels <- max(as.numeric(ggplot_build(plot)$data[[1]]$PANEL))
    if (panels > 3) {
      nCol <- ceiling(sqrt(panels))
      ceiling(panels / nCol)
    } else {
      1
    }
  })
  plot.heights <- nRows / sum(nRows)
  tax.plot.height <- ifelse(sum(nRows) / 4 < length(nRows),
                            length(nRows) * 6,
                            (sum(nRows) / 4) * 8
  )
}


tax.network.cap <- ""

tax.vars.rels <- readRDS(file.path(saveDir, "dt_sigVarTAXrels_condsInModel.rds"))
plot.data0 <- tax.vars.rels[
  ,
  .(taxComps = length(Genus), Family, Genus, Slope),
  #  .(taxComps = length(Species), Species, Slope),
  by = Variable
]
plot.data <- plot.data0[
  ,
  .(varComps = length(Variable), taxComps, Family, Variable, Slope),
  #  by = Species
  by = Genus
]
# Combine the Family and Genus columns into 1
plot.data$Genus <- paste (plot.data$Family, plot.data$Genus)
plot.data$Type <- c(rep("Bug", nrow(plot.data)))
# plot.data$Type <- ifelse(grepl("GTT_AUC", plot.data$Variable), "HOMA_IR", "Insulin")
# plot.data$Species <- gsub("_", " ", plot.data$Species)
ordered.taxa <- levels(with(plot.data, reorder(Genus, varComps)))
# ordered.taxa <- levels(with(plot.data, reorder(Species, varComps)))
ordered.chb <- levels(with(plot.data[Type != "Bug"], reorder(Variable, taxComps)))
ordered.esa <- levels(with(plot.data[Type == "Bug"], reorder(Variable, taxComps)))
range <- c(-5, 5)
l <- mean(range) - 1.5
r <- mean(range) + 1.5
seg.data <- within(plot.data, {
  x0 <- ifelse(Type != "Bug", range[1], range[2])
  y0 <- ifelse(Type != "Bug",
               (match(Variable, ordered.chb) - 0.5) / length(ordered.chb),
               (match(Variable, ordered.esa) - 0.5) / length(ordered.esa)
  )
  x1 <- ifelse(Type != "Bug", l, r)
  y1 <- (match(Genus, ordered.taxa) - 0.5) / length(ordered.taxa)
  #  y1 <- (match(Species, ordered.taxa) - 0.5) / length(ordered.taxa)
  mid <- mean(range)
})
# Switch the variables around

seg.data <- seg.data[order(abs(Slope))]
lab.font.pt <- 6
lab.size <- lab.font.pt * 0.35
taxa.line <- ggplot(data = seg.data) +
  geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1, color = Slope, size = abs(Slope), alpha = abs(Slope))) +
  scale_color_gradientn(name = "Reg. Slope", colors = blue_to_red) +
  scale_size_continuous(name = "abs(Reg. Slope)", range = c(0.05, 1.5)) +
  scale_alpha_continuous(name = "abs(Reg. Slope)", range = c(0.6, 0.8)) +
  geom_tile(aes(x = mid, y = y1), width = r - l, fill = "white", color = "black") +
  geom_text(
    #   data = unique(seg.data[, .(Family, mid, y1)]), aes(x = mid, y = y1, label = Family),
    data = unique(seg.data[, .(Genus, mid, y1)]), aes(x = mid, y = y1, label = Genus),
    #   data = unique(seg.data[, .(Species, mid, y1)]), aes(x = mid, y = y1, label = Species),
    size = lab.size, fontface = "italic"
  ) +
  geom_tile(
    data = seg.data[Type != "Bug"], aes(x = x0 + l, y = y0),
    width = r - l, height = 0.05, fill = "white", color = "black"
  ) +
  geom_text(
    data = unique(seg.data[Type != "Bug", .(Variable, x0, y0)]),
    aes(x = x0 + l, y = y0, label = Variable), size = lab.size
  ) +
  geom_tile(
    data = seg.data[Type == "Bug"], aes(x = x0 + r, y = y0),
    width = r - l, height = 0.05, fill = "white", color = "black"
  ) +
  geom_text(
    data = unique(seg.data[Type == "Bug", .(Variable, x0, y0)]),
    aes(x = x0 + r, y = y0, label = Variable), size = lab.size
  ) +
  scale_x_continuous(
    breaks = c(min(seg.data$x0) + l, 0, max(seg.data$x0) + r),
    labels = c("Link", "Taxa", "Phenotype"),
    position = "top"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(11, "pt"),
    legend.box.margin = unit(c(0, 0, 0, 0), "pt"),
    legend.box.spacing = unit(5, "pt")
  )
taxa.line

######################################
# save the data for the next section
tax.section <- seg.data
######################################

fig.wd <- pnas.2col / 0.8887
figDir = "/Users/isabelle/Downloads/flannery_stagaman_analysis-master/Figure_plots"
ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork.pdf"), width = fig.wd, height = 6)

ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork_HFDHXN.pdf"), width = fig.wd, height = 6)# system("open -a /Applications/Inkscape.app/Figure_plots/mainFig_chb_esa_tax_regNetwork.pdf")

ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork_Gene.pdf"), width = fig.wd, height = 6)# system("open -a /Applications/Inkscape.app/Figure_plots/mainFig_chb_esa_tax_regNetwork.pdf")

ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork_HFDHXN_Gene.pdf"), width = fig.wd, height = 6)# system("open -a /Applications/Inkscape.app/Figure_plots/mainFig_chb_esa_tax_regNetwork.pdf")

ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork_20210212.pdf"), width = fig.wd, height = 6)

ggsave(taxa.line, file = file.path(figDir, "mainFig_chb_esa_tax_regNetwork_HXN_Only.pdf"), width = fig.wd, height = 6)
# Find the most important bacteria per phenotype

path.list.bugs <- as_tibble (seg.data)
path.list.bugs.sorted <- path.list.bugs %>%
  arrange (Variable, desc(Slope))
View(path.list.bugs.sorted)

# Save these pathways
writexl::write_xlsx(path.list.bugs.sorted, "Figure7_Bacteria.xlsx")

# Save for HFD - HXN
writexl::write_xlsx(path.list.bugs.sorted, "HFD_HXN_Bacteria.xlsx")

# Save for HXN only
writexl::write_xlsx(path.list.bugs.sorted, "HXN_Bacteria.xlsx")





### Pathway-Covariate Relationships -->


dfsDir <- "CondsInModel_kos_DFs"
dfsTar <- paste(dfsDir, "tgz", sep = ".")
varPlotDir <- file.path(saveDir, "cpglm_condsInModel_kos_plots")
if (!dir.exists(varPlotDir)) {
  dir.create(varPlotDir)
}
unitPlotDir <- file.path(varPlotDir, "unit_plots")
if (!dir.exists(unitPlotDir)) {
  dir.create(unitPlotDir)
}
kos.unadj.pvals <- readRDS(
  file = file.path(
    saveDir, "list_kos_cpglms_pathways.rds"
  )
)
kos.unadj.pvals.df <- data.frame(kos.unadj.pvals, stringsAsFactors = F)
# colnames (kos.unadj.pvals.df) <- c("GTT_AUC", "GTT_HOMA_IR", "plasma_insulin")
kos.going.in <- comm.tbl # Need to run section 5 to get this info
kos.unadj.pvals.df$unit <- colnames (kos.going.in) # Need to run section 5 to get this info
View(kos.unadj.pvals.df)
colnames (kos.unadj.pvals.df) <- c("GTT_AUC", "GTT_HOMA_IR", "plasma_insulin", "Weight_gain", "Weight_Gain_Daily",
                                   "plasma_leptin",# "plasma_HDL", "Diet_Total_eaten",
                                   "Weight_Food_Efficiency", "plasma_LDL", "plasma_TAG", "LV_TAG", "LV_CHO", "unit")
#     %>%
#  as.matrix() %>%
#  melt(as.is = T)
# Adjust p-values
for (i in 1:(ncol(kos.unadj.pvals.df)-1)) {
  kos.unadj.pvals.df[,i] <- p.adjust(kos.unadj.pvals.df[,i], method = "fdr")
}

# For pathways
kos.adj.pvals <- data.frame(
  #  "var" = c(rep ("GTT_AUC", 237), rep ("GTT_HOMA_IR", 237), rep("plasma_insulin", 237),
  #            rep ("Weight_gain", 237), rep ("Weight_Gain_Daily", 237), rep("plasma_leptin", 237),
  #            rep ("Weight_Food_Efficiency", 237), rep ("plasma_LDL", 237), rep("plasma_TAG", 237),
  #            rep ("LV_TAG", 237), rep ("LV_CHO", 237)),
  #  "unit" = c(rep (kos.unadj.pvals.df$unit, 11)),
  #  "adj.pval" = c(kos.unadj.pvals.df$GTT_AUC, kos.unadj.pvals.df$GTT_HOMA_IR, kos.unadj.pvals.df$plasma_insulin,
  #                 kos.unadj.pvals.df$Weight_gain, kos.unadj.pvals.df$Weight_Gain_Daily, kos.unadj.pvals.df$plasma_leptin,
  #                 kos.unadj.pvals.df$Weight_Food_Efficiency, kos.unadj.pvals.df$plasma_LDL, kos.unadj.pvals.df$plasma_TAG,
  #                 kos.unadj.pvals.df$LV_TAG, kos.unadj.pvals.df$LV_CHO),
  "var" = c(rep ("GTT_AUC", nrow(kos.unadj.pvals.df)), rep ("GTT_HOMA_IR", nrow(kos.unadj.pvals.df)), rep("plasma_insulin", nrow(kos.unadj.pvals.df)), rep("Weight_gain", nrow(kos.unadj.pvals.df)),
            rep ("Weight_Gain_Daily", nrow(kos.unadj.pvals.df)),rep ("plasma_leptin", nrow(kos.unadj.pvals.df)),
            # rep ("plasma_HDL", nrow(kos.unadj.pvals.df)), rep ("Diet_Total_eaten" ,nrow(kos.unadj.pvals.df)),
            rep ("Weight_Food_Efficiency", nrow(kos.unadj.pvals.df)), rep ("plasma_LDL", nrow(kos.unadj.pvals.df)),rep ("plasma_TAG", nrow(kos.unadj.pvals.df)),
            rep ("LV_TAG", nrow(kos.unadj.pvals.df)), rep ("LV_CHO", nrow(kos.unadj.pvals.df))),
  "unit" = c(rep (kos.unadj.pvals.df$unit, (ncol (kos.unadj.pvals.df) - 1 ))),
  "adj.pval" = c(kos.unadj.pvals.df$GTT_AUC, kos.unadj.pvals.df$GTT_HOMA_IR, kos.unadj.pvals.df$plasma_insulin,
                 kos.unadj.pvals.df$Weight_gain, kos.unadj.pvals.df$Weight_Gain_Daily,
                 kos.unadj.pvals.df$plasma_leptin,
                 # kos.unadj.pvals.df$plasma_HDL, kos.unadj.pvals.df$Diet_Total_eaten,
                 kos.unadj.pvals.df$Weight_Food_Efficiency, kos.unadj.pvals.df$plasma_LDL, kos.unadj.pvals.df$plasma_TAG,
                 kos.unadj.pvals.df$LV_TAG, kos.unadj.pvals.df$LV_CHO),
  stringsAsFactors = F
)

#####################################################################################
# For genes
#####################################################################################
kos.adj.pvals <- data.frame(
  "var" = c(rep ("GTT_AUC", 2667), rep ("GTT_HOMA_IR", 2667), rep("plasma_insulin", 2667),
            rep ("Weight_gain", 2667), rep ("Weight_Gain_Daily", 2667), rep("plasma_leptin", 2667),
            rep ("Weight_Food_Efficiency", 2667), rep ("plasma_LDL", 2667), rep("plasma_TAG", 2667),
            rep ("LV_TAG", 2667), rep ("LV_CHO", 2667)),
  "unit" = c(rep (kos.unadj.pvals.df$unit, 11)),
  "adj.pval" = c(kos.unadj.pvals.df$GTT_AUC, kos.unadj.pvals.df$GTT_HOMA_IR, kos.unadj.pvals.df$plasma_insulin,
                 kos.unadj.pvals.df$Weight_gain, kos.unadj.pvals.df$Weight_Gain_Daily, kos.unadj.pvals.df$plasma_leptin,
                 kos.unadj.pvals.df$Weight_Food_Efficiency, kos.unadj.pvals.df$plasma_LDL, kos.unadj.pvals.df$plasma_TAG,
                 kos.unadj.pvals.df$LV_TAG, kos.unadj.pvals.df$LV_CHO),
  stringsAsFactors = F
)
#####################################################################################


kos.sig.pvals <- droplevels(subset(kos.adj.pvals, adj.pval < 0.05))
if (nrow(kos.sig.pvals) == 0) {
  print("No KOs have a significant relationship with any co-variates")
} else {
  metadata <- as(sample_data(stool.kos.phy), "data.frame")
  comm.tbl <- as(otu_table(stool.kos.phy), "matrix")
  comm.tbl <- comm.tbl[row.names(metadata), ]
  kos.sig.pvals$spearman.r <- apply(
    kos.sig.pvals[, 1:2], 1,
    function(row) {
      cor(as.numeric(metadata[, row["var"]]),
          comm.tbl[, row["unit"]],
          method = "spearman",
          use = "complete.obs"
      )
    }
  )
  library (readxl)
  # For pathways
  intermediate <- read_excel ("Pathway_Description.xlsx")
  saveRDS (intermediate, file = "Pathway_Description.rds")
  my.kos.mods.clpsd <- readRDS("Pathway_Description.rds")




  #########################################################################
  # For genes
  #intermediate <- read_excel ("Gene_Descriptions.xlsx")
  #saveRDS (intermediate, file = "Gene_Description.rds")
  #my.kos.mods.clpsd <- readRDS("Gene_Description.rds")
  #########################################################################





  #  my.kos.mods.clpsd <- readRDS(kegg.map.files[2])
  colnames(my.kos.mods.clpsd) <- c("modules", "modules.names")
  # my.kos.mods.clpsd$modules[my.kos.mods.clpsd$modules == "NA"] <- NA
  # my.kos.mods.clpsd$modules.names[my.kos.mods.clpsd$modules.names == "NA"] <- NA
  # write.table(my.kos.mods.clpsd, file="~/Desktop/kos_mods_tbl.txt", sep="\t", quote=F, row.names=F)
  my.kos.mods.clpsd$modules.names[!is.na(my.kos.mods.clpsd$modules.names)] <- paste(
    "MODS:", my.kos.mods.clpsd$modules.names[!is.na(my.kos.mods.clpsd$modules.names)]
  )
  #  my.kos.mods.clpsd$kos.name[!is.na(my.kos.mods.clpsd$kos.name)] <- paste(
  #    "KO:", my.kos.mods.clpsd$kos.name[!is.na(my.kos.mods.clpsd$kos.name)]
  #  )
  colnames (my.kos.mods.clpsd) <- c("unit", "modules.names")
  kos.sig.pvals.modIDs <- merge(kos.sig.pvals,
                                my.kos.mods.clpsd,
                                by.x = "unit",
                                #by.y = "modules",
                                #by.y = "kos",
                                all.x = TRUE
  )
  #  kos.pfam.clpsd <- readRDS(file = "Static_data/kos_pfam_top5matches_clpsd.rds")
  #  kos.pfam.clpsd$pfams[!is.na(kos.pfam.clpsd$pfams)] <- paste(
  #    "PFAMS:", kos.pfam.clpsd$pfams[!is.na(kos.pfam.clpsd$pfams)]
  #  )
  #  kos.sig.pvals.modIDs.pfams <- merge(kos.sig.pvals.modIDs,
  #                                      kos.pfam.clpsd,
  #                                      by.x = "unit",
  #                                      by.y = "kos",
  #                                      all.x = TRUE
  #  )
  kos.sig.pvals.modIDs.pfams <- kos.sig.pvals.modIDs
  kos.plot.list <- list()
  unit.label.list <- list()
  # v <- "child_behav_CBCL_depression_t1"
  for (v in sort(unique(kos.sig.pvals$var))) {
    print(v)
    var.sig.pvals <- subset(kos.sig.pvals.modIDs.pfams, var == v) %>% droplevels()
    if (nrow(var.sig.pvals) == 0) {
      NULL
    } else {
      var.data <- cbind(
        metadata[, v, drop = F],
        comm.tbl[, as.character(var.sig.pvals$unit), drop = F]
      ) %>%
        melt(id.var = v, variable.name = "unit", value.name = "rel.abund")
      plot.data <- merge(var.data, var.sig.pvals, by = "unit", all.x = T)
      plot.data$log.rel.abund <- log(plot.data$rel.abund)
      unit.med.relAbund <- doBy::summaryBy(rel.abund ~ unit,
                                           data = plot.data,
                                           FUN = median
      )
      good.units <- subset(unit.med.relAbund, rel.abund.median > 0)
      if (nrow(good.units) == 0) {
        NULL
      } else {
        plot.data <- subset(plot.data, unit %in% good.units$unit)
        label.base <- plot.data[, c(
          "unit",
          "spearman.r",
          #"kos.name",
          "modules.names" #,
          #"pfams"
        )] %>% unique()
        labs <- NULL
        for (row in c(1:nrow(label.base))) {
          row.lab <- ifelse(
            !is.na(label.base$modules.names[row]),
            label.base$modules.names[row],
            ifelse(
              !is.na(label.base$kos.name[row]),
              paste(label.base$kos.name[row], label.base$pfams[row], sep = " :: "),
              label.base$pfams[row]
            )
          )
          labs <- c(labs, row.lab)
        }
        label.base$plot.label <- labs
        label.data <- merge(
          label.base,
          doBy::summaryBy(log.rel.abund ~ unit, plot.data, FUN = max),
          by = "unit"
        )
        unit.label.list[[v]] <- data.frame(
          "var" = v,
          "unit" = label.data$unit,
          "relationship" = ifelse(label.data$spearman.r < 0, "-",
                                  ifelse(label.data$spearman.r > 0, "+", "flat")
          ),
          "label" = label.data$plot.label
        )
        reg.lines <- data.frame(NULL)
        u <- unique(as.character(plot.data$unit))[1]
        reg.slopes <- data.frame(
          matrix(ncol = 2, nrow = length(unique(as.character(plot.data$unit))))
        ) %>% setNames(c("unit", "reg.slope"))
        row.names(reg.slopes) <- unique(as.character(plot.data$unit))
        for (u in unique(as.character(plot.data$unit))) {
          print(u)
          model.data <- subset(plot.data, unit == u)[, c(2, 3, 8)] # Originally 11 instead of 8
          # unit.df <- readRDS(
          #     file.path(saveDir, dfsDir, paste("df_", v, "-", u, ".rds", sep=""))
          # )
          file.vu <- paste(dfsDir, "/df_", v, "-", u, ".rds", sep = "")
          #untar(file.path(saveDir, dfsTar),
          #      files = file.vu,
          #      compressed = "gzip"
          #)
          setwd ("/Users/isabelle/Downloads/flannery_stagaman_analysis-master/Saved_objs_only_data/")
          unit.df <- readRDS(file.vu)
          #file.remove(file.vu)
          #file.remove(dfsDir)
          unit.df$unit <- as.character(unit.df$unit)
          reg.slopes[u, c("unit", "reg.slope")] <- unit.df[1, c("unit", "m")]
          reg.lines <- rbind(reg.lines, unit.df)
          unit.plot <- ggplot(data = model.data, aes_string(x = v, y = "log.rel.abund")) +
            geom_abline(intercept = unit.df$b, slope = unit.df$m, color = "blue", size = 0.5) +
            geom_point(aes(shape = ifelse(rel.abund == 0, "zero", "not zero")),
                       alpha = 0.6
            ) +
            geom_text(
              data = subset(label.data, unit == u),
              aes(
                label = ifelse(is.na(plot.label),
                               "",
                               mgsub(c(", ", "; "),
                                     c("\n", "\n"),
                                     plot.label,
                                     fixed = T
                               )
                ),
                x = ifelse(spearman.r < 0,
                           max(plot.data[, v], na.rm = T),
                           min(plot.data[, v], na.rm = T)
                ),
                y = log.rel.abund.max,
                hjust = ifelse(spearman.r < 0, 1, 0)
              ),
              vjust = 1, size = 2.5, color = "red"
            ) +
            scale_shape_manual(name = "Relative abundance", values = c(16, 1)) +
            labs(
              y = paste("log(", u, " Rel. Abund.)", sep = ""),
              x = mgsub(
                c(vars.legend$var.name),
                c(vars.legend$description),
                v
              ),
              title = "Linear models on KO abundance (CPGLM)"
            ) +
            theme(
              legend.position = "top",
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 10),
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16)
            )
          setwd ("/Users/isabelle/Downloads/flannery_stagaman_analysis-master/")
          ggsave(unit.plot,
                 file = file.path(
                   unitPlotDir, paste(gsub(".", "_", v, fixed = T), "-", u, ".png", sep = "")
                 ),
                 height = 6,
                 width = 8,
                 dpi = 300,
                 limitsize = FALSE
          )
        }
        unit.label.list[[v]]$reg.slope <- reg.slopes$reg.slope
        var.plot <- ggplot(data = plot.data, aes_string(x = v, y = "log.rel.abund")) +
          geom_abline(
            data = reg.lines, aes(intercept = b, slope = m),
            color = "blue", size = 0.5
          ) +
          geom_point(aes(shape = ifelse(rel.abund == 0, "zero", "not zero")),
                     alpha = 0.6
          ) +
          geom_text(
            data = label.data,
            aes(
              label = ifelse(is.na(plot.label),
                             "",
                             mgsub(c(", ", "; "),
                                   c("\n", "\n"),
                                   plot.label,
                                   fixed = T
                             )
              ),
              x = ifelse(spearman.r < 0,
                         max(plot.data[, v], na.rm = T),
                         min(plot.data[, v], na.rm = T)
              ),
              y = log.rel.abund.max,
              hjust = ifelse(spearman.r < 0, 1, 0)
            ),
            vjust = 1, size = 2.5, color = "red"
          ) +
          scale_shape_manual(name = "Relative abundance", values = c(16, 1)) +
          labs(
            y = "log(Unit Relative Abundance)",
            x = mgsub(
              c(vars.legend$var.name, ".resid"),
              c(vars.legend$description, " Residuals"),
              v
            ),
            title = "Linear models on KO abundance (CPGLM)"
          ) +
          facet_wrap(~unit, scales = "free_y") +
          theme(
            legend.position = "top",
            axis.text.y = element_blank()
          )
        kos.plot.list[[v]] <- var.plot
        nCol <- ceiling(sqrt(length(unique(plot.data$unit))))
        nRow <- ceiling(length(unique(plot.data$unit)) / nCol)
        ggsave(
          var.plot,
          file = file.path(varPlotDir, paste(gsub(".", "_", v, fixed = T), "png", sep = ".")),
          height = nRow * 2,
          width = 8,
          dpi = 300,
          limitsize = FALSE
        )
      }
    }
  }
  # kos.plot.list <- kos.plot.list[!sapply(kos.plot.list, is.null)]
  # unit.label.df <- do.call("rbind", e$unit.label.list)
  unit.label.df <- do.call("rbind", unit.label.list)
  unit.label.df[, 1:4] <- apply(unit.label.df[, 1:4], 2, as.character)
  row.names(unit.label.df) <- NULL
  unit.label.df$var <- mgsub(
    c(vars.legend$var.name, ".resid"),
    c(vars.legend$description, " Residuals"),
    unit.label.df$var
  )
  names(unit.label.df) <- c("Variable", "Unit", "Relationship", "Label", "Slope")
  write.table(
    unit.label.df,
    file = file.path(saveDir, "df_sigVarKOrels_condsInModel.txt"),
    sep = "\t", quote = F, row.names = F
  )
  saveRDS(unit.label.df, file = file.path(saveDir, "df_sigVarKOrels_condsInModel.rds"))
  unit.label.dt <- as.data.table(unit.label.df)
  mods.label.dt <- unit.label.dt[grepl("MODS", Label, ignore.case = F)]
  write.table(mods.label.dt,
              file = file.path(saveDir, "df_sigVarMODrels_condsInModel.txt"),
              sep = "\t", quote = F, row.names = F
  )
  saveRDS(mods.label.dt, file = file.path(saveDir, "dt_sigVarMODrels_condsInModel.rds"))

  nRows <- sapply(kos.plot.list, function(plot) {
    panels <- max(as.numeric(ggplot_build(plot)$data[[1]]$PANEL))
    nCol <- ceiling(sqrt(panels))
    ceiling(panels / nCol)
  })
  plot.heights <- nRows / sum(nRows)
  kos.plot.height <- ifelse(sum(nRows) / 4 < length(nRows),
                            length(nRows) * 8,
                            (sum(nRows) / 4) * 6
  )
}
```


kos.network.cap <- ""



fig.wd <- pnas.2col #/ 0.8887
fig.ht <- pnas.max.ht / 0.965
lab.font.pt <- 6
mods.vars.rels <- mods.label.dt

# Remove the MODS: heading
mods.vars.rels$Label <- gsub('MODS: ', '', mods.vars.rels$Label)
#
plot.data0 <- mods.vars.rels[, .(modComps = length(Label), Label, Slope), by = Variable]
plot.data <- plot.data0[
  ,
  .(varComps = length(Variable), modComps, Variable, Slope),
  by = Label
]
# plot.data$Type <- ifelse(grepl("Insulin", plot.data$Variable), "GTT_AUC", "HOMA_IR")
# plot.data$Label <- strtrim(plot.data$Label, median(nchar(plot.data$Label)))
# grep("manno-heptose", unique(plot.data$Label), value=F)
pnas.max.ht <- 9 # inches
fig.wd <- pnas.2col / 0.8887
fig.ht <- pnas.max.ht / 0.965
new.par <- par(ps = lab.font.pt, family = "Helvetica", fin = c(fig.wd, fig.ht))
to.replace <- data.frame(
  "word" = c(
    "pathway", "; ", "beta", "biosynthesis",
    "system", "polysaccharide", "regulatory", "transport",
    "two", "component", "degradation", "pathogenicity",
    "signature", "deoxyribonuleotide", "zinc", "manganese",
    "; Cell cycle - G2/M transition", "antimicrobial peptide", "sulfate", "Putative",
    "phosphate", "hydroxymethylpyrimidine", "oligosaccharyltransferase", "Oxidation",
    "semialdehyde", "Entner-Doudoroff", "ADP-L-glycero-D-manno-heptose",
    "BRCA1-associated genome surveillance complex ", "; MRX complex"
  ),
  "replace" = c(
    "pwy", ";", "β", "biosyn",
    "sys", "polysac", "reg", "transp",
    "2", "comp", "degrad", "pathogenic",
    "signtr", "dNMP", "Zn", "Mn",
    "", "AMP", "SO4", "Put.",
    "PO4", "HMP", "OST", "oxid'n",
    "semialdhd", "ED", "ADP-LDHep",
    "", ""
  ),
  stringsAsFactors = F
)

# The following fixes the labels (I'm not using this anymore)
max.wd <- 1.26
for (l in 1:length(plot.data$Label)) {
  lab <- plot.data$Label[l]
  lab <- sub("\\(.+\\)[ ;] *", "", lab)
  wd0 <- strwidth(lab, units = "inches")
  if (wd0 > max.wd) {
    lab <- mgsub(to.replace$word, to.replace$replace, lab)
    wd1 <- strwidth(lab, units = "inches")
    if (wd1 > max.wd) {
      lab <- sub(";.+", "", lab)
      wd2 <- strwidth(lab, units = "inches")
      if (wd2 > max.wd) {
        print(lab, quote = F)
      }
    }
  }
  plot.data$Label[l] <- lab
}
par(.pardefault)

plot.data$Type <- c(rep ("Path", nrow(plot.data)))

## Add the bacterial data
tax.section.new <- tax.section [,1:7]
tax.section.new$Type <- c(rep ("Bug", nrow(tax.section)))
tax.section.new <- tax.section.new [,-4]
colnames (tax.section.new) <- colnames (plot.data)
# Combine the tables
plot.data <- rbind (plot.data, tax.section.new, fill = TRUE)
##

## Change plot data columns around
plot.data$Label2 <- plot.data$Label
plot.data$Label <- plot.data$Variable
plot.data$Variable <- plot.data$Label2
plot.data$varComps2 <- plot.data$varComps
plot.data$varComps <- plot.data$modComps
plot.data$modComps <- plot.data$varComps2
plot.data <- plot.data[,1:6]
##

ordered.mods <- levels(with(plot.data, reorder(Label, varComps)))
ordered.chb <- levels(with(plot.data[Type != "Path"], reorder(Variable, modComps)))
ordered.esa <- levels(with(plot.data[Type == "Path"], reorder(Variable, modComps)))
range <- c(-5, 5)
l <- mean(range) - 1.5
r <- mean(range) + 1.5
# r <- 5
#
seg.data <- within(plot.data, {
  x0 <- ifelse(Type != "Path", range[1], range[2])
  y0 <- ifelse(Type != "Path",
               (match(Variable, ordered.chb) - 0.5) / length(ordered.chb),
               (match(Variable, ordered.esa) - 0.5) / length(ordered.esa)
  )
  x1 <- ifelse(Type != "Path", l, r)
  y1 <- (match(Label, ordered.mods) - 0.5) / length(ordered.mods)
  mid <- mean(range)
})
seg.data$Variable <- sub("/", "/\n", seg.data$Variable)
seg.data <- seg.data[order(abs(Slope))]
lab.size <- lab.font.pt * 0.35
# We need to swap the labels with the variables
#seg.data$Label2 <- seg.data$Label
#seg.data$Label <- seg.data$Variable
#seg.data$Variable <- seg.data$Label2
#seg.data <- seg.data[,1:11]
#seg.data$varComps2 <- seg.data$varComps
#seg.data$varComps <- seg.data$modComps
#seg.data$modComps <- seg.data$varComps2
#seg.data <- seg.data[,1:11]

## Wrap the text
#wrapit <- function(text) {
#  wtext <- paste(strwrap(text,width=40),collapse=" \n ")
#  return(wtext)
#}

#library (ggfittext)

#test <- seg.data$Label
#data$wrapped_text <- lapply(test, wrapit)
#seg.data$wrapped_text <- unlist(seg.data$wrapped_text)

mods.line <- ggplot(data = seg.data) +
  geom_segment(
    aes(
      x = x0, xend = x1,
      y = y0, yend = y1,
      color = Slope,
      size = abs(Slope),
      alpha = abs(Slope)
    )
  ) +
  scale_color_gradientn(name = "Reg. Slope", colors = blue_to_red) +
  scale_size_continuous(name = "abs(Reg. Slope)", range = c(0.05, 2.5)) +
  scale_alpha_continuous(name = "abs(Reg. Slope)", range = c(0.6, 2.5)) +
  geom_tile(aes(x = mid, y = y1), width = r - l, height = 0.05, fill = "white", color = "black") +
  geom_text(
    data = unique(seg.data[, .(Label, mid, y1)]),
    aes(x = mid, y = y1, label = Label),
    size = lab.size
  ) +
  geom_tile(
    data = seg.data[Type != "Path"], aes(x = x0 + l, y = y0),
    width = r - l, height = 0.05, fill = "white", color = "black"
  ) +
  geom_text(
    data = unique(seg.data[Type != "Path", .(Variable, x0, y0)]),
    aes(x = x0 + l, y = y0, label = Variable),
    size = lab.size - 0.55
  ) +
  geom_tile(
    data = seg.data[Type == "Path"], aes(x = x0 + r, y = y0),
    width = r - l, fill = "white", color = "black"
  ) +
  geom_text(
    data = unique(seg.data[Type == "Path", .(Variable, x0, y0)]),
    aes(x = x0 + r, y = y0, label = Variable),
    size = lab.size - 0.95
  ) +
  scale_x_continuous(
    # breaks = c(min(seg.data$x0) + l, -1.5 + r),
    breaks = c(min(seg.data$x0) + l,0, max(seg.data$x0) + r),
    # labels = c("Phenotypic Outcomes", "Functional Pathways"),
    labels = c("Taxa (Family and Genus)", "Phenotypic Outcomes", "Functional Pathways"),
    position = "top"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  # geom_fit_text(inherit.aes = TRUE, label = Label, min.size = 1, xmin = 0, xmax = 10) +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(11, "pt"),
    legend.box.margin = unit(c(0, 0, 0, 0), "pt"),
    legend.box.spacing = unit(5, "pt")
  )
mods.line

ggsave(mods.line, file=file.path(saveDir, "mainFig_chb_esa_kos_regNetwork_fdr_best.pdf"), width=fig.wd, height=fig.ht)


ggsave(mods.line, file=file.path(saveDir, "mainFig_chb_esa_kos_regNetwork_all.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file=file.path(saveDir, "mainFig_chb_esa_kos_regNetwork.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file=file.path(saveDir, "mainFig_chb_esa_kos_regNetwork_hfd_hxn.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file=file.path(saveDir, "mainFig_chb_esa_kos_regNetwork_genes.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file = file.path(saveDir, "mainFig_chb_esa_kos_regNetwork.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file = file.path(saveDir, "mainFig_chb_esa_kos_regGenes.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.leg, file = file.path(saveDir, "mainFig_chb_esa_kos_regNetwork_legend.pdf"))
system("open -a /Applications/Inkscape.app/ Figure_plots/mainFig_chb_esa_kos_regNetwork.pdf")

ggsave(mods.line, file = file.path(saveDir, "mainFig_chb_esa_kos_new.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file = file.path(saveDir, "mainFig_XN_20210214.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file = file.path(saveDir, "mainFig_XN_20210217.pdf"), width=fig.wd, height=fig.ht)

ggsave(mods.line, file = file.path(saveDir, "mainFig_XNonly_20210218_all_pathways.pdf"), width=fig.wd, height=fig.ht)

# Find the most important pathways per phenotype

path.list <- as_tibble (mods.vars.rels)
path.list.sorted <- path.list %>%
  arrange (Variable, desc(Slope))
View(path.list.sorted)

# Save these pathways
writexl::write_xlsx(path.list.sorted, "Figure7_Pathways.xlsx") # This includes pathways only in all samples
writexl::write_xlsx(path.list.sorted, "Figure6_Pathways.xlsx") # this includes ALL pathways.

writexl::write_xlsx(path.list.sorted, "XN_Pathways.xlsx") # this is HFD - HXN only
writexl::write_xlsx(path.list.sorted, "XN_Pathways_20200218.xlsx") # this is HFD - HXN only with ALL pathways!


# system("open -a /Applications/Inkscape.app/ Figure_plots/mainFig_chb_esa_kos_regNetwork.pdf")
###  ISSC Talk FIG
# fig.wd <- pnas.2col/0.8887
mods.leg <- cowplot::get_legend(mods.line)
mods.line <- mods.line + theme(legend.position = "none")
ggsave(mods.line, file = file.path(figDir, "mainFig_chb_esa_kos_regNetwork.pdf"), width=fig.wd, height=fig.ht)
ggsave(mods.leg, file = file.path(figDir, "mainFig_chb_esa_kos_regNetwork_legend.pdf"))
system("open -a /Applications/Inkscape.app/ Figure_plots/mainFig_chb_esa_kos_regNetwork.pdf")
# ggsave(mods.line,
#        file = "~/Google Drive/OSU/OSU_Home_Sync/ISSC_talk/mainFig_chb_esa_kos_regNetwork.pdf",
#        width = 190, height = 190 * (fig.ht / fig.wd), units = "mm"
# )
# # system("open '/Users/stagamak/Google Drive/OSU/OSU_Home_Sync/ISSC_talk/mainFig_chb_esa_kos_regNetwork.pdf'")
# system("open -a /Applications/Inkscape.app/ '/Users/stagamak/Google Drive/OSU/OSU_Home_Sync/ISSC_talk/mainFig_chb_esa_kos_regNetwork.pdf'")
```

<!-- ### T6SS KOs/Modules and Taxa Relationships -->

  ```{r t6ss-taxa, fig.height=pnas.2col, eval=FALSE}
set.seed(main.seed)
kos.df <- readRDS(file = file.path(saveDir, "df_sigVarKOrels_condsInModel.rds"))
kos.dt <- as.data.table(kos.df)
t6ss.dt <- kos.dt[grepl("MODS:", Label, ignore.case = T)]
t6ss.kos <- sort(unique(t6ss.dt$Unit))
kos.tbl <- as(otu_table(stool.kos.phy), "matrix")
t6ss.tbl <- kos.tbl[, t6ss.kos]
t6ss.m.dt <- melt(t6ss.tbl, varnames = c("smpl", "kos"), value.name = "rel.kos.abund") %>% as.data.table()
tax.tbl <- as(otu_table(stool.tax.phy), "matrix")
tax.m.dt <- melt(tax.tbl, varnames = c("smpl", "taxon"), value.name = "rel.tax.abund") %>% as.data.table()
t6ss.tax.data <- merge(tax.m.dt, t6ss.m.dt, by = "smpl", all = T, allow.cartesian = T)
cpglm.res <- data.frame()
# txn <- colnames(tax.tbl)[1]
# t6.ko <- t6ss.kos[1]
for (txn in colnames(tax.tbl)) {
  for (t6.ko in t6ss.kos) {
    comp.data <- t6ss.tax.data[taxon == txn & kos == t6.ko]
    if (median(comp.data$rel.tax.abund) > 0 & median(comp.data$rel.kos.abund) > 0) {
      comp.data$log.kos.abund <- log(comp.data$rel.kos.abund)
      tax.t6.cpglm <- try(
        cpglm(rel.tax.abund ~ log.kos.abund, data = comp.data),
        silent = T
      )
      if (class(tax.t6.cpglm) != "try-error") {
        sink(file = "tmp.txt")
        model.sum <- summary(tax.t6.cpglm)
        sink()
        pval <- model.sum$coefficients[2, "Pr(>|t|)"]
        slope <- unname(tax.t6.cpglm$coefficients[2])
        intercept <- unname(tax.t6.cpglm$coefficients[1])
        cpglm.df <- data.frame(
          "taxon" = txn,
          "kos" = t6.ko,
          "m" = slope,
          "b" = intercept,
          "unadj.p" = pval
        )
        cpglm.res <- rbind(cpglm.res, cpglm.df)
      }
    }
  }
}
# View(cpglm.res)
cpglm.res$adj.p <- p.adjust(cpglm.res$unadj.p, method = "fdr")
cpglm.res <- as.data.table(cpglm.res)
sig.assoc <- cpglm.res[adj.p <= 0.05]
# View(sig.assoc)
sig.assoc$n <- LETTERS[1:nrow(sig.assoc)]
sig.assoc <- within(sig.assoc, {
  # We want genus
  # tax.species <- str_extract(taxon, "s__[a-zA-Z0-9 _]+")
  tax.species <- str_extract(taxon, "g__[a-zA-Z0-9 _]+")
  tax.species <- mgsub(c("g__", "_bacterium"), c("", "", "715"), tax.species)
  # plot.lab <- paste(n, " ", tax.species, "\n", kos, sep="")
})
sig.assoc$plot.lab <- factor(sig.assoc$plot.lab, levels = unique(sig.assoc$plot.lab))
# # sig.cors$taxon
sig.t6ss.tax <- t6ss.tax.data[sig.assoc, on = c(taxon = "taxon", kos = "kos")]
sig.t6ss.tax <- within(sig.t6ss.tax, {
  shp <- ifelse(rel.kos.abund == 0 & rel.tax.abund == 0, "zero-zero",
                ifelse(rel.kos.abund == 0, "ko-zero",
                       ifelse(rel.tax.abund == 0, "tax-zero", "no-zero")
                )
  )
})
t6ss.tax.plot <- ggplot(data = sig.t6ss.tax) +
  geom_abline(
    data = sig.assoc, aes(intercept = b, slope = m),
    color = "blue", size = 0.5
  ) +
  geom_point(aes(x = log(rel.kos.abund), y = log(rel.tax.abund), shape = shp), alpha = 0.6) +
  scale_shape_manual(name = "Relative abundance", values = c(16, 1), labels = c("not zero", "zero")) +
  labs(
    y = "log(Taxon Relative Abundance)",
    x = "log(KO Relative Abundance)"
  ) +
  facet_wrap(~n, scales = "free") +
  theme(
    legend.position = "top",
    strip.text = element_text(hjust = 0)
  )
plot.grob <- arrangeGrob(t6ss.tax.plot)
t6ss.dt$Tbl.Label <- mgsub(c(":: ", "; "), c("\n", "\n    "), t6ss.dt$Label)
tax.kos <- sig.assoc[, .(Panel = n, Taxon = gsub("_", " ", tax.species), KO = kos)]
ff <- matrix(c(1, 3, 1), ncol = ncol(tax.kos), nrow = nrow(tax.kos), byrow = T)
tt1 <- ttheme_default(base_size = 8, core = list(fg_params = list(fontface = as.vector(ff))))
tax.ko.tbl <- tableGrob(tax.kos,
                        rows = NULL, theme = tt1
)
kos.labs <- unique(t6ss.dt[order(Unit), .(KO = Unit, Label = Tbl.Label)])
hj <- matrix(c(0.5, 0), ncol = ncol(kos.labs), nrow = nrow(kos.labs), byrow = T)
x <- matrix(c(0.5, 0.01), ncol = ncol(kos.labs), nrow = nrow(kos.labs), byrow = T)
tt2 <- ttheme_default(
  core = list(fg_params = list(hjust = as.vector(hj), x = as.vector(x))),
  colhead = list(fg_params = list(hjust = c(0.5, 0), x = c(0.5, 0.01))),
  base_size = 8
)
ko.lab.tbl <- tableGrob(kos.labs, rows = NULL, theme = tt2)
table.grob <- arrangeGrob(tax.ko.tbl, ko.lab.tbl, nrow = 1)
t6ss.fig <- arrangeGrob(plot.grob, table.grob,
                        nrow = 2,
                        heights = c(0.6, 0.4)
)
ggsave(t6ss.fig,
       file = file.path(figDir, "suppFig_t6ss_tax_regs.pdf"),
       width = pnas.2col + 15, height = pnas.max.ht + 30
       # 2.5 and 2 were original
)
system(paste("open", file.path(figDir, "suppFig_t6ss_tax_regs.pdf")))
# system(paste("open -a /Applications/Inkscape.app", file.path(figDir, "suppFig_t6ss_tax_regs.pdf")))
```

```{r addressing-some-reviewer-comments, eval=FALSE}
smpl.data <- as(sample_data(stool.kos.phy), "data.frame")
ggplot(smpl.data, aes(x = parent_PSI_PC_dysfunc_t1, y = esa_LEC_Turmoil_Combined)) +
  geom_point()
cor.test(smpl.data$parent_PSI_PC_dysfunc_t1, smpl.data$esa_LEC_Turmoil_Combined)
cor.test(
  smpl.data$parent_PSI_PC_dysfunc_t1,
  smpl.data$esa_LEC_Turmoil_Combined,
  method = "spearman"
)
ggplot(smpl.data, aes(x = parent_PSI_PC_dysfunc_t1, y = child_behav_CBCL_depression_t1)) +
  geom_point()
ggsave(file = "Figure_plots/parent-child_dys_vs_CBCL-depression_plot.pdf")
cor.test(smpl.data$parent_PSI_PC_dysfunc_t1, smpl.data$child_behav_CBCL_depression_t1)
cor.test(
  smpl.data$parent_PSI_PC_dysfunc_t1,
  smpl.data$child_behav_CBCL_depression_t1,
  method = "spearman"
)
```




## Generate an output table for the gene associations
slopes <- path.list %>%
  arrange (Variable, desc(Slope))
p.corr <- kos.sig.pvals %>%
  arrange (var, adj.pval)
p.corr$var <- as.factor(p.corr$var)
levels (p.corr$var) <- c("GTT_AUC", "Liver Total Cholesterol", "Liver TAG",
                         "Plasma TAG", "Feed Efficiency")
names (p.corr) <- c("Variable", "Unit", "adj.pval", "spearman.r")
super.sheet <- dplyr::inner_join(p.corr, slopes, by = "Unit")
View(super.sheet)
# Get rid of non-matching rows
super.new <- super.sheet %>% filter (Variable.x == Variable.y)
super.new <- super.new %>%
  select (-Variable.y) %>%
  select (-Label)
super.new <- super.new %>%
  arrange (Variable.x, desc (Slope))
colnames (super.new) <- c("Phenotype", "Unit", "adj.pval", "spearman.r",
                          "Relationship", "Slope")
writexl::write_xlsx(super.new, "genes_output.xlsx")
