#####################################################################
## This project is about analyzing the 16S sequencing data         ##
## for all of the conventional Swiss Webster mouse feeding studies ##
#####################################################################

## Load the required packages
library(dada2)
library(DECIPHER)
library (tibble)
library(biomaRt) # do not need this one

getwd()

## Set the working directory (it is also part of the R-Studio project file though)
## Change this accordingly on your own computer

setwd ("/Users/isabelle/16S/Bbduk")

# path <- "Raw"

# 1st trial
path <- "End_Set"

# 2nd trial
path <- "Start_Set"

list.files(path)

## Set the full path to the Forward end reads (R1) and sort them, 
## so I can reference these later

forwardReads <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE))

# Was this a success?

head(forwardReads)
# Yay! :)

## Now do the same for the paired end reads

pairedReads <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE))

# Did this work?

head(pairedReads)
## It did!!

# Now, we will extract the names of the groups in our dataset, by using 
# several functions in tandem. This will clip the beginning section (,1) 
# of each forwardRead fastq file (the sample name) and save is as a new object

sample.names <- sapply (strsplit(basename(forwardReads), "_"), '[', 1)

# Did we get the sample names we wanted?

sample.names
# We did! Cool!

#########################################################################
### Assess read quality distribution
### This is basically what we did with FastQC
#########################################################################

plotQualityProfile(forwardReads, aggregate=TRUE)
# The forward reads look fantastic
plotQualityProfile(pairedReads, aggregate=TRUE)
# The paired reads, not so much, it appears that we can only use about 220 
# nucleotides or so.

##########################################################################
### Quality filter and trim reads
##########################################################################

# We will now filter these reads for missing data (such as 'N' nucleotides), 
# maximum expected errors (maxEE), and truncate the forward reads to 250 nt 
# and the paired reads to 200 nt.
#
# Let's create a list of new output names, in a new 'filtered' directory, 
# that we will call our new trimmed and cleaned FASTQ files. We will also 
# have these compressed later on, to save space (convert each to a *.gz file)

forwardReads_filt <- file.path(path, "filtered", paste0(sample.names, 
                                                        "_F_filt.fastq.gz"))
# Which will print out
head(forwardReads_filt)

# do the same thing for the paired end reads
pairedReads_filt <- file.path(path, "filtered", paste0(sample.names,
                                                       "_P_filt.fastq.gz"))
head(pairedReads_filt)
# now we are ready to filter the data

# Truncate the R1 files to 250 nt and the R2 files to 200 nt.
# Set the maxEE to two, this is just the "normal" set in the course.
# Set compress to TRUE, because we want to zip these files

filteredFastq <- filterAndTrim(forwardReads, forwardReads_filt, 
                               pairedReads, pairedReads_filt,
                               truncLen = c(280,220),
                               minLen = 50,
                               maxN = 0,
                               maxEE = c(2,4),
                               truncQ=2,
                               rm.phix=TRUE,
                               compress = TRUE, 
                               multithread = TRUE,
                               verbose = TRUE)

head(filteredFastq)

# How many total reads were filtered out?
sum(filteredFastq[,1])-sum(filteredFastq[,2])
# I filtered out 478540 reads!!

##########################################################################
### Learn error rates
##########################################################################

# Now we will learn the error rates of the filtered forward reads. 
error_forward <- learnErrors(forwardReads_filt, multithread=TRUE)

# Now do the same for the paired end reads
error_paired <- learnErrors(pairedReads_filt, multithread=TRUE)

# Let's just plot the error rates
plotErrors (error_forward, nominalQ = TRUE)

# And for the paired end reads
plotErrors (error_paired, nominalQ = TRUE)

##########################################################################
### De-replicate reads
##########################################################################

# To de-replicate our forward filtered reads, run the following code:

forwardReads_filt_derep <- derepFastq(forwardReads_filt, verbose=TRUE)

# Now we need to do the same for our paired-end reads

pairedReads_filt_derep <- derepFastq(pairedReads_filt, verbose=TRUE)

# Let's look at these:
head(pairedReads_filt_derep)

# If we look at the 'names' slot for this object, 
# we will see that it is listing the entire file name, so let's use shorter names
names(forwardReads_filt_derep) <- sample.names
names(pairedReads_filt_derep) <- sample.names

# so, let's look at these names
names(forwardReads_filt_derep)
names(pairedReads_filt_derep)

# Great thus far!
# let's just look at one sample (048) from one of these 'dereplicated' 
# lists (R2) and you will see the following slots:

pairedReads_filt_derep[["048"]]
# We get the following:

# derep-class: R object describing dereplicated sequencing reads
# $uniques: 40118 reads in 22946 unique sequences
# Sequence lengths: min=220, median=220, max=220
# $quals: Quality matrix dimension:  22946 220
# Consensus quality scores: min=7, median=37, max=38
# $map: Map from reads to unique sequences:  51 13909 16068 6276 24 ...
# head (pairedReads_filt_derep[["048"]][["uniques"]],1)

# the average quality score for each sequence

head (pairedReads_filt_derep[["048"]][["quals"]][,1])
# Which gives the following:

# CCTGTTCGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGAAGGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCCTTCTTCTCGCCCACTCAAGGCCCCCAGTTTCAACGGCCGGACGGGGTTGAGCCCCGAATTTTTACCGCTGACTTAAAAGCCCGCCTACGCACCCTTTAAACCC 
# 33.97761 
# CCTGTTCGATACCCGCGCTTTCGAGCTTCAGCGTCAGTAACGCTGCGGCAGGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCTGCCTTCCGCGCACTCAAGGTCTCCAGTTCGCGCCGCACTCTCATGGTTGAGCCACAAAATTTCACGGCACGCTTAAAGACCGGCCTACGCTCCCTTTAAACC
# 33.94068 
# CCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGAACGCCAGTTTCAACGGCTCGACGGCGTTGAGCACCGCTTTTTTACCGCTGACTTGGCATCCCGCCTACGCACCCTTTAAACCC 
# 33.97649 
# CCTGTTCGATACCCACGCTTTCGTGCTTCAGCGTCAGTTGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGACCGGCAGTTTCAACGGCTCGATGAGGTTGAGCCCCACAATTTTACCGCTGACTTACCAGCCCGCCTACGCACCCTTTAAACCC 
# 33.97094 
# CCTGTTCGATACCCACGCTTTCGTGCTTCAGCGTCAGTTGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGACCGGCAGTTTCAACGGCTCGATGAGGTTGAGCCTCACAATTTTACCGCTGACTTACCAGCCCGCCTACGCACCCTTTAAACCC 
# 33.91643 
# CCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCTAGAAAGACAGTTTCAACGGCGAGACGGGGTTGAGCCCCGTCGTTTGACCGCTGACTTGCACTTTCCGCCTGCGCACCCTTTAAACC 
# 33.92456 

# and how each read is assigned to each unique sequence (or map)
head(pairedReads_filt_derep[["048"]][["map"]])
# Which gives us:
# 51 13909 16068  6276    24 16915
# Looks good so far.

##########################################################################
# First, we will use Dada2 to "denoise" the data
# To run this algorithm on our dataset, we proceed with the forward and 
# then the paired 'de-replicated' reads and specify the 'learned error 
# rate' in our datasets for read correction. Let's start with the forward 
# reads

dadaForward <- dada (forwardReads_filt_derep, err=error_forward, 
                     multithread=TRUE)

# print out information on the screen about the number of total reads 
# and unique reads per group
# This new object will contain a lot of information for each sample, 
# including the de-noised data, the clustering of each read, 
# the original error values, updated error values, p-values, etc.

head (dadaForward)
# This gives us the following:
# $`021`
# dada-class: object describing DADA2 denoising results
# 324 sequence variants were inferred from 5554 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
#

dadaPaired <- dada (pairedReads_filt_derep, err=error_paired,
                     multithread = TRUE)
# and look at it
head (dadaPaired)
# We get the following result:
# $`021`
# dada-class: object describing DADA2 denoising results
# 198 sequence variants were inferred from 18532 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
# 
##########################################################################
##########################################################################
##########################################################################
### Generate Contigs                                                  ####
##########################################################################
##########################################################################
##########################################################################
# Now that we have trimmed our samples, de-replicated them, and de-noised
# them, the next step is to create contigs out of the forward and 
# paired-end reads, or merge them together. This is performed by taking 
# the reverse compliment of the paired end read, and aligning it to the 
# forward read. The default states that they must overlap by at least 12 
# nucleotides and have 100% matching identity at these regions of 
# overlap. This function requires both the 'de-noised' file and 
# 'de-replicated' file for both forward and paired. It will merge each 
# de-noised read pair.

contigs <- mergePairs(dadaForward, forwardReads_filt_derep, 
                      dadaPaired, pairedReads_filt_derep, minOverlap = 20,
                      maxMismatch = 0, returnRejects = FALSE,
                      propagateCol = character(0), justConcatenate = FALSE,
                      trimOverhang = TRUE, verbose = FALSE)

# How many unique contigs were formed for sample 48E?

length(contigs[["048"]][["sequence"]])

# Which gives:
# 218
# You can get a lot more information by how these merged by:

contigs[["048"]]

# This will print out the sequences for sample 048
# However, it also gives you a nice summary table.
#     abundance forward reverse nmatch nmismatch nindel prefer accept
# 1        5266       1       2     76         0      0      1   TRUE
# 2        3588       3       1     76         0      0      1   TRUE
# 4        2356       4       3     76         0      0      1   TRUE
# 5        1824       5       8     76         0      0      1   TRUE
# 6        1801       6       8     76         0      0      1   TRUE

# You will see the number of reads found in each 'unique' sequence how 
# it matched. Each row is listed by an index for easy cross-referencing.
# Great overlap, some samples were even over 90!

#########################################################################
#########################################################################
### Create an Amplicon Sequence Variant Table
#########################################################################
#########################################################################

# This is essentially merging all of the contigs into one table 
# (well, a matrix)
# To create the variant table, we run:

seq_table <- makeSequenceTable(contigs)

# You can 'View(seq_table)' if you wish, but it is a matrix constructed 
# with samples as rows, sequence variant as columns, and abundance 
# (reads) as cells.
# Now, let's examine a histogram of sorts of our amplicon sequence 
# variant length distribution with the following command

table (nchar(getSequences(seq_table)))

# This should show the following length distribution:
# 285 336 365 403 404 405 406 407 408 409 410 411 416 418 420 421 422 423 424 425 426 427 
# 1   1   1  24 388 214 106 220   7   2  11   1   4   1   1  18  13 124 270  96   4   1 
# 428 429 430 456 
# 29 152  17   1 

# Generate a graph of these
distribution <- table (nchar(getSequences(seq_table)))
barplot(distribution)

# We can exclude the singletons 
# Here you can pick which parts of the ASV distribution you want to keep

# Select using an index
index <- c(399, 400, 401, 402, 403, 417, 418, 419, 420, 421, 424)
seq_table_length <- seq_table[,nchar(colnames(seq_table)) %in% index]

# Select a range
# I used three ranges, but only the widest range allowed for taxonomy assignment
seq_table_length <- seq_table[,nchar(colnames(seq_table)) %in% 403:410]
seq_table_length <- seq_table[,nchar(colnames(seq_table)) %in% 421:430]
seq_table_length <- seq_table[,nchar(colnames(seq_table)) %in% 403:430] # <- this range

# Or select all
seq_table_length <- seq_table

# I chose to select a range
# If we look at the dimensions of this matrix, we will see:

dim (seq_table_length)
# 81 1703

###########################################################################
###########################################################################
###########################################################################
### Identify and remove Chimeras
###########################################################################
###########################################################################
###########################################################################

# The next step is to remove chimeric reads. Again, these are reads that 
# were improperly amplified during PCR and are a result of 
# 'mixed-ancestry'. Dada2 identified chimeric sequences by identifying 
# reads that can be reconstructed by two more abundant 'parental' 
# sequences.
# To identify and remove chimeric sequences, run the following command:

seq_table_nochimeri <- removeBimeraDenovo(seq_table_length, 
                                          method = "consensus",
                                          multithread = TRUE,
                                          verbose = TRUE)

# This should identify, remove chimeric sequences, and save it as a new 
# matrix
# This gives us the following result:
# Identified 528 bimeras out of 1703 input sequences.
# which means 31.00% of our variants were bimeric sequences

# The new number of bimeric sequences can be determined by checking the 
# dimensions of the new matrix by running:

dim (seq_table_nochimeri)
# which gives 81 1175

sum (seq_table_nochimeri)/sum(seq_table_length)
# which gives 0.993946

#######################################################################
#######################################################################
#######################################################################
### Sanity-check and Summary Table
#######################################################################
#######################################################################
#######################################################################

# Now we can do a summary of each step conducted thus far, as a sanity 
# check, and to see if we want to increase (or decrease) the stringency 
# of our methods
#
# To keep track of your steps, run the following:

getN <- function(x) sum(getUniques(x))
track <- cbind(filteredFastq, sapply(dadaForward, getN), 
               sapply(dadaPaired, getN), sapply(contigs, getN), 
               rowSums(seq_table_nochimeri))
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# If we print out the beginning of this object, we will see a nice summary 
# of our steps

head(track)

# We get the following:
#     input filtered denoisedF denoisedR merged nonchim
# 021 39955    34280     33995     34068  27607   27505
# 022 38094    32193     31944     31826  26110   26017
# 023 36754    31866     31640     31523  27800   27573
# 024 46070    39970     39652     39735  35278   35135
# 025 38187    33520     33144     33221  30106   30052
# 026 34203    29201     29050     28939  24826   24803
# Let's say you want to plot this? How would you do that? 
# Think about this for a minute!

library(reshape2)
library(ggplot2)

# Now let's manipulate (i.e. melt) this matrix into a different format to make 
# it easier to plot (as there are multiple columns).

plot_data <- melt (track, varnames = c("Sample","Stage"))

# If you look at this new object, you will see there are now only three columns, 
# as the number of reads (Values) has been set to a factor called Stage.

head (plot_data)

# But let's say we only want the first 60 rows for clarity?
# Load this package

library(dplyr)

# Now do the work
plot_data <- plot_data %>% arrange (Sample)
number_of_rows <- 60
plot_data_trunc <- plot_data[1:number_of_rows,]

# Okay, now you could plot this with ggplot with the following code and get 
# something super cool to show your peers!

ggplot(data = plot_data_trunc, aes(x = Sample, y = value, fill = Stage)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_classic() + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=16)) +
  theme(axis.title.x=element_blank())
ggsave('trimming_results_end.png', antialias="none")
ggsave('trimming_results.png', antialias="none")

############################################################################
############################################################################
############################################################################
### Classify Sequences 
############################################################################
############################################################################
############################################################################

# This is where we're going to have to filter the sequence table
# We want more than 769 reads to remain. 
# Should replace anything up to 769 with a zero.

# First try, this did not work
library(tidyverse)
seq_tibble <- as_tibble(seq_table_nochimeri)
for (i in names(seq_tibble))
{
  
}
seq_tibble <- seq_tibble %>% filter (seq_tibble[[]] > 630)

# install.packages('data.table') may need to be run if you don't have the
# package
# This also did not work
library(data.table)
outlierReplace = function(dataframe, cols, rows, newValue = NA) {
  if (any(rows)) {
    set(dataframe, rows, cols, newValue)
  }
}

# Still a no go
outlierReplace(seq_tibble, seq_tibble[[1]], 
               which(seq_tibble$TGAGGAATATTGGTCAATGGCCGCAAGGCTGAACCAGCCAAGTAGCGTGAGGGAAGACTGCCCTATGGGTTGTAAACCTCTTTTATGCGGGGATAAAGGCGTCCACGTGTGGATGTTTGCAGGTACCGCATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGTCTTTAAGCGTGCCGTGAAATTTTGTGGCTCAACCATGAGAGTGCGGCGCGAACTGGAGACCTTGAGTGCGCGGAAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCCTGCCGCAGCGTTACTGACGCTGAAGCTCGAAAGCGCGGGTATCGAACAGG < 630), NA)

# Ugh, doesn't work
for (i in 1:ncol(seq_tibble)) {
  for (j in 1:nrow(seq_tibble)) {
    if (seq_tibble[[i,1]] < 630) {
      seq_tibble[[i,1]] == 0
    }
#    outlierReplace(seq_tibble, seq_tibble[[i]], which (j < 630), 0)
  }
}

# Still nothing
for (i in 1:ncol(seq_tibble)) {
  rows <- c(1:nrow(seq_tibble))
  if (seq_tibble[[i,rows]] < 630) {
    seq_tibble[[i, rows]] == 0
  }
}

# Turns out it is this easy!
seq_table_nochimeri_test <- seq_table_nochimeri
seq_table_nochimeri_test[seq_table_nochimeri_test < 769] <- 0
seq_table_nochimeri[seq_table_nochimeri < 769] <- 0

# Now we are ready to assign the taxonomy of our sequence variants. Similar to many 
# other amplicon sequencing pipelines, Dada2 uses a naive Bayesian classifier method 
# to identify the probability that a given sequence came from a particular taxon. 
# We will be using the function assignTaxonomy which uses the sequence table as input, 
# along with a training dataset of known rRNA sequences, and generates taxonomy for 
# sequences that meet a minimum bootstrap confidence level. 
# You can assign your sequences using two different functions, depending on which 
# trainingset you use or what level of taxonomy you are after! These are the 
# assignTaxonomy and/or the addSpecies functions.
# First, let's try the Silva training set. (this is the training set I used for the paper)

taxa <- assignTaxonomy (seq_table_nochimeri, 
                        "~/16S/training/silva_nr_v132_train_set.fa.gz",multithread=TRUE)
taxa <- addSpecies (taxa, "~/16S/training/silva_species_assignment_v132.fa.gz")

# Then, lets try the RDP training set
# Run this command for the training set (not used)
taxa <- assignTaxonomy (seq_table_nochimeri_test, 
                        "~/16S/training/rdp_train_set_16.fa.gz",multithread=TRUE)
taxa <- addSpecies (taxa, "~/16S/training/rdp_species_assignment_16.fa.gz")

# Run this command for Phyloseq (not used)
taxa <- assignTaxonomy (seq_table_nochimeri, 
                        "~/16S/training/rdp_train_set_16.fa.gz",multithread=TRUE)
taxa <- addSpecies (taxa, "~/16S/training/rdp_species_assignment_16.fa.gz")

# The following code uses DECIPHER, which is faster (but I ended up not using)
dna <- DNAStringSet(getSequences(seq_table_nochimeri)) # Create a DNAStringSet from the ASVs
load("~/16S/Reference/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seq_table_nochimeri)

# simply set taxa <- taxid if you want to use this
taxa <- taxid

# now let's look at the 'taxa' matrix

View (taxa)

# You will see that for each sequence variant identified the classification 
# tree assigned to it, from Kingdom to species.
# Let's inspect the taxonomic assignments:

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(n = 20,taxa.print)

######################################################################################
######################################################################################
######################################################################################
### Use the mock community to evaluate error rate
######################################################################################
######################################################################################
######################################################################################

# Recall that one of the samples that was sequenced in this dataset was a 'mock' 
# community of known identity for 10 strains. Not all of the sequence variants from the 
# table were identified in the mock community. Remember, that we have 231 variants in 
# this table. We only have 10 strains in the mock community, so a lot of the rows in 
# this object are going to have a "0", indicating that a specific variant (i.e. row) 
# was not identified in this sample.
#
# We will now use Dada2 to evaluate the accuracy for classifying these known sequences 
# from our mock community as an index of the inherent error rate in our methodology.
# 
# Extract the "Mock" group from the ASV table, save it as a new object called 
# unique_mock, and only keep ASVs with one or more unique sequences (i.e. get rid of 
# the zeros)
#
# First, let's extract the 'mock' sample from our sequencing table

unique_mock <- seq_table_nochimeri["Zym",]

unique_mock <- seq_table_nochimeri_test["Zym",]

# Now let's only keep the variants with 1 or more unique sequences identified in the 
# mock community table.

unique_mock <- sort (unique_mock[unique_mock > 0], decreasing = TRUE)

# Okay, let's see how many variants were identified in our Mock community

cat("DADA2 inferred", length(unique_mock), "sample sequences present in the Mock community.\n")

# As you can see, 8 variants were detected , good.
# Previously, there were 13.
#
# However, just because we had the same number of variants as the number of bugs 
# in our community doesn't mean that they were classified correctly, right? 
# So let's look to see if the taxonomy matches
#
# Let's pull in the sequences of the known strains from our 'mock community'. 
# This is a reference FASTA file called mock_genome.fasta, which is found in the 
# data directory.
#
# Hmm.... The mock_genome.fasta file looks like a good place to figure out what 
# species should be in the mock community.
#
# Now load this FASTA file, be sure it is located in the "Raw" directory
# File path ~/16S/Raw/

reference_mock <- getSequences(file.path(path, "mock_genome_zymo.fasta"))

# And use the grepl function to match the pattern of our variants to the known 
# sequences in our FASTA file

match.ref <- sum(sapply(names(unique_mock), function(x) any(grepl(x, reference_mock))))

which_match <- sapply(names(unique_mock), function(x) any(grepl(x, reference_mock)))

# If you print out the match.ref you will see that there were 8 matches. 
# Previously, there were 9.
# Oh no! My world has fallen apart. 

match.ref

# Now we're going to see how our data stcks up the the expected mock community

library(readxl)
abundance_table <- read_excel("Mock_abun_R.xlsx")
View (abundance_table)
names (abundance_table)

# Now we will plot this data

abundance_table <- as.data.frame(abundance_table)

# That didn't work.
# Let's enter the values manually...maybe that will work better
Calculated <- c(20.39, 17.69, 16.14, 11.22, 9.61, 9.49, 6.35, 9.11)
Expected <- c(17.4, 14.1, 15.5, 9.9, 10.1, 10.4, 4.2, 18.4)
Bacteria <- c("B_subtilis", "L_monocytogenes", "S_aureus", "E_faecalis",
              "E_coli", "S_enterica", "P_aeruginosa", "L_fermentum")
df1 <- tibble(Calculated, Expected, Bacteria)
Plot_table <- melt(df1, id.vars='Bacteria', measure.vars = c("Calculated", "Expected"))
ggplot(data = Plot_table, aes(x = variable, y = value, fill = Bacteria,
                              label = value)) +
  geom_col () +
  geom_text (size = 3, position = position_stack (vjust = 0.5)) +
  labs (x = "Pool", y = "Percent of total composition (%)") +
  theme_classic () + 
  theme (axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  theme (axis.title.x = element_blank())
ggsave('expected_composition.png', antialias="none")

###################################
###### Second try after blasting
###################################

Calculated <- c(19.5, 16.4, 17.3, 10.5, 10.9, 9.0, 5.9, 10.5)
Expected <- c(17.4, 14.1, 15.5, 9.9, 10.1, 10.4, 4.2, 18.4)
Bacteria <- c("B_subtilis", "L_monocytogenes", "S_aureus", "E_faecalis",
              "E_coli", "S_enterica", "P_aeruginosa", "L_fermentum")
df1 <- tibble(Calculated, Expected, Bacteria)
Plot_table <- melt(df1, id.vars='Bacteria', measure.vars = c("Calculated", "Expected"))
ggplot(data = Plot_table, aes(x = variable, y = value, fill = Bacteria,
                              label = value)) +
  geom_col () +
  geom_text (size = 3, position = position_stack (vjust = 0.5)) +
  labs (x = "Pool", y = "Percent of total composition (%)") +
  theme_classic () + 
  theme (axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  theme (axis.title.x = element_blank())

###################################
###### Third try end set only #####
###################################

# Using a cutoff

Calculated <- c(20.5, 17.7, 16.2, 11.3, 9.7, 9.5, 5.9, 9.2)
Expected <- c(17.4, 14.1, 15.5, 9.9, 10.1, 10.4, 4.2, 18.4)
Bacteria <- c("B_subtilis", "L_monocytogenes", "S_aureus", "E_faecalis",
              "E_coli", "S_enterica", "P_aeruginosa", "L_fermentum")
df1 <- tibble(Calculated, Expected, Bacteria)
Plot_table <- melt(df1, id.vars='Bacteria', measure.vars = c("Calculated", "Expected"))
ggplot(data = Plot_table, aes(x = variable, y = value, fill = Bacteria,
                              label = value)) +
  geom_col () +
  geom_text (size = 3, position = position_stack (vjust = 0.5)) +
  labs (x = "Pool", y = "Percent of total composition (%)") +
  theme_classic () + 
  theme (axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  theme (axis.title.x = element_blank())

# Not using a cutoff

Calculated <- c(20.7, 16.2, 17.4, 10.2, 11.2, 10.7, 5.4, 8.3)
Expected <- c(17.4, 14.1, 15.5, 9.9, 10.1, 10.4, 4.2, 18.4)
Bacteria <- c("B_subtilis", "L_monocytogenes", "S_aureus", "E_faecalis",
              "E_coli", "S_enterica", "P_aeruginosa", "L_fermentum")
df1 <- tibble(Calculated, Expected, Bacteria)
Plot_table <- melt(df1, id.vars='Bacteria', measure.vars = c("Calculated", "Expected"))
ggplot(data = Plot_table, aes(x = variable, y = value, fill = Bacteria,
                              label = value)) +
  geom_col () +
  geom_text (size = 3, position = position_stack (vjust = 0.5)) +
  labs (x = "Pool", y = "Percent of total composition (%)") +
  theme_classic () + 
  theme (axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  theme (axis.title.x = element_blank())

########################
# Handoff to phyloseq  #
########################