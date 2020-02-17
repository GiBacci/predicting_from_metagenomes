# Getting data from CF longitudinal survey.
# In principal, this script must be launched one time
# and then can be removed
library(foreach)
library(tidyverse)
source("~/Dropbox/Scripts/R/script.R")

# Importing data
merg.table <- "~/Dropbox/CF_longitudinal_metagenome/Metaphlan Analysis/merged_abundance_table.txt"
phylo <- metaphlan_to_phyloseq(merg.table)

# Reading and adjusting metadata
meta <- read.table("~/Dropbox/CF_longitudinal_metagenome/data/tabella_compl_FINAL.tsv", sep = "\t",
                   header = T, row.names = sample_names(phylo),
                   stringsAsFactors = F)

meta$genotype <- NA
meta$genotype[grep("F508", meta$CFTR)] <- "heterozygote_F508"
meta$genotype[meta$CFTR == "F508/F508"] <- "homozygote_F508"
meta$genotype[is.na(meta$genotype)] <- "other"

meta$Exacerbation <- gsub("ex", "exacerbation", meta$Exacerbation)
meta$Exacerbation <- gsub("pe", "post_exacerbation", meta$Exacerbation)

for(i in seq_along(meta)){
  if(is.character(meta[[i]]))
    meta[[i]] <- factor(gsub("^\\s|\\s$", "", meta[[i]]))
}

meta$Time <- factor(meta$Time, levels = c(0, 1, 2, 3), labels = paste0("t", 0:3))

sample_data(phylo) <- meta

# Building a phyloseq object from metaphlan abundances
# Removing M34 and M24 time 1 since they do not have
# enough assignments
phylo <- subset_samples(phylo, Sample != "M34")
phylo <- subset_samples(phylo, Name.File != "CF_ABM24SS_t1M16")

sample.names <- str_split(sample_names(phylo), "_|-") %>%
  map_chr(~paste(.[2:3], collapse = "_"))
sample_names(phylo) <- sample.names

# Removing G24 t1 since it has too few reads
phylo <- prune_samples(sample_sums(phylo) > 0, phylo)

# Building matrix ad saving
taxa.ab <- otu_table(phylo)
class(taxa.ab) <- "matrix"
attr(taxa.ab, "taxa_are_rows") <- NULL

taxa.meta <- sample_data(phylo)
class(taxa.meta) <- "data.frame"
saveRDS(taxa.meta, "./data/sample_meta.rds")

rownames(taxa.ab) <- tax_table(phylo) %>% as_tibble() %>% pull(s)
taxa.ab <- taxa.ab[!is.na(rownames(taxa.ab)),]
taxa.ab <- t(prop.table(taxa.ab, 2))
saveRDS(taxa.ab, "./data/taxa_ab.rds")

# Listing directories and extracting only those with gene calling
dir_ass <- list.dirs(path = "~/Dropbox/CF_longitudinal_metagenome/AR_heatmap/plotting_script", 
                     recursive = F, full.names = T)
dir_name <- sapply(dir_ass, file.path, "gene_calling")

# Reading rgi formatted files that have been produced by using 
# the script "rgi_analysis.R"
rgi.files <- sapply(dir_name, file.path, "resistance_gene_mod.txt")
rgi.files <- rgi.files[sapply(rgi.files, file.exists)] # keep only existing files

# Parsing rgi files
res <- foreach(r = rgi.files, .combine = rbind) %dopar% {
  # Reading data
  rgi <- read.table(r, sep = "\t", header = T, fill = T, 
                    quote = "", stringsAsFactors = F)
  
  # If no rows are present return NULL
  if(nrow(rgi) == 0) return(NULL)
  
  # Contig length
  rgi$lenght <- as.numeric(gsub("^.*_length_([0-9]+)_.*$", "\\1", 
                                rgi$node))
  
  # Format semicolon separated value
  format.listed <- function(x){
    sorted.class <- sapply(strsplit(x, "\\s?;\\s?"), sort)
    sapply(sorted.class, paste, collapse = "; ")
  }
  
  # Building a humann2-style table
  sample.rgx <- "CF_(AB[A-Z][0-9]{2}SS_t(?:PE)?[0-9]M[0-9]{2})_.*"
  by.taxa <- rgi %>%
    mutate(Drug.Class = format.listed(Drug.Class),
           AMR.Gene.Family = format.listed(AMR.Gene.Family),
           Resistance.Mechanism = format.listed(Resistance.Mechanism)) %>%
    group_by(Best_Hit_ARO, 
             AMR.Gene.Family, 
             Resistance.Mechanism, 
             Drug.Class,
             taxa) %>%
    summarize(cov = sum(nreads)) %>%
    bind_rows(summarize(., cov = sum(cov)) %>% 
                mutate(taxa = "overall") %>% 
                select(names(.))) %>%
    arrange(Best_Hit_ARO) %>%
    ungroup() %>%
    mutate(sample = str_match(r, sample.rgx)[,2])
}
# building final data frame
res <- res[res$sample %in% rownames(taxa.ab),] %>%
  mutate(sample = factor(sample, rownames(taxa.ab)))
complete <- res %>% filter(cov > 0) %>%
  spread(sample, cov, fill = 0)


# getting only overall counts
ar.counts <- complete %>% 
  filter(taxa == "overall") %>%
  select(starts_with("AB")) %>%
  as.matrix() %>%
  `rownames<-` (unique(complete$Best_Hit_ARO))
ar.counts <- ar.counts[,match(rownames(taxa.ab), colnames(ar.counts))]
colnames(ar.counts) <- rownames(taxa.ab)
ar.counts <- t(ar.counts)
ar.counts[is.na(ar.counts)] <- 0
# save count matrix
saveRDS(ar.counts, "./data/ar_counts.rds")

# getting taxa relative abundance
# along with metadata
ar.meta <- complete %>%
  filter(taxa != "overall") %>%
  gather("sample", "value", starts_with("AB")) %>%
  group_by_at(vars(-sample, -value)) %>%
  summarize(value = sum(value)) %>%
  mutate(value = value/sum(value)) %>%
  spread(taxa, value, fill = 0) %>%
  ungroup()
# Saving data frame
saveRDS(ar.meta, "./data/ar_meta.rds")


# Gene counts
gene.counts <- readRDS("/home/giovannib/Dropbox/CF_longitudinal_metagenome/gene_assignments/counts.rds")
colnames(gene.counts) <- str_remove(colnames(gene.counts), "CF_")
r.name <- gene.counts$best.og
gene.counts <- as.matrix(gene.counts[,-1])
rownames(gene.counts) <- r.name
gene.counts <- t(gene.counts)
gene.counts <- gene.counts[rownames(taxa.ab),]
saveRDS(gene.counts, "./data/gene_counts.rds")

gene.meta <- readRDS("/home/giovannib/Dropbox/CF_longitudinal_metagenome/gene_assignments/annotations.rds")
gene.meta <- as.data.frame(gene.meta)
gene.meta <- gene.meta[match(colnames(gene.counts), gene.meta$best.og),]
saveRDS(gene.meta, "./data/gene_meta.rds")
