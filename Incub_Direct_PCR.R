library(dada2)
library(readxl)
library(pdftools)
library(phyloseq)
library(tidyverse)
library(stats)
library(phangorn)
library(DECIPHER)
library(microeco)
library(file2meco)

path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/raw/"
list.files(path)


fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/fQualityPlotFilt.pdf") 

QualityBeforeFilt <- plotQualityProfile(
  fnFs[1:18],    # Input forward read file(s) (file name or character vector of file names)
  aggregate = F  # Logical, whether to plot aggregated quality profile (TRUE) or individual samples (FALSE)
)

QualityBeforeFilt
dev.off() 


#Your reads must still overlap after truncation in order to merge them later! The tutorial is using 2x250 V4 sequence data, so the forward and reverse reads almost completely overlap and our trimming can be completely guided by the quality scores. If you are using a less-overlapping primer set, like V1-V2 or V3-V4, your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.

#Filter and trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

names(filtFs) <- sample.names


any(duplicated(fnFs))
any(duplicated(filtFs))

out <- filterAndTrim(
  fnFs,          # Input forward read file(s) (file name or character vector of file names)
  filtFs,        # Input filtered file(s) (file name or character vector of file names)
  truncLen = 240,  # Truncate reads to this length (bases)
  maxN = 0,        # Maximum number of expected N (ambiguous) bases in a read
  maxEE = 2,       # Maximum expected error rate (higher values are more permissive)
  truncQ=2,
  rm.phix=TRUE,
  multithread = FALSE,  # Logical, whether to use multiple threads for computation
  compress = TRUE  # Logical, whether to compress output files
)

head(out)


pdf("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/fQualityPlotAfterFilterAndTrim.pdf") 
QualityAfterFilt <- plotQualityProfile(filtFs[1:18], aggregate=F)
QualityAfterFilt
dev.off() 

errF <- learnErrors(
  filtFs,                           # Input filtered file(s) (file name or character vector of file names)
  nbases = 1e+09,                   # Number of bases to use for error learning (default: 1 billion)
  errorEstimationFunction = loessErrfun,  # Error estimation function to use
  multithread = TRUE,               # Logical, whether to use multiple threads for computation
  randomize = TRUE                 # Logical, whether to randomize the order of reads before learning errors
)


plotErrors(errF)

#Error rates
pdf("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/ErrorPlot.pdf") 
# 2. Create a plot
Errors <- plotErrors(errF, nominalQ=TRUE)
Errors
# Close the pdf file
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)

names(derepFs) <- sample.names

dadaFs <- dada(
  derepFs,                    # Input filtered and derelicated file(s) (file name or character vector of file names)
  err = errF,                 # Pre-learned error rates (object obtained from learnErrors function)
  selfConsist = FALSE,        # Logical, whether to perform self-consistency iteration for error learning (not necessary when err is provided)
  multithread = TRUE,         # Logical, whether to use multiple threads for computation
  pool = TRUE                # Logical, whether to pool reads during denoising
)

dadaFs[[1]]
#There is much more to the dada-class return object than this (see help("dada-class") for some info), including multiple diagnostics about the quality of each denoised sequence variant, but that is beyond the scope of an introductory tutorial.

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread=TRUE)


dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)



#Check processing of the sequences
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("filtered", "denoisedF", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.table(track, file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/SeqTracking.txt", sep=" ")

taxa <- assignTaxonomy(seqtab.nochim, "D:/OneDrive - Universität Bayreuth/Work/Database/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE, 
                       tryRC=TRUE,
                       minBoot=80)
taxa <- addSpecies(taxa, "D:/OneDrive - Universität Bayreuth/Work/Database/silva_species_assignment_v138.1.fa.gz")


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
head(taxa)

saveRDS(seqtab.nochim, "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/seqtab.rds")
saveRDS(taxa, "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/taxa.rds")


seqtab.nochim=readRDS("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/seqtab.rds")
taxa=readRDS("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/taxa.rds")

##...Construct phylogenetic tree..............##########
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), 
                       anchor=NA,
                       verbose=FALSE)

The phangorn R package is then used to construct a phylogenetic tree. 
Here we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

##Save tree file
tree_file = phy_tree(fitGTR$tree)
ape::write.tree(tree_file, "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/tree.tree")

# Read your phylogenetic tree (replace "tree_file" with your actual file path)
tree <- ape::read.tree("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/tree.tree")

# Read your taxonomy data (replace "taxonomy_file" with your actual file path)
taxonomy_data =as.data.frame(taxa)
taxonomy_data$Taxonomy=paste(taxonomy_data$Kingdom,
                             taxonomy_data$Phylum,
                             taxonomy_data$Class,
                             taxonomy_data$Order,
                             taxonomy_data$Family,
                             taxonomy_data$Genus,
                             taxonomy_data$Species,sep = "; ")

# Extract ASV IDs and corresponding taxonomic assignments
asv_ids <- row.names(taxonomy_data)
tax_assignments <- taxonomy_data$Taxonomy

# Create a mapping between ASV IDs and taxonomic assignments
asv_tax_map <- setNames(tax_assignments, asv_ids)

# Function to replace tip names in the tree
replace_tip_names <- function(tree, name_map) {
  for (i in 1:length(tree$tip.label)) {
    if (tree$tip.label[i] %in% names(name_map)) {
      tree$tip.label[i] <- name_map[[tree$tip.label[i]]]
    }
  }
  return(tree)
}

# Replace tip names in the tree
new_tree <- replace_tip_names(tree_file, asv_tax_map)

# Save the modified tree (replace "output_tree_file" with desired file path)
ape::write.tree(new_tree, "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/tree_with_taxonomy.tree")

tree_file=ape::read.tree("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/processed/tree.tree")

##.......Combine data into a phyloseq object.....########
samdf <- read.table("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/Mapping.txt",header = T)%>%
  mutate(Treatment=recode(Treatment, "Initial"="Cells"))
rownames(seqtab.nochim) <- samdf$sample_id
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE  

rownames(samdf) <- samdf$sample_id

#The full suite of data for this study - the sample-by-sequence feature table, the sample metadata, the sequence taxonomies, and the phylogenetic tree - can now be combined into a single object.

df_phylo <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa),
                   phy_tree(tree_file))

##rename rownames
dna <- Biostrings::DNAStringSet(taxa_names(df_phylo))
names(dna) <- taxa_names(df_phylo)
df_phylo_merge <- merge_phyloseq(df_phylo, dna)
taxa_names(df_phylo_merge) <- paste0("ASV", seq(ntaxa(df_phylo_merge)))

##write phyloseq
saveRDS(df_phylo_merge, file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Phyloseq_Object.rds")

##.......Read Phyloseq..########
df_phylo_merge=readRDS("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Phyloseq_Object.rds")


##Absolute
df_otu=data.frame(t(otu_table(df_phylo_merge)))
colnames(df_otu)=samdf$sample.id
df_otu$asv_id=rownames(df_otu)

df_taxa=data.frame(tax_table(df_phylo_merge))
df_taxa$asv_id=rownames(df_taxa)

df_abs=inner_join(df_otu,df_taxa, by="asv_id")

rownames(df_abs)=df_abs$asv_id

write.csv(df_abs,file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Otu_table_absolute.csv")

##relative
df_rel = transform_sample_counts(df_phylo_merge, function(x) x/sum(x))

df_otu=data.frame(t(otu_table(df_rel)))
colnames(df_otu)=samdf$sample.id
df_otu$asv_id=rownames(df_otu)

df_taxa=data.frame(tax_table(df_rel))
df_taxa$asv_id=rownames(df_taxa)

df_rel=inner_join(df_otu,df_taxa, by="asv_id")

rownames(df_rel)=df_rel$asv_id

write.csv(df_rel,file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Otu_table_relative.csv")


##convert phyloseq to microeco object
df_phylo <- prune_samples(sample_names(df_phylo_merge) != "sa19", df_phylo)

##....MicroEco package.....#####
dataset=phyloseq2meco(df_phylo_merge)

dataset$tidy_dataset()

dataset$sample_sums() %>% range

dataset$rarefy_samples(sample.size = 7000)
dataset$sample_sums() %>% range
dataset$cal_abund()

dir.create("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/taxa_abund")
dataset$save_abund(dirpath = "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/taxa_abund")

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = T)

# save dataset$alpha_diversity to a directory
dir.create("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/alpha_diversity")
dataset$save_alphadiv(dirpath = "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/alpha_diversity")

# If you do not want to calculate unifrac metrics, use unifrac = FALSE
dataset$cal_betadiv(unifrac = T)

dir.create("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/beta_diversity")
dataset$save_betadiv(dirpath = "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/beta_diversity")


##..Incubation...####
dataset$sample_table$Treatment=factor(dataset$sample_table$Treatment,
                                      level=c("Soil","Initial","After_1_day","Control","Glucose","Glutamine"))
# bar plot
t1 <- trans_abund$new(dataset = dataset, 
                      taxrank = "Species", 
                      ntaxa = 20, 
                      groupmean = "Treatment")

t1$plot_bar(others_color = "grey70", 
            legend_text_italic = FALSE,
            xtext_size = 20)+
  scale_x_discrete(limits=c("Initial","After_1_day","Control","Glucose","Glutamine"),
                   labels=c("After_1_day"="Before Incubation"))+
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_Bar.jpeg",height = 8,width = 6)

# alluvial plot
t1 <- trans_abund$new(dataset = dataset, 
                      taxrank = "Class", 
                      ntaxa = 10,
                      groupmean = "Treatment")

t1$plot_bar(use_alluvium = T, 
            clustering = T, 
            others_color = "black",
            xtext_keep = F, 
            xtext_size = 20)+
  labs(y="Relative proportion (%)")+
  
  theme_classic()+
  
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_alluvial.jpeg",height = 10,width = 4)



##..Soil and Cells only...####
dataset$sample_table <- subset(dataset$sample_table, 
                               Treatment %in%c("Soil","Cells"))
dataset$sample_table$Treatment=factor(dataset$sample_table$Treatment,
                                      level=c("Soil","Cells"))

# bar plot
t1 <- trans_abund$new(dataset = dataset, 
                      taxrank = "Class", 
                      ntaxa = 10, 
                      groupmean = "Treatment")
t1$plot_bar(others_color = "grey70", 
            legend_text_italic = FALSE,
            xtext_size = 20)+
  scale_x_discrete(limits=c("Soil","Cells"))+
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_Bar.jpeg",
       height = 8,
       width = 6)

# alluvial plot
t1 <- trans_abund$new(dataset = dataset, 
                      taxrank = "Class", 
                      ntaxa = 10,
                      groupmean = "Treatment")

t1$plot_bar(use_alluvium = T, 
            clustering = T, 
            others_color = "black",
            xtext_keep = F, 
            xtext_size = 20)+
  labs(y="Relative proportion (%)")+
  
  theme_classic()+
  
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_alluvial.jpeg",height = 10,width = 4)

# venn diagram

# merge samples as one community for each group
dataset_venn=clone(dataset)
dataset_venn$sample_table=subset(dataset_venn$sample_table, 
                             Treatment %in%c("Soil","Cells"))
dataset_venn$tidy_dataset()

dataset_venn <- dataset_venn$merge_samples(use_group = "Treatment")

venn_plot <- trans_venn$new(dataset_venn, ratio = NULL)
venn_plot$plot_venn(color_circle = c("black","blue"),alpha = 0.5)

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Venn_Soil_Vs_Cells.jpeg",
       height = 6, width = 8)

# unique and shared genus

dataset_uniq <- trans_venn$new(dataset_venn)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
dataset_feq <- dataset_uniq$trans_comm(use_frequency = T)

# calculate taxa abundance, that is, the frequency
dataset_feq$cal_abund()

# transform and plot
dataset_abu <- trans_abund$new(dataset = dataset_feq, taxrank = "Genus", ntaxa = 10)
dataset_abu$plot_bar(bar_type = "part", 
                     legend_text_italic = F, 
                     color_values = as.vector(pals::polychrome(10)))+ 
  ylab("Frequency (%)")



# Ordination
t1 <- trans_beta$new(dataset = dataset, group = "Treatment", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
t1$plot_ordination(plot_color = "Treatment", plot_type = c("point"))


##......Tidyverse sequence analysed in R.........##########
library(tidyverse)
library(ggtext)
library(pals)
library(readxl)
library(ggalluvial)
library(stringr)
library(cowplot)
library(glue)

metadata <- read_excel("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Phyloseq_Table.xlsx",sheet = "Mapping")


asv_counts <- read_excel("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Phyloseq_Table.xlsx",sheet = "ASVs") %>%
  select(sample_id, starts_with("ASV")) %>%
  pivot_longer(-sample_id, names_to="asv", values_to = "count")


taxonomy <- read_excel("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Phyloseq_Table.xlsx", sheet = "Taxonomy") %>%
  rename_all(tolower) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")

taxonomy_sel=taxonomy%>%
  select(phylum,class,order,family,genus)


asv_rel_abund <- inner_join(metadata, asv_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="asv") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus","species", "asv"),
               names_to="level",
               values_to="taxon") 

##Phylum
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("After_1_day","Control","Glucose","Glutamine"))%>%
  filter(level=="phylum") %>%
  filter(taxon!="unidentified")%>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") 



taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_phylum=inner_join(taxon_rel_abund,taxon_pool, by="taxon")

df_phylum%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))%>%
  mutate(Treatment=factor(Treatment, levels = rev(c("After_1_day","Control","Glucose","Glutamine"))),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon, 
             alluvium = taxon,
             label = taxon,
             fill = taxon)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Phylum",values = c(as.vector(polychrome(6)),"gray"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits=rev)+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(legend.text = element_text(family = "Times",
                                   size=12,
                                   color="black"),
        legend.title = element_text(face="bold",
                                    family = "Times",
                                    size=16,
                                    color="black"),
        legend.key.size = unit(10, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times",
                                   size=12,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times",
                                 size=16,
                                 color="black"),
        strip.text = element_text(face="bold",
                                  family = "Times",
                                  size=16,
                                  color="black")
  )+theme(panel.spacing.x = unit(1, "lines"))

ggsave("epth_Phylum_Alluvial.png",
       height = 8,width = 12,dpi=300)


##genus
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("After_1_day","Control","Glucose","Glutamine"))%>%
  filter(level=="genus") %>%
  #filter(taxon!="unidentified")%>%
  #filter(taxon!="Pseudomonas")%>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_genus=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$genus),])

df_genus=data.frame(df_genus, stringsAsFactors = F)


df_genus%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.1%)","Other (< 0.1%)",
                            if_else(taxon=="*unidentified*","unidentified",
                            paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2=factor(taxon2, levels = unique(taxon2)))%>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Genus",values = c(as.vector(polychrome(24)),"gray"))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times",
                                       size=14,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times",
                                        size=16,
                                        color="black"),
        legend.key.size = unit(10, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(".png",
       height = 8,width = 12,dpi=300)

##......Pseudomonas separated (sequence analysed in R).......####

# Define a function to create the small plot for Pseudomonas
pseudomonas_small_plot <- df_genus %>%
  filter(taxon == "Pseudomonas") %>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment, y = mean_rel_abund, fill = taxon)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = NULL, y = NULL,title = "Pseudomonas") +
  scale_fill_manual(values = "gray") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# Create the main alluvial plot
main_plot <- df_genus %>%
  filter(taxon != "Pseudomonas") %>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon == "Other (< 0.1%)", "Other (< 0.1%)",
                            if_else(taxon=="*unidentified*","unidentified",
                            glue("{taxon} **(<span style = 'font-size:10pt'>{phylum}</span>)**")))) %>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2=factor(taxon2, levels = unique(taxon2)))%>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment,
             y = mean_rel_abund,
             stratum = taxon2,
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "black") +
  geom_stratum() +
  scale_fill_manual(name = "Genus", values = c(as.vector(polychrome(24)), "gray")) +
  #scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(10, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times",
                                 size=16,
                                 color="black"),
        strip.text = element_text(face="bold",
                                  family = "Times",
                                  size=16,
                                  color="black"),
        strip.background.y = element_rect(fill="red"),
        strip.background.x = element_rect(fill="gray"))+
  theme(panel.spacing.x = unit(3, "lines"))

# Combine the main plot and the small Pseudomonas plot
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(pseudomonas_small_plot, x = 0.25, y = 0.8, width = 0.2, height = 0.2)

# Print the final combined plot
print(final_plot)

ggsave("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation.pdf",device = 'pdf',
       width = 16, height = 8, dpi = 300)
       









##......Tidyverse sequence analysed in Qiime2.........##########
library(tidyverse)
library(ggtext)
library(pals)
library(readxl)
library(ggalluvial)
library(stringr)
library(cowplot)
library(glue)
library(microeco)
library(file2meco)
library(magrittr)

abund_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/dada2OutPut/dada2_table.qza"
sample_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/Mapping.txt"
taxonomy_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/taxClassification/dada2_rep_set_classified.qza"
tree_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/aligned2Tree/rooted_tree.qza"
rep_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/dada2OutPut/dada2_rep_set.qza"

df_qiime <- qiime2meco(abund_file_path, 
                       sample_table = sample_file_path, 
                       taxonomy_table = taxonomy_file_path, 
                       phylo_tree = tree_data, 
                       rep_fasta = rep_data, 
                       auto_tidy = TRUE)

df_qiime$tax_table %<>% tidy_taxonomy
tmp <- df_qiime$tax_table
# search each row to add sth
for(i in 1:nrow(tmp)){
  if(any(grepl("__$", tmp[i, ]))){
    matchall <- grep("__$", tmp[i, ])
    tmp[i, matchall] <- paste0(tmp[i, matchall], "unclassified ", gsub(".__", "", tmp[i, matchall[1] - 1]))
  }
}
# then reassign tmp to the raw table
df_qiime$tax_table <- tmp

df_qiime$tidy_dataset()


metadata=data.frame(df_qiime$sample_table)
metadata$sample_id=row.names(metadata)

asv_counts <- t(df_qiime$otu_table)
sample_ids <- rownames(df_qiime$otu_table)
asv_info <- data.frame(ASV_ID = colnames(asv_counts), stringsAsFactors = FALSE)
asv_info$ASV_Number <- paste("ASV", seq_len(ncol(asv_counts)), sep = "")
asv_counts <- asv_counts[, match(asv_info$ASV_ID, colnames(asv_counts))]
asv_counts <- as.data.frame(asv_counts)
colnames(asv_counts) <- asv_info$ASV_Number
asv_counts$sample_id <- row.names(metadata)

taxonomy=data.frame(df_qiime$tax_table)

row.names(taxonomy)=asv_info$ASV_Number
taxonomy$asv=asv_info$ASV_Number
colnames(taxonomy)=c("kingdom", "phylum", "class", "order", "family", "genus","species", "asv")

taxonomy=taxonomy%>%
  pivot_longer(-asv, 
               names_to="taxa", 
               values_to = "taxon")%>%
  mutate(taxon=str_replace(taxon, "[a-z]__", ""),
         taxon=str_replace_all(taxon, "unclassified.*", "Unclassified"))%>%
  mutate(taxon=str_replace(taxon, "Proteobacteria", "Pseudomonadota"))%>%
  mutate(taxon=str_replace(taxon, "Firmicutes", "Bacillota"))%>%
  pivot_wider(names_from = "taxa",values_from = "taxon")

taxonomy_sel=taxonomy%>%
  select(phylum,class,order,family,genus)

asv_counts <- asv_counts%>%
  select(sample_id, starts_with("ASV")) %>%
  pivot_longer(-sample_id, names_to="asv", values_to = "count")
  


asv_rel_abund <- inner_join(metadata, asv_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="asv") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus","species", "asv"),
               names_to="level",
               values_to="taxon") 

##.......Soil Vs Cell at class level.......#######
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("Soil","Initial"))%>%
  filter(level=="class") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_class=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$class),])

df_class=data.frame(df_class, stringsAsFactors = F)


df_class%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.1%)","Other (< 0.1%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=recode(Treatment, "Initial"="SFCE"),
         Treatment=factor(Treatment, levels = c("Soil", "SFCE")))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  scale_fill_manual(name="Class (Phylum)",values = c(as.vector(polychrome(35)),"gray"))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=40),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=40,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=40,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_Class.tiff",
       height = 12, width = 14, dpi=300)

##.......Soil Vs Cell at genus level.......#######
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("Soil","Initial"))%>%
  filter(level=="genus") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.5, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.5%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_genus=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$genus),])

df_genus=data.frame(df_genus, stringsAsFactors = F)


df_genus%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.5%)","Other (< 0.5%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=recode(Treatment, "Initial"="SFCE"),
         Treatment=factor(Treatment, levels = c("Soil", "SFCE")))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  scale_fill_manual(name="Genus",values = c(as.vector(polychrome(32)),"gray"))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_Genus.tiff",
       height = 12, width = 10, dpi=300)


##......Incubation with initial samples...#####

##class
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("Initial","After_1_day","Control","Glucose","Glutamine"))%>%
  filter(level=="class") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_class=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$class),])

df_class=data.frame(df_class, stringsAsFactors = F)


df_class%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.1%)","Other (< 0.1%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=factor(Treatment, levels = c("Initial","After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Class (Phylum)",values = c(as.vector(polychrome(21)),"gray"))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(10, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_Class.tiff",
       height = 12, width = 14, dpi=300)

##genus
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("Initial","After_1_day","Control","Glucose","Glutamine"))%>%
  filter(level=="genus") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.5, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.5%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_genus=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$genus),])

df_genus=data.frame(df_genus, stringsAsFactors = F)


df_genus%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.5%)","Other (< 0.5%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=factor(Treatment, levels = c("Initial","After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Genus",values = c(as.vector(polychrome(17)),"gray"))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(10, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_Genus.tiff",
       height = 12, width = 14, dpi=300)


# Define a function to create the small plot for Pseudomonas
pseudomonas_small_plot <- df_genus %>%
  filter(taxon == "Pseudomonas") %>%
  mutate(Treatment=factor(Treatment, levels = c("Initial","After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment, y = mean_rel_abund, fill = taxon)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = NULL, y = NULL,title = "Pseudomonas") +
  scale_fill_manual(values = "gray") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# Create the main alluvial plot
main_plot <- df_genus %>%
  filter(taxon != "Pseudomonas") %>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon, "^([^<]*)$", "*\\1*")) %>%
  mutate(`taxon2` = if_else(taxon == "Other (< 0.5%)", "Other (< 0.5%)",
                            if_else(taxon == "*unidentified*", "unidentified",
                                    glue("{taxon} **<span style = 'font-size:16pt'> ({phylum}) </span>**")))) %>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2=factor(taxon2, levels = unique(taxon2)))%>%
  mutate(Treatment=factor(Treatment, levels = c("Initial","After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment,
             y = mean_rel_abund,
             stratum = taxon2,
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "black") +
  geom_stratum() +
  scale_fill_manual(name = "Genus", values = c(as.vector(polychrome(16)), "gray")) +
  #scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  guides(fill = guide_legend(ncol = 1)) +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(1, "lines"))

# Combine the main plot and the small Pseudomonas plot

main_plot+patchwork::inset_element(pseudomonas_small_plot,
                                   left = 0.5, 
                                   bottom = 0.4, 
                                   right = 0.8, 
                                   top = 0.95)


ggsave("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_Pseudo.jpeg",
          width = 18, height = 8, dpi = 300)



##......Incubation Only...#####

##class
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("After_1_day","Glucose","Glutamine"))%>%
  filter(level=="class") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_class=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$class),])

df_class=data.frame(df_class, stringsAsFactors = F)


df_class%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.1%)","Other (< 0.1%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Class",values = c(as.vector(polychrome(7)),"gray"))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_Only_Class.tiff",
       height = 12, width = 14, dpi=300)

##genus
taxon_rel_abund=asv_rel_abund %>%
  filter(Treatment%in%c("After_1_day","Glucose","Glutamine"))%>%
  filter(level=="genus") %>%
  group_by(Treatment, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Treatment, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool=taxon_rel_abund%>%
  group_by(taxon)%>%
  summarize(pool=max(mean_rel_abund)<0.1, .groups = "drop")

df_join=inner_join(taxon_rel_abund,taxon_pool, by="taxon")%>%
  mutate(taxon=if_else(pool, "Other (< 0.1%)",taxon))%>%
  group_by(Treatment,taxon)%>%
  summarize(mean_rel_abund=sum(mean_rel_abund))

df_genus=cbind(df_join,taxonomy_sel[match(df_join$taxon,
                                          taxonomy_sel$genus),])

df_genus=data.frame(df_genus, stringsAsFactors = F)


df_genus%>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"))%>%
  mutate(`taxon2` = if_else(taxon=="Other (< 0.1%)","Other (< 0.1%)",
                            if_else(taxon=="*Unclassified*","Unclassified",
                                    paste0(taxon, " **(", phylum, ")**"))))%>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2 = factor(taxon2, levels = unique(taxon2))) %>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x=Treatment, 
             y=mean_rel_abund, 
             stratum = taxon2, 
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "black")+
  geom_stratum() +
  scale_y_break(c(6,72), scales = 2)+
  #geom_bar(stat = "identity",position = "fill", color="black")+
  scale_fill_manual(name="Genus",values = c(as.vector(polychrome(22)),"gray"))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=20,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=24,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(3, "lines"))

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_only_Genus.tiff",
       height = 6, width = 18, dpi=300)


# Define a function to create the small plot for Pseudomonas
pseudomonas_small_plot <- df_genus %>%
  filter(taxon == "Pseudomonas") %>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment, y = mean_rel_abund, fill = taxon)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = NULL, y = NULL,title = "Pseudomonas") +
  scale_fill_manual(values = "gray") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# Create the main alluvial plot
main_plot <- df_genus %>%
  filter(taxon != "Pseudomonas") %>%
  filter(taxon!="Unclassified")%>%
  mutate(taxon = str_replace(taxon, "^([^<]*)$", "*\\1*")) %>%
  mutate(`taxon2` = if_else(taxon == "Other (< 0.1%)", "Other (< 0.1%)",
                            if_else(taxon == "*unidentified*", "unidentified",
                                    glue("{taxon} **<span style = 'font-size:16pt'> ({phylum}) </span>**")))) %>%
  arrange(phylum, taxon2)%>%
  mutate(taxon2=factor(taxon2, levels = unique(taxon2)))%>%
  mutate(Treatment=factor(Treatment, levels = c("After_1_day","Control","Glucose","Glutamine")),
         Treatment=recode(Treatment, "After_1_day"="Before Incubation"))%>%
  ggplot(aes(x = Treatment,
             y = mean_rel_abund,
             stratum = taxon2,
             alluvium = taxon2,
             label = taxon2,
             fill = taxon2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "black") +
  geom_stratum() +
  scale_fill_manual(name = "Genus", values = c(as.vector(polychrome(21)), "gray")) +
  #scale_y_continuous(limits = c(0, 25)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  guides(fill = guide_legend(ncol = 1)) +
  guides(fill=guide_legend(ncol=1))+
  theme_classic() +
  theme(plot.title = element_markdown(family = "Times New Roman",
                                      size=32,
                                      color="black",hjust = 0.5),
        legend.text = element_markdown(family = "Times New Roman",
                                       size=14,
                                       color="black"),
        legend.title = element_markdown(face="bold",
                                        family = "Times New Roman",
                                        size=16,
                                        color="black"),
        legend.key.size = unit(20, "pt"),
        axis.title = element_text(face="bold",
                                  family = "Times New Roman",
                                  size=16),
        axis.text.y = element_text(face="bold",
                                   family = "Times New Roman",
                                   size=16,color="black"),
        axis.text.x=element_text(face="bold",
                                 family = "Times New Roman",
                                 size=16,
                                 color="black"))+
  theme(panel.spacing.x = unit(1, "lines"))

# Combine the main plot and the small Pseudomonas plot
main_plot + patchwork::inset_element(
  pseudomonas_small_plot,
  left = 0.3,
  bottom = 0.5,
  right = 0.6,
  top = 1)

ggsave("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_only_Pseudo.jpeg",
       width = 18, height = 8, dpi = 300)


##....Alpha diversity...########
abund_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Substrate_4(Glu,GT,CA,Gly)_Chem_Incub/Seq_direct_PCR/Data/qiime/dada2OutPut/dada2_table.qza"
sample_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Substrate_4(Glu,GT,CA,Gly)_Chem_Incub/Seq_direct_PCR/Data/Mapping.txt"
taxonomy_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Substrate_4(Glu,GT,CA,Gly)_Chem_Incub/Seq_direct_PCR/Data/qiime/taxClassification/dada2_rep_set_classified.qza"
tree_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Substrate_4(Glu,GT,CA,Gly)_Chem_Incub/Seq_direct_PCR/Data/qiime/aligned2Tree/rooted_tree.qza"
rep_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Substrate_4(Glu,GT,CA,Gly)_Chem_Incub/Seq_direct_PCR/Data/qiime/dada2OutPut/dada2_rep_set.qza"

df_alpha <- qiime2meco(abund_file_path, 
                       sample_table = sample_file_path, 
                       taxonomy_table = taxonomy_file_path, 
                       phylo_tree = tree_data, 
                       auto_tidy = TRUE)

df_alpha$sample_table %<>% subset(Treatment %in% c("Soil","Initial","After_1_day","Control","Glucose","Glutamine"))
df_alpha$sample_table$Treatment %<>% factor(.,levels= c("Soil","Initial","After_1_day","Control","Glucose","Glutamine"))

df_alpha$sample_table %<>% subset(Treatment %in% c("Soil","Initial"))
df_alpha$sample_table$Treatment %<>% factor(.,levels= c("Soil","Initial"))

df_alpha$tidy_dataset()

df_alpha$cal_alphadiv(PD = T)

df_alpha=trans_alpha$new(dataset = df_alpha, group = "Treatment")

df_alpha$cal_diff(method = "anova")

chao1=df_alpha$plot_alpha(measure = "Chao1", y_increase = 0.1, add_sig_text_size = 6)+
    theme(axis.text.x = element_blank())+
    Axis_manupulation_cleanup

shannon=df_alpha$plot_alpha(measure = "Shannon", y_increase = 0.1, add_sig_text_size = 6)+
  theme(axis.text.x = element_blank())+
  Axis_manupulation_cleanup

simpson=df_alpha$plot_alpha(measure = "Simpson", y_increase = 0.1, add_sig_text_size = 6)+
  scale_x_discrete(limit=c("Soil", "Initial"),labels=c("Soil", "SFCE"))+
    Axis_manupulation_cleanup

faith=df_alpha$plot_alpha(measure = "PD", y_increase = 0.1, add_sig_text_size = 6)+
  scale_x_discrete(limit=c("Soil", "Initial"),labels=c("Soil", "SFCE"))+
  Axis_manupulation_cleanup

ggpubr::ggarrange(chao1, shannon, simpson, faith, nrow = 2, ncol = 2, align = "hv")

ggsave("PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Alpha_Soil_Vs_Cell.jpeg",
       height = 6, width = 6, dpi = 300)



lefse_diff <- trans_diff$new(dataset = df_qiime, method = "lefse", p_adjust_method = "none",group = "Treatment")
lefse_diff$res_diffplot_lefse_cladogram(use_taxa_num = 20, use_feature_num = 40, clade_label_level = 4,alpha = 0.2)



#################.......Data exported from Qiime2......#########
library(microeco)
library(file2meco)
library(tidyverse)
library(magrittr)

abund_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/dada2OutPut/dada2_table.qza"
sample_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/Mapping.txt"
taxonomy_file_path <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/taxClassification/dada2_rep_set_classified.qza"
tree_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/aligned2Tree/rooted_tree.qza"
rep_data <- "PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Data/qiime/dada2OutPut/dada2_rep_set.qza"

df_qiime <- qiime2meco(abund_file_path, 
                    sample_table = sample_file_path, 
                    taxonomy_table = taxonomy_file_path, 
                    phylo_tree = tree_data, 
                    rep_fasta = rep_data, 
                    auto_tidy = TRUE)

df_qiime$tax_table %<>% tidy_taxonomy
tmp <- df_qiime$tax_table
# search each row to add sth
for(i in 1:nrow(tmp)){
  if(any(grepl("__$", tmp[i, ]))){
    matchall <- grep("__$", tmp[i, ])
    tmp[i, matchall] <- paste0(tmp[i, matchall], "unclassified ", gsub(".__", "", tmp[i, matchall[1] - 1]))
  }
}
# then reassign tmp to the raw table
df_qiime$tax_table <- tmp

df_qiime$tidy_dataset()

df_qiime$sample_sums() %>% range

df_qiime$rarefy_samples(sample.size = 7200)
df_qiime$sample_sums() %>% range

df_qiime$cal_abund()

##..Incubation...####
df_qiime$sample_table$Treatment=factor(df_qiime$sample_table$Treatment,
                                      level=c("Soil","Initial","After_1_day","Control","Glucose","Glutamine"))
# bar plot
t1 <- trans_abund$new(dataset = df_qiime, 
                      taxrank = "Genus", 
                      ntaxa = 10, 
                      groupmean = "Treatment")

t1$plot_bar(others_color = "grey70", 
            legend_text_italic = FALSE,
            xtext_size = 20,
            xtext_angle = 45)+
  scale_x_discrete(limits=c("Initial","After_1_day","Control","Glucose","Glutamine"),
                   labels=c("After_1_day"="Before Incubation"))+
  Axis_manupulation_cleanup
ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Incubation_qiime2.jpeg",
       height=8, width = 10, dpi=300)

##..Soil and Cells only...####
df_qiime_sub=clone(df_qiime)

df_qiime_sub$sample_table %<>% subset(Treatment %in% c("Soil", "Initial"))
df_qiime_sub$tidy_dataset()

df_qiime_sub$sample_table$Treatment %<>% factor(., levels = c("Soil","Initial"), labels=c("Soil","Cells"))


# bar plot
t1 <- trans_abund$new(dataset = df_qiime, 
                      taxrank = "Genus", 
                      ntaxa = 10, 
                      groupmean = "Treatment")
t1$plot_bar(others_color = "grey70", 
            legend_text_italic = FALSE,
            xtext_size = 20)+
  scale_x_discrete(limits=c("Soil","Cells"))+
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_Bar.jpeg",
       height = 8,
       width = 6)

# alluvial plot
t1 <- trans_abund$new(dataset = df_qiime_sub, 
                      taxrank = "Class", 
                      ntaxa = 10,
                      groupmean = "Treatment")

t1$plot_bar(use_alluvium = T, 
            clustering = T, 
            others_color = "black",
            xtext_keep = F, 
            xtext_size = 20)+
  labs(y="Relative proportion (%)")+
  
  theme_classic()+
  
  Axis_manupulation_cleanup

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Soil_Vs_Cell_alluvial.jpeg",
       height = 8,
       width = 8)



# venn diagram

# merge samples as one community for each group
df_qiime_venn=clone(df_qiime_sub)
df_qiime_venn$tidy_dataset()

df_qiime_venn <- df_qiime_venn$merge_samples(use_group = "Treatment")

venn_plot <- trans_venn$new(df_qiime_venn, ratio = NULL)
venn_plot$plot_venn(color_circle = c("black","blue"),alpha = 0.5)

ggsave(file="PostDoc/Bayreuth/Data_analysis/WP1/Sub_2_Chemistry_Incubation/Seq/Analysis/Figure/Venn_Soil_Vs_Cells.jpeg",
       height = 6, width = 8)

# unique and shared genus

df_qiime_uniq <- trans_venn$new(df_qiime_venn)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
df_qiime_feq <- df_qiime_uniq$trans_comm(use_frequency = T)

# calculate taxa abundance, that is, the frequency
df_qiime_feq$cal_abund()

# transform and plot
df_qiime_abu <- trans_abund$new(df_qiime = df_qiime_feq, taxrank = "Genus", ntaxa = 10)
df_qiime_abu$plot_bar(bar_type = "part", 
                     legend_text_italic = F, 
                     color_values = as.vector(pals::polychrome(10)))+ 
  ylab("Frequency (%)")



# Ordination
df_qiime_sub$cal_betadiv(unifrac = T)
t1 <- trans_beta$new(dataset = df_qiime_sub, group = "Treatment", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
t1$plot_ordination(plot_color = "Treatment", plot_type = c("point"))


# FAPROTAX
t1 <- trans_func$new(df_qiime)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
t1$plot_spe_func_perc()






###########################################################################################################################################################










dir.create("PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund")
dataset$save_abund(dirpath = "PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund")
Then, we calculate the alpha diversity. The result is also stored in the object microtable automatically. As an example, we do not calculate phylogenetic diversity.

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
df_qiime$cal_alphadiv(PD = T)

# save dataset$alpha_diversity to a directory
dir.create("PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund/alpha_diversity")
dataset$save_alphadiv(dirpath = "PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund/alpha_diversity")

We also calculate the distance matrix of beta diversity using function cal_betadiv(). We provide four most frequently used indexes: Bray-curtis, Jaccard, weighted Unifrac and unweighted unifrac.

# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package
df_qiime$cal_betadiv(unifrac = T)

dir.create("PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund/alpha_diversity/beta_diversity")
dataset$save_betadiv(dirpath = "PostDoc/Bayreuth/Data_analysis/Sub_2_Chemistry_Incub/Seq/Analysis/taxa_abund/alpha_diversity/beta_diversity")

# create trans_abund object
# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = df_qiime, 
                      taxrank = "Class", 
                      ntaxa = 20)
# t1 object now include the transformed abundance data t1$abund_data and other elements for the following plotting
We remove the sample names in x axis and add the facet to show abundance according to groups.

t1$plot_bar(others_color = "grey70", 
            facet = "Treatment", 
            xtext_keep = F, 
            legend_text_italic = FALSE)
# return a ggplot2 object


# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = df_qiime, taxrank = "Genus", ntaxa = 20, groupmean = "Treatment")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)




# require ggnested package; see https://chiliubio.github.io/microeco_tutorial/intro.html#dependence
test1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10, high_level = "Phylum", prefix = "\\|")
test1$plot_bar(ggnested = TRUE, facet = c("Treatment"),others_color = "red", xtext_angle = 30)

# fixed number in each phylum
test1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 4)
test1$plot_bar(ggnested = TRUE, xtext_keep = F, facet = c("Treatment"))

# sum others in each phylum
test1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 3, prefix = "\\|")
test1$plot_bar(ggnested = TRUE, high_level_add_other = TRUE, xtext_angle = 30, facet = c("Treatment"))









Then alluvial plot is implemented in the plot_bar function.

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10,groupmean = "Treatment")
# use_alluvium = TRUE make the alluvial plot, clustering =TRUE can be used to reorder the samples by clustering
t1$plot_bar(use_alluvium = TRUE, clustering = F,  others_color = "black",xtext_size = 6,order_x = c("Soil","Initial","After_1_day","Control","Glucose","Glutamine"))


The box plot is an excellent way to intuitionally show data distribution across groups.

# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = df_qiime_sub, taxrank = "Genus", ntaxa = 10)
t1$plot_box(group = "Treatment", xtext_angle = 30)


Then we show the heatmap with the high abundant genera.

# show 40 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Treatment", xtext_keep = FALSE, withmargin = FALSE)


Then, we show the pie chart.

t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Treatment")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)


trans_venn class
The trans_venn class is used for venn analysis. To analyze the unique and shared OTUs of groups, we first merge samples according to the "Group" column of sample_table.

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Treatment")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
# The integer data is OTU number
# The percentage data is the sequence number/total sequence number


When the groups are too many to show with venn plot, we can use petal plot.

# use "Type" column in sample_table
dataset1 <- dataset$merge_samples(use_group = "Treatment")
t1 <- trans_venn$new(dataset1)
t1$plot_venn(petal_plot = TRUE)


Sometimes, we are interested in the unique and shared species. In another words, the composition of the unique or shared species may account for the different and similar parts of ecological characteristics across groups[9]. For this goal, we first transform the results of venn plot to the traditional species-sample table, that is, another object of microtable class.

dataset1 <- dataset$merge_samples(use_group = "Treatment")
t1 <- trans_venn$new(dataset1)

trans_alpha class
Alpha diversity can be transformed and plotted using trans_alpha class. Creating trans_alpha object can return two data frame: alpha_data and alpha_stat.

t1 <- trans_alpha$new(dataset = dataset, group = "Treatment")
# return t1$alpha_stat
t1$alpha_stat[1:5, ]

Then, we test the differences among groups using the KW rank sum test and anova with multiple comparisons.

t1$cal_diff(method = "KW")
# return t1$res_alpha_diff
t1$res_alpha_diff[1:5, ]

t1$cal_diff(method = "anova")
# return t1$res_alpha_diff
t1$res_diff

Now, let us plot the mean and se of alpha diversity for each group, and add the duncan.test (agricolae package) result.

t1$plot_alpha(add_letter = TRUE, measure = "Observed")


We can also use the boxplot to show the paired comparisons directly.

t1$plot_alpha(pair_compare = TRUE, measure = "Chao1")


trans_beta class
The distance matrix of beta diversity can be transformed and plotted using trans_beta class. The analysis referred to the beta diversity in this class mainly include ordination, group distance, clustering and manova. We first show the ordination using PCoA.

# we first create an object and select PCoA for ordination
t1 <- trans_beta$new(dataset = dataset, group = "Treatment", measure = "bray")
# t1$res_ordination is the ordination result list
t1$plot_ordination(plot_color = "Treatment", plot_type = c("point"))


Then we plot and compare the group distances.

# calculate and plot sample distances within groups
t1$cal_group_distance()
# return t1$res_group_distance
t1$plot_group_distance(distance_pair_stat = TRUE)


# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$plot_group_distance(distance_pair_stat = TRUE)


Clustering plot is also a frequently used method.

# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Treatment")


perMANOVA is often used in the differential test of distances among groups.

# manova for all groups
t1$cal_manova(cal_manova_all = TRUE)
t1$res_manova$aov.tab

# manova for each paired groups
t1$cal_manova(cal_manova_paired = TRUE)
t1$res_manova


Differential abundance test is a very important part in the microbial community data analysis. It can be used to find the significant taxa in determining the community differences across groups. Currently, trans_diff class have three famous approaches to perform this analysis: metastat[10], LEfSe[11] and random forest. Metastat depends on the permutations and t-test and performs well on the sparse data. It is used for the comparisons between two groups.

t1 <- trans_diff$new(dataset = df_qiime_sub, method = "lefse", group = "Treatment", alpha = 0.05, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 2)
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:30, width = 0.8)


Then, we show the cladegram of the differential features in the taxonomic tree. There are too many taxa in this dataset. As an example, we only use the highest 200 abundant taxa in the tree and 50 differential features. We only show the full taxonomic label at Phylum level and use letters at other levels to reduce the text overlap.

# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
t1$plot_lefse_cladogram(use_taxa_num = 200, use_feature_num = 200, clade_label_level = 6)


There may be a problem related with the taxonomic labels in the plot. When the levels used are too many, the taxonomic labels may have too much overlap. However, if we only indicate the Phylum labels, the taxa in the legend with marked letters are too many. At this time, you can select the taxa that you want to show in the plot manually like the following operation.

# choose some taxa according to the positions in the previous picture; those taxa labels have minimum overlap
use_labels <- c("c__Deltaproteobacteria", "c__Actinobacteria", "o__Rhizobiales", "p__Proteobacteria", "p__Bacteroidetes", 
                "o__Micrococcales", "p__Acidobacteria", "p__Verrucomicrobia", "p__Firmicutes", 
                "p__Chloroflexi", "c__Acidobacteria", "c__Gammaproteobacteria", "c__Betaproteobacteria", "c__KD4-96",
                "c__Bacilli", "o__Gemmatimonadales", "f__Gemmatimonadaceae", "o__Bacillales", "o__Rhodobacterales")
# then use parameter select_show_labels to show
t1$plot_lefse_cladogram(use_taxa_num = 200, use_feature_num = 50, select_show_labels = use_labels)
# Now we can see that more taxa names appear in the tree


If you are interested in taxonomic tree, you can also use metacoder package[Foster_Metacoder_2017] to plot the taxonomic tree based on the selected taxa. We do not show the usage here.

The third approach is rf, which depends on the random forest[12, 13] and the non-parametric test. The current method can calculate random forest by bootstrapping like the method in LEfSe and only use the significant features. MeanDecreaseGini is selected as the indicator value in the analysis.

# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", rf_taxa_level = "Genus")
# t1$res_rf is the result stored in the object
# plot the result
t2 <- t1$plot_diff_abund(use_number = 1:20, only_abund_plot = FALSE)
gridExtra::grid.arrange(t2$p1, t2$p2, ncol=2, nrow = 1, widths = c(2,2))
# the middle asterisk represent the significances


trans_env class
The environmental variables are very useful in analyzing microbial community structure and assembly mechanisms. We first show the RDA analysis (db-RDA and RDA).

# add_data is used to add the environmental data
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
# use bray-curtis distance to do dbrda
t1$cal_rda(use_dbrda = TRUE, use_measure = "bray")
# t1$res_rda is the result list stored in the object
t1$trans_rda(adjust_arrow_length = TRUE, max_perc_env = 10)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_rda(plot_color = "Group")


# use Genus
t1$cal_rda(use_dbrda = FALSE, taxa_level = "Genus")
# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_rda(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1500, max_perc_tax = 3000, min_perc_env = 200, min_perc_tax = 300)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_rda(plot_color = "Group")


Mantel test can be used to check whether there is significant correlations between environmental variables and distance matrix.

t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel
variable_name	cor_method	corr_res	p_res	significance
Temperature	pearson	0.452	0.001	***
  Precipitation	pearson	0.2791	0.001	***
  TOC	pearson	0.13	0.003	**
  NH4	pearson	-0.05539	0.922	
NO3	pearson	0.06758	0.049	*
  pH	pearson	0.4085	0.001	***
  Conductivity	pearson	0.2643	0.001	***
  TN	pearson	0.1321	0.002	**
  The correlations between environmental variables and taxa are important in analyzing and inferring the factors affecting community structure. In this example, we first perform the differential abundance test and random forest analysis to obtain the important genera. Then we use those taxa to perform correlation analysis.

# first create trans_diff object
t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", rf_taxa_level = "Genus")
# then create trans_env object
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
# calculate correlation
t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_rf$Taxa[1:40])
# return t1$res_cor 
Then, we can plot the correlation results using ggplot2 or pheatmap.

# default ggplot2 method
t1$plot_corr()


# clustering heatmap; require pheatmap package
t1$plot_corr(pheatmap = TRUE)


Sometimes, it is necessary to study the correlations between environmental variables and taxa for different groups.

# calculate correlations for different groups using parameter by_group
t1$cal_cor(by_group = "Group", use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_rf$Taxa[1:40])
# return t1$res_cor
t1$plot_corr()


If you are concerned with the relationship between environmental factors and alpha diversity, you can also use this function.

t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = dataset$alpha_diversity)
t1$plot_corr()


trans_nullmodel class
In recent decades, the integration of phylogenetic analysis and null model promotes the inference of niche and neutral influences on community assembly more powerfully by adding a phylogeny dimension [14, 15]. The trans_nullmodel class provides an encapsulation, including the calculation of the phylogenetic signal, beta mean pairwise phylogenetic distance (betaMPD), beta mean nearest taxon distance (betaMNTD), beta nearest taxon index (betaNTI), beta net relatedness index (betaNRI) and Bray-Curtis-based Raup-Crick (RCbray). The approach for phylogenetic signal analysis is based on the mantel correlogram [16], in which the change of phylogenetic signal is intuitional and clear compared to other approaches. The algorithms of betaMNTD and betaMPD have been optimized to be faster than those in the picante package [3]. The combinations between RCbray and betaNTI (or betaNRI) can be used to infer the strength of each ecological process dominating the community assembly under the specific hypothesis [15]. This can be achievable by the function cal_process() to parse the percentage of each inferred process. We first check the phylogenetic signal.

# generate trans_nullmodel object; use 1000 OTUs as example
t1 <- trans_nullmodel$new(dataset, taxa_number = 1000, add_data = env_data_16S)
# use pH as the test variable
t1$cal_mantel_corr(use_env = "pH")
# return t1$res_mantel_corr
# plot the mantel correlogram
t1$plot_mantel_corr()


betaNRI(ses.betampd) is used to show the 'basal' phylogenetic turnover[16]. Compared to betaNTI, it can capture more turnover information associated with the deep phylogeny. It is noted that there are many null models with the development in the several decades. In the trans_nullmodel class, we randomized the phylogenetic relatedness of species. This shuffling approach fix the observed levels of species ??-diversity and ??-diversity to explore whether the observed phylogenetic turnover significantly differ from null model that phylogenetic relatedness among species are random.

# null model run 500 times
t1$cal_ses_betampd(runs=500, abundance.weighted = TRUE)
# return t1$res_ses_betampd
If we want to plot the betaNRI, we can use plot_group_distance function in trans_beta class. For example, the results showed that the mean betaNRI of TW is extremely and significantly larger that those in CW and IW, revealing that the basal phylogenetic turnover in TW is high.

# add betaNRI matrix to beta_diversity list
dataset$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd
# create trans_beta class, use measure "betaNRI"
t2 <- trans_beta$new(dataset = dataset, group = "Group", measure = "betaNRI")
# transform the distance for each group
t2$cal_group_distance()
# plot the results
g1 <- t2$plot_group_distance(distance_pair_stat = TRUE)
g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)


Sometimes, if you want to perform null model analysis for each group individually, such as one group as one species pool, you can calculate the results for each group, respectively. We can find that, when we perform betaNRI for each group respectively, mean betaNRI between CW and TW are not significantly different, and they are both significantly higher than that in IW, revealing that the strength of variable selection in CW and TW may be similar under the condition that each area is considered as a specific species pool.

# we create a list to store the trans_nullmodel results.
sesbeta_each <- list()
group_col <- "Group"
all_groups <- unique(dataset$sample_table[, group_col])
# calculate for each group, respectively
for(i in all_groups){
  # like the above operation, but need provide 'group' and 'select_group'
  test <- trans_nullmodel$new(dataset, group = group_col, select_group = i, taxa_number = 1000, add_data = env_data_16S)
  test$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
  sesbeta_each[[i]] <- test$res_ses_betampd
}
# merge and reshape to generate one symmetrical matrix
test <- lapply(sesbeta_each, melt) %>% do.call(rbind, .) %>%
  reshape2::dcast(., Var1~Var2, value.var = "value") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
# like the above operation
dataset$beta_diversity[["betaNRI"]] <- test
t2 <- trans_beta$new(dataset = dataset, group = "Group", measure = "betaNRI")
t2$cal_group_distance()
g1 <- t2$plot_group_distance(distance_pair_stat = TRUE)
g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)


BetaNTI(ses.betamntd) can be used to indicate the phylogenetic terminal turnover [15].

# null model run 500 times
t1$cal_ses_betamntd(runs=500, abundance.weighted = TRUE)
# return t1$res_ses_betamntd
1	2	3	4	5
0	-5.637	-5.701	-5.758	-5.531
-5.637	0	-6.101	-6.275	-6.013
-5.701	-6.101	0	-6.333	-6.197
-5.758	-6.275	-6.333	0	-6.141
-5.531	-6.013	-6.197	-6.141	0
RCbray (Bray-Curtis-based Raup-Crick) can be calculated using function cal_rcbray() to assess whether the compositional turnover was governed primarily by drift [17]. We applied null model to simulate species distribution by randomly sampling individuals from each species pool with preserving species occurrence frequency and sample species richness [16].

# result stored in t1$res_rcbray
t1$cal_rcbray(runs = 1000)
# return t1$res_rcbray
As an example, we also calculate the proportion of the inferred processes on the community assembly as shown in the references [15, 16]. In the example, the fraction of pairwise comparisons with significant betaNTI values (|??NTI| > 2) is the estimated influence of Selection; ??NTI > 2 represents the heterogeneous selection; ??NTI < -2 represents the homogeneous selection. The value of RCbray characterizes the magnitude of deviation between observed Bray-Curtis and Bray-Curtis expected under the randomization; a value of |RCbray| > 0.95 was considered significant. The fraction of all pairwise comparisons with |??NTI| < 2 and RCbray > +0.95 was taken as the influence of Dispersal Limitation combined with Drift. The fraction of all pairwise comparisons with |??NTI| < 2 and RCbray < -0.95 was taken as an estimate for the influence of Homogenizing Dispersal. The fraction of all pairwise comparisons with |??NTI| < 2 and |RCbray| < 0.95 estimates the influence of Drift acting alone.

# use betaNTI and rcbray to evaluate processes
t1$cal_process(use_betamntd = TRUE)
# return t1$res_process
t1$res_process
process	percentage
variable selection	4.419
homogeneous selection	48.71
dispersal limitation	0
homogeneous dispersal	8.739
drift	38.13
trans_network class
Network is a frequently used approach to study the co-occurrence patterns in microbial ecology[18-20]. In this part, we describe all the core contents in the trans_network class. The network construction approaches can be classified into two types: correlation-based and non correlation-based. Several approaches can be used to calculate correlations and significances.

We first introduce the correlation-based network. The parameter cal_cor in trans_network is used for selecting the correlation calculation method.

# Use R base cor.test, slow
t1 <- trans_network$new(dataset = dataset, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
# return t1$res_cor_p list; one table: correlation; another: p value
# SparCC method, require SpiecEasi package
# SparCC is very slow, so consider filtering more species with low abundance
t1 <- trans_network$new(dataset = dataset, cal_cor = "SparCC", taxa_level = "OTU", filter_thres = 0.001, SparCC_simu_num = 100)
# When the OTU number is large, use R WGCNA package to replace R base to calculate correlations
# require WGCNA package
t1 <- trans_network$new(dataset = dataset, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
The parameter COR_cut can be used to select the correlation threshold. Furthermore, COR_optimization = TRUE represent using RMT theory to find the optimized correlation threshold instead of the COR_cut[18].

# construct network; require igraph package
t1$cal_network(p_thres = 0.01, COR_optimization = TRUE)
# return t1$res_network
# use arbitrary coefficient threshold to contruct network
t1$cal_network(p_thres = 0.01, COR_cut = 0.7)
# save network
# open the gexf file using Gephi(https://gephi.org/)
# require rgexf package
t1$save_network(filepath = "network.gexf")
We plot the network and present the node colors according to the calculated modules in Gephi.

Now, we show the node colors with the Phylum information and the edges colors with the positive and negative correlations. All the data used has been stored in the network.gexf file, including modules classifications, Phylum information and edges classifications.



# calculate network attributes
t1$cal_network_attr()
# return t1$res_network_attr
Property	Value
Vertex	407
Edge	1989
Average_degree	9.774
Average_path_length	3.878
Network_diameter	9
Clustering_coefficient	0.4698
Density	0.02407
Heterogeneity	1.194
Centralization	0.09908
# classify the node; return t1$res_node_type
t1$cal_node_type()
# return t1$res_node_type
# we retain the file for the following example in trans_func part
network_node_type <- t1$res_node_type
z	module	p	taxa_roles
OTU_50	-1.305	M2	0	Peripheral nodes
OTU_1	-0.04067	M2	0	Peripheral nodes
OTU_55	-1.239	M2	0	Peripheral nodes
OTU_13824	-0.2403	M2	0	Peripheral nodes
OTU_151	-1.372	M2	0.4444	Peripheral nodes
# plot node roles in terms of the within-module connectivity and among-module connectivity
t1$plot_taxa_roles(use_type = 1)


# plot node roles with phylum information
t1$plot_taxa_roles(use_type = 2)


Now, we show the eigengene analysis of modules. The eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module.

t1$cal_eigen()
# return t1$res_eigen
Then we perform correlation heatmap to show the relationships between eigengenes and environmental factors.

# create trans_env object like the above operation
env_data_16S=data.frame(sample_data(leaf_raw))
t2 <- trans_env$new(dataset = cal_filter_final, add_data = env_data_16S[, 5:10])
# calculate correlations
t2$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap
t2$plot_corr()


The function cal_sum_links() is used to sum the links (edge) number from one taxa to another or in the same taxa. The function plot_sum_links() is used to show the result from the function cal_sum_links(). This is very useful to fast see how many nodes are connected between different taxa or within one taxa. In terms of "Phylum" level in the tutorial, the function cal_sum_links() sum the linkages number from one Phylum to another Phylum or the linkages in the same Phylum. So the numbers along the outside of the circular plot represent how many edges or linkages are related with the Phylum. For example, in terms of Proteobacteria, there are roughly total 900 edges associated with the OTUs in Proteobacteria, in which roughly 200 edges connect both OTUs in Proteobacteria and roughly 150 edges connect the OTUs from Proteobacteria with the OTUs from Chloroflexi.

# calculate the links between or within taxonomic ranks
t1$cal_sum_links(taxa_level = "Phylum")
# return t1$res_sum_links_pos and t1$res_sum_links_neg
# require chorddiag package
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10)


The subset_network() function can be used to extract a part of nodes and edges among these nodes from the network. In this function, you should provide the nodes you need using the node parameter.

# this return a sub network that contains all nodes of module M1
t1$subset_network(node = t1$res_node_type %>% .[.$module == "M1", ] %>% rownames, rm_single = TRUE)
# return a new network with igraph class
Then we show the next implemented network construction approach: SPIEC-EASI (SParse InversE Covariance Estimation for Ecological Association Inference) network in SpiecEasi R package [21].

# cal_cor select NA
t1 <- trans_network$new(dataset = dataset, cal_cor = NA, taxa_level = "OTU", filter_thres = 0.0005)
# require SpiecEasi package  https://github.com/zdk123/SpiecEasi
t1$cal_network(network_method = "SpiecEasi")
# see t1$res_network
We also introduce the third network construction approach: Probabilistic Graphical Models (PGM), which is implemented in julia package FlashWeave[22]. If you want to use this method like the following code, you should first install julia language in your computer and the FlashWeave package, and add the julia in the computer path (see FlashWeave part in https://github.com/ChiLiubio/microeco).

# cal_cor select NA
t1 <- trans_network$new(dataset = dataset, cal_cor = NA, taxa_level = "OTU", filter_thres = 0.0001)
# require Julia in the computer path, and the package FlashWeave
t1$cal_network(network_method = "PGM")
# see t1$res_network
trans_func class
Ecological researchers are usually interested in the the funtional profiles of microbial communities, because functional or metabolic data is powerful to explain the structure and dynamics of microbial communities and to infer the underlying mechanisms. As metagenomic sequencing is complicated and expensive, using amplicon sequencing data to predict functional profiles is a good choice. Several software are often used for this goal, such as PICRUSt[23], Tax4Fun[24] and FAPROTAX[25, 26]. These tools are great to be used for the prediction of functional profiles based on the prokaryotic communities from sequencing results. In addition, it is also important to obtain the functions for each taxa or OTU, not just the whole profile of communities. But it is hard to know exact functions of each OTU. FAPROTAX database is a collection of the traits and characteristics of prokaryotes based on the known research results published in books and literatures. We match the taxonomic information of prokaryotes against this database to identify the traits of prokaryotes on biogeochemical roles. We also implement the FUNGuild database[27] to identify the traits of fungi.

# Identify microbial traits
# create object of trans_func
t2 <- trans_func$new(dataset)
# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
t2$cal_spe_func()
# return t2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
methanotrophy	acetoclastic_methanogenesis
OTU_4272	0	0
OTU_236	0	0
OTU_399	0	0
OTU_1556	0	0
OTU_32	0	0
The percentages of the OTUs having the same trait can reflect the functional redundancy of this function in the community or the module in the network.

# calculate the percentages of OTUs for each trait in each module of network
# use_community = FALSE represent calculating module, not community, node_type_table provide the module information
t2$cal_spe_func_perc(use_community = FALSE, node_type_table = network_node_type)
# return t2$res_spe_func_perc
# we only plot some important traits, so we use the default group list to filter and show the traits.
t2$plot_spe_func_perc(select_samples = paste0("M", 1:10))
# M represents module, ordered by the nodes number from high to low


# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("methanotrophy")
# calculate the percentages for communities
t2$cal_spe_func_perc(use_community = TRUE)
# t2$res_spe_func_perc[1:5, 1:2]
methanotrophy	acetoclastic_methanogenesis
0.39	0.04
0.27	0
0.48	0
0.48	0
0.56	0
# then we try to correlate the res_spe_func_perc of communities to environmental variables
t3 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
t3$plot_corr(pheatmap = TRUE)


Tax4Fun requires a strict input file demand on the taxonomic information. To analyze the trimmed or changed OTU data in R with Tax4Fun, we provide a link to the Tax4Fun functional prediction.

t1 <- trans_func$new(dataset)
# install Tax4Fun package and download SILVA123 ref data from  http://tax4fun.gobics.de/
# decompress SILVA123; provide path in folderReferenceData as you put
t1$cal_tax4fun(folderReferenceData = "./SILVA123")
# return two files: t1$tax4fun_KO: KO file; t1$tax4fun_path: pathway file.
# t1$tax4fun_KO$Tax4FunProfile[1:5, 1:2]
K00001; alcohol dehydrogenase [EC:1.1.1.1]	K00002; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]
0.0004823	5.942e-06
0.0005266	4.017e-06
0.0005054	6.168e-06
0.0005109	5.888e-06
0.0005083	5.547e-06
Now, we use pathway file to analyze the abundance of pathway.

# must transpose to taxa row, sample column
pathway_file <- t1$tax4fun_path$Tax4FunProfile %>% t %>% as.data.frame
# filter rownames, only keep ko+number
rownames(pathway_file) %<>% gsub("(^.*);\\s.*", "\\1", .)
# load the pathway hierarchical metadata
data(ko_map)
# create a microtable object, familiar?
func1 <- microtable$new(otu_table = pathway_file, tax_table = ko_map, sample_table = t1$sample_table)
print(func1)
microtable class: sample_table have 90 rows and 4 columns otu_table have 284 rows and 90 columns tax_table have 341 rows and 3 columns

Now, we need to trim data and calculate abundance.

func1$tidy_dataset()
# calculate abundance automatically at three levels: level_1, level_2, level_3
func1$cal_abund()
print(func1)
microtable class: sample_table have 90 rows and 4 columns otu_table have 284 rows and 90 columns tax_table have 284 rows and 3 columns Taxa abundance: calculated for level_1,level_2,level_3

Then, we can plot the abundance.

# bar plot at level_1
func2 <- trans_abund$new(func1, taxrank = "level_1", groupmean = "Group")
func2$plot_bar(legend_text_italic = FALSE)


We can also do something else. For example, we can use lefse to test the differences of the abundances and find the important enriched pathways across groups.

func2 <- trans_diff$new(dataset = func1, method = "lefse", group = "Group", alpha = 0.05, lefse_subgroup = NULL)
func2$plot_lefse_bar(LDA_score = 3, width = 0.8)


Now, we use the ITS amplicon sequencing dataset as an example to show the use of FUNGuild database[27].

# show the ITS dataset preprocessing, the functional identification of OTUs and functional redundancy of modules
data(sample_info_ITS)
data(otu_table_ITS)
data(taxonomy_table_ITS)
# create microtable object
dataset <- microtable$new(sample_table = sample_info_ITS, otu_table = otu_table_ITS, tax_table = taxonomy_table_ITS)
# remove the taxa not assigned in the Kingdom "k__Fungi"http://127.0.0.1:24307/graphics/plot_zoom_png?width=1026&height=822
dataset$tax_table %<>% base::subset(Kingdom == "k__Fungi")
# use tidy_dataset() to make OTUs and samples information consistent across files
dataset$tidy_dataset()
# create trans_network object
t1 <- trans_network$new(dataset = dataset, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0.000001, cor_method = "spearman")
# create correlation network 
t1$cal_network(p_thres = 0.05, COR_cut = 0.6)
# calculate node topological properties
t1$cal_node_type()
node_type_table <- t1$res_node_type
# create trans_func object
t2 <- trans_func$new(dataset)
# identify species traits, automatically select database for prokaryotes or fungi
t2$cal_spe_func()
# calculate abundance-unweighted functional redundancy of each trait for each network module
t2$cal_spe_func_perc(use_community = FALSE, node_type_table = node_type_table)
# plot the functional redundancy of network modules
t2$plot_spe_func_perc(select_samples = paste0("M", 1:10))


Notes
add layers to plot
Most of the plots are generated by applying the ggplot2 package. The important parameters in the plotting functions are configured according to our experience. If the inner parameters can not meet the use, the user can add the layers to the plot like the following operation or make the plot using the data (generally data.frame class) stored in the object.

# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))


clone
R6 class has a special copy mechanism which is different from S3 and S4. If you want to copy an object completely, you should use function clone instead of direct assignment.

# use clone to copy completely
t1 <- clone(dataset)
t2 <- clone(t1)
t2$sample_table <- NULL
identical(t2, t1)
[1] FALSE

# this operation is usually unuseful, because changing t2 will also affect t1
t2 <- t1
t2$sample_table <- NULL
identical(t2, t1)
[1] TRUE

change object
All the classes are set public, meaning that you can change, add or remove the objects stored in them as you want.

# add a matrix you think useful
dataset$my_matrix <- matrix(1, nrow = 4, ncol = 4)
# change the information
dataset$sample_table %<>% .[, -2]


