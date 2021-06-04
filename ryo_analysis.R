

# Libraries
library(tidyverse)


# Data

# Gene data
gene <- read.csv("TCGA-LUSC.htseq_counts.tsv", sep = "\t")

# Survival data
survival <- read.csv("TCGA-LUSC.survival.tsv", sep = "\t")

# Phenotype data
phenotype <- read.csv("TCGA-LUSC.GDC_phenotype.tsv", sep = "\t")

# Gene annotations
annotations <- read.csv("annotations.csv")

patient_data <- survival %>%
  inner_join(phenotype, by = c("sample" = "submitter_id.samples")) 

# Phenotype:
#   Demographics: Gender, Race, age_at_index
#   Diagnosis: age_at_diagnosis, icd_10_diagnoses, tissue_or_organ_of_origin
#   Smoking: tobacco_smoking_history, pack_year_smoked.exposures, years_of_tobacco_smoking_onset

# Functions
replace_period_with_hyphen <- function(arr){
  gene_pat2 <- c()
  
  for(i in arr){
    temp1 <- str_replace_all(i, '\\.', '-')
    gene_pat2 <- append(gene_pat2, temp1)
    
  }
  gene_pat2
}

sample_to_patient <- function(arr){
  gene_pat2 <- c()
  
  for(i in arr){
    temp1 <- substr(i, 1, 12)
    gene_pat2 <- append(gene_pat2, temp1)
    
  }
  gene_pat2
}
# ----






# Risk of early onset lung cancer?

m1 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history
        , data = phenotype)
summary(m1)
summary(lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + pack_years_smoked.exposures
   , data = phenotype))

phenotype %>%
  group_by(tobacco_smoking_history) %>%
  count()
  select(initial_weight.samples)



# Ensemble ID
ensembl_id <- gene$Ensembl_ID
  
  
# Sample IDs
genes_id <- colnames(gene)

phenotype_id <- phenotype$submitter_id.samples

genes_id2 <- replace_period_with_hyphen(genes_id)


# Figure out which samples are in all 3 dataframes

all_id <- intersect(genes_id2, phenotype_id)

all_patients <- sample_to_patient(all_id)

# Filter the datasets by common samples


genes_index <- genes_id2 %in% all_id
phenotype_index <- phenotype_id %in% all_id

phenotype2 <- phenotype[phenotype_index,]
genes2 <- as.data.frame(t(gene[, genes_index]))
colnames(genes2) <- ensembl_id
rownames(genes2) <- replace_period_with_hyphen(rownames(genes2))

annotations <- annotations[2:3]

ensembl_id2 <- c()
annotation_2 <- c()

for(i in ensembl_id){
  i2 <- substr(i, 1, 15)
  ensembl_id2 <- append(ensembl_id2, i2)
  
  if(i2 %in% annotations$gene_id){
    the_index <- which(annotations$gene_id == i2)
    the_type <- annotations[the_index,2]
    annotation_2 <- append(annotation_2, the_type)
  }
  else{
    annotation_2 <- append(annotation_2, NA)
  }
}

# Get only protein coding genes
pc_index <- annotation_2 == 'protein_coding'
pc_index2 <- which(pc_index %in% c(TRUE))
genes_pc <- genes2[, pc_index2]


genes_pc

# Get top 10 pc's
# Scale & center?

# Remove zero-variance genes
genes_var <- genes_pc %>%
  select(-(names(genes_pc[, sapply(genes_pc, function(v) var(v, na.rm=T)==0)])))

pca1 <- prcomp(genes_pc)



top_10_pc <- as.data.frame(pca1$x[,1:10])



covariates <- c('gender.demographic',
                'race.demographic',
                'ethnicity.demographic',
                'age_at_diagnosis.diagnoses',
                'tobacco_smoking_history')


phenotype3 <- phenotype2 %>%
  select(age_at_diagnosis.diagnoses
         , gender.demographic
         , race.demographic
         , tobacco_smoking_history)

# Changing age at diagnosis in years

# phenotype3$age_at_diagnosis.diagnoses <- phenotype3$age_at_diagnosis.diagnoses/365.25

rownames(phenotype3) <- phenotype2$submitter_id.samples

finaldat <- cbind(phenotype3, top_10_pc)



summary(m1)
m.pc.1 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1
              , data = finaldat)
summary(m.pc.1)
m.pc.2 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2
             , data = finaldat)
summary(m.pc.2)
m.pc.3 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3
             , data = finaldat)
summary(m.pc.3)
m.pc.4 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4
             , data = finaldat)
summary(m.pc.4)
m.pc.5 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5
             , data = finaldat)
summary(m.pc.5)
m.pc.6 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6
             , data = finaldat)
summary(m.pc.6)
m.pc.7 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
             , data = finaldat)
summary(m.pc.7)
m.pc.8 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
             , data = finaldat)
summary(m.pc.8)
m.pc.9 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9
             , data = finaldat)
summary(m.pc.9)
m.pc.10 <- lm(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + tobacco_smoking_history + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
             , data = finaldat)
summary(m.pc.10)

####################
# Mixed models
####################
library(lme4)
library(nlme)
# load("df_train.RData")
# load("df_test.RData")
load("df_train_v2.RData") # imports as df_train
load("df_test_v2.RData") # imports as df_test



pclist <- colnames(df_train)[8:ncol(df_train)]



# df_train$IDBoot <- sample(1:nrow(df_train), nrow(df_train), replace = T)

grps <- rep(1:3, each=141)


variances <- c()
for (i in 1:100) {
  print(i)
  df_train$IDBoot <- as.factor(sample(grps))

  tryCatch({
    lmm.boot <- lme(age_at_diagnosis.diagnoses ~ gender.demographic + tobacco_smoking_history + race.demographic +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11
                    , random = ~ -1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 | IDBoot
                    , method = "REML"
                    , data = df_train)
    
    var.boot <- (summary(lmm.boot)$sigma)**2
    variances <- c(variances, var.boot)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  
}

vardf <- data.frame(variances = variances)


library(ggplot2)

ggplot(data = vardf, aes(variances)) +
  geom_density() +
  xlim(60, 65)





lmm.tob <- lme(age_at_diagnosis.diagnoses ~ gender.demographic + race.demographic + ethnicity.demographic +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11
                , random = ~ 1 + PC1 + PC2 + PC4 | tobacco_smoking_history
                , method = "REML"
                , data = df_train)
summary(lmm.tob)$sigma ** 2

# genes_var: all protein-coding genes used in model

head(genes_var[1:3])




# Manhattan plot code


manhattandat <- phenotype3 %>%
  transmute(age_at_diagnosis.diagnoses
            , ID = row.names(phenotype3)) %>%
  left_join(genes_var %>%
              mutate(
                ID = row.names(genes_var))
            , by = c("ID")
  )



# Create manhattan dataset

p.vals <- c()
for (i in 3:ncol(manhattandat)) {
  pct <- ((i-2) / (ncol(manhattandat) - 2)) * 100
  print(paste(round(pct, 2), "%"))
  g <- colnames(manhattandat)[i]
  call <- paste0("age_at_diagnosis.diagnoses ~ ", g)
  model <- lm(
    formula = call
    , data = manhattandat)
  mod.sum <- summary(model)$coefficients
  if (dim(mod.sum)[1] < 2) {
    p.vals <- c(p.vals, NA)
  } else {
   p <- mod.sum[2, 4]
   p.vals <- c(p.vals, p)
  }
}

manhattan.pvals <- data.frame(gene = colnames(manhattandat)[3:ncol(manhattandat)]
           , geneid = substr(colnames(manhattandat)[3:ncol(manhattandat)], 1, 15)
           , pval = p.vals)

ensg_mapping <- read.csv("mart_export.txt", sep="\t")

library(stringr)
manhattan.pvals <- manhattan.pvals %>%
  left_join(ensg_mapping, by = c("geneid" = "Gene.stable.ID"))  %>%
  group_by(geneid) %>%
  mutate(rank = row_number(Transcript.start..bp.)
         , chr_start = sprintf("%09d", Transcript.start..bp.)
         , chrm = ifelse(str_pad(Chromosome.scaffold.name, width=2, side="left", pad="0") == "0X", "CX"
                         , ifelse(str_pad(Chromosome.scaffold.name, width=2, side="left", pad="0") == "0Y", "CY", str_pad(Chromosome.scaffold.name, width=2, side="left", pad="0")))) %>%
  arrange(geneid, rank) %>%
  filter(rank == 1) %>%
  mutate(geneloc = paste0("chr", chrm, ":", chr_start)) %>%
  arrange(geneloc, pval)

bonf <- log(0.05 / nrow(manhattan.pvals))
ggplot(data = manhattan.pvals
       , aes(x = geneloc
             , y = log(pval)
             , color = chrm)) +
  geom_point(size=0.75
             , alpha=0.7) +
  # facet_grid( ~ chrm, scales = "free_x", switch = "x", space = "free_x") +
  geom_hline(yintercept = bonf
             , color = "darkred"
             , linetype = "dashed") +
  scale_y_reverse() +
  xlab("Gene Chromosome Position, Ordered") +
  ylab("P-value (log)") +
  theme(axis.ticks.x=element_blank()
        , axis.text.x=element_blank()
        # , legend.position="none"
        , axis.text.y = element_text(size=12)
        , axis.title.x = element_text(size=16)
        , axis.title.y = element_text(size=16)
        # , panel.spacing = unit(0.1, "lines")
        # , strip.background = element_blank()
        # , strip.placement = "outside"
        )



ggplot(data = dat) +
  aes(x = subcat, y = value, fill = subcat) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(mapping = aes(label = paste0(value, "%")), vjust = -0.5) +
  facet_wrap( ~ category, strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  xlab("x-axis label")
