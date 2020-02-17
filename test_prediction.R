library(caret)
library(DESeq2)
library(compositions)
library(vegan)

# inporting gene counts
genes <- readRDS("./data/gene_counts.rds")

# inporting proportion of taxa
taxa <- readRDS("./data/taxa_ab.rds")

# sample data
meta <- readRDS("./data/sample_meta.rds")

# CLR transformation
taxa.clr <- clr(taxa)
class(taxa.clr) <- "matrix"

genes.vst <- DESeqDataSetFromMatrix(t(genes), colData = meta, 
                                    design = ~ 1)
genes.vst <- estimateSizeFactors(genes.vst)
genes.vst <- vst(genes.vst)
genes.vst <- t(assay(genes.vst))

taxa.pca <- prcomp(taxa.clr, center = T, scale. = T)
genes.pca <- prcomp(genes.vst, center = T, scale. = T)

taxa.x <- taxa.pca$x[,taxa.pca$sdev^2 >= 1]
genes.x <- genes.pca$x[,genes.pca$sdev^2 >= 1]

colnames(taxa.x) <- gsub("PC", "T", colnames(taxa.x))
colnames(genes.x) <- gsub("PC", "F", colnames(genes.x))

alpha <- diversity(taxa, index = "invsimpson")
alpha <- data.frame(A=(alpha - mean(alpha))/sd(alpha))

y <- meta$genotype
full.data <- cbind(y = y, taxa.x, genes.x, alpha)
full.data <- droplevels(full.data[full.data$y != "other",])

set.seed(1239)
train <- createDataPartition(full.data$y, times = 1, 
                             list = FALSE, p = .8)

model.train <- full.data[train,]
model.test <- full.data[-train,]

cntr <- trainControl(method = "adaptive_cv",
                     number = 10, repeats = 3,
                     search = "random",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     sampling = "up",
                     savePredictions = T,
                     adaptive = list(min = 5, alpha = 0.05, 
                                     method = "gls", 
                                     complete = FALSE))

set.seed(123)
rf <- train(y ~ ., data = model.train, 
            method = "rf",
            trControl = cntr, tuneLength = 15,
            verbose = FALSE, scale = FALSE,
            metric = "ROC")

set.seed(123)
gbm <- train(y ~ ., data = model.train, 
             method = "gbm",
             trControl = cntr, tuneLength = 15,
             verbose = FALSE,
             metric = "ROC")

set.seed(123)
svm <- train(y ~ ., data = model.train, 
             method = "svmRadial",
             trControl = cntr, tuneLength = 15,
             verbose = FALSE, scale = FALSE,
             metric = "ROC")

model.list <- list(RF = rf,
                   GBM = gbm,
                   SVM = svm)
resamps <- resamples(model.list)
bwplot(resamps)

library(pROC)
rocCurve <- function(model, data){
  best <- model$bestTune
  pred <- model$pred
  d <- data.frame(pred[,names(best)])
  extr <- apply(d, 1, function(x) all(x == best))
  pred <- pred[extr,]
  
  p <- predict(model, data, type = "prob")[,1]
  roc_obj <- roc(data$y, p, levels = levels(data$y))
  
  plot.roc(pred$obs, pred$heterozygote_F508)
  plot.roc(roc_obj)
}

old <- par(mfrow = c(3, 2), oma=c(0,3,1,0))
invisible(sapply(model.list, rocCurve, model.test))
mtext( 'Train data', side=3, line=-1, at=.3, outer=TRUE, cex=.8, font = 2)
mtext( 'Test data', side=3, line=-1, at=.8, outer=TRUE, cex=.8, font = 2)
mtext( 'RF', side=3, line=-5.2, at=-.04, outer=TRUE, cex=.8, font = 2)
mtext( 'GBM', side=3, line=-16.2, at=-.03, outer=TRUE, cex=.8, font = 2)
mtext( 'SVM', side=3, line=-27.2, at=-.03, outer=TRUE, cex=.8, font = 2)
par(old)


get_roc <- function(model, data){
  roc_obj <- roc(data$y, 
                 predict(model, data, type = "prob")[,1],
                 levels = levels(model.test$y))
  setNames(as.vector(ci(roc_obj)), 
           c("lower", "ROC", "upper"))
}
roc <- t(sapply(model.list, get_roc, model.test))
accuracy <- sum(model.test$y == predict(rf, model.test)) / nrow(model.test)


v.imp <- varImp(rf, useModel = T, scale = F)$importance
vars <- rownames(v.imp)
v.imp <- data.frame(vars = vars,
                    imp = v.imp[[1]],
                    var.type = substr(vars, 1, 1))

library(ggbeeswarm)
library(tidyverse)
top5 <- quantile(v.imp$imp, .95)
ggplot(v.imp, aes(x = var.type, y = imp)) +
  geom_quasirandom(shape = 1) +
  geom_hline(yintercept = top5,
             linetype = 2,
             color = "red") +
  theme_bw(base_size = 10, base_family = "Helvetica",
           base_line_size = 0.25)

varContrib <- function(data.pca){
  # Get coordinates for variables (loadings X standard deviation)
  var.coord <- t(t(data.pca$rotation) * data.pca$sdev)
  # Get quality of representation on the factor map (coordinates^2)
  var.cos2 <- var.coord^2
  # Get variable contributions to the PCs (relative importance)
  t(t((var.cos2*100)) / colSums(var.cos2))
}

top.var <- as.character(v.imp$vars[v.imp$imp >= top5])

# Genes contribution to PCs
genes.contrib <- varContrib(genes.pca)
colnames(genes.contrib) <- gsub("PC", "F", colnames(genes.contrib))
genes.contrib <- genes.contrib[, grep("F", top.var, value = T)]
genes.contrib <- data.frame(best.og = rownames(genes.contrib),
                            contrib = rowMeans(genes.contrib),
                            stringsAsFactors = F)

meta.genes <- readRDS("data/gene_meta.rds")

meta.genes %>%
  separate_rows(`COG cat`, sep = ", ") %>%
  left_join(genes.contrib, by = "best.og") %>%
  filter(`COG cat` != "S") %>%
  group_by(`COG cat`) %>%
  summarise(contrib = sum(contrib)) %>%
  mutate(`COG cat` = fct_reorder(`COG cat`, contrib)) %>%
  ggplot(aes(x = `COG cat`, y = contrib)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8)) +
  theme_classic(base_size = 10, base_family = "Helvetica",
                base_line_size = 0.25) +
  ylab("Average contribution on selected components (%)")


# Taxa contribution ot PCs
taxa.contrib <- varContrib(taxa.pca)
colnames(taxa.contrib) <- gsub("PC", "T", colnames(taxa.contrib))
taxa.contrib <- taxa.contrib[, grep("T", top.var, value = T)]
taxa.contrib <- data.frame(taxa = rownames(taxa.contrib),
                           contrib = rowMeans(taxa.contrib),
                           stringsAsFactors = F)

taxa.contrib %>%
  mutate(taxa = fct_reorder(taxa, contrib,
                            .fun = sum, 
                            .desc = T)) %>%
ggplot(aes(x = taxa, y = contrib)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8)) +
  theme_classic(base_size = 10, base_family = "Helvetica",
                base_line_size = 0.25) +
  ylab("Average contribution on selected components (%)") +
  xlab("Taxa") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
# Cleaning output
taxa.contrib %>%
  mutate(taxa = fct_reorder(taxa, contrib,
                            .fun = sum, 
                            .desc = T)) %>%
  filter(!grepl("virus", taxa)) %>%
  filter(!grepl("phage", taxa)) %>%
  filter(!grepl("unclassified", taxa)) %>%
  filter(!grepl("DNA", taxa)) %>%
  ggplot(aes(x = taxa, y = contrib)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4)) +
  coord_cartesian(xlim = c(1, 39.9)) +
  theme_classic(base_size = 10, base_family = "Helvetica",
                base_line_size = 0.25) +
  ylab("Average contribution on selected components (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank())
