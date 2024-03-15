library(readr)
library(ggplot2)
#1a 
you <- read_tsv("Q1.data") 
phenotype <- you[,2001]
genotype <- you[,1:2000]

#p-values
p_values <- numeric(length = 2000)
for (i in 1:2000) {
  model <- lm(as.matrix(phenotype) ~ genotype[[i]])
  p_values[i - 1] <- summary(model)$coefficients[2, 4]
}
observed_p_values <- -log10(p_values)
expected_p_values <- -log10(ppoints(length(p_values)))

#QQplot
qqplot(expected_p_values, observed_p_values,
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
       main = "QQ Plot of SNP-Phenotype Association")

#Bonferroni
alpha <- 0.05
bonferroni_threshold <- (alpha / length(p_values))
significant_snps <- which(observed_p_values > bonferroni_threshold)

#read
cat("significantly associated SNPs at α = 0.05:",
    length(significant_snps), "\n")

significant_snps

#1b
result <- prcomp(genotype, scale. = TRUE)
pcs <- result$x
plot(1:length(result$sdev), result$sdev^2,
     type = "b", pch = 16, xlab = "Principal Component",
     ylab = "Variance", main = "Scree Plot")
ggplot(data.frame(PC1 = pcs[, 1], PC2 = pcs[, 2])) +
  geom_point(aes(x = PC1, y = PC2), size = 3) +
  labs(x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)",
       title = "PC1 vs PC2 Plot")
pc1 <- pcs[,1]
length(pc1[pc1 < 0])

#1c
p_values <- numeric(length = ncol(genotype))
for (i in 1:ncol(genotype)) {
  model <- lm(as.matrix(phenotype) ~ genotype[[i]] + pc1)
  p_values[i] <- summary(model)$coefficients[2, 4]  # p-value for the SNP coefficient
}
observed_p_values <- -log10(p_values)
expected_p_values <- -log10(ppoints(length(p_values)))

qqplot(expected_p_values, observed_p_values,
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
       main = "QQ Plot of SNP-Phenotype Association with PC1")

alpha <- 0.05
bonferroni_threshold <- (alpha / length(p_values))

significant_snps <- which(observed_p_values > bonferroni_threshold)

cat("Number of significantly associated SNPs at α = 0.05 (Bonferroni correction):",
    length(significant_snps), "\n")

