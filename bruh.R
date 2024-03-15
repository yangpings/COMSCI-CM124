# Load data
genetic_data <- read.csv("train.genotype.txt", sep = " ")
phenotypes <- read.csv("train.phenotype.txt")
test_data <- read.csv("test.genotype.txt", sep = " ")

combined_data <- cbind(phenotypes, genetic_data)
num_snps <- ncol(combined_data) - 1  
new_column_names <- c("phenotypes", paste0("SNP", 1:num_snps))
coles<-c(paste0("SNP",1:num_snps))
colnames(combined_data) <- new_column_names
colnames(test_data) <- coles


#training
library(glmnet)
X <- as.matrix(combined_data[, -1])
y <- combined_data$phenotypes
lasso_model <- cv.glmnet(X, y, alpha = 0.5)
selected_features <- coef(lasso_model, s = "lambda.min")[, 1]
selected_featuress <- names(selected_features[selected_features != 0])
selected_data <- cbind(phenotypes, combined_data[,2],combined_data[,8],combined_data[,19],combined_data[,24])
selected_data <- cbind(selected_data, combined_data[,33],combined_data[,35],combined_data[,47],combined_data[,51])
selected_data <- cbind(selected_data, combined_data[,55],combined_data[,56],combined_data[,59],combined_data[,61])
selected_data <- cbind(selected_data, combined_data[,71],combined_data[,90],combined_data[,91],combined_data[,103])
selected_data <- cbind(selected_data, combined_data[,105],combined_data[,118],combined_data[,128],combined_data[,133])
selected_data <- cbind(selected_data, combined_data[,134],combined_data[,149],combined_data[,150],combined_data[,153])
selected_data <- cbind(selected_data, combined_data[,157],combined_data[,164],combined_data[,165],combined_data[,179])
selected_data <- cbind(selected_data, combined_data[,188])
linear_model <- lm(X.1.445386 ~ ., data = selected_data)
summary(linear_model)

test_selected_data<- cbind(test_data[,1],test_data[,7],test_data[,18],test_data[,23])
test_selected_data <- cbind(test_selected_data, test_data[,32],test_data[,34],test_data[,46],test_data[,50])
test_selected_data <- cbind(test_selected_data, test_data[,54],test_data[,55],test_data[,58],test_data[,60])
test_selected_data <- cbind(test_selected_data, test_data[,70],test_data[,89],test_data[,90],test_data[,102])
test_selected_data <- cbind(test_selected_data, test_data[,104],test_data[,117],test_data[,127],test_data[,132])
test_selected_data <- cbind(test_selected_data, test_data[,133],test_data[,148],test_data[,149],test_data[,152])
test_selected_data <- cbind(test_selected_data, test_data[,156],test_data[,163],test_data[,164],test_data[,178])
test_selected_data <- cbind(test_selected_data, test_data[,187])

coefficients <- coef(linear_model)
intercept <- coefficients[1]
linear_combination <- as.matrix(test_selected_data) %*% coefficients[-1] + intercept
predicted_phenotypes <- as.vector(linear_combination)
test_results <- data.frame(Phenotypes = predicted_phenotypes)
write.csv(test_results, "test_predictions.csv", row.names = FALSE)
