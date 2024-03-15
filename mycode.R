#problem 3.2
pheno.geno <- read.csv("gwas.pheno", sep = " ")
gwas.geno <- read.csv("gwas.geno", sep = " ")
overall_p <- function(my_model) { #extract p-value
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

A <- list() # 1st column of phenotype vs genotypes
for (i in 1:ncol(gwas.geno)){
  A[i] <- overall_p(lm(pheno.geno[,1] ~ gwas.geno[,i]))
}

B <- list()# 2nd column of phenotype vs genotypes
for (i in 1:ncol(gwas.geno)){
  B[i] <- overall_p(lm(pheno.geno[,2] ~ gwas.geno[,i]))
}

C <- list()# 3rd column of phenotype vs genotypes
for (i in 1:ncol(gwas.geno)){
  C[i] <- overall_p(lm(pheno.geno[,3] ~ gwas.geno[,i]))
}

D <- list()# 4th column of phenotype vs genotypes
for (i in 1:ncol(gwas.geno)){
  D[i] <- overall_p(lm(pheno.geno[,4] ~ gwas.geno[,i]))
}

# determine threshold below FWER
for (j in 1:382){
  if (A[j] <= 0.00013089005){
    print(A[j])
  }
  else{
    print('no')}
}

#result: 
#2nd phenotype 258th SNP = 2.2e-16
#3rd phenotype 119th SNP = 4.687267e-07 
#4th phenotype 43th SNP = 1.120269e-08ㄌ

#problem 3.3
hist(unlist(A), breaks = 25)
hist(unlist(B), breaks = 25)
hist(unlist(C), breaks = 25)
hist(unlist(D), breaks = 25)

#problem 3.4
second<- gwas.geno[,258]
third<- gwas.geno[,119]
fourth<- gwas.geno[,43]

boxplot(pheno.geno[,2]~second, xlab="genotype of 258th SNP", ylab = "2nd phenotype")
boxplot(pheno.geno[,3]~third, xlab="genotype of 119th SNP", ylab = "3rd phenotype")
boxplot(pheno.geno[,4]~fourth, xlab="genotype of 43th SNP", ylab = "4th phenotype")

#problem 4.2
test.geno <- read.csv("ridge.test.geno", sep = " ")
train.geno <- read.csv("ridge.training.geno", sep = " ")
test.pheno <- read.csv("ridge.test.pheno")
train.pheno <- read.csv("ridge.training.pheno")
#naming lambdas
lamb1 <- 0.001
lamb2 <- 2
lamb3 <- 5
lamb4 <- 8
#setting empty lists
A <- list()
B <- list()
C <- list()
D <- list()
#loops
for (j in 1:1000){
  inverse <- solve(t(train.geno[,j]) %*% train.geno[,j] + lamb1 * diag(1))
  A[j] <- inverse %*% t(train.geno[,j]) %*% train.pheno[,1]
}

for (j in 1:1000){
  inverse <- solve(t(train.geno[,j]) %*% train.geno[,j] + lamb2 * diag(1))
  B[j] <- inverse %*% t(train.geno[,j]) %*% train.pheno[,1]
}

for (j in 1:1000){
  inverse <- solve(t(train.geno[,j]) %*% train.geno[,j] + lamb3 * diag(1))
  C[j] <- inverse %*% t(train.geno[,j]) %*% train.pheno[,1]
}

for (j in 1:1000){
  inverse <- solve(t(train.geno[,j]) %*% train.geno[,j] + lamb4 * diag(1))
  D[j] <- inverse %*% t(train.geno[,j]) %*% train.pheno[,1]
}

# MSE = Var+bias^2. From the prompt, data is centered so bias = 0.
g<- as.matrix(test.geno)
a <- unlist(A)
aa <- matrix(a, nrow = 1000, ncol = 1)
AA <- list()
for(g in 1:1000){
  temp <- test.geno[,g] * aa[g,1]
  AA[g] <- temp
}

b <- unlist(B)
bb <- matrix(b, nrow = 1000, ncol = 1)
BB <- list()
for(g in 1:1000){
  temp <- test.geno[,g]*bb[g,1]
  BB[g] <- temp
}

c <- unlist(C)
cc <- matrix(c, nrow = 1000, ncol = 1)
CC <- list()
for(g in 1:1000){
  temp <- test.geno[,g] * cc[g,1]
  CC[g] <- temp
}

d <- unlist(D)
dd <- matrix(d, nrow = 1000, ncol = 1)
DD <- list()
for(g in 1:1000){
  temp <- test.geno[,g] * dd[g,1] 
  DD[g] <- temp
}

MSE1<- list()
  for (g in 1:99){
   temp<- (a[g] - test.pheno[g,1]) ^ 2
   MSE1[g] = temp
  }
RMSE1<- sum(unlist(MSE1)) / 99
MSE2<- list()
for (g in 1:99){
  temp<- (b[g] - test.pheno[g,1])^ 2
  MSE2[g] = temp
}
RMSE2<- sum(unlist(MSE2)) / 99
MSE3<- list()
for (g in 1:99){
  temp<- (c[g] - test.pheno[g,1]) ^ 2
  MSE3[g] = temp 
}
RMSE3<- sum(unlist(MSE3)) / 99

MSE4<- list()
for (g in 1:99){
  temp<- (d[g] - test.pheno[g,1]) ^ 2
  MSE4[g] = temp 
}
RMSE4<- sum(unlist(MSE4)) / 99

MSE<- c(RMSE1,RMSE2,RMSE3,RMSE4)
lambda<- c(lamb1,lamb2,lamb3,lamb4)

plot(lambda,MSE)
  
