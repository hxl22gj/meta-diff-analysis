# Load required libraries
library(nlme)                    
library(lme4)                    
library(lattice)
library("stats")
library(dplyr)
library(ggplot2)
library(broom)
library(survival)
library(survminer)
library ("plyr")
library(purrr)

abundance<-"abundance.txt"
metadata<-"meta.txt"

# Data reading and preparation
data <- read.csv(abundance, sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
data1 <- as.data.frame(t(data))
#data arcsine square-root transformation
#numeric_columns <- sapply(data1, is.numeric)
#data1[numeric_columns] <- asin(sqrt(data1[numeric_columns]))
data1$ID <- rownames(data1)
data1 <- data1[, c("ID", setdiff(names(data1), "ID"))] 

meta <- read.csv(metadata, sep = "\t", check.names = FALSE, header = TRUE)
x <- merge(meta, data1, by = "ID")
x$cohort <- factor(x$cohort)
x$phenotype <- factor(x$phenotype)
x$`age-category` <- factor(x$`age_category`)
x$gender <- factor(x$gender)
x$antibiotics <- factor(x$antibiotics)
table(x$cohort)
table1<-function(x){return(table(x$`age_category`))}  #phenotype age_category gender antibiotics
ddply(x1, .(cohort),table1)
x$phenotype <- relevel(x$phenotype, ref = "HC")



# function for linear regression model  age+sex+antibiotics+CD
linear_regression_group <- function(x) {
  gene_name2<-c()
  gene_name<-c()
  pvalue_CD<-c()
  OR_CD<-c()
  LOW_CD<-c()
  HIGH_CD<-c()
  pvalue_UC<-c()
  OR_UC<-c()
  LOW_UC<-c()
  HIGH_UC<-c()
  pvalue_Male<-c()
  OR_Male<-c()
  LOW_Male<-c()
  HIGH_Male<-c()
  pvalue_Senior<-c()
  OR_Senior<-c()
  LOW_Senior<-c()
  HIGH_Senior<-c()
  number_1<-c()
  number_2<-c()
  for (i in 16:ncol(x))
  {
    g1=colnames(x)[i]
    tryCatch({ 
      gro <- x[,5]
      species<-x[,g1]
      aa <-na.omit(data.frame(gro,species))
      model <- lm(x[,i] ~ phenotype + age_category+gender, data = x) #+  antibiotics
      tidy_model <- tidy(model, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)
      gene_name2=c(gene_name2,g1)
      number<-ddply(aa,.(gro),nrow)
      number_1=c(number_1,number[1,2])
      number_2=c(number_2,number[2,2])
      pvalue_CD=c(pvalue_CD,as.numeric (tidy_model[2,5]))
      OR_CD=c(OR_CD,as.numeric (tidy_model[2,2]))
      LOW_CD=c(LOW_CD,as.numeric (tidy_model[2,6]))
      HIGH_CD=c(HIGH_CD,as.numeric (tidy_model[2,7]))
     }
      ,
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(data.frame(gene_name2,number_1,number_2,OR_CD,LOW_CD,HIGH_CD,pvalue_CD))
}

xIBD<-x[x$cohort!="Stinki",]
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="LSS-PRISM",] #geneder
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="CS-PRISM",] #Age+antibiotics
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="CD-china"|x$cohort=="NL-IBD",] #Age+geneder
xIBD<-xIBD[xIBD$phenotype!="UC",]

xIBD$phenotype <- factor(xIBD$phenotype)
xIBD$phenotype <- relevel(xIBD$phenotype, ref = "HC")

myrsa_cor<-ddply(xIBD, .(cohort),linear_regression_group)
write.table(myrsa_cor,file="all-CN-prism-species.lm.csv",sep = ",",row.names=FALSE)



# function for Linear Mixed-Effects -3 fixed cofactors
linear_mixed_group <- function(x) {
  results <- list()
  
  # Loop through the columns of interest
for (i in 16:ncol(x)) {
  tryCatch({
    col_name <- colnames(x)[i]
    x$y <- x[, i]
    fit.m3 <- lme(y ~ phenotype   + age_category + gender, # + antibiotics  
                  method = "ML", 
                  data = x, 
                  random = list(subject_accession = ~ 1, study_accession_database = ~ 1), # Corrected random effects formula
                  na.action = na.exclude)
    coef.fit.m3 <- summary(fit.m3)$tTable
    
    # Store results for each column in a list
    results[[col_name]] <- cbind(
      tax_name = col_name,
      diagnosisCD_value = coef.fit.m3[2, 1],
      diagnosisCD_t = coef.fit.m3[2, 4],
      diagnosisCD_p = coef.fit.m3[2, 5],
      #diagnosisUC_value = coef.fit.m3[3, 1],
      #diagnosisUC_t = coef.fit.m3[3, 4],
      #diagnosisUC_p = coef.fit.m3[3, 5],
      groupHC_mean = mean(x[x$phenotype == "HC", i], na.rm = TRUE),
      groupCD_mean = mean(x[x$phenotype == "CD", i], na.rm = TRUE) #,
      #groupUC_mean = mean(x[x$phenotype == "UC", i], na.rm = TRUE)
    )
  }, error = function(e){
    # Optionally, print the error message
    message("Error in column: ", col_name, "; Error: ", e$message)
  })
}
return (do.call(rbind, results))
}

xIBD<-x[x$cohort!="Stinki",]
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="LSS-PRISM",] #choose cohort based on cofactors
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="CS-PRISM",] #choose cohort based on cofactors
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD<-x[x$cohort=="CD-china"|x$cohort=="NL-IBD",] #choose cohort based on cofactors
xIBD<-xIBD[xIBD$phenotype!="UC",]

xIBD$phenotype <- factor(xIBD$phenotype)
xIBD$phenotype <- relevel(xIBD$phenotype, ref = "HC")

myrsa_cor<-ddply(xIBD, .(cohort),linear_mixed_group)
write.csv(myrsa_cor,file="linear_mixed_cn_cofactor_species.csv",sep='\t',row.names=FALSE)


# function for Linear Mixed-Effects -no fixed cofactors
linear_mixed_no <- function(x) {
  results <- list()
  
  # Loop through the columns of interest
  for (i in 16:ncol(x)) {
    tryCatch({
      col_name <- colnames(x)[i]
      x$y <- x[, i]
      fit.m3 <- lme(y ~ phenotype,
                    method = "ML", 
                    data = x, 
                    random = list(subject_accession = ~ 1, study_accession_database = ~ 1), # Corrected random effects formula
                    na.action = na.exclude)
      coef.fit.m3 <- summary(fit.m3)$tTable
      
      # Store results for each column in a list
      results[[col_name]] <- cbind(
        tax_name = col_name,
        diagnosisCD_value = coef.fit.m3[2, 1],
        diagnosisCD_t = coef.fit.m3[2, 4],
        diagnosisCD_p = coef.fit.m3[2, 5],
        #diagnosisUC_value = coef.fit.m3[3, 1],
        #diagnosisUC_t = coef.fit.m3[3, 4],
        #diagnosisUC_p = coef.fit.m3[3, 5],
        groupHC_mean = mean(x[x$phenotype == "HC", i], na.rm = TRUE),
        groupCD_mean = mean(x[x$phenotype == "CD", i], na.rm = TRUE) #,
        #groupUC_mean = mean(x[x$phenotype == "UC", i], na.rm = TRUE)
      )
    }, error = function(e){
      # Optionally, print the error message
      message("Error in column: ", col_name, "; Error: ", e$message)
    })
  }
  return (do.call(rbind, results))
}

xIBD<-x[x$cohort!="Stinki",]  
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD$phenotype <- factor(xIBD$phenotype)
xIBD$phenotype <- relevel(xIBD$phenotype, ref = "HC")

myrsa_cor<-ddply(xIBD, .(cohort),linear_mixed_no)
write.csv(myrsa_cor,file="linear_mixed_all_nocofactor_species.csv",sep='\t',row.names=FALSE)

#For test
fit.m3 <- lme(Adlercreutzia_equolifaciens ~ phenotype+ age_category + antibiotics + gender, 
              method = "ML", 
              data = xIBD, 
              random = list(subject_accession = ~ 1, study_accession_database = ~ 1), # Corrected random effects formula
              na.action = na.exclude)
coef.fit.m3 <- summary(fit.m3)$tTable


#K-W diff#
diff<-function(x){
gene_name2<-c()
gene_name<-c()
pvalue<-c()
qvalue<-c()
mean_1<-c()
median_1<-c()
mean_2<-c()
median_2<-c()
number_1<-c()
number_2<-c()


for (i in 16:ncol(x))
{
  g1=colnames(x)[i]
  gene_name=c(gene_name,g1)
  val<-as.numeric(x[,i])
  gro <- x[,5]
  aa <-na.omit(data.frame(gro,val))
  tryCatch({y=kruskal.test(val~gro, data=aa)
  pvalue=c(pvalue,y$p.value)
  gene_name2=c(gene_name2,g1)
  mean <- tapply(aa[,2],aa[,1],mean)
  median <- tapply(aa[,2],aa[,1],median)
  number<-ddply(aa,.(gro),nrow)
  number_1=c(number_1,number[1,2])
  number_2=c(number_2,number[2,2])
  a<-as.data.frame(mean)
  b<-as.data.frame(median)
  mean_1=c(mean_1,a[1,])
  mean_2=c(mean_2,a[2,])
  median_1=c(median_1,b[1,])
  median_2=c(median_2,b[2,])
  }, error=function(e){
    print(i)})
}  
qvalue=p.adjust(pvalue,'BH')
return (data.frame(gene_name2,mean_1, mean_2,number_1,median_1,median_2,number_2, pvalue,qvalue))
}

xIBD<-x[x$cohort!="Stinki",]  
xIBD<-xIBD[xIBD$phenotype!="UC",]
xIBD$phenotype <- factor(xIBD$phenotype)
xIBD$phenotype <- relevel(xIBD$phenotype, ref = "HC")
myrsa_cor<-ddply(xIBD, .(cohort),diff)

write.table(myrsa_cor,file="diff.all.CD.species.csv",sep = ",",row.names=FALSE)


#combine all results
getwd()
setwd("C:/Users/86156/Desktop/IBD-reanalysis/all.diff")
file_list <- list.files(pattern = "\\.csv$")
file_list
data_list <- list()

# Loop through the files
for (i in seq_along(file_list)) {
  # Read each file
  temp_data <- read.csv(file_list[i])
  
  # Create a new ID column by combining column 1 and column 2
  # Assuming column 1 and column 2 are the first two columns
  temp_data$ID <- paste(temp_data[,1], temp_data[,2], sep = "_")
  
  # Add the modified dataframe to the list
  data_list[[i]] <- temp_data
}

# Now, you can merge all data frames using the approach described previously
# For example, using a loop and merge() function in base R
merged_data <- data_list[[1]]  # Start with the first dataframe
for (i in 2:length(data_list)) {
  merged_data <- merge(merged_data, data_list[[i]], by = "ID", all = TRUE)
}

write.csv(merged_data,"merged_diff.csv")

