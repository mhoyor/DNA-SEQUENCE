# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Examination project in R 
# Author: 
# Email: 
# Submission date:
# Version: 1
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
install.packages("qdap")
library(tidyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(qdap)
library(reshape)

#Working directory
getwd()
#setwd("/Users/user/Documents/R_Examproject")Setting working directory

# Data Management ---------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: data import, variable assignment, dataset reorganisation (merge + long format),
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
getwd()
## Importing the datasets
dir_name <- "./Data/Original_data/"
file_name_pk <- "BPI889_PK_33.csv"
file_name_snp <- "BPI889_SNP_33.txt"

pk_path <- paste(dir_name, file_name_pk, sep = "")
snp_path <- paste(dir_name, file_name_snp, sep = "")

data_pk <- read.csv(pk_path, header = T, as.is = T, fill=T, sep=",")
data_snp <- read.table(snp_path, header = T, as.is = T, fill=TRUE, sep="")

#Making the PAT.ID column the first column and adding a index column in data_snp
data_snp <- cbind(ID = rownames(data_snp), data_snp)
rownames(data_snp) <- 1:nrow(data_snp)

#changing variable name in data_pk
colnames(data_pk)[colnames(data_pk) == "X"] = "ID"

#convert from wide to long 
tidy_pk <- gather(data_pk, TIME, Tco, Time.0.15.h:Time.24.h, factor_key=TRUE)

#convert Tco to numeric 
tidy_pk$Tco <- as.numeric(as.character(tidy_pk$Tco))

#extract the string in time to leave only numbers and convert to numeric 
tidy_pk[[6]] <- unlist(genXtract(tidy_pk[[6]], "Time.", ".h"))
tidy_pk$TIME <- as.numeric(as.character(tidy_pk$TIME))

#rename columns
tidy_pk <- rename(tidy_pk, c(WT = Weight..kg.,
                             HT = Height..cm., 
                             AGE = Age..yrs.))

#merge two datasets
data_all <- merge(tidy_pk, data_snp, by= "ID", all = TRUE)

# Variable calculations ---------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: calculation of body size measurement, categorization of body size measurement, 
# PK variable calculation
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  

#TBW calculation 
data_all$TBW[data_all$Sex== "M"]<- (2.447-0.09156*data_all$AGE) + 
                                    (0.1074*data_all$HT) + 
                                      (0.3362*data_all$WT)  #TBW calculation for Male(M)
data_all$TBW[data_all$Sex== "F"]<- (-2.097) + (0.1067*data_all$HT) + (0.2466*data_all$WT) #TBW calculation for Female(F)

#creating a new column 'CTBW' to categorize 'TBW'
data_all$CTBW <- 'NA' #assigning NA vale to the new column

data_all$CTBW[data_all$TBW <= 40]<- 0  #Convert BMI
data_all$CTBW[data_all$TBW > 40]<-1
data_all$CTBW <- factor(data_all$CTBW, levels = c(0, 1),  labels = c("Below", "Above"))


#calculating the PK variables (t1/2, Cmax, AUC)---------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#converting the time columns which are characters to numeric
data_pk[, 2:16] <- lapply(data_pk[, 2:16, drop = FALSE], as.numeric)

#Calculating the constant K from the data_pk data.frame for all individuals
#using timepoint 12hrs and 2hrs
#NOTE: K can't be negative
# k = (log(C) - log(C ))/10   where 6 is the difference in Time(t)

#calculating K
data_pk$K_el <- (log(data_pk$Time.2.h) - log(data_pk$Time.8.h))/6

#for missing values in the time frame 8hrs column
#using the next column after it to replace the missing value
data_pk$Time.8.h[is.na(data_pk$Time.8.h)] <- data_pk$Time.6.h[is.na(data_pk$K_el)]
data_pk$K_el <- (log(data_pk$Time.2.h) - log(data_pk$Time.8.h))/6

#Calculating AUC using the "data_AUC" data.frame
data_AUC <- data_pk[, colnames(data_pk)[c(1, 2:16, 21)]]
data_AUC[is.na(data_AUC)] <- 0
str(data_AUC)
#we have Time (t) in hours 0.15, 0.3, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 12, 24
#AUC = Tr(1) + Tr(2) +...+ Tr(n) +Cn /k
#where Tr(n)  =(C(n) +C(n-1) )/2*(t(n) -t(n-1))
#values in () are subsets
data_AUC$Tr1 <- ((data_AUC$Time.0.3.h + data_AUC$Time.0.15.h)/2) * (0.3 - 0.15)
data_AUC$Tr2 <- ((data_AUC$Time.0.5.h + data_AUC$Time.0.3.h)/2) * (0.5 - 0.3)
data_AUC$Tr3 <- ((data_AUC$Time.0.75.h + data_AUC$Time.0.5.h)/2) * (0.75 - 0.5)
data_AUC$Tr4 <- ((data_AUC$Time.1.h + data_AUC$Time.0.75.h)/2) * (1 - 0.75)
data_AUC$Tr5 <- ((data_AUC$Time.1.5.h + data_AUC$Time.1.h)/2) * (1.5 - 1)
data_AUC$Tr6 <- ((data_AUC$Time.2.h + data_AUC$Time.1.h)/2) * (2 - 1)
data_AUC$Tr7 <- ((data_AUC$Time.2.5.h + data_AUC$Time.2.h)/2) * (2.5 - 2)
data_AUC$Tr8 <- ((data_AUC$Time.3.h + data_AUC$Time.2.5.h)/2) * (3 - 2.5)
data_AUC$Tr9 <- ((data_AUC$Time.4.h + data_AUC$Time.3.h)/2) * (4 - 3)
data_AUC$Tr10 <- ((data_AUC$Time.5.h + data_AUC$Time.4.h)/2) * (5 - 4)
data_AUC$Tr11 <- ((data_AUC$Time.6.h + data_AUC$Time.5.h)/2) * (6 - 5)
data_AUC$Tr12 <- ((data_AUC$Time.8.h + data_AUC$Time.6.h)/2) * (8 - 6)
data_AUC$Tr13 <- ((data_AUC$Time.12.h + data_AUC$Time.8.h)/2) * (12 - 8)
data_AUC$Tr.last <- data_AUC$Time.12.h / data_AUC$K_el                 

data_AUC$AUC <- rowSums(data_AUC[ , c(18:31)], na.rm=TRUE)
toExclude <- colnames(data_AUC)[18:32]
tidy_AUC <- gather(data_AUC, Tmax, C, Time.0.15.h:Time.24.h, factor_key=TRUE, -ID)

#data_all <- gather(data_all, SNP, SNPcount, T134A:A990C, factor_key=TRUE)


## Re-arrange the data
tidy_AUC <- tidy_AUC[-c(3:16)]
### Reorder variables
tidy_AUC <- tidy_AUC[, c("ID", "C", "Tmax", "AUC", "K_el")]

### Reorder rows
tidy_AUC <- tidy_AUC[with(tidy_AUC, order(ID, -C)),]

pk_first <- tidy_AUC[!duplicated(data_all$ID), ]

pk_first$Tmax <- unlist(genXtract(pk_first$Tmax, "Time.", ".h"))
pk_first$Tmax <- as.numeric(as.character(pk_first$Tmax))

#calculate t1/2 as t0.5
pk_first$t0.5 <- log(2)/pk_first$K_el
tidy_AUC$t0.5 <- log(2)/tidy_AUC$K_el

#rename C to Cmax
colnames(pk_first)[colnames(pk_first) == "C"] = "Cmax"

#merge dataframes together

#calculating 
# Data Exploration --------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: numerical summary of PK variables, graphical assessment of 1) PK profiles,
# 2) PK variable correlations, 3)PK variable-SNP correlations, 
# 4) PK variable-body size measurement correlation with linear regression
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
#Display summary of PK varibles and Standard Deviation (SD)
cols_1 <- c("Cmax", "t0.5", "AUC")
summary(pk_first[cols_1])
sapply(pk_first[cols_1], sd)

#1) Creating a spaghetti plot
ggplot(data_all,aes(x=Tco,y=TIME,
              group=CTBW))+geom_line()

ggplot(data=data_all, aes(x=TIME, y=Tco, group=ID)) +
  geom_line() +geom_hline(yintercept = 0, linetype="dashed", color = "red")

#2) Scatterplots 
pairs(~ Cmax + t0.5 + AUC, data=pk_first)

pk_first <- pk_first[, c(1, 2, 4, 6, 3, 5)]

ggpairs(pk_first[, 2:4])

#3) Box whisker plots using SNP versus TIME
#Gather the SNP column using data_all, after dropping all duplicate data in data_all
data_snp <- merge(pk_first, data_snp, by= "ID", all = TRUE)
#data_snp$t0.5AUC = interaction(data_snp$AUC, data_snp$t0.5)

ggplot(aes(y=t0.5, x=AUC), data=data_snp) +
  geom_boxplot(aes(fill=SNP))


data_snp <- gather(data_snp, SNP, SNPcount, T134A:A990C, factor_key=TRUE, -ID)

ggplot(data = data_snp, aes(y = AUC, x = SNP, fill = SNP)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 16,
               outlier.size = 2, notch = FALSE) +
  labs(title = 'Box plot of T0.5 vs SNPs') +
  scale_fill_brewer(palette = "Blues") +
  theme_classic()

#3) Merge 
data_all <- data_all[with(data_all, order(ID, -Tco)),]
data_all <- data_all[!duplicated(data_all$ID), ]
data_all <- merge(data_all, pk_first, by= "ID", all = TRUE)

data_all <- data_all[ -c(6:7) ]
### Reorder variables
data_all <- data_all[, c("ID", "Tmax", "Cmax", "AUC", 
                         "t0.5", "K_el", "CTBW", "TBW", 
                         "WT", "HT", "AGE", "Sex", "T134A", 
                         "A443G", "G769C", "G955C", "A990C" )]

ggplot(data=data_all, aes(x=t0.5, y=TBW, shape = CTBW, color=CTBW))
                      + geom_point() + geom_abline(intercept = 40, slope = -5, color="red", 
                                                                   linetype="solid", size=1.5)

ylim = # Statistical testing -----------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: ANOVA of PK variables for SNPs, t-test of PK variable for body size measurement groups 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
#)AUC and SNPcount
data_snp$SNPcount <- factor(data_snp$SNPcount)
data_snp$AUC <- as.numeric(data_snp$AUC)
levels(data_snp$AUC)
anova_one_way_AUC <- aov(AUC~SNPcount, data = data_snp)
summary(anova_one_way_AUC)

TukeyHSD(anova_one_way_AUC) 

#) CMAX and SNPcount
anova_one_way_CMAX <- aov(Cmax~SNPcount, data = data_snp)
summary(anova_one_way_CMAX)

TukeyHSD(anova_one_way_CMAX)



#) T test for t0.5 on CTBW
t.test(t0.5 ~ CTBW, paired = FALSE, data = data_all)
