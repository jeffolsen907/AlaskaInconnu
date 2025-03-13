#Data analysis for Inconnu from Alaska. The study examines genetic diversity and population structure.
#Use base R and tidyverse functions

# Open packages--------------------------------------------------------------------------------------------------------------
library(adegenet)
library(allelematch)
library(ape)
library(cowplot)
library(genepop)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(gtools)
library(hierfstat)
#library(pegas) #the pegas function Fst masks the genepop function Fst so did not open pegas.  The genepop Fst function has more options.
#library(PopGenKit) #PopGenKit has not been updated for R
library(poppr)
library(RColorBrewer)
library(readr)
library(reshape2)
library(tidyverse)
library(viridis)



# Set print to 200 lines-----------------------------------------------------------------------------------------------------
options(tibble.print_max = 200, tibble.print_min = 200)

# Open functions-------------------------------------------------------------------------------------------------------------

# Tibble to Genepop 
# Create a dataframe (tibble) with Individual and Population in first two columns followed by columns for each locus in two or three digit format ("0103", "001003").  
# The dataframe (tibble) must be organize as follows: 
#     First column contains individual identifiers and is named IndList.
#     Second column contains population identifiers and is named PopList.
#     The remaining columns contain the genotypes for each locus in either two (e.g., 0102) or three (e.g., 001002) digit code with no symbol separating alleles.
#     The locus columns are named for each locus
# The individual and population columns are both used here and one will be selected for use in GenePop (see subset step below) depending upon analysis.
Tib2Genepop <- function(x, y, z) {
  inputGP <- x
  # Creates blank dataframe with same number of columns as inputGP and no. columns - 1 rows.
  NewGP <-  data.frame(matrix(nrow = length(inputGP) - 1, ncol = length(inputGP))) 
  # Assign columns one and two of NewGP the column names of inputGP column 2 through end.  
  NewGP[,1:2] <-  names(inputGP[,2:length(inputGP)]) 
  # Assign the column names of inputGP to NewGP
  names(NewGP) <-  names(inputGP) 
  # Assign the first row, column 1 of NewGP the string "Genepop Input'.
  NewGP[1,1:2] <- 'Genepop Input' 
  # Create LociCount vector with one element that is the number of columns of NewGP minus 2.
  LociCount <-  length(NewGP)-2 
  # Add a comma to each character in the column IndList in the dataframe inputGP.  The comma is needed for genepop. 
  inputGP$IndList <- paste(inputGP$IndList,",", sep = "")
  # Add a comma to each character in the column PopList in the dataframe inputGP. The comma is needed for genepop. 
  inputGP$PopList <- paste(inputGP$PopList,",", sep = "")
  # Create a one row empty dataframe 
  Pop <-  as.data.frame(t(c(rep(NA, length(inputGP))))) 
  # Assign the column names of inputGP to Pop
  names(Pop) <-  names(inputGP) 
  # Assign the first row, column 1 and 2 of Pop the string 'Pop'
  Pop[1:2] <-  'Pop' 
  # Create a PopN vector with the names of each population from the the column inputGP$PopList in inputGP.
  PopN <-  unique(inputGP$PopList) 
  # for loop appends the datafram NewGP with the Pop dataframe and each population, in turn, from the dataframe inputGP.
  for(i in 1:length(PopN)) {
    tempGP <-  subset(inputGP, inputGP$PopList == PopN[i])
    NewGP <-  rbind(NewGP, Pop, tempGP)  
  }
  # Remove either column titled IndList or PopList from NewGP
  #NewGP <-  subset(NewGP, select = -c(PopList)) #Retaining the Indlist column for intra-watersded assignment test.
  NewGP_Ind <-  subset(NewGP, select = -c(PopList)) #Retaining the Poplist column for the TNWRRBT study
  NewGP_Pop <-  subset(NewGP, select = -c(IndList)) #Retaining the Indlist column for the TNWRRBT study.  
  # Write NewGP_Ind and NewGP_Pop to file.
  write.table(NewGP_Ind,y,quote=FALSE,row.names=F,col.names=F,sep="\t",na="")   
  write.table(NewGP_Pop,z,quote=FALSE,row.names=F,col.names=F,sep="\t",na="")    
}

#============================================================================================================================


# 1=PREPARE DATA=============================================================================================================

# 1.1-Named vector linking short Aggregation code to long Aggregation code.-------------------------------------------------
Agg <- c("EKOB0805", "ESEL0803", "ETAG0701", "EYUK0801", "EALA1401", "ESUL0702", "ESUL0802", "ESUL0806", "EINN1104")
AggCode <- c("A1", "B1", "B2", "C1", "C2", "C3", "C3", "C3", "C4")
AggCon <- setNames(AggCode, Agg)

# 1.2-Named vector linking new Inconnu locus names to old locus names.------------------------------------------------------
LocusRename <- read_csv(file = "InconnuLocusNameChange.csv")
LocusRenameCon <- setNames(LocusRename$New, LocusRename$Former)

# 1.3-Named vector linking Kuskokwim individuals to short Aggregation code.-----------------------------------------------
KuskoAggCodes <- c(TonzonaRiver = "D1", MiddleFork = "D2", BigRiver = "D3")
KuskoID <-  read_csv(file = "InconnuKuskokwimSamples.csv", skip_empty_rows = TRUE) %>% 
  mutate(Ind = str_replace_all(Ind, "[.]|-", "_")) %>% 
  mutate(Location = str_replace_all(Location, " ", "")) %>% 
  mutate(AggCode = recode(Location, !!!KuskoAggCodes))
KuskoAggCodesInd <- set_names(KuskoID$AggCode,KuskoID$Ind)

# 1.4-read csv input files for NWAK and YK and combine.---------------------------------------------------------------------
# NWAK
input_NWAK <- read_csv(file = "InconnuGenotypes_NWAK.csv", skip_empty_rows = TRUE) %>% 
  mutate(Aggregation = str_sub(Ind, start = 1, end = 8)) %>% 
  select(Ind, Aggregation, everything()) %>% 
  arrange(Ind) %>% 
  drop_na(Ind) %>% 
  filter(str_detect(Aggregation, "^EKOB08|^ESEL08|^ETAG07")) %>% 
  mutate(Ind = str_replace_all(Ind, "[.]|-", "_")) %>% 
  pivot_longer(cols = !c("Ind", "Aggregation"), names_to = "Allele", values_to = "Alleletype") %>% 
  mutate(Locus = str_sub(Allele, end = -3)) %>% 
  mutate(Locus = str_replace_all(Locus, "[.]|-", "")) %>% 
  mutate(Locus = recode(Locus, !!!LocusRenameCon)) %>% 
  filter(Locus %in% LocusRenameCon) %>% 
  mutate(Allele = paste("a", str_sub(Allele, start = -1), sep = "")) %>% 
  mutate(AggCode = recode(Aggregation, !!!AggCon, .default = Aggregation)) %>%
  select(Ind, Aggregation, AggCode, Locus, Allele, Alleletype) %>% 
  arrange(AggCode, Ind, Locus)
# YK
input_YK <- read_csv(file = "InconnuGenotypes_YK.csv", skip_empty_rows = TRUE) %>% 
  mutate(Aggregation = str_sub(Ind, start = 1, end = 8)) %>% 
  select(Ind, Aggregation, everything()) %>% 
  arrange(Ind) %>% 
  drop_na(Ind) %>% 
  filter(str_detect(Aggregation, "^EALA|^EBIG|^EINN|^EKUS|^ESUL|^EYUK")) %>% 
  mutate(Ind = str_replace_all(Ind, "[.]|-", "_")) %>% 
  pivot_longer(cols = !c("Ind", "Aggregation"), names_to = "Allele", values_to = "Alleletype") %>% 
  mutate(Locus = str_sub(Allele, end = -3)) %>% 
  mutate(Locus = str_replace_all(Locus, "[.]|-", "")) %>% 
  mutate(Locus = recode(Locus, !!!LocusRenameCon)) %>% 
  filter(Locus %in% LocusRenameCon) %>% 
  mutate(Allele = paste("a", str_sub(Allele, start = -1), sep = "")) %>% 
  mutate(AggCode = recode(Ind, !!!KuskoAggCodesInd, .default = Aggregation)) %>% 
  mutate(AggCode = recode(Aggregation, !!!AggCon, .default = AggCode)) %>% 
  select(Ind, Aggregation, AggCode, Locus, Allele, Alleletype) %>% 
  arrange(AggCode, Ind, Locus) 
# Combine NWAK and YK
input <- bind_rows(input_NWAK, input_YK) %>% 
  mutate_at("Alleletype", as.character)

# 1.5-Allele count.---------------------------------------------------------------------
# Identify rare alleles and possible scoring errors for correction.
AlleleCount <- input %>% 
  group_by(Alleletype) %>% 
  tally() %>% 
  arrange(n)
 
# 1.6-Create input file with 6 digit genotypes.----------------------------------------
inputLong <-  input %>% 
  mutate_at("Alleletype", ~str_pad(., width = 3, side = "left", pad = "0")) %>% 
  pivot_wider(names_from = Allele, values_from = Alleletype) %>% 
  mutate(Genotype = paste(a1, a2, sep = "/")) %>% 
  select(Ind, AggCode, Locus, Genotype) 

# 1.7-Confirm genotypes have 7 characters (e.g., 243/257)
countChar <- mutate(inputLong, char = nchar(Genotype))
distinct(countChar, char)

# 1.8-Identify individuals in input that are missing genotypes at more than 25% of loci.-------------------------------------
DropIndividuals <- group_by(inputLong,Ind,Genotype,Locus) %>% 
  tally() %>% 
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(Genotype == "000/000", freq > 0.25)

saveRDS(DropIndividuals, "DropIndividuals.rds")

# 1.9-Identify Loci in input that are missing genotypes at more than 25% of individuals.-------------------------------------
DropLoci_missingGT <- group_by(inputLong, Locus, Genotype, Ind) %>% 
  tally() %>% 
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(Genotype == "000/000", freq > 0.25)

saveRDS(DropLoci_missingGT, "DropLoci_missingGT.rds")

# 1.10-Remove Loci and Individuals from input that are missing genotypes at 25% or more individuals or loci----------------
inputLong <- inputLong %>% 
  filter(!(Ind %in% DropIndividuals$Ind)) %>% 
  filter(!(Locus %in% DropLoci_missingGT$Locus)) 

# 1.11-Subsample 100 from aggregations with more than 100 samples.---------------------------------------
inputLongSub <- inputLong %>% 
  pivot_wider(names_from = Locus, values_from = Genotype) %>% 
  group_by(AggCode) %>% slice_sample(n = 100) %>%   
  pivot_longer(cols = !c("Ind", "AggCode"), names_to = "Locus", values_to = "Genotype") %>% 
  arrange(AggCode, Ind) %>% 
  ungroup()


#============================================================================================================================


# 2=IDENTIFY AND REMOVE MATCHING GENOTYPES===================================================================================

# Use package allelematch to identify matching individuals (geneotypes) presumably sampled multiple times

# 2.1-Use the separate and unite functions to create separate columns in input for each allele at each locus.--------------------------
input_Allele <- inputLongSub %>% 
  separate(Genotype, c("A1", "A2"), sep = "/") %>% 
  pivot_longer(c("A1", "A2"), names_to = "Allele", values_to = "Alleletype") %>% 
  unite(Locus_Allele, Locus:Allele, sep = "_") %>% 
  pivot_wider(names_from = Locus_Allele, values_from = Alleletype)
  
# 2.2-Use the amDataset function in allelematch to produce an input dataset for allelematch from input_Allele.-------------
AM_input <- amDataset(input_Allele, missingCode = "000", indexColumn = 1, ignoreColumn = c(2))

# 2.3-Use the amUniqueProfile function to determine the optimum setting------------------------------------------------------ 
amUniqueProfile(AM_input, doPlot = TRUE)

# 2.4-Use the amUnique function to identify unique genotypes (individuals) and matching individuals.------------------------- 
# (used alleleMismatch = 2 from amUniqueProfile, allowing for ~ 5% allele mismatch)
PW_unique <- amUnique(AM_input, alleleMismatch = 2)
# Send PW_unique results to web browser as html (see allelematch documentation for amUnique function in R).
summary(PW_unique, html = TRUE)
# Send PW_unique results to local directory as csv file (see allelematch documentation for amUnique function in R).
summary(PW_unique, csv = "myUnique.csv")
# Define UniqueOut
UniqueOut <- read_csv(file = "myUnique.csv")

# 2.5-Further analyze individuals that were unclassified or with multiple matches.-------------------------------------------
# (see tutorial Example 3 in allelematch supplementary documentation).
# Preview individuals that were unclassified.
select(filter(UniqueOut, rowType == "UNCLASSIFIED"), c(uniqueGroup:matchIndex,Psib))
# If there are unclassified samples then use amPairwise to evaluate those samples (see tutorial Example 3 in allelematch supplementary documentation).
# The following code is used for unclassified samples but first see the tutorial Example 3 for details.
# This code was not used because no samples were unclassified
#Not used.  PW_unclassified <- amPairwise(PW_unique$unclassified, PW_unique$unique, alleleMismatch = 7)
# Send PW_unclassified results to web browser as html to inspect unclassified individual.
#Not used.  summary(PW_unclassified, html = TRUE)
# Preview individuals with multiple matches.
select(filter(UniqueOut, rowType == "MULTIPLE_MATCH"), c(uniqueGroup:matchIndex,Psib))
# Because individuals with multiple matches were present, use amPairwise to evaluate each match (see tutorial Example 3 in allelematch supplementary documentation).
# This code was not used because there were no multiple matches
#Not used.  PW_mismatch <- amPairwise(PW_unique$multipleMatches, PW_unique$unique, alleleMismatch = 6)
#Send PW_mismatch results to web browser as html to inspect individuals with multiple matches.
#Not used.  summary(PW_mismatch, html = TRUE)

# 2.6-Create vector of samples to remove because of matches.-----------------------------------------------------------------
# Use the pull function to pull out a column.
SampRemove <- UniqueOut %>% 
  filter(nUniqueGroup > 1) %>% 
  select(uniqueGroup:score) %>% 
  filter(!(uniqueIndex == matchIndex)) %>% 
  filter(!(rowType == "MULTIPLE_MATCH")) %>% 
  filter(!(Psib == "NA")) %>% 
  pull(matchIndex)

saveRDS(SampRemove, "SampRemove.rds")
SampRemove <- readRDS("SampRemove.rds")

# 2.7-Remove matching samples from input and input_pivot.--------------------------------------------------------------------
#Use filter function to remove samples in sampRemove from inputv5_pivot.
inputLongSub <- inputLongSub %>%
  filter(!(Ind %in% SampRemove))

# 2.8-TABLE S1: Genotype data for analysis--------------------------------------
TableS1 <- inputLongSub %>%
  mutate(Ind = factor(Ind, levels = unique(Ind))) %>% 
  group_by(Ind) %>% 
  mutate(ID = cur_group_id()) %>% 
  select(ID, everything())

# 2.8.1-Write TABLE S1 to directory-----------------------------------
write.table(TableS1,"TableS1_AlaskaInconnu.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",") 

# 2.9-Pivot input so rows are individual genotypes and columns are loci.----------------------------------------------------
inputWideSub <- inputLongSub %>% 
  pivot_wider(names_from = Locus, values_from = Genotype) 

# 2.10-Save long and wide input as rds files--------------------------------
saveRDS(inputLongSub, "inputLongSub.rds")
saveRDS(inputWideSub, "inputWideSub.rds")

#============================================================================================================================


# 3=PREPARE INPUT FOR FIRST LEVEL ANALYSIS USING HIERFSTAT, GENEPOP==========================================================

# 3.1-START HERE to avoid resubsampling step 1.5 above.---------------------------------------------------------------------- 
inputWideSub <- readRDS("inputWideSub.rds")
inputLongSub <- readRDS("inputLongSub.rds")

# 3.2-Determine sample size by Collection code and River code.--------------------------------------------------------------- 
group_by(inputWideSub, AggCode) %>% tally()

# 3.3-Simplify locus names in input_pivot from original (long) to short code (e.g., Loc1, Loc2...)-------------------------
# Define vector of original (long) locus names
LocNameLong <- names(select(inputWideSub,-(Ind:AggCode)))
# Define vector of short vector names
LocNameCode <- paste("Loc",rep(c(1:length(LocNameLong))),sep="")
# Define vector associating long locus names to codes. This object is only for reference, not used in subsequent code.
LocNameCon <- setNames(LocNameCode,LocNameLong)
saveRDS(LocNameCon, "LocNameCon.rds")
# Use rename_at function to rename multiple columns (all locus names). The tilda identifies function
inputWideSub <- inputWideSub %>% 
  rename_at(vars(all_of(LocNameLong)), ~ LocNameCode)

# 3.4-Prepare input data for use in hierfsat packages------------------------------------------------------------------------
# Vector of populations aggregation code (AggList) for adegenet input
PopList <- inputWideSub$AggCode
# Vector of individuals (for adegenet input)
IndList <- inputWideSub$Ind
# Vector of loci
LociList <- LocNameCode
# Use 'select' to remove columns without genotypes.
input4AD <- select(inputWideSub, -(Ind:AggCode)) 
# Convert to adegenet input
AD <- df2genind(input4AD,sep = '/', ncode=4, ind.names=IndList, pop=PopList, NA.char="000/000", type="codom") 
# Convert to hierfsat input
HF <- genind2hierfstat(AD) 

# 3.5-Prepare input data for use in genepop package--------------------------------------------------------------------------  
# Create a dataframe (tibble) with Individual and Population in first two columns followed by columns for each locus in three digit format ('152160').  
# The individual and population columns are both used here and one will be selected for use in GenePop depending upon analysis.
input2GP <- pivot_longer(inputWideSub, c("Loc1":"Loc20"), names_to = "Loci", values_to = "Genotype") %>% 
  #remove slash "/" separating alleles
  mutate(Genotype = str_replace_all(Genotype, "/", "")) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  #two columns (Ind and AggCode) and all loci
  select(Ind, AggCode, Loc1:Loc20) %>% 
  rename(IndList = Ind, PopList = AggCode)

# Use function Tib2Genepop to convert input2GP tibble to genepop input files for individual and populations.
Tib2Genepop(input2GP, "AKInconnu_GP_Ind_Final.txt", "AKInconnu_GP_Pop_Final.txt")


#============================================================================================================================


# 4-FIRST LEVEL ANALYSIS=====================================================================================================
inputWideSub <- readRDS("inputWideSub.rds")

# 4.1-Use hierFstat to compute Ho, Hs, Fis, Fst for all loci and populations-------------------------------------------------
B <-  basic.stats(HF) 
PopHs_mean <- round(colMeans(B$Hs,na.rm=TRUE),3)
WC <- wc(HF)

# 4.2-Use Poppr to estimate the Shannon-Weiner index for each locus----------------------------------------------------------
# Used to select most informative locus of locus pair if linked
SW_Nuc <- locus_table(AD, index = 'shannon')  
# Coerce SW_Nuc to tibble
SW_Nuc <- as_tibble(rownames_to_column(as.data.frame(SW_Nuc), var = "Locus")) 
# Use "filter" in dplyr package to remove last row (says "mean" in Locus column) showing mean SW estimates for all loci
SW_Nuc <- filter(SW_Nuc, Locus != "mean") 

#4.3-Use Genepop to test HWP and GD----------------------------------------------------------------------------------------
locinfile <-  "AKInconnu_GP_Pop_Final.txt"
basic_info(locinfile, outputFile = "AKInconnu_GP_Pop_Out.txt")
# Test HWP (Hardy_Weinberg proportions) for all loci and collections
test_HW(locinfile, outputFile = "AKInconnu_GP_HW.txt", enumeration = TRUE, dememorization = 10000, batches = 100, iterations = 5000)
# Test GD (gametic disequilibrium) for all pairs of loci across collections
test_LD(locinfile, outputFile = "AKInconnu_GP_GD.txt", dememorization = 10000, batches = 100, iterations = 5000)

# 4.4-Evaluate results of HWP tests from genepop-----------------------------------------------------------------------------

# First, test assumption of global HWP across all loci and collections
# use cumulative binomial probability distribution to derive 95% CI to the number of table-wide significant tests (p < 0.05). See Waples 2015. 
# Modify the test_HW output using python script (outside of R) to create a locus x collection matrix of p-values. 
# Read data file matrix (created using python script) of HWP p-values 24 pops. The col_types argument explicitly specifies the column data type (c=character, d=double).
inputHWP <- read_csv("TNWRRBT_GP_HW_Pvalues.txt", col_types = paste("c",strrep("d",24), sep=""), na = "-") 
# Create tibble from inputHWP that sums the number p-values < 0.05 for each locus.
LocLTalfa <- tibble('Loc' = inputHWP$Loc, 'NLocLTalfa' = rowSums(inputHWP[-1] < 0.05, na.rm = TRUE)) 
# Create tibble from inputHWP that sums the number p-values < 0.05 for each population.
PopLTalfa <- tibble('Pop' = names(inputHWP[-1]), 'NPopLTalfa' = colSums(inputHWP[-1] < 0.05, na.rm=TRUE)) 
# Table-wide number of HWP tests counts number of tests (p-values) in inputHWP (without counting unmeaningful tests denoted by NA) 
NoTestsHWP <- inputHWP %>% 
  select(-(Loc)) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  mutate(Tot = rowSums(.)) %>% 
  select(Tot)
NoTestsSigHWP <- sum(LocLTalfa$NLocLTalfa)
binomtestHWP <- dbinom(1:200, size = NoTestsHWP$Tot, prob = 0.05)
binomtestHWP <- round(binomtestHWP,4)
# cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumHWP <- cumsum(binomtestHWP) 

# Second, test assumption of HWP for each collection across loci and each locus across collections

# Population-level test of HWP across loci. Genepop data file summarized using python script.
PopSigTest <-  read_csv("TNWRRBTGenepopHWP_PopSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
PopBinomTest <- dbinom(0:length(LocLTalfa$Loc), size = length(LocLTalfa$Loc), prob = 0.05)
PopBinomTestTable <- tibble('NumSigLoci' = 0:length(LocLTalfa$Loc), 'PopBinomExp' = round(PopBinomTest * length(PopLTalfa$Pop), 2), 'PopBinomObs' = as.vector(table(factor(PopSigTest$ST, levels = 0:length(LocLTalfa$Loc)))))
PopBinomTestTableTidy <- gather(PopBinomTestTable, 'PopBinomExp', 'PopBinomObs', key = 'ObsOrExp', value = 'NumAggregations')
# Plot results for collections across loci
PopBinomTestPlot <- ggplot(PopBinomTestTableTidy) + 
  geom_bar(aes(NumSigLoci, NumAggregations, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Aggregations") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,55,5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2, 0.8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5))

# Locus-level test of HWP across populations. Genepop data file summarized using python script.
LocSigTest <-  read_csv("TNWRRBTGenepopHWP_LocSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
LocBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
LocBinomTestTable <- tibble('NumSigAggregations' = 0:length(PopLTalfa$Pop), 'LocBinomExp' = round(LocBinomTest * length(LocLTalfa$Loc), 2), 'LocBinomObs' = as.vector(table(factor(LocSigTest$ST, levels = 0:length(PopLTalfa$Pop)))))
LocBinomTestTableTidy <- gather(LocBinomTestTable, 'LocBinomExp', 'LocBinomObs', key = 'ObsOrExp', value = 'NumLoci')
# Plot results for loci across collections
LocBinomTestPlot <- ggplot(LocBinomTestTableTidy) + 
  geom_bar(aes(NumSigAggregations, NumLoci, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Loci") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,35,5), limits = c(0,35)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2,0.8), 
        axis.title = element_text(size=12), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5)) 

# 4.5-FIGURE S1: Binomial distribution of HWP results for loci and populations-----------------------------------------------
# Print Figure S1 to directory
tiff('FigS1_TNWRRBTPaper.tiff', height = 6.5, width = 8.5, units = 'in', compression = 'none', res = 300) 
ggarrange(PopBinomTestPlot, LocBinomTestPlot, ncol=1, nrow=2, align = 'v')
dev.off()

# 4.6-Evaluate results of GD tests from genepop------------------------------------------------------------------------------

# read data file matrix of GD p-values for all loci pairs for all 24 populations.  Data from genepop summarized using python script.
inputGD <- read_csv("TNWRRBTGenepopGD_Pvalues.txt", col_types = paste('c', strrep('d',24), strrep('i',2), sep=''), na='NA') 
NoTestsGD <- length(inputGD$LocXLoc)*length(PopLTalfa$Pop)
# Table-wide number of GD tests counts number of tests (p-values) in inputGD (without counting unmeaningful tests denoted by NA) 
NoTestsGD <- inputGD %>% 
  select(-c(LocXLoc,ST,BST)) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  mutate(Tot = rowSums(.)) %>% 
  select(Tot)
NoTestsSigGD <- sum(inputGD$ST)
binomtestGD <- dbinom(1:10000, size = NoTestsGD$Tot, prob = 0.05)
binomtestGD <- round(binomtestGD,4)
# cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumGD <- cumsum(binomtestGD) 
arrange(select(inputGD, LocXLoc, ST), desc(ST))

GDBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
GDBinomTestTable <- tibble('NumSigAggregations' = 0:length(PopLTalfa$Pop), " Expected" = round(GDBinomTest * length(inputGD$LocXLoc), 2), " Observed" = as.vector(table(factor(inputGD$ST, levels = 0:length(PopLTalfa$Pop)))))
GDBinomTestTableTidy <- gather(GDBinomTestTable, " Expected", " Observed", key = 'ObsOrExp', value = 'NumLocPairs')

# plot of count of locus pairs versus number of significant collections
GDBinomTestPlot <- ggplot(GDBinomTestTableTidy) + 
  ggtitle("Locus Pairs") +
  geom_bar(aes(NumSigAggregations, NumLocPairs, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,850)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.3,0.8),
        axis.title = element_text(size=12), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=10),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title = element_text(size = 14, vjust = -10, hjust = 0.5)) +
  annotate("text", x = c(13), y = 75, label = c("OMS00154 x\nOmy_arp_630"), size = 4) +
  annotate("segment", x = 13, xend = 13, y = 40, yend = 5, colour = "black", size=1, alpha=1, arrow=arrow(length = unit(0.3, "cm")))

# 4.7-FIGURE S2: GD results showing number of locus pairs versus number of significant collections---------------------------
# Print Figure S2 to directory.
tiff('FigS2_TNWRRBTPaper.tiff', height = 5.5, width = 8.5, units = 'in', compression = 'none', res = 300) 
print(GDBinomTestPlot)
dev.off()


# 4.8-Compute Mean MAF per locus---------------------------------------------------------------------------------------------
input <- readRDS("input.rds")
input_pivot <- readRDS("input_pivot.rds")

MAF <- pivot_longer(input_pivot, -(Ind:RiverNo), names_to = "Locus", names_ptypes = "d", values_to = "Genotype") %>% 
  separate(Genotype, into=(c("A1","A2"))) %>% 
  select(Locus, A1, A2) %>% 
  pivot_longer(c("A1","A2"), names_to = "Allele", values_to = "Score") %>% 
  filter(!(Score == "00")) %>% 
  group_by(Locus, Score) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(freq == min(freq))

# 4.9-TABLE S3: descriptive statistics of loci tested in first level analysis-----------------------------------------------
# Tibble of locus short and long names.
LocusName <- tibble(Locus = names(select(input_pivot, -(Ind:RiverNo))), LocusLong = pull(distinct(input,Locus)))
# Create dataframe of descriptive stats for all loci from hierFstat and add mean MAF computed above
ByLoc <-  B$perloc %>% 
  rownames_to_column(var = "Locus") %>% 
  # coerce ByLoc to tibble
  as_tibble() %>% 
  # create columns of statistics for each locus
  mutate(Ref = c(rep("a", 55)), FisWC = WC$per.loc$FIS, FstWC = WC$per.loc$FST, SWH = SW_Nuc$H) %>% 
  # move Loc17 (OMS00154) in row 18 to bottom. This locus is linked to Loc30 (Omy_arp_630) and was removed for final analysis.
  slice(c(1:16,18:55,17)) %>% 
  inner_join(MAF, by = "Locus") %>% 
  inner_join(LocusName, by = "Locus") %>% 
  rename("MAF" = "freq") %>% 
  select(LocusLong, Ref, MAF, Ho, Hs, FisWC, FstWC, SWH) %>% 
  rename("Locus" = "LocusLong") %>% 
  # add a rwo of column means for diversity stats 
  bind_rows(summarise_all(., ~ifelse(is.numeric(.), mean(.), "mean"))) %>% 
  # use "mutate_if" in dplyr package to round number to 3 decimals
  mutate_if(is.double, round, 3) %>% 
  mutate_at(vars(MAF:SWH), as.character) %>% 
  mutate_at(vars(MAF:SWH), ~str_pad(.,5,"right","0")) %>% 
  # add loci that were excluded from the data set because mMAF was less than 0.01
  full_join(DropLoci_lowMAF, by = "Locus") %>%  
  # add locus that was exluded from the set because the assay failed in more than 25% of individuals.
  full_join(DropLoci_missingGT, by = "Locus") %>% 
  select(Locus, MAF, Ho, Hs, FisWC, FstWC, SWH)

# write Table S3 to directory.
write_csv(ByLoc, "TableS3_TNWRRBT.txt", col_names = TRUE) 

# 4.10-TABLE S4: descriptive statistics of aggregations tested in first level analysis-----------------------------------------------
ByPopHo <- round(colMeans(B$Ho,na.rm=TRUE),3) 
ByPopHs <- round(colMeans(B$Hs,na.rm=TRUE),3) 
ByPopFis <- round(colMeans(B$Fis,na.rm=TRUE),3) 
# Create dataframe of descriptive stats for all aggregations from hierFstat
ByPop <- tibble(Ho = ByPopHo, Hs = ByPopHs, FisWC = ByPopFis) %>% 
  rbind(c(mean(ByPopHo), mean(ByPopHs), mean(ByPopFis))) %>% 
  mutate_if(is.double, round, 3) %>% 
  mutate(Acode = c(names(ByPopHo), "mean")) %>% 
  select(Acode, everything()) %>% 
  mutate_at(vars(Ho:FisWC), as.character) %>% 
  mutate_at(vars(Ho:FisWC), ~str_pad(.,5,"right","0"))

# write Table S4 to directory.
write_csv(ByPop, "TableS4_TNWRRBT.txt", col_names = TRUE) 

# 4.11-Remove Locus 17 from data set and sort by River Number.----------------------------------------------------------------
# GD test indicated Loci 17 and 29 (OMS00154 and Omy_arp_630) are linked.
# Removed locus 17 because it has the lower overall mMAF.
inputFinal <- input_pivot %>% 
  select(-(Loc17))
input <- input %>%
  filter(!(Locus == "OMS00154"))

# 4.12-Save inputFinal and input as rds files for use in next steps.---------------------------------------------------
saveRDS(inputFinal, "inputFinal.rds")
saveRDS(input, "input.rds")


#============================================================================================================================

# Structure input
write.struct(HF, ilab = IndList, pop = PopList, MARKERNAMES = TRUE, MISSING = -9, fname = "AKSheefish.str")





input<-read_csv("AlaskaInconnu.csv") %>%
  pivot_longer(c("BWF2_1":"Sle014_2"), names_to = "Locus_Allele", values_to = "Alleletype") %>% 
  # use the mutate function in dplyr to replace periods and dashes with underscore in strings in the columns "FishN" and "Locus"
  mutate(FishN = str_replace_all(FishN, "[.]|-", "_")) %>% 
  mutate(Locus_Allele = str_replace_all(Locus_Allele, "[.]|-", "")) %>% 
  mutate(Allele = paste("A", str_sub(Locus_Allele, start = -1), sep = "")) %>% 
  mutate(Locus = str_sub(Locus_Allele, end = -3)) %>% 
  # remove YRap collection because it is a possible mixture and not a confirmed spawning area
  filter(PopN != "YRap") %>% 
  # remove YTan collection because the sample size (9) is too small.
  filter(PopN != "YTan") %>% 
  # add column PopCode for each collection using named vector PopCon
  mutate(AggCode = recode(PopN, !!!AggCon)) %>% 
  select(FishN, PopN, AggCode, Locus, Allele, Alleletype) %>% 
  rename(Ind = FishN) %>% 
  arrange(AggCode) %>% 
  pivot_wider(names_from = Allele, values_from = Alleletype) %>% 
  mutate(Genotype = paste(A1, A2, sep = "/")) %>% 
  select(Ind, AggCode, Locus, Genotype) %>% 
  mutate(Genotype = str_replace_all(Genotype, c("0/0" = "000/000")))



# FIGURE 7 for Brown et al. --MDS of pairwise Fst----------------------------------------

inputM = read.table("AKSheeFishpwwc_poster.txt", header = T)
inputM[is.na(inputM)] = 0
Pops = c("NWK", "NWS", "YIN", "YSL", "YAL", "YTN", "YFL", "KBI", "KMF", "KTO")
rownames(inputM)<-Pops
colnames(inputM)<-Pops
FstMDS<-cmdscale(inputM,eig=FALSE,k=2) |> 
  as_tibble(.name_repair = "unique") |> 
  rename(Axis1 = 1, Axis2 = 2) |> 
  mutate(Pop = Pops) |> 
  select(Pop, everything())

Fig7Plot <- ggplot(data = FstMDS, aes(x = Axis1, y = Axis2)) +
  #Plot type  
  geom_point(stroke = 1, aes(shape = Pop)) +
  geom_text_repel(aes(label = Pop), force = 40, force_pull = 50, size = 3.0) +
  #layout
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=11), axis.text = element_text(size=9), legend.position = "none") +
  #overrides
  scale_shape_manual(values = c(rep(17,3), rep(19,2), rep(15,5))) +
  #scale_shape_manual(values = c(17, 15, 18, 2, 0, 5, 6, 2, 5, 5)) +
  scale_y_continuous(limits = c(-0.03,0.03), breaks = seq(-0.03,0.03,0.01), name = "Axis 2") +
  scale_x_continuous(limits = c(-0.03,0.03), breaks = seq(-0.03,0.03,0.01), name = "Axis 1") +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash')


# Print FIGURE 7 to directory
ggsave(filename = "Fig7_Brown_etal.tiff", plot = Fig7Plot, device = "tiff",
       height=2.5, width = 3.0, units="in", dpi=600)


