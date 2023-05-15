library(stringr)
library(data.table)
library(purrr)
library(tibble)
library(dplyr)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(forcats)
registerDoParallel(cores=28)
#Takes in pileup base (pileup4) and pileup quality (pileup5)
#and returns one column with just high quality bases
remove_low_quality_bases_from_Row = function(pileup4, pileup5){
    qualities = as.integer(charToRaw(as.character(pileup5)))
    highqualitybases = ""
    for(j in 1:length(qualities)){
      if(qualities[j] > 69){
        highqualitybases = paste0(highqualitybases,substr(pileup4, j, j))
      }
    }
    return(highqualitybases)
}
#Takes in refBase (pileup3), pileup base (pileup4), and pileup quality (pileup5)
#and returns the number of correct and total bases read at this position
remove_low_quality_bases_from_Row_and_make_coverage_and_correct = function(pileup3,pileup4, pileup5){
  qualities = as.integer(charToRaw(as.character(pileup5)))
  bases = as.character(str_split_fixed(as.character(pileup4),"",Inf))
  pileup54 = cbind(qualities,bases)
  pileup54 = pileup54[pileup54[,1]>69,2]
  empty_df = data.frame(c("A","G","C","T"))
  colnames(empty_df) = "pileup54"
  base_table = as.data.frame(merge(as.data.frame(table(pileup54)),empty_df, all = TRUE))
  rm(pileup54)
  base_table$Freq[is.na(base_table$Freq)] = 0
  coverage = sum(base_table$Freq)
  correct = base_table[base_table$pileup54==pileup3,2]
  cvacr = c(coverage, correct)
  return(c(coverage, correct))
}


# add_mutation_rate = function(pileup){
#   pileup$correctBaseCount = str_count(pileup$HQBases,pattern = pileup$RefBase)
#   pileup$ATCGCount = str_count(pileup$HQBases,pattern = 'A|T|C|G')
#   pileup$mutationRate = 1-pileup$correctBaseCount/pileup$ATCGCount
#   return(pileup[,c("Coordinate","ATCGCount","correctBaseCount","mutationRate")])
# }

# getMutationRates = function(pileup){
#   pileup = pileup[between(pileup$Coordinate,1101,1150),]
#   pileup$HQBases = mcmapply(remove_low_quality_bases_from_Row, pileup[,4], pileup[,5])
#   pileup = add_mutation_rate(pileup)
#   return(pileup)
# }


convertOutputListsToDataFrames = function(outputList){
  coverage_df = do.call(cbind,lapply(outputList,"[", ,"Coverage"))
  correct_df = do.call(cbind,lapply(outputList,"[", ,"Correct"))
  mutation_df = do.call(cbind,lapply(outputList,"[", ,"MutationRate"))
  rownames(coverage_df) = rownames(correct_df) = rownames(mutation_df) = rownames(outputList[[1]])
  colnames(coverage_df) = colnames(correct_df) = colnames(mutation_df) = 
    as.character(do.call(cbind,lapply(outputList,"[",1,"Sample")))
  return(list("Coverage" = coverage_df,
              "Correct" = correct_df,
              "Mutation_Rates" = mutation_df))
}
fileSetAndFileNames = function(path_to_pileups){
  fileSet = list.files(path = path_to_pileups,
                       pattern = "_pileup.txt", full.names = T)
  fileNames = str_split_fixed(fileSet,"/",Inf)
  fileNames = str_split_fixed(fileNames[,ncol(fileNames)],"_",Inf)[,c(2:4,7)]
  fileNames = apply(fileNames, 1, paste, collapse="_")
  return(cbind(fileSet,fileNames))
}
generateOutputMatrices = function(path_to_pileups, startingCoordinate, endingCoordinate){
  fSF = fileSetAndFileNames(path_to_pileups)
  fileSet = fSF[,1]
  fileNames = fSF[,2]
  correctBaseMatrix = coverageMatrix = mutationRateMatrix = data.frame(matrix(ncol = length(fileSet), nrow = 50))
  colnames(correctBaseMatrix) = colnames(coverageMatrix) = colnames(mutationRateMatrix) = fileNames
  outputList = foreach(i = 1:length(fileSet)) %dopar%{
    pileup = read.delim(fileSet[i], sep = " ", header = F,
                        col.names = c("SeqName","Coordinate","RefBase","ReadBases","BaseQuals","AlignQuals"))
    pileup = pileup[between(pileup$Coordinate,startingCoordinate,endingCoordinate),]
    cvacr = as.data.frame(
      t(mcmapply(remove_low_quality_bases_from_Row_and_make_coverage_and_correct, pileup[,3], pileup[,4], pileup[,5])))
    colnames(cvacr) = c("Coverage", "Correct")
    cvacr$MutationRate = 1-cvacr$Correct/cvacr$Coverage
    cvacr$Sample = fileNames[i]
    rownames(cvacr) = seq(startingCoordinate, endingCoordinate)
    
    cvacr
  }
  convertedList = convertOutputListsToDataFrames(outputList)
  coverage_df = as.data.frame(convertedList$Coverage)
  df_rows = row.names(coverage_df)
  rows.to.keep = !(df_rows %in% c("1108"))
  #rows.to.keep = rep(TRUE,40)
  coverage_df = coverage_df[rows.to.keep,]
  coverage_df["Number_of_Total_Bases_in_Sample",] = colSums(coverage_df)
  correct_df = as.data.frame(convertedList$Correct)
  correct_df = correct_df[rows.to.keep,]
  correct_df["Number_of_Correct_Bases_in_Sample",] = colSums(correct_df)
  mutation_df = as.data.frame(convertedList$Mutation_Rates)
  mutation_df = mutation_df[rows.to.keep,]
  mutation_df["Average_Mutation_Rate_in_Sample",] = 
    1-correct_df[nrow(correct_df),]/coverage_df[nrow(coverage_df),]
  return(list("Coverage" = coverage_df,
              "Correct" = correct_df,
              "Mutation_Rates" = mutation_df))
}
mutations_per_clone = function(){
  fSF = fileSetAndFileNames(path_to_pileups)
  fileSet = fSF[,1]
  fSF[,2] = fileNames
  correctBaseMatrix = coverageMatrix = mutationRateMatrix = data.frame(matrix(ncol = length(fileSet), nrow = 50))
  colnames(correctBaseMatrix) = colnames(coverageMatrix) = colnames(mutationRateMatrix) = fileNames
  outputList = foreach(i = 1:length(fileSet)) %dopar%{
    pileup = read.delim(fileSet[i], sep = " ", header = F,
                        col.names = c("SeqName","Coordinate","RefBase","ReadBases","BaseQuals","AlignQuals"))
    pileup = pileup[between(pileup$Coordinate,startingCoordinate,endingCoordinate),]
    cvacr = as.data.frame(
      t(mcmapply(remove_low_quality_bases_from_Row_and_make_coverage_and_correct, pileup[,3], pileup[,4], pileup[,5])))
    colnames(cvacr) = c("Coverage", "Correct")
    cvacr$MutationRate = 1-cvacr$Correct/cvacr$Coverage
    cvacr$Sample = fileNames[i]
    rownames(cvacr) = seq(startingCoordinate, endingCoordinate)
    
    cvacr
  }
  convertedList = convertOutputListsToDataFrames(outputList)
  coverage_df = as.data.frame(convertedList$Coverage)
  df_rows = row.names(coverage_df)
  rows.to.keep = !(df_rows %in% c("1108"))
  #rows.to.keep = rep(TRUE,40)
  coverage_df = coverage_df[rows.to.keep,]
  coverage_df["Number_of_Total_Bases_in_Sample",] = colSums(coverage_df)
  correct_df = as.data.frame(convertedList$Correct)
  correct_df = correct_df[rows.to.keep,]
  correct_df["Number_of_Correct_Bases_in_Sample",] = colSums(correct_df)
  mutation_df = as.data.frame(convertedList$Mutation_Rates)
  mutation_df = mutation_df[rows.to.keep,]
  mutation_df["Average_Mutation_Rate_in_Sample",] = 
    1-correct_df[nrow(correct_df),]/coverage_df[nrow(coverage_df),]
  return(list("Coverage" = coverage_df,
              "Correct" = correct_df,
              "Mutation_Rates" = mutation_df))
}
generateAndSaveOutput = function(output_location = "~/Documents/Programming/Rahul_Sequencing/JH4Seq/outputs/output_data/"){
  path_to_pileups = file.path(output_location,"pileups")
  outputMatrices = generateOutputMatrices(path_to_pileups, startingCoordinate = 1100, endingCoordinate = 1139)
  odl = file.path(output_location,"output_dfs")
  dir.create(odl, showWarnings = FALSE)
  write.csv(outputMatrices$Coverage, file.path(odl,"coverage.csv"))
  write.csv(outputMatrices$Correct, file.path(odl,"correct.csv"))
  write.csv(outputMatrices$Mutation_Rates, file.path(odl,"mutation_rates.csv"))
  return(outputMatrices)
}

outputMatrices = generateAndSaveOutput(output_location = "~/Documents/Programming/Rahul_Sequencing/JH4Seq_2/outputs/output_data")
coverage_df = outputMatrices$Coverage
correct_df = outputMatrices$Correct
mutation_df = outputMatrices$Mutation_Rates
colnames(mutation_df) = str_split_fixed(colnames(mutation_df),"_",Inf)[,2]

mutation_df_reordered <- as_tibble(mutation_df,rownames = "Position") %>% 
                                     dplyr::select("Position",c(ends_with("A"),ends_with("B"))) %>%
  dplyr::select("Position","1A":"9A","10A","1B":"8B")
mutation_rates = mutation_df_reordered[nrow(mutation_df_reordered),c(2:ncol(mutation_df_reordered))]
mutation_rates = as_tibble(cbind.data.frame(Names = names(mutation_rates), MutationRate = t(mutation_rates),
                           SampleType = c(rep("WT",4),rep("LDHA_KO",4),rep("AIDKO",2),
                                          rep("WT",3),rep("LDHA_KO",3),rep("AIDKO",2))))
mutation_rates$SampleType = as_factor(mutation_rates$SampleType)
mplot = ggplot(mutation_rates, aes(x = MutationRate, y = SampleType)) + coord_flip() + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(trans='log10') +
  geom_boxplot(aes(fill=SampleType)) + scale_fill_manual(values = c("#B3B3B3","#118040","skyblue")) + 
  geom_point(size=2) +
  stat_compare_means(comparisons = list(c("WT","AIDKO")), label = "p.signif")
pdf("/Volumes/ChaudhuriLab/Ryan/Rahul Sequencing/220705_JH4Seq_MutationRate.pdf")
mplot
dev.off()

t.test(as.numeric(mutation_rates[1,1:8]),as.numeric(mutation_rates[1,9:10]))
coverageMatrix$AIDKO_avg = rowMeans(coverageMatrix[grep("KO", names(coverageMatrix))])
coverageMatrix$WT_avg = rowMeans(coverageMatrix[grep("WT", names(coverageMatrix))])
coverageMatrix$WToverAIDKO = coverageMatrix$WT_avg/coverageMatrix$AIDKO_avg
View(coverageMatrix)
mutationRateMatrix$AIDKO_avg = rowMeans(mutationRateMatrix[grep("AID-KO", names(mutationRateMatrix))])
mutationRateMatrix$WT_avg = rowMeans(mutationRateMatrix[grep("WT", names(mutationRateMatrix))])
mutationRateMatrix$WToverAIDKO = mutationRateMatrix$WT_avg/mutationRateMatrix$AIDKO_avg
View(mutationRateMatrix)
