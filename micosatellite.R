#!/usr/bin/env Rscript

#Code to install packages if necessary, and read them with library function
required_packages <- c("ggplot2","readxl","tidyverse","adegenet","genepop","hierfstat","here","mapdata",
                       "mapplots","data.table","colorspace","poppr","knitr")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}


# --- Obtain data --- #

#set your session in current directory of R file
here::set_here()
getwd()


# import the data within specified sheet
data <- read_excel("SexAsex4Ammar.xlsx", 
                  sheet = "Ammar")

#covert data to data frame
data <- as.data.frame(data)


# --- Data map visualization --- #
# summarize data by population and location
pop_data <- data %>%
  group_by(Pop, Long, Lat,Profil) %>%
  summarize(n = n()) %>%
  mutate(percent = n / sum(n)) %>%
  select(Pop, Long, Lat,Profil, n, percent)

#draw map Representative the popular individual across France map within reproductive mode
# The area of the France Region;
xlim <- c(-4.3,9.1)
ylim <- c(41, 51)

# Creates an xyz object for use with the function draw.pie
xyz <- make.xyz(pop_data$Long, pop_data$Lat, pop_data$n, pop_data$Profil)

# Colors used
col <- c("#003366", "#CCCC66", "#CC3366")

# The plot of the pie chart above the map
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),
    mgp = c(2.5, 0.5, 0), family = "Arial")
basemap(xlim =  c(-4.3,9.1), ylim =  c(41, 51), bg = "white",
        main = "Distribution of population of popular rust in France")
map("france", fill=FALSE, col="light blue", xlim = xlim, ylim = ylim, add = TRUE)
draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.3, col = col)
legend.pie(-5, 42, labels = c("Asex", "NA", "Sex"), 
           radius = 0.2, bty = "n", col = col, cex = 0.8, label.dist = 1.5)




#draw map Representative the popular individual across France map within asexual linage

# summarize data by population and location
pop_data_mll <- data %>%
  group_by(Pop, Long, Lat,Profil,Mll) %>%
  summarize(n = n()) %>%
  mutate(percent = n / sum(n)) %>%
  select(Pop, Long, Lat,Profil,Mll, n, percent)

# Create new column based on condition make all sex mode as 1000 linage
pop_data_mll$Mll_new <- ifelse(pop_data_mll$Mll >= 1000, 1000, pop_data_mll$Mll)

# Get unique values in Mll_new column
mll_unique <- unique(pop_data_mll$Mll_new)

# Generate color palette with number of colors equal to number of unique values in Mll_new
colors <- qualitative_hcl(length(mll_unique))

xyz_new <- make.xyz(pop_data_mll$Long, pop_data_mll$Lat, pop_data_mll$n, pop_data_mll$Mll_new)

# The plot of the pie chart above the map
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),
    mgp = c(2.5, 0.5, 0), family = "Arial")
basemap(xlim =  c(-4.3,9.1), ylim =  c(41, 51), bg = "white",
        main = "Distribution of population of popular rust in France")
map("france", fill=FALSE, col="light blue", xlim = xlim, ylim = ylim, add = TRUE)
draw.pie(xyz_new$x, xyz_new$y, xyz_new$z, radius = 0.3, col = colors)
legend.pie(-4.3, 41, labels = mll_unique , 
           radius = 0.2, bty = "n", col = colors, cex = 0.8, label.dist = 1.3)


# --- Genetic diversity analysis --- #
#create table of genotype

# Select columns with "Mlp" in the name and the first column as isolate id
genotype_cols <- c("Isolate", grep("Mlp", names(data), value = TRUE))
genotype_data <- data[, genotype_cols]

#make isolate id as column names 
rownames(genotype_data) <- genotype_data[,1]
#delete the first column
genotype_data$Isolate = NULL


# Convert to genind object
data_GenInd <- df2genind(
  X = genotype_data,   #data.frame containing allele data only 
  sep = NULL,
  ncode = 3,           #an optional integer giving the number of characters used for coding one genotype at one locus.
  ind.names = rownames(genotype_data),  # individuals names
  loc.names = colnames(genotype_data),  # markers names
  pop = data$Pop,                       # giving the population of each individual
  NA.char = "999",                      # string corresponding to missing allele 999 or 999999
  ploidy = 2,
  type = "codom",                       #codom' stands for 'codominant' (e.g. microstallites, allozymes)
  strata = NULL,
  hierarchy = NULL,
  check.ploidy = getOption("adegenet.check.ploidy")
)


# ca marche !!
summary(data_GenInd)
# additional information
data_GenInd@loc.n.all
data_GenInd@all.names
levels(data_GenInd@pop) #110 pops

Nb_Pop = length(levels(data_GenInd@pop)) # number of population
is.numeric(Nb_Pop)




# Computation of some population genetics indices
# different vectors that we want to calculate
N <- vector(mode ="integer", length = Nb_Pop)
GsurN<- vector(mode ="numeric", length = Nb_Pop)
Ho <- vector(mode ="numeric", length = Nb_Pop)
Hs <- vector(mode ="numeric", length = Nb_Pop)
Fis <- vector(mode ="numeric", length = Nb_Pop)
rbarD<- vector(mode ="numeric", length = Nb_Pop)



Pop <- levels(data_GenInd$pop)

# Using Poppr to retrieve the first indices N, G/N et rbarD
Table_PPR <- poppr(data_GenInd) # marche bien

rbarD <- Table_PPR$rbarD[1:Nb_Pop]
N <- Table_PPR$N[1:Nb_Pop]
GsurN <- (Table_PPR$MLG[1:Nb_Pop]-1)/(Table_PPR$N[1:Nb_Pop]-1)


#calcul de la proba de DL
ProbaLD <- vector(mode ="numeric", length = Nb_Pop)
for (i in 1:Nb_Pop) {
  i=1
  Temp_Sample <- popsub(data_GenInd, Table_PPR$Pop[i])
  ProbaLD[i] = ia(Temp_Sample, sample = 999, plot = F)[4] 
}

#Use of Fstat to calculate the other indices HO, HE, Fis and Ar

data_Fstat <-genind2hierfstat(data_GenInd)

# Boucles qui permet de calculer les indices souhaités pour chaque population (et pas la globalité)
a =0
for (i in levels(data_Fstat$pop) ){
  Poptmp <- data_GenInd[which(data_GenInd$pop==i),]
  #Poptmp  <-genind2hierfstat(Poptmp)
  fstat_basic_Temporel <- basic.stats(Poptmp)
  a = a+1
  Ho[a]<-fstat_basic_Temporel$overall["Ho"]  
  Hs[a]<-fstat_basic_Temporel$overall["Hs"]  
  Fis[a]<-fstat_basic_Temporel$overall["Fis"] 
}



# calculation of the allelic richness
Obj_Ar <- allelic.richness(data_Fstat[which(data_Fstat[,"pop"]!="NA"),])
Ar_per_loc <- Obj_Ar$Ar
Ar<- vector(mode ="numeric", length = Nb_Pop)
for (i in 1:Nb_Pop){ 
  Ar[i] = mean(Ar_per_loc[,i])
}




# Formatting of the final table
Tab_Indices_per_pop <- rbind(N, GsurN, Ar, Ho, Hs, Fis, rbarD)
colnames(Tab_Indices_per_pop) <- Pop
Tab_Indices_per_pop <- t(Tab_Indices_per_pop)
kable(Tab_Indices_per_pop, digits = 3)









