---
title: "population genetics analysis Mlp"
output: html_notebook
---

## Author: Ammar Abdalrahem

------------------------------------------------------------------------

# 1. Dependencies

Remember to re-run this code every time you re-open this R Notebook.

```{r}
#Code to install packages if necessary, and read them with library function

required_packages <- c("ggplot2","readxl","tidyverse","genepop","hierfstat","here","mapdata",
                       "mapplots","data.table","grDevices","colorspace","adegenet","poppr","pegas","ape","ade4","remotes","ggtree")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

```

------------------------------------------------------------------------

# 2. Obtain data


```{r}
# get the path of the current R script
path <- dirname(rstudioapi::getSourceEditorContext()$path)

# set the working directory to the path of the current R script
setwd(path)

# check the current working directory
getwd()


# import the data within specified sheet
data <- read_excel("SexAsex4Ammar.xlsx", 
                  sheet = "Ammar")

#covert data to data frame
data <- as.data.frame(data)


# Remove rows containing "sex" in "Profil" except for those sampled at "Prelles"
data <- data[!(data$Profil == "sex" & !grepl("Prelles", data$Site, ignore.case = TRUE)), ]

# remove unknown Profil & population
data <- data[!(data$Profil == "NA" | data$Pop == "NA"), ]

#take a look to data

head(data)

```
------------------------------------------------------------------------


# 3. Data map visualization

just to look at where each population is located 

```{r}
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
xyz <- make.xyz(pop_data$Long, pop_data$Lat, pop_data$percent, pop_data$Profil)

# Colors used
my_colors <- c("#EDF060", "#F34213")
transparent_colors <- adjustcolor(my_colors, alpha.f = 0.5)

# The plot of the pie chart above the map
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),
    mgp = c(2.5, 0.5, 0), family = "Arial")
basemap(xlim =  c(-4.3,9.1), ylim =  c(41, 51), bg = "white",
        main = "Distribution of population of poplar rust in France")
map("france", fill=FALSE, col="light blue", xlim = xlim, ylim = ylim, add = TRUE)
draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.3, col = transparent_colors)
legend.pie(-4.5, 41.4, labels = c("Asex","Sex"), 
           radius = 0.3, bty = "n", col = my_colors, cex = 1, label.dist = 1.5)
```
------------------------------------------------------------------------

# 4. Microsatellite marker data


```{r}
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

Nb_Pop = length(levels(data_GenInd@pop)) # number of population
is.numeric(Nb_Pop)

#convert Genind to Genepop format
data_Genpop <- genind2genpop(data_GenInd, process.other=TRUE)
data_Genpop # See result
```

------------------------------------------------------------------------

# 5. Calculate the Genetic distance as Euclidean for NJ tree
 
```{r}
#Edwards' distance (Euclidean)
data.dist.edwards <- dist.genpop(x=data_Genpop, method=2, diag=T, upper=T)
is.euclid(data.dist.edwards, plot=TRUE, print=TRUE, tol=1e-10) # FALSE = Yes

#Calculate and test NJ tree for Euclidean 
data.nj <- nj(data.dist.edwards) # Calculates the tree 
# Test tree quality - plot original vs. reconstructed distance
plot(as.vector(data.dist.edwards), as.vector(as.dist(cophenetic(data.nj))),
     xlab="Original distance", ylab="Reconstructed distance")
abline(lm(as.vector(data.dist.edwards) ~ as.vector(as.dist(cophenetic(data.nj)))), col="red")
cor.test(x=as.vector(data.dist.edwards), y=as.vector(as.dist(cophenetic (data.nj))), alternative="two.sided") # Testing the correlation
# Linear model for above graph
summary(lm(as.vector(data.dist.edwards) ~ as.vector(as.dist(cophenetic(data.nj))))) # Prints summary text
# Plot a basic tree - see ?plot.phylo for details
plot.phylo(x=data.nj, type="phylogram")



```

------------------------------------------------------------------------

# 6. Calculate the Genetic distance as non-Euclidean for NJ tree
the final tree for me don't make any sense 
 
```{r}
# Standard Nei's genetic distance (D) 1972 (not Euclidean)
data.dist.D <- dist.genpop(x=data_Genpop, method=1, diag=T, upper=T)
is.euclid(data.dist.D, plot=TRUE, print=TRUE, tol=1e-10) # FALSE 


data.dist.D.m <- cailliez(distmat=data.dist.D, print=FALSE, tol=1e-07, cor.zero=TRUE)
is.euclid(data.dist.D.m, plot=TRUE, print=TRUE, tol=1e-10) # TRUE = OK

#Calculate and test NJ tree for Euclidean 
data.nj <- nj(data.dist.D.m) # Calculates the tree 
# Test tree quality - plot original vs. reconstructed distance
plot(as.vector(data.dist.D.m), as.vector(as.dist(cophenetic(data.nj))),
     xlab="Original distance", ylab="Reconstructed distance")
abline(lm(as.vector(data.dist.D.m) ~ as.vector(as.dist(cophenetic(data.nj)))), col="red")
cor.test(x=as.vector(data.dist.D.m), y=as.vector(as.dist(cophenetic (data.nj))), alternative="two.sided") # Testing the correlation
# Linear model for above graph
summary(lm(as.vector(data.dist.D.m) ~ as.vector(as.dist(cophenetic(data.nj))))) # Prints summary text
# Plot a basic tree - see ?plot.phylo for details
plot.phylo(x=data.nj, type="phylogram",edge.width=1, main="NJ")

```
------------------------------------------------------------------------

# 7. Microsatellite marker data as Mll

```{r}
# Convert to genind object
data_GenInd_mll <- df2genind(
  X = genotype_data,   #data.frame containing allele data only 
  sep = NULL,
  ncode = 3,           #an optional integer giving the number of characters used for coding one genotype at one locus.
  ind.names = rownames(genotype_data),  # individuals names
  loc.names = colnames(genotype_data),  # markers names
  pop = data$Mll,                       # giving the population of each individual
  NA.char = "999",                      # string corresponding to missing allele 999 or 999999
  ploidy = 2,
  type = "codom",                       #codom' stands for 'codominant' (e.g. microstallites, allozymes)
  strata = NULL,
  hierarchy = NULL,
  check.ploidy = getOption("adegenet.check.ploidy")
)


# ca marche !!
summary(data_GenInd_mll)

Nb_Pop_mll = length(levels(data_GenInd_mll@pop)) # number of mll
is.numeric(Nb_Pop_mll)

#Do you think it is important to dd coordinates to genotype data?

#convert Genind to Genepop format
data_Genpop_mll <- genind2genpop(data_GenInd_mll, process.other=TRUE)
data_Genpop_mll # See result
```

------------------------------------------------------------------------
# 8. Calculate the Genetic distance as non-Euclidean for NJ tree for Mll
the final tree for me don't make any sense 
 
```{r}
# Standard Nei's genetic distance (D) 1972 (not Euclidean)
data.dist.D_mll <- dist.genpop(x=data_Genpop_mll, method=1, diag=T, upper=T)
is.euclid(data.dist.D_mll, plot=TRUE, print=TRUE, tol=1e-10) # FALSE 


data.dist.D_mll.m <- cailliez(distmat=data.dist.D_mll, print=FALSE, tol=1e-07, cor.zero=TRUE)
is.euclid(data.dist.D.m, plot=TRUE, print=TRUE, tol=1e-10) # TRUE = OK

#Calculate and test NJ tree for Euclidean 
data.nj.mll <- nj(data.dist.D_mll.m) # Calculates the tree 
# Test tree quality - plot original vs. reconstructed distance
plot(as.vector(data.dist.D_mll.m), as.vector(as.dist(cophenetic(data.nj.mll))),
     xlab="Original distance", ylab="Reconstructed distance")
abline(lm(as.vector(data.dist.D_mll.m) ~ as.vector(as.dist(cophenetic(data.nj.mll)))), col="red")
cor.test(x=as.vector(data.dist.D_mll.m), y=as.vector(as.dist(cophenetic (data.nj.mll))), alternative="two.sided") # Testing the correlation
# Linear model for above graph
summary(lm(as.vector(data.dist.D_mll.m) ~ as.vector(as.dist(cophenetic(data.nj.mll))))) # Prints summary text
# Plot a basic tree - see ?plot.phylo for details
plot.phylo(x=data.nj.mll, type="phylogram",edge.width=1, main="NJ")
```

------------------------------------------------------------------------

# 9. try manual analysis

```{r}
# Assuming your genotype data is stored in a data frame called genotype_data

# Create an empty data frame to store the separated values
separated_data <- list()

# Iterate over each column in the genotype data
for (col in colnames(genotype_data)) {
  # Split each three-digit value into two columns
  separated_values <- cbind(substr(genotype_data[, col], 1, 3), substr(genotype_data[, col], 4, 6))
  
  # Add the separated values to the new data frame
  separated_data <- cbind(separated_data, separated_values)
  

}

# Convert the values in the matrix to numeric
separated_data <- apply(separated_data, 2, as.numeric)

# Assign column names to the separated data
colnames(separated_data) <- rep(colnames(genotype_data), each = 2)

# Assign row names to the separated data
rownames(separated_data) <- rownames(genotype_data)

# check the separated data
head (separated_data)
class(separated_data[1,1])

dist_v <- dist(separated_data, method = "euclidean", diag = TRUE, upper = TRUE)

#Neighbor-Joining Tree Estimation
#This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
tree <- nj(dist_v)

# Create group information based on 'Profil' column
groupInfo <- split(data$Isolate, data$Profil)

# Group the tree labels based on the group information
tree <- groupOTU(tree, groupInfo)

# Plot the circular tree with grouped labels based on 'Profil' information
ggtree(tree, aes(color = group)) +
    geom_tiplab(size = 1)

# Generate bootstrap support values for the tree
boot_tree <- boot.phylo(tree, dist_v, FUN = nj, B = 100)

```

------------------------------------------------------------------------
