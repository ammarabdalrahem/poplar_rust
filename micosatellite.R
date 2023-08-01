---
title: "population genetics analysis Mlp"
output: html_notebook
---

## Author: Ammar Abdalrahem

------------------------------------------------------------------------

## 1. Dependencies

Remember to re-run this code every time you re-open this R Notebook.

```{r, eval=TRUE}
#Code to install packages if necessary, and read them with library function

required_packages <- c("ggplot2","readxl","tidyverse","genepop","hierfstat","here","mapdata",
                       "mapplots","data.table","grDevices","colorspace","adegenet","poppr","pegas","ape","ade4","remotes","ggtree","ggtreeExtra")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

```



## 2. Obtain data


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
#data <- data[!(data$Profil == "sex" & !grepl("Prelles", data$Site, ignore.case = TRUE)), ]

# remove unknown Profil & population
data <- data[!(data$Profil == "NA" | data$Pop == "NA"), ]

#take a look to data

head(data)

```


## 3. Define MLG and MLL

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

data_Genclone <- as.genclone(data_GenInd)

#Define MLG according threshold = 0
mlg_assignments <- mlg.filter(data_Genclone, threshold = 0, distance = "diss.dist", threads = 1L,missing = "asis") #he "asis" option is used to indicate that missing data should be treated "as is," meaning that missing values will be retained as NA in the distance matrix.

#add MLG result to the table
genotype_data$MLG <- mlg_assignments

#Define the MLL 
#Choosing a threshold
#After you have chosen a genetic distance and a filtering algorithm
#choose threshold to represent the minimum genetic distance at which two individuals would be considered from different clonal lineages.

data_filtered <-filter_stats(data_Genclone, distance = diss.dist, plot = TRUE, missing = "asis")
# “farthest neighbor” algorithm.
#Arnaud-Haond et al. 2007, @bailleul2016rclone
print(farthest_thresh <- cutoff_predictor(data_filtered$farthest$THRESHOLDS))

# “UPGMA ” algorithm.
print(average_thresh  <- cutoff_predictor(data_filtered$average$THRESHOLDS))

# “nearest neighbor” algorithm.
print(nearest_thresh  <- cutoff_predictor(data_filtered$nearest$THRESHOLDS))


#Define the MLL threshold 0.5, algorithm	Farthest neighbor

mll_assignments<- mlg.filter(data_Genclone, threshold = farthest_thresh, distance = "diss.dist", threads = 1L, missing = "asis") 


genotype_data$MLL <- mll_assignments

#write.csv2(genotype_data, "data_new.csv")

```


## Identification the individuals profile(sex/asex)

```{r}
# First identification by cluster
grp <- find.clusters(data_GenInd, method = "kmeans", stat = "BIC", n.pca= 90 , n.clust= 2, n.iter=100000, n.start=100)
cluster_assignments <- grp$grp
genotype_data$cluster <- cluster_assignments

dapc1 <- dapc(data_GenInd, grp$grp, n.pca= 90 , n.clust=2, n.da = 100)


# results visualization  
scatter(dapc1)
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:2), lab="", ncol=1, xlab="individuals")

#found special individual 
# Convert the posterior probabilities to a data frame
posterior_data <- as.data.frame(dapc1$posterior)

# Add individual names as a column in the data frame
posterior_data$Individual <- rownames(posterior_data)

# Melt the data frame for visualization
melted_data <- reshape2::melt(posterior_data, id.vars = "Individual", variable.name = "Cluster", value.name = "Probability")

# Create a scatter plot with individual names as labels
ggplot(melted_data, aes(x = Cluster, y = Probability, color = Cluster, label = Individual)) +
  geom_point() +
  xlab("Cluster") +
  ylab("Probability") +
  labs(color = "Cluster") +
  theme_minimal() +
  geom_text(nudge_y = 0.02)  # Add labels slightly above the data points


#reomve uncertin cluster
# Assuming you have stored the posterior probabilities in a variable named 'posterior'
# Subset the main table to include only individuals with probability 1.00 in either cluster
new_data <- data[posterior_data[, 1] == 1 | posterior_data[, 2] == 1, ]

# View the filtered data
print(filtered_data)

write.csv2(filtered_data, "data_new.csv")
#how to check now ?
```




## 4. Data map visualization

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




----------------------------------------------------------------

#  try NJ it with popper
 
```{r}
# Calculate distance matrix and build NJ tree
dist <- diss.dist(data_GenInd)
tree <- nj(dist)

# Create group information based on 'Profil' column
groupInfo <- split(data$Isolate, data$Profil)

# Group the tree labels based on the group information
tree <- groupOTU(tree, groupInfo)

# Create a dataframe for annotation
dat1 <- data.frame(
  ID = data$Isolate,
  lat = data$Lat,
  long = data$Long,
  Location = data$Site,
  Group = data$Profil,
  Year = data$Year
)

# Create the ggtree plot with circular layout
options(ignore.negative.edge=TRUE)
p <- ggtree(tree, aes(color = group))
p

# Use %<+% of ggtree to add annotation dataset to the tree
p1 <- p %<+% dat1

# Initialize fill scale using new_scale_fill() from ggnewscale package
p2 <- p +
  geom_fruit(
    data = dat1,
    geom = geom_col,
    mapping = aes(y = ID, x = Year, fill = Location),  # Map 'Location' to fill
    pwidth = 0.4,
    offset = 0.01,
    axis.params = list(
      axis = "x",  # Add x-axis text
      text.angle = -45,  # Adjust text angle
      hjust = 0,  # Adjust horizontal position
      text.size = 2.5,
      line.size = 0.4,
      line.color = "black"
    ),
    grid.params = list(color = "black", linetype = 5, size = 0.4, alpha = 0.8)  # Add grid lines
  ) +
  scale_shape_manual(
    values = 1:length(unique(dat1$Location))  # Set shape values
  ) +
  theme(
    #legend.position = c(1.15, 0.5),  # Adjust legend position
    legend.background = element_rect(fill = NA),  # Set legend background
    legend.title = element_text(size = 9),  # Adjust legend title size
    legend.text = element_text(size = 7),  # Adjust legend text size
    legend.spacing.y = unit(0.3, "cm")  # Adjust legend spacing (y orientation)
  ) 
   #+ geom_treescale(fontsize=2, linesize=0.3, x=-10, y=1000)

p2
# Modify legend titles
p2 <- p2 +
  labs(fill = "Location", color = "Profile") 
  

p2
#ggsave("tree_plot2.png", p2, width = 17, height = 8, dpi = 600)  # Save the plot with desired dimensions

```
----------------------------------------------------------------






