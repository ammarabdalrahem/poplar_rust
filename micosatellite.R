#!/usr/bin/env Rscript

#Code to install packages if necessary, and read them with library function
required_packages <- c("ggplot2","readxl","tidyverse","poppr","adegenet","genepop","hierfstat","here","ggforce","ggspatial")
for (package in required_packages) {
  if (package %in% row.names(installed.packages())) {
    library(package, character.only = TRUE)
  } else {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}


# --- Obtain data --- #

#set your data directory
#anyway to make the directory what i use automatic?
#set your in current directory 
here::set_here()
getwd()



# import the data within specified sheet
data <- read_excel("SexAsex4Ammar.xlsx", 
                  sheet = "Ammar")

data <- as.data.frame(data)

# Calculate number of non-missing values for each individual
#non_missing <- rowSums(!is.na(data[,19:42]))

# Remove individuals with no scored loci
#data_clean <- data[non_missing > 0,]


# --- work on data --- #
# summarize data by population and location
pop_data <- data %>%
  group_by(Pop, Long, Lat,Profil) %>%
  summarize(n = n()) %>%
  mutate(percent = n / sum(n)) %>%
  select(Pop, Long, Lat,Profil, n, percent)



# example data

# The area of the France Region;
xlim <- c(-4.3,9.1)
ylim <- c(41, 51)

# Creates an xyz object for use with the function draw.pie
xyz <- make.xyz(pop_data$Long, pop_data$Lat, pop_data$percent, pop_data$Profil)

# Colors used
col <- c("#003366", "#CCCC66", "#CC3366")

# The plot of the pie chart above the map
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),
    mgp = c(2.5, 0.5, 0), family = "Arial")
basemap(xlim =  c(-4.3,9.1), ylim =  c(41, 51), bg = "white",
        main = "Distribution of population of popular rust in France")
map('world2Hires', xlim = xlim, ylim = ylim, add = TRUE)
draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.3, col = col)
legend.pie(-4.3, 41, labels = c("sex", "asex", "NA"), 
           radius = 0.2, bty = "n", col = col, cex = 0.8, label.dist = 1.3)



#text(121.5, 12.1, "Tuna Species:", cex = 0.8, font = 2)
#text(123.4, 14.3, "Camarines\nNorte", cex = 0.8, font = 2)
#text(124.8, 13.75, "Catanduanes", cex = 0.8, font = 2)
#text(122.9, 13.5, "Camarines\nSur", cex = 0.8, font = 2)
#text(124.1, 13.25, "Albay", cex = 0.8, font = 2)
#text(124.3, 12.8, "Sorsogon", cex = 0.8, font = 2)
#text(123.5, 12, "Masbate", cex = 0.8, font = 2)
mtext("Data Source: FAb Statistics Authority", side = 1, outer = TRUE,
      adj = 1, cex = 0.8, font = 3)








# Create a base map
world <- map_data("world")

# Plot the base map and add points
ggplot(data_clean, aes(x = Long, y = Lat)) +
  geom_circle(aes(r = sqrt(Value/pi)), alpha = 0.5, fill = "blue", color = NA) +
  geom_text(aes(label = ID), size = 2, alpha = 0.5) +
  theme_void() +
  coord_sf(xlim = c(-130, -65), ylim = c(24, 50), expand = FALSE)

# Table of genotypes
tail(data_clean[,19:42 ]) 
names(data_clean[,19:42 ]) 


# Convert to genind 

data_GenInd <- df2genind(
  X = data_clean[,19:42 ],
  sep = NULL,
  ncode = 3,
  ind.names = data_clean[,1]  ,
  pop = data_clean[,3],
  NA.char = "999",
  ploidy = 2,
  type = "codom", 
  strata = NULL,
  hierarchy = NULL,
  check.ploidy = getOption("adegenet.check.ploidy")
)


# ca marche !!
summary(data_GenInd)



## Identification des individus appartenant aux lignées asexuées
# Set the initial value of the flag to FALSE
flag <- FALSE

# Repeat until the flag is TRUE
while(!flag) {
  
  # Find clusters with k=2 on the data
  grp2 <- find.clusters(data_GenInd, max.n=50, n.pca=60, n.clust=2) 
  plot (grp2)
  # Check if the second cluster contains the expected number of individuals (206)
  if (summary(grp2$grp)[2] == 206) { 
    flag <- TRUE  # Set the flag to TRUE to exit the loop
  }
  
}



# visualisation des résultats
sex_asex = grp2$grp
summary(grp2$grp)
table.value(table(pop(data_GenInd), sex_asex), col.lab=paste("  cluster", 1:2), row.lab=c(levels(data_GenInd@pop)))









#marker data
marker_data <- data %>% select(contains("Mlp"))



#1.	Multidimensional scaling (MDS) plots 
#unction to compute pairwise distances between all 
#rows and then convert the resulting object to a square matrix using the as.matrix() function.

d <- dist(marker_data, method = "euclidean")
square_matrix <- as.matrix(d)

# Calculate MDS coordinates
mds_coordinates <- cmdscale(d, k = 2)

# Plot the MDS coordinates
plot(mds_coordinates, type = "n")
text(mds_coordinates, labels = row.names(marker_data), cex = 0.7, pos = 3)


# Perform PCA on the data


# create a leaflet map centered on the France
map <- leaflet() %>%
  addTiles() %>%
  setView(lng = 46.2276, lat = 2.2137, zoom = 4)

# add markers to the map for each location
map_markers <- addMarkers(map, 
                          lng = data$Long, 
                          lat = data$Lat, 
                          popup = data$Isolate)



















# summarize data by population
pop_data <- data %>%
  group_by(Pop) %>%
  summarize(n = n()) %>%
  mutate(percent = n / sum(n))

# create markers with pie charts for each population
markers <- pop_data %>%
  # use pie charts to visualize population proportions
  with(pie(percent, labels = Pop)) %>%
  # add circle markers
  addCircleMarkers(
    lng = data$Long, 
    lat = data$Lat, 
    radius = 10, # adjust the radius size
    fillColor = "#FFFFFF", # fill color
    color = "#000000", # border color
    weight = 1, # border weight
    opacity = 1 # opacity
  )

# add the markers to the map
map %>% 
  addLegend(
    position = "topright",
    pal = colorNumeric(palette = "Set1", domain = pop_data$Pop),
    values = pop_data$Pop,
    title = "Population"
  ) %>% 
  addLayersControl(
    overlayGroups = "Markers",
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  addCircleMarkers(data = markers, group = "Markers")


