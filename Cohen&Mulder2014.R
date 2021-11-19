library(rglobi)
library(tidyverse)

food = read.table("H:/Literature/Cohen&Mulder2014_135FoodWebs.txt", header = T,
                  sep = "\t")
https://figshare.com/ndownloader/files/5629329

food = read.table("https://figshare.com/ndownloader/files/5629329", header = T,
                  sep = "\t")
doof = read.table("https://figshare.com/ndownloader/files/5629326", header = T,
                  sep = "\t")
unique(food$Feeding.Preference)
unique(food$Genus.Morphon)

taxa.info = food %>% group_by(Genus.Morphon) %>% 
  summarise(Feeding.Preference = first(Feeding.Preference)) 
trythis = str_split(taxa.info$Feeding.Preference, " ", n = 2, simplify = T)

trythat = str_split(trythis[,2], " ", n = 2, simplify = T)

mat = get_interaction_matrix(food[1:39,"Genus.Morphon"], 
                             food[1:39,"Genus.Morphon"],
                             "eats")
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = as.matrix(mat)