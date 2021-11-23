library(rglobi)
library(tidyverse)
library(stringr)

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
trythis = as.data.frame(str_split(taxa.info$Feeding.Preference, " ", n = 2, simplify = T))


trythat = str_split(trythis[,2], " ", n = 2, simplify = T)

taxainfo = cbind(taxa.info, trythis)
taxainfo$is = taxainfo$V2
taxainfo$is[which(str_detect(taxainfo$V2, "mite"))] = "mite"
taxainfo$is[which(str_detect(taxainfo$V2, "nematode"))] = "nematode"

taxainfo$eats = NA

mat = matrix(0, nrow = 262, ncol = 262,
             dimnames = list(c(taxainfo$Genus.Morphon,"bacteria","fungi","plants","detritus"),
                             c(taxainfo$Genus.Morphon,"bacteria","fungi","plants","detritus")))

mat["fungi", which(taxainfo$V1=="Fungivore" | 
                   taxainfo$V1=="Microphytophage" |
                   taxainfo$V1=="Omnivore")] = 1

mat["bacteria", which(taxainfo$V1=="Bacterivore" |
                      taxainfo$V1=="Incidental"|
                      taxainfo$Feeding.Preference == "Omnivore nematode")] = 1

mat["plants", which(taxainfo$V1=="Macrophytophage" |
                    taxainfo$V1=="Plant-feeding"|
                    taxainfo$V1=="Omnivore")] = 1

mat["detritus", which(taxainfo$V1=="Substrate-ingesting"|
                      taxainfo$V1=="Substrate-inhabiting"|
                      taxainfo$V1=="Omnivore")] = 1

mat[which(taxainfo$is=="nematode"), which(taxainfo$Feeding.Preference == "Predatory mite (attacking nematodes)" |
                                          taxainfo$Feeding.Preference == "Predating nematode (consuming nematodes)"|
                                          taxainfo$Feeding.Preference == "Generalist mite"|
                                          taxainfo$Feeding.Preference == "Omnivore nematode"|
                                          taxainfo$Feeding.Preference == "Omnivore mite" |
                                          taxainfo$Feeding.Preference == "Parasitizing mite (hosts are mites or nematodes)")] = 1

mat[which(taxainfo$is=="mite"), which(taxainfo$Feeding.Preference == "Predatory mite (attacking arthropods)"|
                                      taxainfo$Feeding.Preference == "Generalist mite" |
                                      taxainfo$Feeding.Preference == "Omnivore mite" |
                                      taxainfo$Feeding.Preference == "Parasitizing mite (hosts are mites or nematodes)"|
                                      taxainfo$Feeding.Preference == "Omnivore insect")] = 1

mat[which(taxainfo$Feeding.Preference == "Fungivore insects and pauropods"|
          taxainfo$Feeding.Preference == "Plant-feeding collembolans and symphylids"),
    which(taxainfo$Feeding.Preference == "Predatory mite (attacking arthropods)"|
          taxainfo$Feeding.Preference == "Generalist mite" |
          taxainfo$Feeding.Preference == "Omnivore mite" |
          taxainfo$Feeding.Preference == "Omnivore insect")] = 1

#I think that covers everything?


which(str_detect(taxainfo$V2, "mite"))
View(taxainfo[taxainfo$Feeding.Preference == "Plant-feeding collembolans and symphylids",])

mat = get_interaction_matrix(food[1:39,"Genus.Morphon"], 
                             food[1:39,"Genus.Morphon"],
                             "eats")
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = as.matrix(mat)