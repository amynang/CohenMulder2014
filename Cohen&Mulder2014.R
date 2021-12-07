library(rglobi)
library(tidyverse)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(NetIndices)

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
g = graph_from_adjacency_matrix(mat, mode = "directed")
e <- get.edgelist(g)
df <- as.data.frame(cbind(e,E(g)$weight))
colnames(df) = c("Resource","Consumer")
df$type = ifelse(df$Resource %in% taxainfo$Genus.Morphon, "predation", 
                 ifelse(df$Resource == "plants", "herbivory", 
                        "detritivory"))

trlomn = TrophInd(mat)

df$TL = trlomn$TL[match(df$Consumer, rownames(trlomn))]

hairball <- as_tbl_graph(df) %>% 
  activate(edges) %>% 
  mutate(type = as.character(type),
         TL = as.numeric(TL))

ggraph(hairball, layout = 'linear', sort.by = TL, use.numeric = T) + 
  geom_edge_arc(aes(colour = type), fold = T)


tkplot(g)

# individual attribute dataframes for each location
att = split(food, with(food, Web.ID)) 
# add basal resources
for (i in 1:length(att)) {
  att[[i]] = att[[i]] %>% add_row(Genus.Morphon = c("bacteria","fungi","plants","detritus"))
}

# individual matrices for each location
web = vector(mode = "list", length=length(att))
for (i in 1:length(att)) {
  web[[i]] = mat[att[[i]]$Genus.Morphon,
                 att[[i]]$Genus.Morphon]
}

# g = graph_from_adjacency_matrix(web[[1]], mode = "directed")
# tkplot(g)
# 
# 
# lay=matrix(nrow=dim(web[[1]])[1],ncol=2) # create a matrix with one column as runif, the other as trophic level
# lay[,1]=runif(nrow(web[[1]]), min = 1, max = 90)
# lay[,2]=TrophInd(web[[1]])$TL
# 
# for (i in 1:length(lay[,2])) { 
# lay[i,2] = lay[i,2] + rnorm(1, mean = 0,sd = .2)
# }
# plot(g, 
#      layout=lay,
#      vertex.size=8,
#      edge.arrow.size=.3,
#      edge.width=.5)



mat = get_interaction_matrix(food[1:39,"Genus.Morphon"], 
                             food[1:39,"Genus.Morphon"],
                             "eats")
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = as.matrix(mat)