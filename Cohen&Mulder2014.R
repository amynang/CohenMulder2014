library(rglobi)
library(tidyverse)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(NetIndices)
library(fluxweb)
library(soilfoodwebs)

# get the data
food = read.table("https://figshare.com/ndownloader/files/5629329", header = T,
                  sep = "\t")
# what's wrong here?
doof = read.table("https://figshare.com/ndownloader/files/5629326", header = T,
                  sep = "\t")

# who's there
unique(food$Genus.Morphon)
# are you ready to order
unique(food$Feeding.Preference)

# who's there and what do they like
taxa.info = food %>% group_by(Genus.Morphon) %>% 
  summarise(Feeding.Preference = first(Feeding.Preference)) 

# this mostly results in V1=blabla-eating, V2=broad taxon
trythis = as.data.frame(str_split(taxa.info$Feeding.Preference, " ", n = 2, simplify = T))
#trythat = str_split(trythis[,2], " ", n = 2, simplify = T)

# a bit more juggling
taxainfo = cbind(taxa.info, trythis)
taxainfo$is = taxainfo$V2
taxainfo$is[which(str_detect(taxainfo$V2, "mite"))] = "mite"
taxainfo$is[which(str_detect(taxainfo$V2, "nematode"))] = "nematode"

#taxainfo$eats = NA

############################ The Meta-foodweb ##################################

# assembling the meta-foodweb (all taxa across plots, adding missing basal resources)
mat = matrix(0, nrow = 262, ncol = 262,
             dimnames = list(c(taxainfo$Genus.Morphon,"bacteria","fungi","plants","detritus"),
                             c(taxainfo$Genus.Morphon,"bacteria","fungi","plants","detritus")))

# next, we use the info that was packed into "Feeding.Preference" to create the interaction matrix:

# fungi are consumed by...
mat["fungi", which(taxainfo$V1=="Fungivore" | 
                   taxainfo$V1=="Microphytophage" |
                   taxainfo$V1=="Omnivore")] = 1
# bacteria are consumed by...
mat["bacteria", which(taxainfo$V1=="Bacterivore" |
                      taxainfo$V1=="Incidental"|
                      taxainfo$Feeding.Preference == "Omnivore nematode")] = 1
# plants are consumed by...
mat["plants", which(taxainfo$V1=="Macrophytophage" |
                    taxainfo$V1=="Plant-feeding"|
                    taxainfo$V1=="Omnivore")] = 1
# detritus is consumed by...
mat["detritus", which(taxainfo$V1=="Substrate-ingesting"|
                      taxainfo$V1=="Substrate-inhabiting"|
                      taxainfo$V1=="Omnivore")] = 1
# nematodes are consumed by...
mat[which(taxainfo$is=="nematode"), which(taxainfo$Feeding.Preference == "Predatory mite (attacking nematodes)" |
                                          taxainfo$Feeding.Preference == "Predating nematode (consuming nematodes)"|
                                          taxainfo$Feeding.Preference == "Generalist mite"|
                                          taxainfo$Feeding.Preference == "Omnivore nematode"|
                                          taxainfo$Feeding.Preference == "Omnivore mite" |
                                          taxainfo$Feeding.Preference == "Parasitizing mite (hosts are mites or nematodes)")] = 1
# mites are consumed by...
mat[which(taxainfo$is=="mite"), which(taxainfo$Feeding.Preference == "Predatory mite (attacking arthropods)"|
                                      taxainfo$Feeding.Preference == "Generalist mite" |
                                      taxainfo$Feeding.Preference == "Omnivore mite" |
                                      taxainfo$Feeding.Preference == "Parasitizing mite (hosts are mites or nematodes)"|
                                      taxainfo$Feeding.Preference == "Omnivore insect")] = 1
# insects, pauropoda, collembola, symphyla are consumed by...
mat[which(taxainfo$Feeding.Preference == "Fungivore insects and pauropods"|
          taxainfo$Feeding.Preference == "Plant-feeding collembolans and symphylids"),
    which(taxainfo$Feeding.Preference == "Predatory mite (attacking arthropods)"|
          taxainfo$Feeding.Preference == "Generalist mite" |
          taxainfo$Feeding.Preference == "Omnivore mite" |
          taxainfo$Feeding.Preference == "Omnivore insect")] = 1

# I think that covers everything?

########################## Site specific foodwebs ##############################

# individual attribute dataframes for each site
att = split(food[,c(1:9,29)], with(food, Web.ID))
for (i in 1:length(att)) { 
  #att[[i]]$Abundance = exp(att[[i]]$Log.Abundance.)
  att[[i]]$Abundance = 10^att[[i]]$Log.Abundance.
  att[[i]]$averageMass.mg = (10^att[[i]]$Log.averageMass.)*1e-3 # from log(micrograms/m^2) to milligrams/m^2)
  att[[i]]$Biomass.mg = (10^att[[i]]$Log.Biomass.)*1e-3         # from log(micrograms/m^2) to milligrams/m^2) 
  
  # broad-group specific conversion of dry mass (micrograms) to fresh mass (mg),
  # based on Mercer (2001) for annelids and insects, 
  # van den Hoogen (2019) for nematodes (20% of fresh weight)
  # and Yaninek & Gnanvossou (1993) for mites (30% of fresh weight)
  # Mercer's equations work with base 10 logarithms and grams!
  # so dry mass is multiplied by 1e-3 to go from milligrams to grams
  # then fresh mass (grams) is multiplied by 1e3 to go down to milligrams
  # this is unnecessary for mites and nematodes since we are working with percentages
  # the "taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)]" bit 
  # is an unintuitive way to figure out which group each taxon belongs to
  att[[i]]$freshMass.mg = case_when(taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)] %in% 
                                      c("enchytraeids",
                                        "lumbricids") ~ (10^(0.9282 +1.0899*log10(att[[i]]$averageMass*1e-3)))*1e3,
                                    taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)] %in% 
                                      c("insects and pauropods",
                                        "insect",
                                        "collembolans and symphylids") ~ (10^(0.6111 +1.0213*log10(att[[i]]$averageMass*1e-3)))*1e3,
                                    taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)] == "mite" ~ (att[[i]]$averageMass)/.3,
                                    taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)] == "nematode" ~ (att[[i]]$averageMass)/.2)
  
  # metabolic rate (J/h) as a function of bodymass (mg), using the "linear model" from Ehnes 2011 
  att[[i]]$ind.met.rate = exp(23.055335 + 0.695071*log(att[[i]]$freshMass.mg) - 0.68642*(1/(8.62*1e-5*(att[[i]]$Average.T+273.15))))
  # population level metabolism (J/h)
  att[[i]]$pop.met.rate = att[[i]]$ind.met.rate * att[[i]]$Abundance
  
  # add missing basals 
  att[[i]] = att[[i]] %>% add_row(Genus.Morphon = c("bacteria","fungi","plants","detritus"), 
                                  Biomass = 1, # dummy !!!will mess up preferences if that is left to be done by fluxing()!!!
                                  pop.met.rate = 0,
                                  Average.T = att[[i]][1,"Average.T"])
  
  # resource based assimilation efficiency, from Lang et al. 2017
  temp.kT <- ((273.15+att[[i]]$Average.T)-293.15)/(8.62*1e-5*(273.15+att[[i]]$Average.T)*293.15) # explain 
  att[[i]]$efficiency = case_when(att[[i]]$Genus.Morphon == "detritus" ~ exp(-1.670)*exp(.164*temp.kT) / (1 + exp(-1.670)*exp(.164*temp.kT)),
                                  att[[i]]$Genus.Morphon ==   "plants" ~ exp(0.179) *exp(.164*temp.kT) / (1 + exp(0.179) *exp(.164*temp.kT)),
                                                                  # everything else (including fungi & bacteria???) gets predation efficiency
                                                                  TRUE ~ exp(2.266) *exp(.164*temp.kT) / (1 + exp(2.266) *exp(.164*temp.kT)))
  
  #att[[i]] = att[[i]] %>% relocate(Sepal.Width:Petal.Length, .before = Sepal.Length)

}

# individual matrices for each site
web = vector(mode = "list", length=length(att)) # an empty list...
for (i in 1:length(att)) { 
  # ...filled with subsets of the meta-foodweb, adding the missing basal resources
  web[[i]] = mat[c(att[[i]]$Genus.Morphon,"bacteria","fungi","plants","detritus"),
                 c(att[[i]]$Genus.Morphon,"bacteria","fungi","plants","detritus")]
}


################### Calculating energy fluxes using {fluxweb} ##################
# i.e. the method by Barnes et al (2018) ####

just1flux = fluxing(web[[8]], #This needs to be replaced by a preference matrix, controlling omnivores' diet
                    att[[8]]$Biomass,
                    att[[8]]$pop.met.rate,
                    att[[8]]$efficiency,
                    bioms.prefs = TRUE, # This needs to change
                    bioms.losses = FALSE, 
                    ef.level = "prey")

#### Try {soilfoodwebs} ####



#### graph nonsense ####
taxainfo$is[match(att[[i]]$Genus.Morphon, taxainfo$Genus.Morphon)]

exp(28.019958 + .558099*log(metabolic.rates$mg.fresh.mass) - .803007*(1/(8.62*1e-5*(Temp+273.15))))

View(att[[8]])

# individual matrices for each plot
web = vector(mode = "list", length=length(att)) # an empty list...
for (i in 1:length(att)) { 
  # ...filled with subsets of the meta-foodweb, adding the missing basal resources
  web[[i]] = mat[c(att[[i]]$Genus.Morphon,"bacteria","fungi","plants","detritus"),
                 c(att[[i]]$Genus.Morphon,"bacteria","fungi","plants","detritus")]
}




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