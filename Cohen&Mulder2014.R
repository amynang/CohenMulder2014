library(rglobi)
library(tidyverse)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(NetIndices)
library(fluxweb)
library(soilfoodwebs)
library(vegan)

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
# taxa.info.2 = food %>% group_by(Genus.Morphon) %>% 
#   summarise(Feeding.Preference = first(Feeding.Preference),
#             AvgMass = mean(10^Log.averageMass.)) 


# this mostly results in V1=blabla-eating, V2=broad taxon
trythis = as.data.frame(str_split(taxa.info$Feeding.Preference, " ", n = 2, simplify = T))
#trythat = str_split(trythis[,2], " ", n = 2, simplify = T)

# a bit more juggling
taxainfo = cbind(taxa.info, trythis)
taxainfo$is = taxainfo$V2
taxainfo$is[which(str_detect(taxainfo$V2, "mite"))] = "mite"
taxainfo$is[which(str_detect(taxainfo$V2, "nematode ") | 
                    endsWith(taxainfo$V2,"nematode"))] = "nematode"

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
  
#  # metabolic rate (J/h) as a function of bodymass (mg), using the "linear model" from Ehnes 2011 
#  att[[i]]$ind.met.rate = exp(23.055335 + 0.695071*log(att[[i]]$freshMass.mg) - 0.68642*(1/(8.62*1e-5*(att[[i]]$Average.T+273.15))))
#  # population level metabolism (J/h)
#  att[[i]]$pop.met.rate = att[[i]]$ind.met.rate * att[[i]]$Abundance
  
  # population level metabolism (J/h)
  # the function will take as arguments the nth element of vectors Abundance, 
  # freshMass.mg, freshMass.mg and return a list of vectors
  # every vector contains sampled bodymasses for the nth taxon
  att[[i]]$pop.met.rate = pmap(list(ceiling(att[[i]]$Abundance), # how many draws
                                    att[[i]]$freshMass.mg,       # mean
                                    att[[i]]$freshMass.mg),      # sd
                               rlnormtrunc.intuitive) %>%  # then is passed on to another mapping
    map(.,  # the sum of metabolic energy loss of all individuals
            # metabolic rate (J/h) as a function of bodymass (mg), using the "linear model" from Ehnes 2011 
        ~ sum(exp(23.055335 + 0.695071*log(.) - 0.68642*(1/(8.62*1e-5*(att[[i]][1,]$Average.T+273.15)))))
    ) %>% do.call(rbind, .) %>% as.vector()  #...as vector
  
  
  # add missing basals 
  att[[i]] = att[[i]] %>% add_row(Genus.Morphon = c("bacteria","fungi","plants","detritus"), 
                                  Biomass.mg = sum(att[[i]]$Biomass.mg), # makes omnivores eat equally from everything
                                  pop.met.rate = 0,
                                  Average.T = att[[i]][1,"Average.T"])
  
  # resource based assimilation efficiency, from Lang et al. 2017
  temp.kT <- ((273.15+att[[i]]$Average.T)-293.15)/(8.62*1e-5*(273.15+att[[i]]$Average.T)*293.15) # explain 
  att[[i]]$efficiency = case_when(att[[i]]$Genus.Morphon == "detritus" ~ exp(-1.670)*exp(.164*temp.kT) / (1 + exp(-1.670)*exp(.164*temp.kT)),
                                  att[[i]]$Genus.Morphon ==   "plants" ~ exp(0.179) *exp(.164*temp.kT) / (1 + exp(0.179) *exp(.164*temp.kT)),
                                                                  # everything else (including fungi & bacteria???) gets predation efficiency
                                                                  TRUE ~ exp(2.266) *exp(.164*temp.kT) / (1 + exp(2.266) *exp(.164*temp.kT)))
  
  #att[[i]] = as.matrix(att[[i]]) 

}

# individual matrices for each site
web = vector(mode = "list", length=length(att)) # an empty list...
for (i in 1:length(att)) { 
  # ...filled with subsets of the meta-foodweb, adding the missing basal resources
  web[[i]] = mat[c(att[[i]]$Genus.Morphon),
                 c(att[[i]]$Genus.Morphon)]
}


library(cubature)

# If I have 1 bodymass distribution pre taxon across all locations I can do this once
# Repeating 135 times is... crazy


int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dlnormtrunc.intuitive(x, m=mu1, s=sd1, p=1)
  f2 <- dlnormtrunc.intuitive(x, m=mu2, s=sd2, p=1)
  pmin(f1, f2)
}

for (i in 1:length(att)) { 
  # number of taxa (removing basals)
  n = length(att[[i]]$freshMass.mg)-4
  # this is not wrong but probably not the simplest way?
  body.mat = replicate(n, att[[i]]$freshMass.mg)
  bodymat = body.mat[1:n,]
  bodymat[,] = 0
  
  # which cells in the interaction matrix are non zero
  ind = which(web[[i]][1:n,1:n] != 0, arr.ind = T)
  
  # vector of bodymasses
  bodymasses = att[[i]][ ,"freshMass.mg"]
  
  for (j in 1:nrow(ind)) { # for every predator-prey pair
    # we calculate prey suitability as the integral of the overlap of that prey's 
    # bodymass distribution and the optimal prey distribution for that predator
    # assuming OPPMR = 10^.6 (optimal prey ~4 times smaller than predator cf Brose 2006)
    overlap = cubintegrate(int_f, 0, Inf,
                           mu1=bodymasses[ind[j,][2]]/10^.6, # predator/10^.6
                           sd1=bodymasses[ind[j,][2]]/10^.6,
                           mu2=bodymasses[ind[j,][1]],       # prey
                           sd2=bodymasses[ind[j,][1]])$integral
    bodymat[ind[j,][1],ind[j,][2]] <- overlap
  }
  # this is so much slower!!!
  # for (j in 1:n) {
  #   for (k in 1:n) { 
  #     overlap <- cubintegrate(int_f, 0, Inf, 
  #                             mu1=bodymasses[j]/10^.6, 
  #                             sd1=bodymasses[j]/10^.6,
  #                             mu2=bodymasses[k],
  #                             sd2=bodymasses[k]
  #     )$integral
  #     bodymat[k,j] <- overlap
  #   }
  #   
  #   
  # }
  checkthat = web[[i]][1:n,1:n]*bodymat
  checkthat = vegan::decostand(checkthat,"total", 2)
  web[[i]][1:n,1:n] = checkthat
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(i, '/', length(att)))
  #Sys.sleep(.05)
  if (i == length(att)) cat('- Done!')
}

# https://rpsychologist.com/calculating-the-overlap-of-two-normal-distributions-using-monte-carlo-integration
int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dlnormtrunc.intuitive(x, m=mu1, s=sd1, p=1)
  f2 <- dlnormtrunc.intuitive(x, m=mu2, s=sd2, p=1)
  pmin(f1, f2)
}
integrate(int_f, -Inf, Inf, mu1=5.740768e-04/10^.6, mu2=2.877200e-04, 
                            sd1=5.740768e-04/10^.6, sd2=2.877200e-04)$value
#integrate(int_f, -Inf, Inf, mu1=55, mu2=65, sd1=5, sd2=6)

dlnormtrunc.intuitive = function(x, m, s, p=.9) {
  trnc <- EnvStats::dlnormTrunc(x, 
                                meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                sdlog = sqrt(log(1 + (s^2 / m^2))), 
                                min = qlnorm((1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))), 
                                max = qlnorm(1-(1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))))
  return(trnc)
}

rlnormtrunc.intuitive = function(n, m, s, p=.9) {
  trnc <- EnvStats::rlnormTrunc(n, 
                                meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                sdlog = sqrt(log(1 + (s^2 / m^2))), 
                                min = qlnorm((1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))), 
                                max = qlnorm(1-(1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))))
  return(trnc)
}

df <- data.frame(
  Data=factor(rep(c("D1", "D2"), each=2000)),
  weight=c(rlnormtrunc.intuitive(2000, m=5.740768e-04/10, s=5.740768e-04/10, p=1),
           rlnormtrunc.intuitive(2000, m=2.877200e-04,  s=2.877200e-04, p=1))
)
#df$weight = log10(df$weight)
d1dens <- with(df, density(weight[Data == "D1"], 
                           from = min(weight), 
                           to = max(weight)))
d2dens <- with(df, density(weight[Data == "D2"], 
                           from = min(weight),
                           to = max(weight)))
joint <- pmin(d1dens$y, d2dens$y)

df2 <- data.frame(x = rep(d1dens$x, 3), 
                  y = c(d1dens$y, d2dens$y, joint),
                  Data = rep(c("D1", "D2", "overlap"), each = length(d1dens$x)))

ggplot(df2, aes(x, y, fill = Data)) + 
  geom_area(position = position_identity(), color = "black") +
  scale_fill_brewer(palette = "Pastel2") +
  theme_bw()



################### Calculating energy fluxes using {fluxweb} ##################
# i.e. the method by Barnes et al (2018) ####

just1flux = fluxing(web[[8]], #This needs to be replaced by a preference matrix, controlling omnivores' diet
                    att[[8]]$Biomass.mg,
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




g = graph_from_adjacency_matrix(web[[1]], mode = "directed")
e <- get.edgelist(g)
df <- as.data.frame(cbind(e,E(g)$weight))
colnames(df) = c("Resource","Consumer")
df$type = ifelse(df$Resource %in% taxainfo$Genus.Morphon, "predation", 
                 ifelse(df$Resource == "plants", "herbivory", 
                        "detritivory"))

trlomn = TrophInd(web[[1]])

df$TL = trlomn$TL[match(df$Consumer, rownames(trlomn))]

hairball <- as_tbl_graph(df) %>% 
  activate(edges) %>% 
  mutate(type = as.character(type),
         TL = as.numeric(TL))

ggraph(hairball, layout = 'linear' #, sort.by = type
       , use.numeric = F) + 
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

g = graph_from_adjacency_matrix(web[[1]], mode = "directed")
tkplot(g)
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
