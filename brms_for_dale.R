
library(brms)

##########################################################################
#PREPARE DATA
##########################################################################

# Calculate mean richness for each species in the ATTRIBUTE_SPECIES column
mean_richness_by_species <- aggregate(richness ~ ATTRIBUTE_SPECIES, data = meta_adult, FUN = mean)

# Rename the 'richness' column from the aggregated result to 'mean_richness'
colnames(mean_richness_by_species)[colnames(mean_richness_by_species) == "richness"] <- "mean_richness"

# Merge the mean_richness back into the original dataset
meta_adult <- merge(meta_adult, mean_richness_by_species, by = "ATTRIBUTE_SPECIES", all.x = TRUE)

# Select specific columns without repeating column names
meta_adult = data.frame(meta_adult)
data_repeat <- meta_adult[, c( "richness", "ATTRIBUTE_SPECIES","ATTRIBUTE_SPECIES", "mean_richness.x","ATTRIBUTE_PC1score","ATTRIBUTE_PC2score")]
colnames(data_repeat) = c("phen","species","phylo","spec_mean_cf", "spec_mean_pc1", "spec_mean_pc2")
data_repeat = data_repeat[data_repeat[,"species"]!= "BLANK",]

for(i in 1:nrow(data_repeat)) {
  sp = data_repeat[i,"species"]
 data_repeat[i,"species"]  = meta_seedling[meta_seedling[,"ATTRIBUTE_SPECIES"] == sp,"genu_species"][1]
  data_repeat[i,"phylo"]  = meta_seedling[meta_seedling[,"ATTRIBUTE_SPECIES"] == sp,"genu_species"][1]
  data_repeat[i,"spec_mean_pc1"]  = meta_seedling[meta_seedling[,"ATTRIBUTE_SPECIES"] == sp,"ATTRIBUTE_PC1score"][1]
  data_repeat[i,"spec_mean_pc2"]  = meta_seedling[meta_seedling[,"ATTRIBUTE_SPECIES"] == sp,"ATTRIBUTE_PC2score"][1]
}

data_repeat = data_repeat[data_repeat[,"species"] %in% mytree22$tip.label,]

data_repeat$spec_mean_pc1 = as.numeric(data_repeat$spec_mean_pc1)
data_repeat$spec_mean_pc2 = as.numeric(data_repeat$spec_mean_pc2)
is.numeric(data_repeat$spec_mean_pc1)

#write.csv(data_repeat, "data_repeat_for_dale")
##########################################################################
#SET UP THE PHYLOGENY
##########################################################################

#phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
#A <- ape::vcv.phylo(phylo)

Tree<-ape::read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")
selected_labels = unique(meta_seedling$genu_species[meta_seedling$genu_species%in%Tree$tip.label])

sub_tree <- ape::drop.tip(Tree, tip = Tree$tip.label[!(Tree$tip.label %in% selected_labels),drop = T])

A <- ape::vcv.phylo(sub_tree, corr=T) # the covariance must be positive definite
diag(A) <- diag(A)+0.000001 # we should add a small value to the diagonal as an offset to avoid an error

###Now run phylotrends analyses on the thing 
mytree22<- multi2di(sub_tree)
mytree22$edge.length <- pmax(mytree22$edge.length,1/365) # use a day as minimum branch length

quartz()
plot(mytree22, cex = 0.8, tip.color = colors)

#write.csv(A, "A_for_dale")
##########################################################################
#FIT MODEL
##########################################################################

model_repeat1 <- brm(
  #phen ~ spec_mean_pc1 + spec_mean_pc2 + (1|gr(phylo, cov = A)) + (1|species),
   phen ~ scale(spec_mean_pc1) + scale(spec_mean_pc2) + (1|gr(phylo, cov = A)),
  data = data_repeat,
  family = gaussian(),
  data2 = list(A = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

summary(model_repeat1)

post <- posterior_samples(model_repeat1)

# Probability that the effect of spec_mean_pc1 is positive
prob_pc1 <- mean(post$b_spec_mean_pc1 > 0)

# Probability that the effect of spec_mean_pc2 is positive
prob_pc2 <- mean(post$b_spec_mean_pc2 > 0)

prob_pc1
prob_pc2

