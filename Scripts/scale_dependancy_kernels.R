#
# Scale dependancy analysis
#
# A. Boyé 20/08/2021
#
#------------

## -------------
# Packages
## -------------

# Data manipulation
library(tidyverse)

# Analysis
library(ade4)
library(BAT)
library(vegan)

# Parallelisation
library(furrr)
plan(multiprocess)

set.seed(64)

## -------------
# Data
## -------------

load("data/inv.mat.Rdata")
load("data/traits_v.Rdata")
load("data/spatial_scales.Rdata")

## -------------
# Trait space
## -------------

#Prepare the fuzzy coded dataframe
blocks <- c(5, 4, 3, 2, 5, 3, 4, 3, 3, 4, 3, 4, 6, 3, 4, 3)

names(blocks) <- c("Max_size","Max_desc","Max_cycles","Nb_repro_cycles","Life_duration",
                   "Repro_tech","Ovopo_site","Egg/Egg mass","Diss_pot","Attach_substr",
                   "Body_flex","Body form","Feed_hab","Diet_pref","Resp_aqua_stages","Aqua_stages")  

w <- prep.fuzzy.var(as.data.frame(traits_v), blocks)

#Run the FCA
fca1 <- dudi.fca(w, scannf = FALSE, nf = ncol(w))

# Plot the FCA
biplot(fca1)

# Retrieve species and modality loadings
sp_loadings <- fca1$li %>% rownames_to_column("Species")
modality_loadings <- fca1$co %>% rownames_to_column("Modality")

# Format the trait_space
inv_trait <- sp_loadings %>%
  column_to_rownames("Species")

## -------------
# Community matrix
## -------------

inv_com <- inv.mat %>%
  gather(Species, Abundance, -Site, -Years, -Attribute) %>%
  filter(Species %in% rownames(inv_trait)) %>%
  right_join(Spatial_scales, .) %>%
  mutate(Country = "New-Zealand")

## -------------
# Function to apply to each scale
## -------------

compute_kernel_indic <- function(scale, n.dim = 5, abundance = FALSE){
  
  # Create directory to save
  #-------------------------
  
  if(abundance){
    label <- paste(scale, "weighted", sep="_")
  }else{
    label <- paste(scale, "unweighted", sep="_")
  }
  
  directory_name <- paste("outputs",label, sep="/")
  
  if(!label %in% dir("outputs")){
    dir.create(directory_name)
  }
  
  # Aggregate data to the defined scale
  #------------------------------------
  
  inv_com_mat <- inv_com %>%
    group_by(across(c({{ scale }}, "Years", "Species"))) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    spread(Species, Abundance, fill = 0) %>%
    unite(Row_name, c(scale, Years)) %>%
    column_to_rownames("Row_name")

   inv_com_mat_pa <- inv_com %>%
      group_by(across(c({{ scale }}, "Years", "Species"))) %>%
      summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>%
      ungroup() %>%
	  mutate(Abundance = if_else(Abundance >0, 1, 0)) %>%
      spread(Species, Abundance, fill = 0) %>%
      unite(Row_name, c(scale, Years)) %>%
      column_to_rownames("Row_name")

  # Remove samples with less species than the dimensions of the trait space
  #------------------------------------------------------------------------
  
  inv_com_mat <- inv_com_mat[rowSums(inv_com_mat_pa)>n.dim,colSums(inv_com_mat)>0]
  
  # Compute the kernels and kernels hotspots
  #-----------------------------------------
  
  trait_kernel <- kernel.build(comm =inv_com_mat,
                               trait = inv_trait[,1:n.dim], 
                               method = "gaussian",
                               abund = abundance,
							   cores = 1)
  
  trait_kernel_hotspots <- kernel.hotspots(trait_kernel, prop = 0.5)

  # Save it
  object_name <- c(paste(label, "trait_kernel", sep ="_"),paste(label, "trait_kernel_hotspots", sep ="_"))
  assign(object_name[1],trait_kernel)
  assign(object_name[2],trait_kernel_hotspots)
  save(list = object_name, file = paste(directory_name,paste(label, "trait_kernels.Rdata", sep = "_"),sep="/"))
  
  # Compute alpha diversity
  #------------------------
  
  functional_richness <- future_map_dbl(trait_kernel@HVList, kernel.alpha)
  functional_richness_hotspots <- future_map_dbl(trait_kernel_hotspots@HVList, kernel.alpha)
  functional_evenness <- future_map_dbl(trait_kernel@HVList, kernel.evenness)
  functional_evenness_hotspots <- future_map_dbl(trait_kernel_hotspots@HVList, kernel.evenness)
  
  trait_kernel_alpha <- data.frame(functional_richness,
                                   functional_richness_hotspots,
                                   functional_evenness,
                                   functional_evenness_hotspots)
  
  # Save it
  object_name <- paste(label, "trait_kernel_alpha", sep ="_")
  assign(object_name,trait_kernel_alpha)
  save(list = object_name, file = paste(directory_name,paste(label, "trait_kernel_alpha.Rdata", sep = "_"),sep="/"))
  
  # Compute redundancy
  #-------------------
  
  # Compute raw originality 
  trait_kernel_origin <- future_map(trait_kernel@HVList, kernel.originality, relative = FALSE) %>%
    future_map_dfr(.,~as.data.frame(t(.)))
  
  rownames(trait_kernel_origin) <- rownames(inv_com_mat)
  
  # Save it
  object_name <- paste(label, "trait_kernel_origin", sep ="_")
  assign(object_name,trait_kernel_origin)
  
  save(list = object_name, file = paste(directory_name,paste(label, "trait_kernel_origin.Rdata", sep = "_"),sep="/"))
  
}

param_grid <- expand.grid(scale = rev(c("Country", "Catchment", "Ecoregions", "Island", "Site")),
            abundance = rev(c(TRUE, FALSE)))

scale <- as.list(param_grid$scale)
abundance <- as.list(param_grid$abundance)

future_pmap(list(x=scale,y=abundance), ~compute_kernel_indic(scale = .x, abundance = .y))
