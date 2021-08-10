# fun_fg_spp - Function to generate species attributions to functional units
# fun_met - Function to calculate functional metrics: UTC, functional redundancy

## --------------------------------------------------------------
## Name: fun-fg-spp-function.R
## Description: Function to generate species attributions to functional units
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: October 2016
## Outputs: Function named 'fun_fg_spp'
## Args:
## tdf = Dataframe of trait data in bins
## em = Species x Site matrix
## e = Number of ecoregions assessed
## name = Name of the method to attach to objects being saved
## --------------------------------------------------------------

fun_fg_spp <- function(tdf, em, e, name) {
  
  ## gp ##
  
  # dataframe of all unique trait combinations for the trait data
  gp <- tdf %>% 
    distinct %>%
    # unique trait combinations
    mutate(group = 1:nrow(.))
  # add an identifying column called group to the dataframe
  
  ## fg_spp ##
  
  tdf <- tdf %>%  
    rownames_to_column(var = "species")
  # convert rownames to explicit column for use in join
  
  # dataframe of species attributions to functional units
  fg_spp <- inner_join(tdf, gp)
  
  # dataframe of number of species per functional unit
  count <- table(fg_spp$group) %>% 
    # count species per functional unit
    as.data.frame(., stringsAsFactors = FALSE) %>% 
    # convert to dataframe
    mutate(Var1 = as.numeric(Var1))
  # convert to numeric to allow join
  
  # join number of species per functional unit to species attributions
  fg_spp <- fg_spp %>% 
    left_join(count, by = c("group" = "Var1")) %>% 
    # join count to dataframe
    dplyr::select(species, group, count = Freq, everything())
  # reorder columns
  
  # save: fg_spp
  saveRDS(fg_spp, paste0("FctTradeoffs/df_fg_spp_", name, ".rds"))
  # species attributions to functional groups
  
  return(fg_spp)
  
}

## --------------------------------------------------------------
## Name: fun-met-function.R
## Description: Function to calculate functional metrics: UTC, functional redundancy
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: October 2017
## Outputs: Function named 'fun_met'
## Args:
## tdf = Dataframe of trait data in bins
## em = Matrix of Species x Site data
## fd = Dataframe of functional dispersion results
## e = Number of ecoregions assessed
## name = Name of the method to attach to objects being saved
## --------------------------------------------------------------

fun_met <- function(tdf, em, fd, e, name) {
  
  #### Generate species attributions to functional units ####
  
  ## fg_spp ##
  
  # run fun_fg_spp function
  fg_spp <- fun_fg_spp(tdf = tdf, em = em, e = e, name = name)
  
  #### Species richness per ecoregion ####
  
  eco_spp <- table(eco$eco) %>% 
    # species richness per ecoregion
    as.data.frame %>% 
    # convert table to dataframe
    dplyr::rename(eco = Var1, spp = Freq) %>% 
    # rename columns
    mutate(eco = as.character(eco))
  
  #### Unique Trait Combinations ####
  
  ## utc_pool ##
  
  # calculate utc for the species pool and use this for scaling purposes
  utc_pool <- mvfd(as.matrix(tdf), em, 0, calc.ovr = 0) 
  # 0 is for resolution
  # doesn't run with overlap on
  pool_traitspace <- utc_pool$utc[1]
  # full trait space
  
  ## utc_eco ##
  
  # calculate sutc (scaled utc) for every ecoregion
  utc_eco <- mvfd(as.matrix(tdf), em, 0, traitspace = pool_traitspace, calc.ovr = 0) 
  
  ## utc_df ##
  
  # create dataframe of ecoregions, sutc, utc and spp
  utc_df <- as.data.frame(cbind(rownames(em), utc_eco$utc), stringsAsFactors = FALSE) %>% 
    # extract data from mvfd function into dataframe
    dplyr::rename(eco = V1, utc = V2) %>% 
    # rename columns of dataframe
    mutate(across(.cols=utc, .fns=as.numeric))
  #vars=vars(utc), funs=as.numeric)
  
  # convert variables to numeric
  
  ## fr ##
  fr <- left_join(eco_spp, utc_df, by = "eco") %>%
    # join species richness to UTC data
    mutate(fred = spp/utc) %>%
    # calculate functional redundancy %>% 
    left_join(fd, by = "eco")
  # join functional dispersion data
  
  # save: fr
  saveRDS(fr, paste0("FctTradeoffs/df_fr_", name, ".rds"))
  # dataframe of functional metrics standardised between ecoregions (data.frame)
  
  return(fr)
  
}
