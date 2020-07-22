### This is a short script to make the subsetted data from  Roswell et al. 2019 PLoS One https://doi.org/10.5061/dryad.c3rr6q1 and used for figures in Roswell et al. 2020 Oikos.



#read in data from .csv on Dryad
XerDate<-read.csv("https://datadryad.org/stash/downloads/file_stream/93833") %>%
  dplyr::filter(keep == "Y" &
           bee_species != "sp" &
           bee_species != "species" &
           bee_genus != "sandwasp" &
           bee_genus != "sand wasp" &
           bee_genus != "Anacrabro") %>% #note non-bee accidentally uploaded
  tidyr::unite(col="sr", site, sampling_round, remove=F) %>%
  tidyr::unite(col="siteday",site, sday, remove=F) %>% droplevels()

XerDay<-XerDate %>%
  dplyr::group_by(sr, sampling_round, site) %>%
  dplyr::mutate(sday1=dplyr::dense_rank(sday))



set.seed(5062, sample.kind = "Rounding") #due to turnover in set.seed() with R 3.6.0

makedat1 <- XerDay %>%
  dplyr::group_by(sr, sday) %>%
  dplyr::do(filter(.,length(unique(boutstart)) > 4)) %>% #drop dates with unsufficient sampling effort
  dplyr::do(dplyr::filter(., boutstart %in% sample(unique(boutstart)
                                                   , size=5 #select 5 bouts from each day of sampling at each site (equal effort)
                                                   , replace=F))) %>%
  dplyr::ungroup()



beeObs<-makedat1 %>%
  dplyr::select(uniqueID, bee, start=boutstart, sr, sday1, site, sampling_round) %>%
  # FOUR SITE-ROUNDS SELECTED TO DEMONSTRATE PITFALLS IN DIVERSITY COMPARISON
  dplyr::filter(sr %in% c(
    "IAS_3"
    , "Cold Soil_5"
    , "Lord Stirling_4"
    , "Fox Hill_5"
    ))


beeAbunds <- beeObs %>% #summary table used in MS
  dplyr::group_by(sr, bee) %>%
  dplyr::summarize(abund=n()) %>%
  tidyr::spread(sr, abund, fill=0)

usethis::use_data(beeObs, overwrite = TRUE)
usethis::use_data(beeAbunds, overwrite = TRUE)
