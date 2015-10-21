
dplyrConvert <- function(){

  data(stockholm,package = "prostateCancerStockholm")
  library(tidyr)
  library(dplyr)
  library(Biobase)
  geoData <- camcap
  
  zscore <- function(x) t(scale(t(x)))
  
  iqrs <- apply(log2(exprs(geoData)), 1, IQR,na.rm=TRUE)
  
  camcap <- tbl_df(data.frame(Probe=featureNames(geoData),IQR=iqrs,zscore(exprs(geoData)))) %>%
    gather(geo_accession,Expression,- c(Probe,IQR))
  
  fd <- tbl_df(fData(geoData)) %>% rename(Probe = ID,Gene=Symbol) %>% 
    select(Probe, Gene, Entrez_Gene_ID)
  
  pd <- tbl_df(pData(geoData))
  

  camcap <- camcap %>% full_join(pd) %>% full_join(fd)

  ##selecting most variable probe for each gene

    varProbes <- camcap %>% group_by(Gene) %>% 
    summarise(Probe = Probe[which.max(IQR)])
  
  camcap <- inner_join(camcap, varProbes,by="Probe") %>% rename(Gene = Gene.x) %>% select(-c(Gene.y))
  camcap
}
