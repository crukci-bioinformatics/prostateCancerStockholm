
dplyrConvert <- function(){

  data(stockholm,package = "prostateCancerStockholm")
  library(tidyr)
  library(dplyr)
  library(Biobase)
  geoData <- stockholm
  
  zscore <- function(x) t(scale(t(x)))
  
  iqrs <- apply(log2(exprs(geoData)), 1, IQR,na.rm=TRUE)
  
  mu <-   apply(log2(exprs(geoData)), 1, mean,na.rm=TRUE)
  sd <- apply(log2(exprs(geoData)), 1, sd,na.rm=TRUE)
  
  camcap <- tbl_df(data.frame(Probe=featureNames(geoData),IQR=iqrs,mu,sd,log2(exprs(geoData)))) %>%
    gather(geo_accession,Expression,- c(Probe,IQR))
  
  fd <- tbl_df(fData(geoData)) %>% rename(Probe = ID,Gene=Symbol) %>% 
    select(Probe, Gene, Entrez_Gene_ID)
  
  pd <- tbl_df(pData(geoData))
  
  
  camcap <- camcap %>% full_join(pd) %>% full_join(fd)

  ##selecting most variable probe for each gene

}
