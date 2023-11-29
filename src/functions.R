#' Calculate rarefied and extrapolated diversity metrics
#'
#' @name extrar_divmetrics
#' @description Function to calculate rarefied and extrapolated diversity metrics- Richness, Shannon, Simpson
#' according to Chao et al. 2014 in the context of Hill numbers from a table with species abundance data
#' @param datafile dataframe with information of species abundance per site and year. It must contain at least
#' four columns: SiteID, year, species and abundance
#' @return a dataframe with diversity metrics per site and year
#' @export
#' @import dplyr
#' @import vegan
#' @import readr
#' @import tidyr
#' @import iNEXT

#data.spp <- read_csv("data/fia/fia_species.csv")
extrar_divmetrics <- function (datafile){
  
  #datafile_year <- datafile %>%
  #mutate(sitesyearID = paste0(siteID,year))
  
  list.sites <- datafile %>%
    group_by(siteID) %>%
    summarise(list(abundance))
  
  list.sites2 <- list.sites[[2]]
  
  names(list.sites2) <- list.sites[[1]]
  
  x <- iNEXT::iNEXT(list.sites2, q=0, datatype="abundance",size=NULL)
  #ggiNEXT(x, type=1, se=TRUE, facet.var="none", color.var="none", grey=FALSE)
  
  fia.inextparams <- as.data.frame(x[3]) %>%
    dplyr::select(siteID = AsyEst.Site,
                  divmetric = AsyEst.Diversity,
                  obsmetric = AsyEst.Observed,
                  estmetric = AsyEst.Estimator,
                  se.metric = AsyEst.s.e.)
  
  rich.df <- fia.inextparams %>%
    filter(divmetric=="Species richness") %>%
    dplyr::select(siteID,
                  q0.obs = obsmetric,
                  q0.est = estmetric,
                  q0.se = se.metric)
  
  shannon.df <- fia.inextparams %>%
    filter(divmetric=="Shannon diversity") %>%
    dplyr::select(siteID,
                  q1.obs = obsmetric,
                  q1.est = estmetric,
                  q1.se = se.metric)
  
  simpson.df <- fia.inextparams %>%
    filter(divmetric=="Simpson diversity") %>%
    dplyr::select(siteID,
                  q2.obs = obsmetric,
                  q2.est = estmetric,
                  q2.se = se.metric)
  
  div.params <- rich.df %>%
    left_join(.,shannon.df, by="siteID") %>%
    left_join(.,simpson.df,by="siteID")
  
  dataname <-  deparse(substitute(datafile))
  write_csv(div.params,paste0("data-raw/",dataname,"_inext_params.csv"))
}

################## MULTIDATE##########

multidate <- function(data, formats){
  a<-list()
  for(i in 1:length(formats)){
    a[[i]]<- as.Date(data,format=formats[i])
    a[[i]][a[[i]]>Sys.Date() | a[[i]]<as.Date("1000-01-01")]<-NA
    a[[1]][!is.na(a[[i]])]<-a[[i]][!is.na(a[[i]])]
  }
  a[[1]]
}


