bestmod <- function(data, f, threshval = NULL) { 
  
  filter(validtab, threshval == threshval) %>% 
              group_by(lineID) %>% 
              slice_max(!!as.name(f), with_ties = F, n = 1) %>% 
              select(lineID, kclass) %>% 
    mutate(bestmod = factor(kclass),
           bestmodtype = factor(ifelse(bestmod == 'nospw', 'nospw', 'spw'))) %>% 
    select(-kclass) 

}
