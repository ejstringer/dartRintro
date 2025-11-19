library(dartRverse)
library(tidyverse)

tympo <- gl.read.dart("Report_DTym25-13579_SNP_2.csv",
                      ind.metafile = "Tympo_metadata.csv")


# filter ------------------------------------------------------------------
gl.smearplot(tympo)

tympo_filtered <- gl.filter.callrate(tympo, threshold = 0.95,method = "loc")
tympo_filtered <- gl.filter.rdepth(tympo_filtered,lower = 5, upper=50)
tympo_filtered <- gl.filter.callrate(tympo_filtered, threshold = 0.95, 
                                     method = "ind")
tympo_filtered <- gl.filter.reproducibility(tympo_filtered)


nInd(tympo_filtered)
nLoc(tympo_filtered)

gl.report.monomorphs(tympo_filtered)
gl.report.secondaries(tympo_filtered)
gl.report.maf(tympo_filtered)

tympo_filtered <- gl.filter.secondaries(tympo_filtered)

## pcoa --------------------------------------------------------------------
pop(tympo_filtered) <- tympo_filtered@other$ind.metrics$pop
pc <- gl.pcoa(tympo_filtered)
gl.pcoa.plot(pc, tympo_filtered, yaxis = 1, xaxis = 2)



## diversity ---------------------------------------------------------------

pop(tympo_filtered) <- paste(tympo_filtered@other$ind.metrics$pop,
                                tympo_filtered@other$ind.metrics$year,
                                sep = '-')




ar <- gl.report.allelerich(tympo_filtered)
ar$`Richness per population` %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
  #  filter(popsize>1) %>% 
  ggplot(aes(year, mean_richness, colour = pop))+
  geom_point(aes(size = popsize))+
  geom_smooth(method = 'lm', aes(group=pop))+
  theme_classic()

h <- gl.report.heterozygosity(tympo_filtered)
h %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
  ggplot(aes(year, He, colour = pop))+
  geom_point()+
  geom_smooth(method = 'lm', aes(group = pop))+
  theme_classic()


