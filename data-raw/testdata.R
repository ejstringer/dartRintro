library(dartRverse)
library(tidyverse)

# data --------------------------------------------------------------------
prjdir <- getwd()
setwd('../dartR.intro/inst/extdata/')

tympo.gl <- gl.read.dart('Report_DTym25-13579_SNP.csv',
                         ind.metafile = 'Tympo_metadata.csv')
tympo.gl@other$history
usethis::use_data(tympo.gl, overwrite = TRUE)


# test silico

tympo.silico <- gl.read.silicodart('Report_DTym25-13579_SilicoDArT.csv',
                                   ind.metafile = 'Tympo_metadata.csv')

setwd(prjdir)


# CHECKS ------------------------------------------------------------------
#gl.subsample.ind(tympo.gl, n = 20)
tympo.gl
gl.report.monomorphs(tympo.gl)
## smearplot ---------------------------------------------------------------
gl.smearplot(tympo.gl)

## filter ------------------------------------------------------------------

tympo.gl_filtered <- gl.filter.callrate(tympo.gl, 
                                        threshold = 0.95,method = "loc")
tympo.gl_filtered <- gl.filter.rdepth(tympo.gl_filtered,lower = 5, upper=50)
tympo.gl_filtered <- gl.filter.callrate(tympo.gl_filtered, threshold = 0.95,method = "ind")
tympo.gl_filtered <- gl.filter.reproducibility(tympo.gl_filtered)


nInd(tympo.gl_filtered)
nLoc(tympo.gl_filtered)

gl.report.monomorphs(tympo.gl_filtered)
gl.report.secondaries(tympo.gl_filtered)

tympo.gl_filtered <- gl.filter.secondaries(tympo.gl_filtered)

## pcoa --------------------------------------------------------------------
pop(tympo.gl_filtered) <- tympo.gl_filtered@other$ind.metrics$pop
popNames(tympo.gl_filtered)
#levels(pop(tympo.gl_filtered)) <- c('Cookanalla', 'North Canberra', 'Jerra East', 'Jerra West', 'Unknown')
pc <- gl.pcoa(tympo.gl_filtered)
p<-gl.pcoa.plot(pc, tympo.gl_filtered, yaxis = 2, xaxis = 1)



# drop and reassign -------------------------------------------------------


pcdf <- pc$scores %>% as.data.frame() %>% 
  mutate(id = rownames(.)) %>% 
  left_join(tympo.gl_filtered@other$ind.metrics[1:2])

### find inds to delete and inds to reassign. 

## to delete: 

pcdf %>% 
  filter(PC1 < -1, PC2 <1,
         pop != 'Tuggeranong') #AA24149

pcdf %>% 
  filter(PC1 > -1, PC2 <1, 
         pop != 'Googong', pop != 'Royalla') #AA24002 AA24001

## reassign

pcdf %>% 
  filter(PC1 > 1, PC2 >1, 
         pop != 'Kowen') # AA24117 AA24155 AA24250 AA24545





## diversity ---------------------------------------------------------------

pop(tympo.gl_filtered) <- paste(tympo.gl_filtered@other$ind.metrics$pop,
                                tympo.gl_filtered@other$ind.metrics$year,
                                sep = '-')




ar <- gl.report.allelerich(tympo.gl_filtered)
ar$`Richness per population` %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
  #  filter(popsize>1) %>% 
  ggplot(aes(year, mean_richness, colour = pop))+
  geom_point(aes(size = popsize))+
  geom_smooth(method = 'lm', aes(group=pop))+
  theme_classic()

h <- gl.report.heterozygosity(tympo.gl_filtered)
h %>% 
  tidyr::separate(pop, into = c('pop', 'year'), sep = '-') %>% 
  ggplot(aes(year, He, colour = pop))+
  geom_point()+
  geom_smooth(method = 'lm', aes(group = pop))+
  theme_classic()

###


# silico ------------------------------------------------------------------



silico_locind_filt <- gl.filter.callrate(tympo.silico, threshold = 0.95,method = "loc")
silico_locind_filt <- gl.filter.rdepth(silico_locind_filt,lower = 5, upper=50)
silico_locind_filt <- gl.filter.callrate(silico_locind_filt, threshold = 0.95,method = "ind")
silico_locind_filt <- gl.filter.reproducibility(silico_locind_filt)


silico_locind_filt@pop %>% table


pcsil <- gl.pcoa(silico_locind_filt)
gl.pcoa.plot(pcsil, silico_locind_filt, yaxis = 2, xaxis = 1)

