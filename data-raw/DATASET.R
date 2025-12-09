## code to prepare `DATASET` dataset goes here


library(dartRverse)
library(tidyverse)
#load all the data
# https://r-pkgs.org/data.html#sec-data-data-raw

# load --------------------------------------------------------------------
tympos1 <- gl.read.dart("./data-raw/SNPs_original.csv", 
                        ind.metafile = "./data-raw/GED2025_corr.csv")
silico <- gl.read.silicodart(filename ="./data-raw/Report_DTym25-10518_10_moreOrders_SilicoDArT_1.csv",
                             ind.metafile = './data-raw/GED2025_corr.csv')

tym <- gl.filter.secondaries(tympos1)


1-(nLoc(tym)/nLoc(tympos1))

seconds <- which(!tympos1@loc.names %in% tym@loc.names)

table(tympos1@other$loc.metrics$clone %in% tympos1@other$loc.metrics$clone[seconds])

distsec <- round(table(table(tympos1@other$loc.metrics$clone))*0.12)
clonetb <- table(tympos1@other$loc.metrics$clone)
clonelist <- list()
for (i in 1:5) {
  clonelist[[i]]<-names(clonetb)[clonetb == i][1:distsec[i]]
}

indexsec <- which(tympos1@other$loc.metrics$clone %in% unlist(clonelist))
length(indexsec)
# downsample loci ---------------------------------------------------------

nLoc(tympos1)
systematic_sample <- seq(1, nLoc(tympos1), length.out = 3210)

systematic_sample <- c(1:5132,seconds[1:5000])

systematic_sample <- indexsec
tympos <- gl.keep.loc(tympos1, loc.list = tympos1@loc.names[systematic_sample])
tympos

gl.report.secondaries(tympos)

# wild pops only ----------------------------------------------------------

index <- tympos@other$ind.metrics$Colony=="Wild" 
index <- ifelse(is.na(index), FALSE, index)
tw <- tympos[index,]
index2 <- tw@other$ind.metrics$group !="Monaro"
tw2 <- tw[index2,]
#tw2_5 <- gl.drop.ind(tw2, ind.list =  "AA61626")
tw2_5 <- gl.drop.ind(tw2, ind.list =  "CK1 hatchling")
tw_final <- tw2_5
nInd(tw_final)
tw_final <- gl.filter.monomorphs(tw_final)
nLoc(tw_final)
# generic pop names -------------------------------------------------------

pop(tw_final) <- factor(ifelse(grepl('West', tw_final@other$ind.metrics$pop),
                               'W_Can',as.character(tw_final@other$ind.metrics$group)))

pop(tw_final) <- factor(ifelse(grepl('Cook', tw_final@other$ind.metrics$pop),
                               'E_Can', as.character(tw_final@pop)))
table(tw_final@pop)
tw_final@other$ind.metrics[tw_final@pop == 'unknown',]


# random lat lon ----------------------------------------------------------

gl.map.interactive(possums.gl)
poplocation <- cbind(id = possums.gl@ind.names, 
                     pop = possums.gl@pop, 
                     possums.gl@other$latlon) %>% 
  filter(pop %in% LETTERS[c(2:4,5)]) %>% 
  mutate(pop = case_when(
   # pop == 'F' ~ 'D',
    pop == 'C' ~ 'B',
    .default = pop
  )) %>% 
  mutate(pop = case_when(
    pop == 'D' ~ 'N_Can',
    pop == 'E' ~ 'S_Can',
    pop == 'B' ~ 'W_Can'
  )) %>% 
  group_by(pop) %>% 
  summarise(min_lat = min(lat),
            max_lat = max(lat),
            min_lon = min(lon),
            max_lon = max(lon)) 

idpop <- data.frame(id = tw_final@ind.names, pop = tw_final@pop)

poplocationE <- poplocation %>%
  rbind(data.frame(pop = 'E_Can', 
                   min_lat = poplocation$min_lat[3],
                   max_lat = poplocation$max_lat[3],
                   min_lon = poplocation$max_lon[2]-0.03,
                   max_lon = poplocation$min_lon[1]-0.03))

individuals<-left_join(idpop, poplocationE)
table(tw_final@ind.names == individuals$id)
individuals<-individuals %>% 
  #filter(complete.cases(min_lat)) %>% 
  rowwise() %>% 
  mutate(#lat =runif(1, min_lat, max_lat),
    #lon = runif(1, min_lon, max_lon),
    lat = rnorm(1, mean = mean(c(min_lat, max_lat)), 0.01),
    lon = rnorm(1, mean = mean(c(min_lon, max_lon)), 0.01),
    lat = ifelse(is.nan(lat), NA, lat),
    lon = ifelse(is.nan(lon), NA, lon))

individuals %>% names

tw_final@other$latlon$lat <- individuals$lat
tw_final@other$latlon$lon <- individuals$lon

gl.map.interactive(tw_final)

# new pop with unique alleles ---------------------------------------------

gl.alf(tw_final) %>% filter(alf2>0, alf2 < 1) %>%  arrange(alf2) %>% head

mat <- as.matrix(tw_final)

which(mat[, '13662966-58-C/T'] > 0)

countHet <- apply(mat, 2, function(x)  sum(x== 1, na.rm = T))
diverseInd <- names(sort(rowSums(mat[,names(countHet[countHet==1])]==1, na.rm =T),decreasing = T)[4])

tw_final@other$ind.metrics[tw_final@other$ind.metrics$id==diverseInd,]


gl.drop.ind(tw_final,diverseInd, mono.rm=TRUE)


# tidy metadata -----------------------------------------------------------

tw_final@other$ind.metrics %>% head

tw_final@other$ind.metrics$pop <- tw_final@pop
tw_final@other$ind.metrics$lat <- tw_final@other$latlon$lat
tw_final@other$ind.metrics$lon <- tw_final@other$latlon$lon
tw_final@other$ind.metrics$age <- tw_final@other$ind.metrics$Age

tw_final@other$ind.metrics[,c('id', 'pop', 'lat', 'lon', 'year', 'sex', 'age')] %>% head
tw_final@other$ind.metrics <- tw_final@other$ind.metrics[,c('id', 'pop', 'lat', 'lon', 'year', 'sex', 'age')]

tw_final@other$ind.metrics %>% head

tw_final@other$ind.metrics[tw_final@other$ind.metrics$id == diverseInd,]
tw_final@other$ind.metrics[tw_final@other$ind.metrics$id == diverseInd,-1] <- NA
tw_final@other$ind.metrics$pop[tw_final@other$ind.metrics$id == diverseInd] <- "unknown"
tw_final@other$ind.metrics$year[tw_final@other$ind.metrics$id == diverseInd] <- "unknown"
tw_final@other$ind.metrics[tw_final@other$ind.metrics$id == diverseInd,]


# silico file ------------------------------------------------------------
keepsil <- which(silico@loc.names %in% tw_final@other$loc.metrics$CloneID)
length(keepsil)
length(unique(c(1:4010, keepsil)))

silico_loc <- gl.keep.loc(silico, silico@loc.names[unique(c(1:4000, keepsil))])
silico_locind <- gl.keep.ind(silico_loc, ind.list = tw_final@ind.names)

popid <- tw_final@other$ind.metrics[,c('id', 'pop')] %>% 
  rename(popSNP= pop)

silico_locind@other$ind.metrics <-silico_locind@other$ind.metrics %>% 
  left_join(popid)

silico_locind <- gl.report.monomorphs(silico_locind)

# csv keep files ----------------------------------------------------------


loc_data <- read.csv('./data-raw/SNPs_original.csv',
                     header = F)
table(loc_data$V2 %in% tw_final@other$loc.metrics$clone)/2
## ids are in row 7
loc_data[1:9, 30:40]
row6 <- loc_data[7,]

## find which column ids start
locmet <- which(row6 == 'RepAvg')
indkeep <- which(row6 %in% tw_final@ind.names) 

loc_data_keep <- loc_data[,c(1:locmet, indkeep)]

# dataframe to get new ids in same order as SNP csv
newids <- data.frame(id = tw_final@other$ind.metrics$id[order(tw_final@other$ind.metrics$year)],
                     id2 = paste0('AA', 24001:(24000+nInd(tw_final)))) 



new_loc_ids <- data.frame(id = unlist(loc_data_keep[7,-c(1:24)]))%>% 
  left_join(newids)

new_loc_ids[new_loc_ids$id %in% new_loc_ids$id[duplicated(new_loc_ids$id)],]
# assign new ids
table(loc_data_keep[7,-c(1:24)]==new_loc_ids$id)
loc_data_keep[7,-c(1:24)] <- new_loc_ids$id2

tw_final@other$ind.metrics$id[order(tw_final@other$ind.metrics$year)] <- newids$id2
tw_final@ind.names <- tw_final@other$ind.metrics$id

new_loc_ids$id %>% duplicated %>% table

## Remove loci
#keeploc <- which(loc_data$V2 %in% tw_final@other$loc.metrics$clone)
locdataV1 <- data.frame(V1 = loc_data$V1) %>%  
  mutate(v = V1,
         v = sub('\\|F\\|0--', '_', v),
         v = sub('\\|F\\|0-', '_', v),
         v = gsub('>', '_', v),
         v = gsub(':', '_', v),
         v = sub('-', '_', v)) %>% 
  separate(v, sep = '_', into = LETTERS[1:4]) %>% 
  mutate(nucliotide = paste(C, D, sep = '/'),
         locnames = paste(A, B, nucliotide, sep = '-')) %>% 
  select(V1, locnames)

keeploc <- which(locdataV1$locnames %in% tw_final@loc.names)

loc_data_keep_s <-loc_data_keep[c(1:7,keeploc),]

## remove duplicated ids
dups_in_keep <- duplicated(as.vector(loc_data_keep_s[7,]))

loc_data_keep_nodups <- loc_data_keep_s[,!dups_in_keep]
# save new data

(nrow(loc_data_keep_nodups)-7)/2
nLoc(tw_final)

# write.csv(loc_data_keep_nodups,
#           './data-raw/Report_DTym25-13579_SNP.csv',
#           row.names = F) 
# 
# ## DELETE top row of file and move to extdata

write.table(loc_data_keep_nodups,
          './inst/extdata/Report_DTym25-13579_SNP.csv', sep = ',',
          row.names = F, col.names = F) 


# csv keep silico  --------------------------------------------------------

loc_data <- read.csv('./data-raw/Report_DTym25-10518_10_moreOrders_SilicoDArT_1.csv',
                     header = F)
table(loc_data$V1 %in% silico_locind@other$loc.metrics$CloneID)
## ids are in row 7
loc_data[1:9, 30:40]
row6 <- loc_data[7,]

## find which column ids start
locmet <- which(row6 == 'Reproducibility')
indkeep <- which(row6 %in% silico_locind@ind.names) 

loc_data_keep <- loc_data[,c(1:locmet, indkeep)]

# dataframe to get new ids in same order as SNP csv

new_loc_ids

new_loc_ids[new_loc_ids$id %in% new_loc_ids$id[duplicated(new_loc_ids$id)],]
# assign new ids
table(loc_data_keep[7,-c(1:locmet)]==new_loc_ids$id)
loc_data_keep[7,-c(1:locmet)] <- new_loc_ids$id2


new_loc_ids$id %>% duplicated %>% table

## Remove loci
keeploc <- which(loc_data$V1 %in% silico_locind@loc.names)


#keeploc <- which(locdataV1$locnames %in% tw_final@loc.names)

loc_data_keep_s <-loc_data_keep[c(1:7,keeploc),]

## remove duplicated ids
dups_in_keep <- duplicated(as.vector(loc_data_keep_s[7,]))

loc_data_keep_nodups <- loc_data_keep_s[,!dups_in_keep]
# save new data

#(nrow(loc_data_keep_nodups)-7)/2
nLoc(silico_locind)


write.table(loc_data_keep_nodups,
            './inst/extdata/Report_DTym25-13579_SilicoDArT.csv', sep = ',',
            row.names = F, col.names = F) 


# metadata ----------------------------------------------------------------
metaweights <- tw_final@other$ind.metrics %>% 
  mutate(svlmin = case_when(
    age == 'A' ~ 48,
    age == 'SA' ~ 38,
    age == 'J' ~ 20
  ),
  svlmax = case_when(
    age == 'A' ~ 85,
    age == 'SA' ~ 48,
    age == 'J' ~ 38
  ),
  mweight = case_when(
    age == 'A' ~ 5.5,
    age == 'SA' ~ 4,
    age == 'J' ~ 3
  )) %>% 
  rowwise() %>% 
  mutate(svl = runif(1, min = svlmin, max = svlmax),
         weight = rnorm(1, mweight, 0.5)) %>% 
  mutate(svl = ifelse(is.nan(svl),NA, svl),
         weight = ifelse(weight < 1, 1.13, weight),
         weight = (svl*0.01+1)*weight)

tw_final@other$ind.metrics$svl = round(metaweights$svl,1)
tw_final@other$ind.metrics$weight = round(metaweights$weight,2)
tw_final@other$ind.metrics$age <- factor(tw_final@other$ind.metrics$age)
tw_final@other$ind.metrics$species <- 'Tympanocryptis lineata'

boxplot(metaweights$svl ~ metaweights$age)
boxplot(metaweights$weight ~ metaweights$age)

plot(metaweights$svl, metaweights$weight, col = metaweights$age)
summary(lm(metaweights$weight~metaweights$svl+metaweights$age))



levels(tw_final@other$ind.metrics$pop) <- c('Googong', 'Kowen', 'Royalla', 'Unknown', 'Tuggeranong')

pop(tw_final) <- tw_final@other$ind.metrics$pop
tw_final@pop
head(tw_final@other$ind.metrics)
table(tw_final@pop)

write.csv(tw_final@other$ind.metrics, './inst/extdata/Tympo_metadata.csv',
          row.names = F)

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

