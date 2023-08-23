#Analyses reported in the review article 'Broadening applications of stochastic patch occupancy models over three decades'

#Code authors: Crone E, Gutierrez-Arellano C, Hodgson J
#Date: August 2023

#INDEX
# Species Numeralia                                               (line 29) 
# Species IUCN status                                             (line 58) 
#H1. Conservation and management studies have increased with time (line 87)
#H2. a. Taxonomic diversity has increased with time               (line 111)
#    b. There is a taxonomic bias towards insects                 (line 153)
#H3. a. Biome diversity has increased with time                   (line 162)
#    b. There is a bias towards temperate biomes                  (line 204)
#H4. Propotion of suitable habitat has varied with time           (line 216)
#    a. Landscape statistics                                      (line 258)
#    b. Temporal trends                                           (line 309)
#H5. Study length varies with interests                           (line 321)
#H6. Model complexity varies depending on interest                (line 338)

library(tidyverse)
library(rredlist)
library(MASS)
library(nnet)
library (vegan)

#Import database
db = read.csv("SPOM_database.csv", na.strings = c("NA"))

#Species numeralia --------

#No. of species by class
db_species<-db%>%
  distinct(class, species)%>%
  group_by(class)%>%
  summarise(total_cases = n(), 
            percentage = (n() / nrow(db)) * 100)

sum(db_species$total_cases)
#106 species


#Percentage of cases by class
db_class<-db%>%
  group_by(class)%>%
  summarise(total_cases = n(), 
            percentage = (n() / nrow(db)) * 100)
#Insects constitute a ~quarter of all study cases


#Percentage of cases by group
db_group<-db%>%
  group_by(Group)%>%
  summarise(total_cases = n(), 
            percentage = (n() / nrow(db)) * 100)
#Butterflies constitute 18.5% of all study cases


#Species IUCN Red List Categories-------

# NB: requires token to run (request from IUCN Red List API)

species<-as.data.frame(unique(db$species))
colnames(species)<-'species'

category<-data.frame()
for (i in species$species){
  err<-try({
    cat<-as.data.frame(rl_search(i, key= token)[["result"]][["category"]])
    names(cat)<-'category'
    cat$species<-i
  })
  if(!inherits(err,"try-error")){
    category<-rbind(category,cat)
  }
}
#Ignore errors, these are for species not in the IUCN Red List

ref<- merge(species, category, by = "species", all.x = TRUE)

#count number of species per category ignoring 'least concern' species
ref_count<-ref%>%
  drop_na(category)%>%
  filter(category!='LC')%>%
  count(category)


# H1 ----------------------------------------------------------------------
## Conservation and management studies have increased with time

#main interest per case
#A—description/application/assessment of SPOM; 
#B—SPOMs comparison; 
#C —insight to metapopulation dynamic of species; 
#D— modelling for conservation or management purposes

db$meth1 = as.numeric(grepl("A", db$main_interest))
db$meth2 = as.numeric(grepl("B", db$main_interest))
db$methods = db$meth1+db$meth2
db$methods = as.factor(db$methods> 0)
db$basic = as.factor(grepl("C", db$main_interest))
db$conserv = as.factor(grepl("D", db$main_interest))

#Study unit: paper 

#Remove duplicates of the papers with more than one species/landscape
dbuniquepaper<-db[!duplicated(db$title,),]

summary(glm(conserv ~ year,family = binomial, data = dbuniquepaper))
#increase Estimate=0.112, p=0.0047 **

# H2 ----------------------------------------------------------------------
##a. Taxonomic diversity has increased with time

#Study unit: species by paper

#Select relevant information from database
db_div<-db %>% dplyr::select(ID, class, year, biome,methods,basic,conserv)


#Divide chronologically ordered data in 9 bins with 15 cases per bin
db_div<-db_div[order(db_div$year),]

db_div$bin<- NA
db_div$bin[1:15] <- 1 
db_div$bin[16:30] <- 2
db_div$bin[31:45] <- 3
db_div$bin[46:60] <- 4
db_div$bin[61:75] <- 5
db_div$bin[76:90] <- 6
db_div$bin[91:105] <- 7
db_div$bin[106:120] <- 8
db_div$bin[121:135] <- 9

# Median years
medyear<-db_div%>%
  group_by(bin)%>% 
  summarise(medyear=median(year))

#Shannon diversity index
db_div_mat<-as.data.frame.matrix(table(db_div[, c("bin", "class")])) 
tax_div<-as.data.frame(diversity(db_div_mat[-1], index="shannon"))
colnames(tax_div) <- c("diversity")
tax_div<- tibble::rownames_to_column(tax_div, "bin")
tax_div$bin<-as.numeric(tax_div$bin)
tax_div<-merge(tax_div, medyear, by="bin")

summary(lm(diversity~medyear,data=tax_div))
#increase Estimate:0.03164 p=0.0180 *
with(tax_div,plot(diversity~medyear))

table(db_div$class)

##b. There is a taxonomic bias towards insects (class with the highest occurrence)

#Cases divided in 'insect' and 'non-insect'
levels(db_div$class)
db_div$insect<- ifelse(db_div$class=="Insecta", TRUE, FALSE)
summary(glm(insect~conserv, data = db_div, family = binomial))
#NS, non-significant bias towards insects in methods papers


# H3 ----------------------------------------------------------------------
## a. Biome diversity has increased with time 

# Study unit: landscape by paper

#Select relevant information from database
db_div_biome<-db %>% dplyr::select(ID,title, year, study_site, no_patches, year, biome, methods,basic,conserv)

#Remove duplicates of the papers with more than one species in a unique landscape
dups<- db_div_biome[c("title","study_site", "no_patches")]
db_div_biome<-db_div_biome[!duplicated(dups),]

#Divide chronologically ordered data in 9 bins with 8/9 cases per bin
db_div_biome<-db_div_biome[order(db_div_biome$year),]

db_div_biome$bin<- NA
db_div_biome$bin[1:8] <- 1
db_div_biome$bin[9:16] <- 2
db_div_biome$bin[17:24] <- 3
db_div_biome$bin[25:32] <- 4
db_div_biome$bin[33:40] <- 5
db_div_biome$bin[41:48] <- 6
db_div_biome$bin[49:56] <- 7
db_div_biome$bin[57:65] <- 8
db_div_biome$bin[66:74] <- 9

# Median years
medyear<-db_div_biome%>%
  group_by(bin)%>% 
  summarise(medyear=median(year))

db_div_biome_mat<-as.data.frame.matrix(table(db_div_biome[, c("bin", "biome")])) 
biome_div<-as.data.frame(diversity(db_div_biome_mat[-1], index="shannon"))
colnames(biome_div) <- c("diversity")
biome_div<- tibble::rownames_to_column(biome_div, "bin")
biome_div$bin<-as.numeric(biome_div$bin)
biome_div<-merge(biome_div, medyear, by="bin",all=F)
summary(lm(diversity~medyear,data=biome_div))
#NS
with(biome_div,plot(diversity~medyear))


##b. There is a bias towards temperate biomes 

#Cases divided in 'non-temperate' and 'temperate'

levels(db_div_biome$biome)
db_div_biome$temperate<- ifelse(db_div_biome$biome=="Temperate broadleaf and mixed forests"|db_div_biome$biome=="Temperate Conifer Forests"|db_div_biome$biome=="Temperate Grasslands; Savannas and Shrublands", TRUE, FALSE)

summary(glm(temperate~conserv, data = db_div_biome, family = binomial))
#Estimate= -1.2432,p=0.019* there is a bias towards non-temperate biomes in conservation papers



# H4 ----------------------------------------------------------------------
## Landscape fragmentation has varied with time

# Study unit: landscape by paper

###Prepare 'area' variables

#infer missing total patch area
db$tot_patch_area[is.na(db$tot_patch_area)]<-  with(db[is.na(db$tot_patch_area),],mean_patch_area*no_patches)

#convert km^2 and m^2 to hectares
nafalse <- function(x) {
  x & !is.na(x)
}

 #total extent
db$totext_ha<-db$tot_extent  
db$totext_ha[db$area_unit %in% c("km","m")]<-NA
db$totext_ha[nafalse(db$area_unit == "km^2")]<- 100*db$tot_extent [nafalse(db$area_unit == "km^2")]
db$totext_ha[nafalse(db$area_unit == "m^2")]<- (db$tot_extent[nafalse(db$area_unit == "m^2")]/10000)
summary(db$totext_ha)

 #total patch area
db$totpat_ha<-db$tot_patch_area  
db$totpat_ha[db$area_unit %in% c("km","m")]<-NA
db$totpat_ha[nafalse(db$area_unit == "km^2")]<- 100*db$tot_patch_area [nafalse(db$area_unit == "km^2")]
db$totpat_ha[nafalse(db$area_unit == "m^2")]<- (db$tot_patch_area[nafalse(db$area_unit == "m^2")]/10000)

#Calculate fragmentation —proportion of patch/extent areas 
plot(db$totext_ha,db$totpat_ha)
plot(log10(db$totext_ha),log10(db$totpat_ha))
db$frag<-(db$totpat_ha/db$totext_ha)


#Select relevant information from database
dbvar<-db %>% dplyr::select(ID, title, year, study_site, no_patches, totext_ha, totpat_ha,frag)

#Remove duplicates of the papers with more than one species in a unique landscape
dups<- dbvar[c("title","study_site", "no_patches")]
dbvar<-dbvar[!duplicated(dups),]


# a. Landscape statistics--------------------------------------------------------------

##total extent in ha 
min(dbvar$totext_ha,na.rm=T)
max(dbvar$totext_ha,na.rm=T)
median(dbvar$totext_ha,na.rm=T)
summary(log10(dbvar$totext_ha))
sd(dbvar$totext_ha, na.rm=TRUE)
hist(dbvar$totext_ha)
hist(log10(dbvar$totext_ha),breaks=20)
#Total extent IQR
totext<-dbvar[!is.na(dbvar$totext_ha),c("totext_ha"),]
IQR(totext)
#total extent ranges from 0.1 to 4x10^8 ha (median 8180, IQR 63310)

##total patch area in ha 
options(digits=12)
summary(dbvar$totpat_ha,)
summary(log10(dbvar$totpat_ha))
sd(dbvar$totpat_ha, na.rm=TRUE)
hist(dbvar$totpat_ha)
hist(log10(dbvar$totpat_ha))
#Total patch area IQR
totpat<-dbvar[!is.na(dbvar$totpat_ha),c("totpat_ha"),]
IQR(totpat)
#total patch area ranges from 0.07 to 3.5x10^7 ha (median 162, IQR 426)

##number of patches
summary(dbvar$no_patches)
sd(dbvar$no_patches, na.rm=TRUE)
#Number of patches IQR
nopatch<-dbvar[!is.na(dbvar$no_patches),c("no_patches"),]
IQR(nopatch)
#number of patches ranges from 6 to 4000 (median 74, IQR 180.5 )


#Proportion of habitat 
options(digits=12)
summary(dbvar$frag)
summary(log10(dbvar$frag))
dbvar[which.max(dbvar$frag),] #ID 69
dbvar[which.min(dbvar$frag),] #ID 7
sd(dbvar$frag, na.rm=TRUE)
hist(dbvar$frag)
hist(log10(dbvar$frag))
#Fragmentation IQR
Frag<-dbvar[!is.na(dbvar$frag),c("frag"),]
IQR(Frag)

#proportion of habitat ranges from 0.000016 to 0.997 (median 0.025, IQR 0.0738)

# b. Temporal trends ----------------------------------------------------------------
summary(lm(log10(frag)~year,data=dbvar))
#Estimate= 0.059, p= 0.013*, proportion of habitat increases with time 
summary(lm(log10(totext_ha)~year,data=dbvar))
#NS no significant trend in total extent (w and w/o outliers)

#plot habitat coverage and extent values vs year
par(mfcol=c(2,2),mex=0.6)
with(dbvar,plot(frag~year,log="y",pch=4, ylab="fraction habitat"))
with(dbvar,plot(totext_ha~year,log="y",pch=4, ylab="study area extent, ha"))


# H5 ----------------------------------------------------------------------
## Study length varies with interests

##Study unit: study length by paper 

#Keep cases with different study lengths 
dbtime<-db[!duplicated(db$title,db$total_time_study,),]
dim(dbtime)

#Study length stats
summary(dbtime$total_time_study)
par(mfcol=c(1,1),mex=0.6)
hist(dbtime$total_time_study)

summary(glm.nb(total_time_study~conserv, dbtime))
#NS, length of conservation studies is not significantly longer than other studies 

# H6 ----------------------------------------------------------------------
## Model complexity varies depending on interest

#NS relationship between model complexity and study interest. Simple models have been applied with conservation and management aims

#Study unit: model(s) by paper

#Remove duplicates of species/landscapes
db_model<-db[!duplicated(db[, c("title", "model")]),]
dim(db_model)

levels(db_model$model)

db_model$model <- sub('[+]','modif',db_model$model)

db_model$mod1 = as.factor(grepl("IFM|Levins", db_model$model))
db_model$mod2 = as.factor(grepl("modif", db_model$model))
db_model$mod3 = as.factor(grepl("Author's|SRLM|PRM|MANAGE", db_model$model))

#Cases of model within papers
db_mod_1<- db_model[db_model["mod1"] == "TRUE",]
db_mod_1$mod_complex=1 #1: IFM/Levins
db_mod_2<- db_model[db_model["mod2"] == "TRUE",]
db_mod_2$mod_complex=2 #2: modif (IFM modified)
db_mod_3<- db_model[db_model["mod3"] == "TRUE",]
db_mod_3$mod_complex=3 #3: Author's/SRLM/PRM/MANAGE

db_mod_comp<- rbind(db_mod_1, db_mod_2, db_mod_3)

db_mod_comp$mod_complex<-as.factor(db_mod_comp$mod_complex)

c<-multinom(mod_complex~conserv, data =db_mod_comp)
summary(c)
cz <- summary(c)$coefficients/summary(c)$standard.errors
cp <- (1 - pnorm(abs(cz), 0, 1)) * 2
cp
#NS, the models used in conservation papers are not significantly more complex




  