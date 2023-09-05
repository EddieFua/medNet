rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/CHD_pos_apLCMS_624.bin')
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/neg_624.bin')
dat = read.delim('/Users/fuyinghao/Documents/Mediator analysis/Second data/CHDWB_first500_sampleid_mapping_pos.txt')
head(dat)
dim(g)
dim(dat)
dim(x)
###clean the colname of g and x
colnames(g) = gsub('.mzML','',colnames(g))
colnames(x) = gsub('.mzML','',colnames(x))

###there are some mismatch between two datasets
dat = dat[grep('CHD',dat$Sample.ID),]
dim(dat)
which(is.na(match(dat$File.Name,colnames(x))))
dim(x)
dat[c(441,442),]
dat = dat[-c(439:444),]
dim(dat)
which(is.na(match(colnames(x),dat$File.Name)))[-(1:4)]
x = x[,-which(is.na(match(colnames(x),dat$File.Name)))[-(1:4)]]
dim(x)
dim(dat)
length(unique(dat$Sample.ID))
###some sample were experimented for twice
positive_id = matrix(,ncol = 3)[-1,]
colnames(positive_id) = colnames(dat)
for (i in unique(substr(dat$Sample.ID,1,10))) {
  if(nrow(dat[substr(dat$Sample.ID,1,10) == i,])>3){
    positive_id = rbind(positive_id,dat[substr(dat$Sample.ID,1,10) == i,][1:3,])
  }else{
    positive_id = rbind(positive_id,dat[substr(dat$Sample.ID,1,10) == i,])
  }
}
dim(positive_id)
x = x[,-which(is.na(match(colnames(x),positive_id$File.Name)))[-(1:4)]]
dim(x)

File_Name_negative = paste0(substr(positive_id$File.Name,1,10),stringr::str_pad(as.numeric(substr(positive_id$File.Name,11,13))+1, 3, side = "left", "0"))
Sample_ID_nagative = paste0(substr(positive_id$Sample.ID,1,11),as.numeric(substr(positive_id$Sample.ID,12,12))+1)
negative_id = data.frame(File_Name_negative,Sample_ID_nagative,positive_id$Comment)
colnames(negative_id) = colnames(positive_id)
g = g[,-which(is.na(match(colnames(g),negative_id$File.Name)))[-(1:4)]]
dim(g)
###normalize data
normalize_data = function(data){
  data_normalize = log(data+1)
  return(data_normalize)
}
g_normalize = cbind(g[,1:4],normalize_data(g[,5:ncol(g)]))
x_normalize = cbind(x[,1:4],normalize_data(x[,5:ncol(x)]))
combine_triplet = function(dat){
  dat[dat==0]=NA
  mean_dat = rowMeans(dat,na.rm = T)
  return(mean_dat)
}
positive = cbind(x[,1:4],sapply(1:((ncol(x)-4)/3),function(i)combine_triplet(x_normalize[,(5+(i-1)*3):(4+(i)*3)]),simplify = T))
negative = cbind(g[,1:4],sapply(1:((ncol(g)-4)/3),function(i)combine_triplet(g_normalize[,(5+(i-1)*3):(4+(i)*3)]),simplify = T))
positive[which(is.nan(positive),arr.ind = T)]=0
negative[which(is.nan(negative),arr.ind = T)]=0
colnames(negative) = c(colnames(g)[1:4],unique(substr(negative_id$Sample.ID,1,10)))
colnames(positive) = c(colnames(g)[1:4],unique(substr(positive_id$Sample.ID,1,10)))
###remove batch effect
library(sva)
library(bladderbatch)
library(preprocessCore)
batch = positive_id$Comment[match(unique(substr(negative_id$Sample.ID,1,10)),substr(negative_id$Sample.ID,1,10))]
data = rbind(positive,negative)
dim(data)
edata = data[,5:dim(data)[2]]
combat_edata = ComBat(dat = edata, batch = batch) #remove batch effect
combat_edata <- normalize.quantiles(combat_edata)
data[,5:dim(data)[2]] = combat_edata
dim(data)

nrow(positive)
nrow(negative)
library(metapone)
pos<- data[1:8807,]
neg<- data[8808:nrow(data),]
match.tol.ppm=5
data(hmdbCompMZ)
pos.adductlist = c("M+H", "M+NH4", "M+Na", "M+ACN+H","M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H")
neg.adductlist = c("M-H","M-2H","M-2H+Na","M-2H+K", "M-2H+NH4","M-H2O-H","M-H+Cl", "M+Cl", "M+2Cl")
dat<-list(pos,neg)
type=c("pos","neg")

#####metapone annotation
concate<-function(a)
{
  b<-a[1]
  if(length(a)>1)
  {
    for(i in 2:length(a)) b<-paste(b, a[i], sep=",")
  }
  return(b)
}

if (length(dat)!=length(type)){
  message("The lengths of dat and type do not match")
  stop()
}

adductlist <- list()
for (i in 1:length(type)) {
  if (type[i]=='pos') {adductlist[[i]] <- pos.adductlist}
  else if (type[i]=='neg') {adductlist[[i]] <- neg.adductlist}
  else {message('Element in type should be pos or neg')
    stop()}
}

comp.mz<-new("list")
for(i in seq_len(length(adductlist))) comp.mz[[i]]<-hmdbCompMZ[hmdbCompMZ[,3] %in% adductlist[[i]],]

to.remove<-rep(FALSE, 2)
for(i in seq_len(length(dat)))
{
  if(is.null(dat[[i]])) to.remove[i]<-TRUE
}

if(sum(to.remove)>0)
{
  dat<-dat[-which(to.remove)]
  comp.mz<-comp.mz[-which(to.remove)]
  adductlist<-adductlist[-which(to.remove)]
}
nnn<-1

for(n in seq_len(length(dat)))
{
  this.dat<-dat[[n]]
  this.adductlist<-adductlist[[n]]
  this.comp.mz<-comp.mz[[n]]
  
  this.comp.mz[,1]<-as.character(this.comp.mz[,1])
  this.comp.mz[,3]<-as.character(this.comp.mz[,3])
  this.comp.mz[,2]<-as.numeric(this.comp.mz[,2])
  
  #this.matched<-foreach(nnn=1:length(this.adductlist), .combine=rbind) %dopar%
  
  
  match.fun<-function(ion, comp.mz, this.dat, match.tol.ppm)
  {
    mass.match<-function(x, known.mz, match.tol.ppm=5)
    {
      mass.matched.pos<-rep(0, length(x))
      for(i in seq_len(length(x)))
      {
        this.d<-abs((x[i]-known.mz)/x[i])
        if(sum(!is.na(this.d))>0)
        {					
          if(min(this.d,na.rm=TRUE) < match.tol.ppm/1e6) mass.matched.pos[i]<-1
        }
      }
      return(mass.matched.pos)
    }
    this.kegg<-comp.mz[comp.mz[,3] %in% ion,]
    m<-mass.match(this.dat[,1], this.kegg[,2], match.tol.ppm=match.tol.ppm)
    
    this.matched<-NULL
    if(sum(m)>0)
    {
      this.matched<-matrix(0, ncol=2+ncol(this.dat)+ncol(comp.mz), nrow=1)
      
      for(i in which(m==1))
      {
        d<-abs(this.dat[i,1]-this.kegg[,2])/this.dat[i,1]
        sel<-which(d <= (match.tol.ppm*1e-6))
        for(j in seq_len(length(sel))) this.matched<-rbind(this.matched, c(i, sel[j],this.dat[i,],this.kegg[sel[j],]))
      }
      
      this.matched<-this.matched[-1,]
      if(is.null(nrow(this.matched))) this.matched<-matrix(this.matched, ncol=2)
      this.matched<-this.matched[,seq(-1,-2)]
    }
    this.matched
  }
  z<-bplapply(this.adductlist, match.fun, comp.mz=this.comp.mz, this.dat=this.dat, match.tol.ppm=match.tol.ppm)
  this.matched<-z[[1]]
  if(length(z)>1)
  {
    for(m in 2:length(z))
    {
      this.matched<-rbind(this.matched, z[[m]])
    }
  }
  
  if(n == 1)
  {
    matched<-this.matched
  }else{
    matched<-rbind(matched, this.matched)
  }
}
dim(matched)
matched = as.data.frame(matched)
matched = t(apply(matched, 1, unlist))
makenumcols<-function(df){
  df<-as.data.frame(df)
  cond <- apply(df, 2, function(x) {
    x <- x[!is.na(x)]
    all(suppressWarnings(!is.na(as.numeric(x))))
  })
  numeric_cols <- names(df)[cond]
  df[,numeric_cols] <- apply(df[,numeric_cols],2, as.character) # deals with factors
  df[,numeric_cols] <- sapply(df[,numeric_cols], as.numeric)
  return(df)
}
matched = makenumcols(matched)

# save(matched,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/matched.RData')
###delete element of matched which does not have C name
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/matched.RData')

load('/Users/fuyinghao/Documents/Mediator analysis/Second data/HMDB4.0_database.bin')
metabolite.table[1:5,1:5]
metabolite.table$KEGG_compound_ID[metabolite.table$KEGG_compound_ID==''] = NA
metabolite.table$KEGG_compound_ID = gsub("\t", "", metabolite.table$KEGG_compound_ID)
temp_id = metabolite.table$KEGG_compound_ID[match(matched$HMDB_ID,metabolite.table$HMDB_ID)]
matched = matched[!is.na(temp_id),]
matched = as.matrix(matched)
rownames(matched) = na.omit(temp_id)
dim(matched)
matched = matched[match(unique(rownames(matched)),rownames(matched)),]
dim(matched) ####matched 3317 426
matched[1:5,1:5]
matched = makenumcols(matched)
###intake
library(readxl)
exposure_dat = read_xlsx('/Users/fuyinghao/Documents/Mediator analysis/Second data/URC_Aim1a_NWO_FFQ_demographics092016_grams.xlsx',sheet = 'URC_alldata_V1')
dim(exposure_dat)
exposure_dat$CHD_ID_metabolomics = toupper(exposure_dat$CHD_ID_metabolomics)
idx = match(exposure_dat$CHD_ID_metabolomics,colnames(matched))
which(is.na(idx))
idx = na.omit(idx)
exposure_dat = exposure_dat[-86,]
matched = matched[,idx]

exposure_variable = c('Gender','Age','Race','ModerateActivity','EDUCATION','INCOME','Hx_Dz','TBALC_TOBUSE','RegionPfat_TotalBody') 
outcome_variable = 'BMI'
exposure_variable_a = read_xlsx('/Users/fuyinghao/Documents/Mediator analysis/Second data/URC_Aim1a_NWO_FFQ_demographics092016_grams.xlsx',sheet = 'FFQ_DataDictionary')
exposure_variable_a = exposure_variable_a[1:87,]
exposure_variable = c(exposure_variable,exposure_variable_a$`Report Column Name`)
exposure_variable = toupper(exposure_variable)
colnames(exposure_dat) = toupper(colnames(exposure_dat))
outcom_dat = exposure_dat$BMI
temp_id = match(exposure_variable,colnames(exposure_dat))
exposure_variable[which(is.na(temp_id))]
exposure_variable[is.na(temp_id)] = c('FF_GROUP_VEGETBLFRUITFIBER_TOT','FF_GROUP_SUGARYBEVG_TOT_GRAMS','FF_GROUP_SUGARYBEVG_TOTAL_KCAL')
temp_id = match(exposure_variable,colnames(exposure_dat))
exposure_variable[which(is.na(temp_id))]
exposure_dat = exposure_dat[,temp_id]


######final dat, 179 samples, 3317 mediators, 96 exposures
exposure = exposure_dat
dim(exposure)
outcome = outcom_dat
mediator = t(matched)
dim(mediator)
save(exposure,outcome,mediator,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/without_network.RData')
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/without_network.RData')

#####network
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/human.graph.degree.less.than.20.bin')
library(igraph)
g <- as.undirected(g) # change to undirected
g_name <- V(g)$name
r_name <- g_name[substring(g_name,1,1)=='R']

add_edge_func <- function(g, name){
  g_name <- V(g)$name
  idx <- which(g_name==name)
  nb <- unique(neighbors(g, idx)$name)
  if (length(nb)>1){
    for (i in 1:(length(nb)-1)){
      for (j in (i+1):length(nb)){
        g <- add_edges(g, c(which(V(g)$name==nb[i]), which(V(g)$name==nb[j])), color = 'red')
      }
    }
  }
  g <- delete.vertices(g, which(V(g)$name==name))
  print(length(V(g)$name))
  return(g)
}

for (r in r_name){g <- add_edge_func(g, r)}


all_compound <- V(g)$name
sel_compound <- intersect(all_compound, colnames(mediator))
length(sel_compound) #1035

idx_g <- c()
for (i in sel_compound) {
  idx_g <- c(idx_g, which(all_compound==i))
}

adj_new <- as_adjacency_matrix(g)[idx_g,][,idx_g]

dim(adj_new) # 1035 * 1035

mediator <- mediator[,which(colnames(mediator) %in% sel_compound)]
dim(mediator) #179*1035
# adj_new = as.matrix(adj_new)
tri_mediator <- Matrix::triu(adj_new)
network <- Matrix::which(tri_mediator!=0,arr.ind = T)
mediator_network = network
mediator_network[,1] = colnames(adj_new)[network[,1]]
mediator_network[,2] = colnames(adj_new)[network[,2]]
mediator_graph = graph_from_edgelist(mediator_network)
vertices_d = degree(mediator_graph)[which(degree(mediator_graph)>10&degree(mediator_graph)<20)]

library(KEGGREST)
listDatabases()
b = matrix(NA,nrow = length(vertices_d),ncol = 2)
for (i in 1:length(vertices_d)) {
  b[i,1] = names(keggFind("compound", names(vertices_d)[i]))
  b[i,2] = keggFind("compound", names(vertices_d)[i])
}
b = cbind(b,as.numeric(vertices_d))
colnames(b) = c('Entry','Compound','Degree')
write.csv(b,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/10<vertice<20.csv')

mediator =  mediator[,which(colnames(mediator) %in% V(mediator_graph)$name)]
exposure = exposure[,-6]
exposure = as.matrix(exposure)
exposure = cbind(GENDER = as.factor(exposure[,1]),apply(exposure[,2:ncol(exposure)],2,as.numeric))
response = ifelse(outcome>=25,1,0)
mediator_network = unique(mediator_network)
mediator_network[,1] = match(mediator_network[,1],colnames(mediator))
mediator_network[,2] = match(mediator_network[,2],colnames(mediator))
mediator_network = makenumcols(mediator_network)
save(mediator,response,exposure,mediator_network,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')

###856 mediators,179 samples, 95 exposures, 918 edges
dim(mediator)
dim(network)
dim(exposure)



