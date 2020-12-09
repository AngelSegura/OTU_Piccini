## OTU ANALISIS 
# Gradient from Salto Grande to Punta del Este- Uruguay South America
# Piccini et al 
#Corresponding author:  Claudia Piccini <claudia.piccini@gmail.com> ; Angel Segura <asegura@cure.edu.uy>

R
## Required data files 
# Conectivity_matrix.csv
# otu_table_from_biom.txt
#output_FASTA_Hamming.csv

## Required packages
require(vegan)
require(sads)
require(meteR)
require(fBasics)
require(ecolottery)
require(e1071)

# From Supporting material del paper de Baselga https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1466-8238.2009.00490.x
source("DATA/geb_490_sm_appendix_s1.r")
source("DATA/geb_490_sm_appendix_s2.r")


# Pseudo R squared for SAD curves
# R2m according to Nature Ecology-Evolution methods pg 6
r2m<-function(SpeciesAbund, modelSAD){
obs<-rad(SpeciesAbund)[,2]
pred<-radpred(modelSAD)[,2]

denom<-sum((log(obs)-log(pred))^2) 
numer<-sum((log(obs)-mean(log(obs)))^2) 
r2m<- 1- (denom/numer)
r2m
}

### Load data
data<-read.table("DATA/otu_table_from_biom.txt", header=TRUE, sep=",", dec=".")
summary(data)
colnames(data)

##  DISTANCE PAIRS AMONG SITES
dist_mat<-read.table("DATA/Conectivity_matrix.csv", header=TRUE, sep=",", dec=".",na.strings = "")

Dist_km<-t(dist_mat)[lower.tri(t(dist_mat))]
dist_mat_1<- as.dist(t(dist_mat))

# DISTANCE FROM HEADWATERS TO THE OUTLET
dist_fromOutlet<-as.numeric(max(dist_mat[1,])-dist_mat[1,])


# Data manipulation and exploration
## ABUNDANCE
# Extraigo la data
data_abund1<-data[,c(12,8,2,6,7,9,5,3)] # SAco las BC  que no sirven % OJO QUE ESTOY DEJANDO LAS DE 2014
dim(data_abund1)
colnames(data_abund1)

length(which(rowSums(data_abund1)<=1)) # numero de singletons

selSingle <-which(rowSums(data_abund1)<=1)
data_abund<-data_abund1[-selSingle,1:6] # Saco los singletons o los que no habia aqui

data_names<-data[-selSingle,13:19]

data_ID<-data.frame(ID=data[-selSingle,1], data_abund)

## Data 2014
data_abund14<-data_abund1[-selSingle,c(7,8)] # Remove unuseful BC

all_data<-data.frame(data_abund, data_abund14)

## Descriptive analysis
Ntot<-colSums(all_data) # Ntotal
sum(Ntot)

# Richness
S<-colSums(all_data!=0)
S
Ssingle<-colSums(data_abund1!=0) # Exclude singletons

plot(S, Ssingle, xlim=c(0,10000), ylim=c(0, 10000))
abline(a=1000, b=1) # 

specnumber(t(all_data)) # observed number of species
raremax <- min(rowSums(t(all_data)))
Srare <- rarefy(t(all_data), raremax)
#jpeg("Figures/RarefactionCurveSpecies.jpeg", height=5, width=5, res=300, units="in")
rarecurve(t(all_data), step = 20, sample = raremax, col = "blue", cex = 0.6, ylab="OTU richness")
dev.off()

# Equiability (Simpson Index)
SimpsonIndex<-diversity(all_data, index="simpson", MARGIN=2)
SimpsonExp<- 0.79*Ntot^-0.24
log10(SimpsonExp)
log10(SimpsonIndex)

# Average Total OTU abundance

mean(Ntot)
sd(Ntot)

mean(S)
sd(S)

mean(SimpsonIndex)
sd(SimpsonIndex)

## FIGURE 1

#jpeg("Figures/FIGURE_1_TotalAbundance_andS_distance_EvenessSrarefied.jpeg", height=5, width=3.5, res=300, units="in")
x11(height=5,width=3.5)
par(mfrow=c(2,1), mar=c(2,4,1,1))
options(scipen=0)
plot(dist_fromOutlet, Ntot[1:6], axes=FALSE,ylab="", xlim=rev(range(dist_fromOutlet)),log="y", type="b", pch=19, ylim=c(2e4,6e5), xlab="", cex.lab=0.9);axis(2, cex.axis=0.9);box(); axis(1, at=c(0,200,400,600))
mtext("Total OTU abundance", 2, line=2,cex=0.9)
points(dist_fromOutlet[c(1,6)], Ntot[7:8], col="gray", pch=19)
legend("topleft", legend="A", bty="n")
legend("bottomleft", legend=c("2013", "2014"), pch=19, col=c("black", "grey"), bty="n",cex=0.9)

par(mar=c(4,4,0,1)); 
options(scipen=999)
plot(dist_fromOutlet, Srare[1:6],xlim=rev(range(dist_fromOutlet)),xlab="",type="b", axes=TRUE,ylab="", ylim=c(1e2,8e3),cex.lab=0.9, pch=19);
mtext("OTU Richness", 2, line=2,cex=0.9)
points(dist_fromOutlet[c(1,6)], S[7:8], col="grey", pch=19)
legend("topleft", legend="B", bty="n")
text(770,0, "Reservoir", pos=4, cex=0.8)
text(580,0, "Riverine", pos=4, cex=0.8)
text(200,0, "Estuarine", pos=4, cex=0.8)
mtext("Distance from outlet (km)", 1, line=2.5)
dev.off()


################################################################################################
## Ecolottery- Coalescent model theoretical predictions
################################################################################################

# Function to construct SAD from ecolotery object 
SAD_coales<-function(coales_out, local=TRUE){ # either local or pool community
if(local==TRUE) out_data<-coales_out$com$sp else out_data<-coales_out$pool$sp
y<-sort(table(out_data),decreasing=TRUE)
x<-1:length(y)
cbind(x,y)
}

#Data matrix in Coalesc format and vice versa
data2coalesc<-function(mat){ # first row names, second row abundances
sp<-rep(mat[,1], times=mat[,2])
ind<-1:length(sp)
data.frame(ind, sp)
}

coalesc2data<-function(coalesc_output){
spp_pool<-names(table(coalesc_output$pool$sp))
abund_pool<-table(coalesc_output$pool$sp)
#trait_pool<- unique(coalesc_output$pool$sp) # debo seleccionar una spp de cada y poner su trait
spp_com<-names(table(coalesc_output$com$sp))
abund_com<-table(coalesc_output$com$sp)
list(com=data.frame(spp_com, abund_com=as.vector(abund_com)), pool=data.frame(spp_pool,pool_abund=as.vector(abund_pool)))
}


## Function for Gaussian environmental filtering
filt_gaussian <- function(t, x) exp(-(x - t)^2/(2*sigma^2))

# Define simulation parameters
J=2e5 # FIXED Similar to mean(Ntot)
m=0.05 # FIXED- Intensity of migration from REGIONAL Community ## THEN SHOULD BE VARIABLE Tucker et al 2015


# Observed communities
data_in<-data.frame(paste("OTU", 1:nrow(data_abund[data_abund[,1]>0,]), sep=""), data_abund[data_abund[,1]>0,1])
a<-data.frame(name=paste("a",c(1:30)),abund=matrix(round(runif(30,0,10))), nrow=30)
Names<- paste("OTU",1:nrow(data_abund),sep="")
OBS_comm1<-data2coalesc(data.frame(Names,rowSums(data_abund))) # Community assuming all sites
OBS_comm<- data.frame(OBS_comm1, tra1=seq(0,1,length=nrow(OBS_comm1)))
OBS_comm[1:100,] # First hundred cases from regional community 

NroOTUs<-length(rowSums(data_abund))
AbundOTUs<-sum(data_abund)
Abund_promOTUs<-round(AbundOTUs/NroOTUs)

# Function for richness and biodiversiy index
f.sumstats <- function(com) array(dimnames=list(c("S", "Es")),
                                       c(vegan::specnumber(table(com[,2])),
                                         vegan::diversity(table(com[,2]))))

# Generate list to save coalescent output
overallSIM<-list()

MU_VECT<-seq(0.2,0.9, length=6) # Equal position in the niche axis

y1<- 1+MU_VECT*0 # Homogeneous sigma across mu
y2<-0.1+MU_VECT*0 # Homogeneous sigma, low sigma 
y3<- 0.25- 0.2*MU_VECT # Sigma changes with mu, environmental gradient

SIGMA_VEC<-list(y1,y2,y3)

j=2
for(j in 1:3){  

## Changes environmental filter in each simulation
mu_vect<-MU_VECT # Siempre fijo
sigma_vect<-SIGMA_VEC[[j]]

SIM_COM<-list()
for(site in 1:6){
mu<-mu_vect[site] # Mu is a function of site
sigma<-sigma_vect[site] # the niche width varies with position towards the outer end
comm_aux <- coalesc(J, m, filt = function(x) filt_gaussian(mu, x),pool=OBS_comm) #esto funciona con un pool similar al de toda la composicion del Gradiente
SimComm_aux<-coalesc2data(comm_aux)$com
SIM_COM[[site]]<-SimComm_aux
print(site)
}

overallSIM[[j]]<-list(SIM_COM,OBS_comm)
print(paste("we are in simulation",j))
}

## Function to analize ecolotery results and estimate biodiversity metrics (beta, richness, etc)

EcoloteryBiodivIndex<-function(SIM_LIST, A){
# SIM_LIST- a List with output from ecoloterry transformed to a list of OTUS
# A number of simulation to perform calculations

DF_j_1<-Reduce(function(x, y) merge(x, y, all = TRUE, by="spp_com"), SIM_LIST[[A]][[1]], accumulate=FALSE)
DF_j_1<-replace(DF_j_1, is.na(DF_j_1), 0)
colnames(DF_j_1)<-c("spp_com", "C1","C2","C3","C4", "C5","C6")

# Extract average trait values from OTUs
OTUtra1<-aggregate(SIM_LIST[[A]][[2]]$tra1, list(SIM_LIST[[A]][[2]]$sp), mean)
colnames(OTUtra1)<- c("spp_com","trait")

# Merge both matrix
OTU_1_all<-merge(DF_j_1, OTUtra1, by="spp_com") # matriz con OTUS y valores de abundancia y valores de los traits

# Separate vector of traits
Traits<- OTU_1_all$trait

# Matrix of Sites X Otus
OTU_SITE<- OTU_1_all[,2:7]


## Species Richness per site
S<- colSums(OTU_SITE>0)

## Eveness
SimpsonIndexSim<-diversity(OTU_SITE, index="simpson", MARGIN=2)

## Beta diversity
betaSor_SIM<-as.matrix(as.dist(beta.sor(t(OTU_SITE))))
betaSim_SIM<-as.matrix(as.dist(beta.sim(t(OTU_SITE))))
betaNes_SIM<-as.matrix(as.dist(beta.nes(t(OTU_SITE))))

# Distancia entre los mu (Optimos de los traits). Asume un mach perfecto entre el trait filter optima y la posicion espacial
DistSim<-as.matrix(vegdist(data.frame(y=seq(0.2,0.9, length=6),x=rep(0,6)), method="euclidean"))
diag(DistSim)<-NaN # Excluyo la comparacion del sitio consigo mismo.

# VALORES DE LA CORRELACION ENTRE DISTANCIA y COMPONENTES de BETA
P_sor<-cor.test(as.numeric(betaSor_SIM), as.numeric(DistSim), na.rm=TRUE)$p.value
rho_sor<-cor.test(as.numeric(betaSor_SIM), as.numeric(DistSim), na.rm=TRUE)$estimate

P_sim<-cor.test(as.numeric(betaSim_SIM), as.numeric(DistSim), na.rm=TRUE)$p.value
rho_sim<-cor.test(as.numeric(betaSim_SIM), as.numeric(DistSim), na.rm=TRUE)$estimate

P_nes<-cor.test(as.numeric(betaNes_SIM), as.numeric(DistSim), na.rm=TRUE)$p.value
rho_nes<-cor.test(as.numeric(betaNes_SIM), as.numeric(DistSim), na.rm=TRUE)$estimate

BetaDistmat<-matrix(c(P_sor,P_sim,P_nes, rho_sor,rho_sim,rho_nes), nrow=3)
colnames(BetaDistmat)<- c("Pvalue", "rho Estimate")
rownames(BetaDistmat)<- c("SOR", "SIM","NES")

BetaINFO<-list(DISTANCE=DistSim, SOR=betaSor_SIM,SIM=betaSim_SIM,NES=betaNes_SIM)

## SADs
# Lognormal SAD
R2ln<-unlist(lapply(OTU_SITE, function(X){ XX<-X[X>0]; fln<-fitlnorm(XX); r2m(XX, fln)})) # Ajusto una SAD, y calculo el R2 estimado usando la funcion r2m
R2zm<-unlist(lapply(OTU_SITE, function(X){ XX<-X[X>0]; fln<-fitmand(XX); r2m(XX, fln)})) # Ajusto una SAD, y calculo el R2 estimado usando la funcion r2m
R2ls<-unlist(lapply(OTU_SITE, function(X){ XX<-X[X>0]; fln<-fitls(XX); r2m(XX, fln)})) # Log series METE and UNTB when sample is large
R2mz<-unlist(lapply(OTU_SITE, function(X){ XX<-X[X>0]; fln<-fitmzsm(XX); r2m(XX, fln)})) # The mZSM describes the SAD of a sample taken from a neutral metacommunity under random drift.
R2mat<-rbind(R2ln,R2zm,R2ls,R2mz)


return(list(S=S, SimpsonIndex=SimpsonIndexSim, BetaDist=BetaDistmat, BetaInfo=BetaINFO, SAD_R2mat=R2mat))

}

A1<-EcoloteryBiodivIndex(SIM_LIST=overallSIM, A=1)
A1
A2<-EcoloteryBiodivIndex(SIM_LIST=overallSIM, A=2)
A2
A3<-EcoloteryBiodivIndex(SIM_LIST=overallSIM, A=3)
A3

## FIGURE 2
    x11()
    par(mfcol=c(5,3), mar=c(5,4,1,1))
    A=1
    plot(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT,SIGMA_VEC[[1]]), type="n", xlim=c(0,1),ylim=c(0,0.7), ylab="", xlab="")
    for(i in 1:length(MU_VECT)) lines(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT[i],SIGMA_VEC[[A]][i]),col=i)
    legend("topleft",legend="No selection", bty="n")
    mtext("Environmental" ,2,line=2.7, cex=0.7)
    mtext("filter" ,2,line=1.8, cex=0.7)

    plot(seq(0.2,0.9,length=6), A1$S, ylim=c(1500,4000), xlab="",ylab="",xlim=c(0,1), pch=19)
    mtext("Richness (S)" ,2,line=2, cex=0.7)

    plot(A1$BetaInfo$DISTANCE, A1$BetaInfo$SOR, ylim=c(0,1), pch=19, xlab="", ylab="")
    mtext(expression(beta[sor]), 2, line=2, cex=1)

    plot(A1$BetaInfo$DISTANCE, A1$BetaInfo$SIM, ylim=c(0,1), pch=19, xlab="", ylab="")
    mtext(expression(beta[sim]), 2, line=2, cex=1)

    plot(A1$BetaInfo$DISTANCE, A1$BetaInfo$NES, ylim=c(0,1), pch=19, xlab="", ylab="")
    mtext(expression(beta[nes]), 2, line=2, cex=1)

    A=2
    plot(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT,SIGMA_VEC[[2]]), type="n", xlim=c(0,1),ylim=c(0,7), ylab="", xlab="")
    for(i in 1:length(MU_VECT)) lines(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT[i],SIGMA_VEC[[A]][i]),col=i)
    legend("topleft",legend="Homogeneous selection", bty="n")

    plot(seq(0.2,0.9,length=6), A2$S, ylim=c(1500,3200), xlab="",ylab="",xlim=c(0,1), pch=19)

    plot(A2$BetaInfo$DISTANCE, A2$BetaInfo$SOR, ylim=c(0,1), pch=19, xlab="",ylab="")
    plot(A2$BetaInfo$DISTANCE, A2$BetaInfo$SIM, ylim=c(0,1), pch=19, xlab="",ylab="")
    plot(A2$BetaInfo$DISTANCE, A2$BetaInfo$NES, ylim=c(0,1), pch=19, xlab="",ylab="")

    A=3
    plot(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT,SIGMA_VEC[[3]]), type="n", xlim=c(0,1),ylim=c(0,7), ylab="", xlab="")
    for(i in 1:length(MU_VECT)) lines(seq(0,1, length=100), dnorm(seq(0,1, length=100), MU_VECT[i],SIGMA_VEC[[A]][i]),col=i)
    legend("topleft",legend="Heterogeneous selection", bty="n")

    plot(seq(0.2,0.9,length=6), A3$S, ylim=c(1500,3200),xlim=c(0,1), pch=19, ylab="", xlab="")

    plot(A3$BetaInfo$DISTANCE, A3$BetaInfo$SOR, ylim=c(0,1), pch=19, xlab="",ylab="")

    plot(A3$BetaInfo$DISTANCE, A3$BetaInfo$SIM, ylim=c(0,1), pch=19, xlab="", ylab="")
    plot(A3$BetaInfo$DISTANCE, A3$BetaInfo$NES, ylim=c(0,1), pch=19, xlab="", ylab="")
    mtext("Relative distance", 1, line=2.5, at=-0.5)
      
dev.off()

## Richness
x11()
par(mfrow=c(3,1))
plot(seq(0.2,0.9,length=6), A1$S, ylim=c(1000,3200))
plot(seq(0.2,0.9,length=6), A2$S, ylim=c(1000,3200))
plot(seq(0.2,0.9,length=6), A3$S, ylim=c(1000,3200))

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

### SPECIES ABUNDANCE DISTRIBUTION SADS
## Utilzare el paquete sads (que incluye todos los modelos de vegan, y meteR)
SA<-data_abund[data_abund[,1]>0,1]
FB<-data_abund[data_abund[,2]>0,2]
CA<-data_abund[data_abund[,3]>0,3]
CO<-data_abund[data_abund[,4]>0,4]
MO<-data_abund[data_abund[,5]>0,5]
PE<-data_abund[data_abund[,6]>0,6]
SA14<-all_data[all_data[,7]>0,7]
PE14<-all_data[all_data[,8]>0,8]

ALL_SAD<-list(SA, FB,CA,CO,MO,PE)

## GAMMA DIVERSITY SAD
nrow(data_abund)

GAMMA_SAD<-rowSums(data_abund)[rowSums(data_abund)>0]
length(GAMMA_SAD)

lnSAD_gamma<-fitlnorm(GAMMA_SAD) # lognormal
gmSAD_gamma<-fitgamma(GAMMA_SAD)

zoSAD_gamma<-fitzipf(GAMMA_SAD)  # Solo zipf
zmSAD_gamma<-fitmand(GAMMA_SAD)  # Zipf Mandelbrot- Modelo generalizado

gsSAD_gamma<-fitgs(GAMMA_SAD)    # Motomura's Geometric Series Preemption
bsSAD_gamma<-fitbs(GAMMA_SAD)    # Broken-stick distribution (MacArthur 1960)

lsSAD_gamma<-fitls(GAMMA_SAD)    # log series (METE) UNTB when sample is large
mzSAD_gamma<-fitmzsm(GAMMA_SAD)  # metacommunity Zero-sum multinomial distribution. The mZSM describes the SAD of a sample taken from a neutral metacommunity under random drift.

R2_gammaSAD<-c(
r2m(GAMMA_SAD, lnSAD_gamma),
r2m(GAMMA_SAD, gmSAD_gamma),
r2m(GAMMA_SAD, zoSAD_gamma),
r2m(GAMMA_SAD, zmSAD_gamma),
r2m(GAMMA_SAD, gsSAD_gamma),
r2m(GAMMA_SAD, bsSAD_gamma),
r2m(GAMMA_SAD, lsSAD_gamma),
r2m(GAMMA_SAD, mzSAD_gamma)
)
round(R2_gammaSAD,2)


## Rarity patters 
length(which(100*cumsum(sort(GAMMA_SAD,decreasing=TRUE))/sum(GAMMA_SAD) < 50)) # cuantas especies acumulan el 50% de la abundancia

length(which(100*cumsum(sort(GAMMA_SAD,decreasing=TRUE))/sum(GAMMA_SAD)<80)) # cuantas especies acumulan el 80% de la abundancia

length(GAMMA_SAD) # Total number of OTUS

## LOCAL SAD
#### SAD indivudual

x11(height=4, width=12)
par(mfrow=c(1,6), mar=c(5,2,1,1))
plot(sort(SA, decreasing=TRUE), log="y", xlim=c(1,7000)) # equivalent to plot(rad(SA))
#points(sort(SA14, decreasing=TRUE),col="grey")
plot(sort(FB, decreasing=TRUE), log="y", xlim=c(1,7000))
plot(sort(CA, decreasing=TRUE), log="y", xlim=c(1,7000))
plot(sort(CO, decreasing=TRUE), log="y", xlim=c(1,7000))
plot(sort(MO, decreasing=TRUE), log="y", xlim=c(1,7000))
plot(sort(PE, decreasing=TRUE), log="y", xlim=c(1,7000))
#points(sort(PE14, decreasing=TRUE),col="grey")

## AJUSTO MODELOS CON Sads # SEE BELOW SAVED OBJECT
SADmodels<-list()
for(i in 1:6){
SAD_aux<-ALL_SAD[[i]]

lnSAD<-fitlnorm(SAD_aux) # lognormal
gmSAD<-fitgamma(SAD_aux) # Geometric mean

zoSAD<-fitzipf(SAD_aux)  # Zipf
zmSAD<-fitmand(SAD_aux)  # Zipf Mandelbrot-

gsSAD<-fitgs(SAD_aux)    # Motomura's Geometric Series Preemption
bsSAD<-fitbs(SAD_aux)    # Broken-stick distribution (MacArthur 1960)

lsSAD<-fitls(SAD_aux)    # log series (METE) UNTB when sample is large
mzSAD<-fitmzsm(SAD_aux)  # metacommunity Zero-sum multinomial distribution. The mZSM describes the SAD of a sample taken from a neutral metacommunity under random drift.

SADmodels[[i]]<-list(ln=lnSAD,gm=gmSAD,zo=zoSAD,zm=zmSAD,gs=gsSAD,bs=bsSAD,ls=lsSAD,mz=mzSAD)
}

## SADs 2014
SADmodels14<-list()
for(i in 1:2){
SAD_aux<-list(SA14,PE14)[[i]]

lnSAD<-fitlnorm(SAD_aux) # lognormal
gmSAD<-fitgamma(SAD_aux)

zoSAD<-fitzipf(SAD_aux)  # Solo zipf
zmSAD<-fitmand(SAD_aux)  # Zipf Mandelbrot- Modelo generalizado

gsSAD<-fitgs(SAD_aux)    # Motomura's Geometric Series Preemption
bsSAD<-fitbs(SAD_aux)    # Broken-stick distribution (MacArthur 1960)

lsSAD<-fitls(SAD_aux)    # log series (METE) UNTB when sample is large
mzSAD<-fitmzsm(SAD_aux)  # metacommunity Zero-sum multinomial distribution. The mZSM describes the SAD of a sample taken from a neutral metacommunity under random drift.

SADmodels14[[i]]<-list(ln=lnSAD,gm=gmSAD,zo=zoSAD,zm=zmSAD,gs=gsSAD,bs=bsSAD,ls=lsSAD,mz=mzSAD)
}


## Adjusted R squared SADs
## CALCULO r2m para cada modelo y sitio
R2mat<- matrix(NaN, ncol=8, nrow=8) # Columnas SITIOS
rownames(R2mat)<-rownames(modelAIC)
colnames(R2mat)<-c(colnames(modelAIC),colnames(modelAIC14))
for(i in 1:8){
r2_SA<-r2m(ALL_SAD[[1]], SADmodels[[1]][[i]])
r2_FB<-r2m(ALL_SAD[[2]], SADmodels[[2]][[i]])
r2_CA<-r2m(ALL_SAD[[3]], SADmodels[[3]][[i]])
r2_CO<-r2m(ALL_SAD[[4]], SADmodels[[4]][[i]])
r2_MO<-r2m(ALL_SAD[[5]], SADmodels[[5]][[i]])
r2_PE<-r2m(ALL_SAD[[6]], SADmodels[[6]][[i]])
r2_SA14<-r2m(SA14,       SADmodels14[[1]][[i]]) # del 2014
r2_PE14<-r2m(PE14,       SADmodels14[[2]][[i]]) # del 2014
R2mat[i,]<-c(r2_SA,r2_FB,r2_CA,r2_CO,r2_MO,r2_PE, r2_SA14, r2_PE14)
}

round(R2mat,2)

# ZM
mean(R2mat[4,])
sd(R2mat[4,])
# LN
mean(R2mat[1,])
sd(R2mat[1,])

## FIGURE 4 SADs
x11(height=3.5, width=12)
par(mfrow=c(1,6), mar=c(4.5,4,1,0))
plot(sort(SA, decreasing=TRUE), log="y", xlim=c(1,7000),col="grey" ,xlab="",ylab="") # equivalent to plot(rad(SA))
mtext("OTU abundance",2, line=2)
legend("topright", legend="A", bty="n")
legend(3000, 9000, legend=c("2013","2014"), lty=c(1,2),bty="n")
lines(radpred(SADmodels[[1]]$zm), col="black")
lines(radpred(SADmodels14[[1]]$zm),col="black",lty=2)
par(mar=c(4.5,3,1,1))
plot(sort(FB, decreasing=TRUE), log="y", xlim=c(1,7000), col="grey",xlab="",ylab="")
lines(radpred(SADmodels[[2]]$zm), col="black")
legend("topright", legend="B", bty="n")
plot(sort(CA, decreasing=TRUE), log="y", xlim=c(1,7000), col="grey",xlab="",ylab="")
lines(radpred(SADmodels[[3]]$zm), col="black")
legend("topright", legend="C", bty="n")
plot(sort(CO, decreasing=TRUE), log="y", xlim=c(1,7000), col="grey",xlab="",ylab="")
lines(radpred(SADmodels[[4]]$zm), col="black")
legend("topright", legend="D", bty="n")
plot(sort(MO, decreasing=TRUE), log="y", xlim=c(1,7000), col="grey",xlab="",ylab="")
lines(radpred(SADmodels[[5]]$zm), col="black")
legend("topright", legend="E", bty="n")
plot(sort(PE, decreasing=TRUE), log="y", xlim=c(1,7000), col="grey",xlab="",ylab="")
lines(radpred(SADmodels[[6]]$zm), col="black")
legend("topright", legend="F", bty="n")
lines(radpred(SADmodels14[[2]]$zm),col="black",lty=2)
mtext("Species rank", 1, line=2.5, at=-2e4)
dev.off()


# Parameters from the ZM SAD
zm_pars<-matrix(unlist(lapply(SADmodels, function(y) coef(y$zm) )),nrow=6, byrow=TRUE)
colnames(zm_pars)<-c("N","s","v")

zm_pars14<-matrix(unlist(lapply(SADmodels14, function(y) coef(y$zm) )),nrow=2, byrow=TRUE)
colnames(zm_pars14)<-c("N","s","v")

ZM_pars<-round(rbind(zm_pars,zm_pars14),2)
summary( ZM_pars)
colMeans(ZM_pars)
apply(ZM_pars, 2, sd)


## BETA DIVERSITY PATTERNS
########################################################################################################################################################
## TUCKER BETA DIVERSITY NULL TEST
########################################################################################################################################################
########################################################################################################################################################
## This takes a long time to run
{
DATA_BETA<-t(data_abund) # species x site matrix

## Load required source files and libraries
library(vegan)
source("/home/user/Escritorio/TIOJORGE/PaperSyCosasINprogress/Piccini_OTUs/Tucker/MetacommunityDynamicsFctsOikos.R")
source("/home/user/Escritorio/TIOJORGE/PaperSyCosasINprogress/Piccini_OTUs/Tucker/PANullDevFctsOikos.R")


dune<-DATA_BETA# PResence absence
patches<-nrow(dune)
################

##Calculate beta-diversity for metacommunity
dune_bin <- siteXsp_prep(dune, plot_names_in_col1=FALSE)

beta_comm <- vegdist(dune_bin, "jaccard")  	# Calculate beta-diversity

res_beta_comm <- as.matrix(as.dist(beta_comm))
diag(res_beta_comm) <- NA

beta_div_stoch <- apply(res_beta_comm, 2, FUN="mean", na.rm=TRUE) 	# Calculate patch/site mean value

## Calculate PA null deviation
nullobj <- null_distributions(dune_bin, test_func = vegdist, method = "jaccard", reps=999)		# generate null distributions
nulldev <- null_deviations(dune_bin, nullobj, vegdist, method="jaccard")		# Calculate null deviations using the expected value method

### Store results, where each value is the PA beta-null deviation value for a pairwise comparison (distance matrix) between patches
res_nulldev <- as.matrix(as.dist(nulldev$null_devs_eval))
diag(res_nulldev) <- NA
PA_null_dev <- apply(res_nulldev, 2, FUN="mean", na.rm=TRUE) #Calculate patch mean value
PA_null_devSD <- apply(res_nulldev, 2, FUN="sd", na.rm=TRUE) #Calculate patch sd value

plot(PA_null_dev , ylim=c(-0.2,0.2)): abline(h=0, col="gray")
segments(1:6, PA_null_dev-PA_null_devSD, 1:6, PA_null_dev+PA_null_devSD)


### NULL DEVIATIONS FOR ABUNDANCE

### Prepare and calculate abundance beta-null deviation metric
## Adjusted from Stegen et al 2012 GEB
bbs.sp.site <- dune
rand <- 999
null.alphas <- matrix(NA, ncol(dune), rand)
null.alpha <- matrix(NA, ncol(dune), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
gamma <- ncol(bbs.sp.site) #gamma
obs_beta <- 1-mean.alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

##Generate null patches
for (randomize in 1:rand) {  
	null.dist = dune
	for (species in 1:ncol(null.dist)) {
		tot.abund = sum(null.dist[,species])
		null.dist[,species] = 0
		for (individual in 1:tot.abund) {
			sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
			null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
		}
	}
	
	##Calculate null deviation for null patches and store
	null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
	null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
	expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
	null.alpha <- mean(null.alphas[,randomize])
	null.alpha.comp <- c(null.alpha.comp, null.alpha)
	
	bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
	diag(bucket_bray) <- NA
	bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
	
} ## end randomize loop

## Calculate beta-diversity for obs metacommunity
beta_comm_abund <- vegdist(dune, "bray")
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
diag(res_beta_comm_abund) <- NA
# output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
beta_div_abund_stochSD <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)

# output abundance beta-null deviation
abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)

plot(abund_null_dev,ylim=c(0,1)) # Suggest niche dynamics


dataBetaNull<-rbind(abund_null_dev,PA_null_dev,PA_null_devSD)
#write.csv(dataBetaNull, "ResultsBetaNull.csv")
}


dataBetaNULL<-read.table("ResultsBetaNull.csv", header=TRUE, sep=",", dec=".") ### ESTO ES DE CORRER LOS METODOS DE ARRIBA. TARDAN BASTANTE
summary(as.numeric(dataBetaNULL[1,-1]))

#jpeg("Figures/FigureX_BetaNullDeviation.jpeg", height=4, width=4, res=300, units="in")
x11(height=4,width=4)
par(mar=c(5,4,1,1))
plot(dist_fromOutlet,dataBetaNULL[1,-1], xlim=c(800,0),ylim=c(0,1), pch=19, xlab="Distance from outlet (Km)", ylab="", cex.lab=1)
mtext(expression(beta["Null deviation"]), 2,line=2, cex=1.5)
dev.off()

### COMPONENTS OF BETA DIVERSITY
betaSor<-beta.sor(t(data_abund))  # total beta
betaNes<-beta.nes(t(data_abund))# betanestedness
betaSim<-beta.sim(t(data_abund))# beta turnover

betaSor14<-beta.sor(t(data_abund1[-selSingle,]))
betaSim14<-beta.sim(t(data_abund1[-selSingle,]))
betaNes14<-beta.nes(t(data_abund1[-selSingle,]))
dist_mat
Dist_km2014<-c(Dist_km[1:5], 0, 724, 100,174, 353,483, 241,483,74,253,383,341,383,179,309, 415,309,130,594,130,724,0,724) # PARA USAR CON LOS DATOS 2014

dist2014col<-c(6,7,12,13,17,18,21,22,24,25,27,28) # Para seleccionar las comparaciones que son con 2014


## FIGURE 5 Distance decay in beta diversity 
#jpeg("Figures/Figure3_BetadiversityDistanceDecay.jpeg", height=5, width=2, res=300, units="in")
x11(height=5, width=2)
par(mfrow=c(3,1),mar=c(3,5,1,1))
plot(Dist_km,matrix(betaSor), type="p",ylim=c(0,1),ylab=expression(beta[sor]),cex.lab=1.5,xlab="",cex.axis=1.1, xlim=c(0,800), pch=19) #
points(Dist_km2014[dist2014col], matrix(betaSor14)[dist2014col], col="grey",pch=19)
legend("topleft", legend="A", bty="n")
legend("bottomright", legend=c("2013","2014"), pch=19, col=c("black", "grey"), bty="n")
par(mar=c(3,5,1,1))
plot(Dist_km,matrix(betaSim), pch=19,type="p", log="", ylab=expression(beta[sim]),ylim=c(0,1),cex.lab=1.5,xlab="",cex.axis=1.1, xlim=c(0,800)) # Variability is caused by replacement or turnover of species
points(Dist_km2014[dist2014col], matrix(betaSim14)[dist2014col], col="grey",pch=19)
legend("topleft", legend="B", bty="n")
par(mar=c(4,5,0,1))
plot(Dist_km,matrix(betaNes), pch=19, type="p",ylim=c(0,1), ylab=expression(beta[nest]),cex.lab=1.5,cex.axis=1.1, xlab="Distance (km)", xlim=c(0,800)) # not by nestedness (o sea que se ve un recambio de especies a medida que las estaciones estan mas lejos)
points(Dist_km2014[dist2014col], matrix(betaNes14)[dist2014col], col="grey",pch=19)
legend("topleft", legend="C", bty="n")

dev.off()

## SIGNIFICANT DISTANCE DECAY- not significan in nestedness
cor.test(Dist_km,matrix(betaSor)) 
cor.test(Dist_km,matrix(betaSim))
cor.test(Dist_km,matrix(betaNes))



########################################################################################################################################################
## READ SEQUENCES FROM FASTA FILE- Jeraldos Method
########################################################################################################################################################
data_fasta<-read.table("DATA/output_FASTA_Hamming.csv", header=TRUE, dec=".")# Sale de out.uc que envio Claudia Piccini siguiendo los scripts de Jeraldo- Borre los N
summary(data_fasta)

data_id_fasta<-data.frame(Q=data_fasta[,"LabelQuery"], T=data_fasta[,"LabelTarget"], Dist=data_fasta[,4])

#
data_abund<-data_abund1[-selSingle,1:6] # Saco los singletons o los que no habia aqui
data_names<-data[-selSingle,13:19]
data_ID<-data.frame(ID=data[-selSingle,1], data_abund)

data4fasta<-data.frame(data_ID, data_abund)

plot(sort(data4fasta[,2], decreasing=TRUE), log="y")

ID_dom<-list()
ID_rare<-list()
percent=0.01 # Las dominantes se definen como el 1% de las spp following Jeraldo et al 2012 PNAS  Quantification of the relative roles of niche and neutral
# Los resultados de los histogramas no cambian usando el 1%, 5% o 7% como dominantes
#processes in structuring gastrointestinal microbiomes www.pnas.org/cgi/doi/10.1073/pnas.1206721109
for(i in 1:6){

data_aux<-data_ID[,i+1]

S<-sum(data_aux>1)
DOM<-round(S*percent) # number of dominant species
selDOM<-which(data_aux>=sort(data_aux, decreasing=TRUE)[DOM])
data_select<-data_ID[selDOM,1]

selRARE<-which(data_aux<sort(data_aux, decreasing=TRUE)[DOM] & data_aux>1)
data_selectRARE<-data_ID[selRARE,1] ## Aqui seleccionar tambien la abundancia del sitio donde el OTU era dominante

ID_dom[[i]]<-data_select # Selecciono las dominantes por sitio.
ID_rare[[i]]<-data_selectRARE # seleeciono las raras 
}

# Calculo de distancia
DIST_MAT_site<-list()
for(j in 1:6){ # Salto a Punta del Este

ID_domSite<-ID_dom[[j]]
ID_rareSite<- ID_rare[[j]]
data_dist<-numeric()


for(k in 1:length(ID_dom[[j]])){
    
    DomOTU<-ID_domSite[k] # selecciono la primera OTU dominante
    
    OtherDom_2remove<-ID_domSite[-k] # Genero la lista de las demás OTUS dominantes
    
    selDistQ<-which(data_id_fasta[,1]%in%DomOTU) # selecciono las OTU dominantes en las secuencias de Distancias
    selDistT<-which(data_id_fasta[,c(2)]%in%DomOTU)

selDom2remove<-which(OtherDom_2remove%in%data_id_fasta[c(selDistQ,selDistT),c(1,2)]) # Reviso que las distancias no sean con otras secuencias dominantes
if(length(selDom2remove)>0) stop("Hay dominantes en la serie")

data_dist0<-data_id_fasta[c(selDistQ,selDistT),3]
data_dist<-c(data_dist, data_dist0)
    } # end loop K (dominant species)
DIST_MAT_site[[j]]<-(100-data_dist)/100

} # end loop j (SITES)

## NULL MODEL- From Jeraldo et al 2012
colnames(data_fasta)
summary(data_fasta[,3]) # Average sequence length =207
# HAMMING Function
# Input is a matriz with each OTU sequence
hammingAMS <- function(X) {
        uniqs <- unique(as.vector(X))
        U <- X == uniqs[1]
        H <- t(U) %*% U
        for ( uniq in uniqs[-1] ) {
            U <- X == uniq
            H <- H + t(U) %*% U
        }
    nrow(X) - H
}


alpha<- 0.01 # Fraction of species undergoing niche dynamics (1-alpha) species following neutral dynamics
N=5000 # Number of OTUs
L=200 # Length of genome sequence

TotAbund<-1e5 # Total abundance of the system


Neutral_SEQ_mat<- matrix(sample(c(1,2,3,4), L*(1-alpha)*N, replace=TRUE), ncol=(1-alpha)*N, nrow=L)

alpha*N # number of OTU following niche dynamics

Niche_center<-sample(c(1,2,3,4), L, replace=TRUE) # Random niche center

Niche_SEQ_mat<-matrix(Niche_center, nrow=L, ncol=alpha*N) # Mutations with respect to niche center

# Number of mutation in niche species is sampled from exponential distribution
NroMutations<-round(rexp((N*alpha)-1, rate=0.1))

for(i in 1:((N*alpha)-1)){
indice<-sample(1:200, NroMutations[i], replace=FALSE)

oldSeq<-Niche_SEQ_mat[indice,i+1]
NewBase<-sample(1:4, NroMutations[i],replace=TRUE)
while(sum(NewBase==oldSeq)>0) NewBase<-replace(NewBase, oldSeq==NewBase, sample(1:4, sum(oldSeq==NewBase),replace=TRUE)) # Secure there is mutation

Niche_SEQ_mat[indice,i+1]<-NewBase # Random mutations of species
}

dim(Niche_SEQ_mat)

DistNiche<-NroMutations
RelAbundNiche<- rexp(-DistNiche)
RelAbundNeutral<- 0.1*rexp((1-alpha)*N)

Conv_fac<-TotAbund/ sum(c(RelAbundNeutral,RelAbundNiche))
Abundance<-Conv_fac*c(RelAbundNiche,RelAbundNeutral)

SEQ_MAT<-cbind(Niche_SEQ_mat, Neutral_SEQ_mat)


## REVISAR ESTO QUE NO E VE BIEN! ESTAN LAS DOS MATRICES IGUALES
DIST_NULL_MODEL_niche<-hammingAMS(Niche_SEQ_mat)/L
DIST_NULL_MODEL_neutral<-hammingAMS(Neutral_SEQ_mat)/L


## FIGURE 6

x11(height=6, width=8)
par(mfrow=c(2,3),mar=c(5,4,1,1))
hist(DIST_MAT_site[[1]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02),freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="A", bty="n",cex=1.5)

hist(DIST_MAT_site[[2]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02), freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
#legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="B", bty="n",cex=1.5)

hist(DIST_MAT_site[[3]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02), freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
#legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="C", bty="n",cex=1.5)

hist(DIST_MAT_site[[4]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02), freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
#legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="D", bty="n",cex=1.5)

hist(DIST_MAT_site[[5]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02), freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
#legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="E", bty="n",cex=1.5)

hist(DIST_MAT_site[[6]], xlim=c(0,1),breaks=seq(0,0.36, by=0.02), freq=FALSE, main="",xlab="Distance to nearest OTU");box()
lines(density(DIST_NULL_MODEL_niche), col="darkgrey")
lines(density(DIST_NULL_MODEL_neutral), col="lightgrey",lty=2)
#legend("top", legend=c("Observed","Null Model-Niche", "Null Model- Neutral"), pch=c(22,NaN, NaN), lty=c(NaN, 1,2),bty="n", col=c("black","grey","grey"))
legend("topright", legend="F", bty="n",cex=1.5)

dev.off()


















