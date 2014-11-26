#!/usr/bin/env Rscript
# FONCTIONS ########################################################################
creation_matrice_vide <- function(data,matPB)
{
  aa<- unique(data[,"AA"])
  pb <- rownames(matPB)
  aa <- as.character(aa)
  matriceVide <- matrix(0, length(aa), length(pb))
  rownames(matriceVide) <- aa
  colnames(matriceVide) <- pb
  return(matriceVide)
}

calcul_pourcentage <- function(case,vecteur)
{
  proportion <- case/sum(vecteur)*100
  return(proportion)
}

proportion_table <- function(data,matPB,tableAModifier)
{
  tableModifiee <- creation_matrice_vide(data,matPB)
  for (i in 1:dim(tableAModifier)[1])
  {
    for (j in 1:dim(tableAModifier)[2])
    {
      tableModifiee[i,j]<-round(calcul_pourcentage(tableAModifier[i,j],tableAModifier[i,]),2)
    }
  }
  return(tableModifiee)
}

mesure_distance_angle <- function(v,w)
{
  somme <- 0
  for (i in 1:(8/2))
  {
    somme <- somme + ((v[2*i-1]-w[2*i-1])^2+(v[2*i]-w[2*i])^2)
  }
  return((1/8)*somme)
}

recherche_distance_minimale <- function(matPB,vAA)
{
  PB <- c()
  distMinimale <- 10000000000
  dist <- c()
  for (i in 1:dim(matPB)[1])
  {
    dist <- mesure_distance_angle(matPB[i,],vAA)
    if(dist < distMinimale)
    {
      distMinimale <- dist
      PB <- c(rownames(matPB)[i])
    }
  }
  return(PB)
}

# MAIN #############################################################################
#_________________________Importation des données__________________________________#

## mutual.csv
data <- read.table("results/mutual.csv", sep=";", dec=".", header=FALSE)
label <- c("PDB","Numero","résidu.chaine","AA","STRUCTURE","?1","?2","?3","?4","?5","?6","?7","BP1","BP2","x","ACC","N.H..>O","O..>H.N","N.H..>O","O..>H.N","TCO","KAPPA","ALPHA","PHI","PSI","X.CA","Y.CA","Z.CA")
colnames(data) <- label
data[,"AA"] <- toupper(data[,"AA"])
data <- data.frame(data)
data[,"PDB"] <- as.character(data[,"PDB"])
data[,5] <- as.character(data[,5])
data[,6] <- as.character(data[,6])
data[,7] <- as.character(data[,7])
data[,8] <- as.character(data[,8])
data[,9] <- as.character(data[,9])
data[,10] <- as.character(data[,10])
data[,11] <- as.character(data[,11])
data[,12] <- as.character(data[,12])

data[data[,6]=="",6] <- " "
data[data[,7]=="",7] <- " "
data[data[,8]=="",8] <- " "
data[data[,9]=="",9] <- " "
data[data[,10]=="",10] <- " "
data[data[,11]=="",11] <- " "
data[data[,12]=="",12] <- " "

structure <- paste(data[,5],data[,6],data[,7],data[,8],data[,9],data[,10],data[,11],data[,12],sep="")
data <- cbind(data[,1:4],"STRUCTURE"=structure,data[,13:28])

## matrice-pentapeptide.csv
matPB <- read.table("results/Kohonen/matrice-pentapeptide.csv",sep=";",dec=".",header=FALSE)
rownames(matPB) <- matPB[,1]
matPB <- as.matrix(matPB[,2:9])

#______________________________Initialisation______________________________________#

## Matrice AA / PB
matAAPB <- creation_matrice_vide(data,matPB)

## Tableau avec les PB
tableauAvecPB <- c()
tableauAvecPB <- cbind(data[,1:4],"PB"="",data[,5:21])
tableauAvecPB[,"PB"]<-as.character(tableauAvecPB[,"PB"])
tableauAvecPB[,"Numero"]<-as.numeric(tableauAvecPB[,"Numero"])

cheminFichierPb <- paste("results/Kohonen/","listePB.txt",sep="")
file.create(cheminFichierPb)
con <- file(cheminFichierPb, open = "w")
#_______________________________Comptage___________________________________________#
for(prot in unique(data[,"PDB"]))
{
  cat(paste("XD ",prot,"\n",sep=""),file=con)
  ligneFichier <- ""
  table <- data[data[,"PDB"]==prot,c("PDB","Numero","AA","STRUCTURE","PHI","PSI")]  
  compteur80caracteres <- 0
  for(i in 3:(dim(table)[1]-2))
  {
    compteur80caracteres <- compteur80caracteres+1
    AA <- table[i,"AA"]
    angles <- c(table[i-2,"PSI"],table[i-1,"PHI"],table[i-1,"PSI"],table[i,"PHI"],table[i,"PSI"],table[i+1,"PHI"],table[i+1,"PSI"],table[i+2,"PHI"])
    PB <- recherche_distance_minimale(matPB,angles)
    matAAPB[AA,PB] <- matAAPB[AA,PB]+1
    tableauAvecPB[tableauAvecPB[,"PDB"]==prot & tableauAvecPB[,"Numero"]==i,"PB"] <- PB 
    ligneFichier <- paste(ligneFichier,PB,sep="")
    if( compteur80caracteres == 79 )
    {
      compteur80caracteres <- 0
      ligneFichier <- paste(ligneFichier,"\n",sep="")
      cat(ligneFichier,file=con)
      ligneFichier <- ""
    }
  }
  ligneFichier <- paste(ligneFichier,"\n",sep="")
  cat(ligneFichier,file=con)
  cat("-------------------------------------------------------------------------------\n",file=con)
}
close(con)

for(prot in unique(tableauAvecPB[,"PDB"]))
{
  table <- tableauAvecPB[tableauAvecPB[,"PDB"]==prot,c(2:dim(data)[2])] 
  colnames(table)[1]<-"#"
  write.table(table,file=paste("results/Kohonen/newDssp/",prot,".dssp",sep=""),sep=";",dec=",",row.names=TRUE,col.names=NA)
}

matAAPBPourcentage <- proportion_table(data,matPB,matAAPB)
write.table(matAAPBPourcentage,file="results/Kohonen/matrice-AA-PB.csv",sep=";",dec=".",row.names=TRUE,col.names=NA)
