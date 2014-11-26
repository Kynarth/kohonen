#!/usr/bin/env Rscript
## FONCTIONS####################################################################

#___________________________Initialisations____________________________________#
init_matrice_vide <- function(hauteur,largeur,tailleVecteur)
{
  myArray <- array(dim=c(hauteur,largeur,tailleVecteur)) 
  return(myArray)
}

initialisation_matrice <- function(hauteur,largeur,tailleVecteur)
{
  myArray <- array(dim=c(hauteur,largeur,tailleVecteur)) 
  for (i in 1:dim(myArray)[1])
  {
    for (j in 1:dim(myArray)[2])
    {
      myArray[i,j,] <- sample(-180:180,tailleVecteur)
    }
  }
  return(myArray)
}

#_________________________Apprentissage________________________________________#
calcul_rayon <- function(r0,t,to)
{
  return(r0/(1+(t/to)))
}

calcul_taux_apprentissage <- function(n0,t,to)
{
  return(n0/(1+(t/to)))
}

mesure_distance_angle <- function(v,w,tailleVecteur)
{
  somme <- 0
  for (i in 1:(tailleVecteur/2))
  {
    somme <- somme + ((v[2*i-1]-w[2*i-1])^2+(v[2*i]-w[2*i])^2)
  }
  return((1/tailleVecteur)*somme)
}

recherche_distance_minimale <- function(v,myArray,tailleVecteur)
{
  coord <- c(0,0)
  distMinimale <- 10000000000
  dist <- c()
  for (i in 1:dim(myArray)[1])
  {
    for (j in 1:dim(myArray)[2])
    {
      dist <- mesure_distance_angle(v,myArray[i,j,],tailleVecteur)
      if(dist < distMinimale)
      {
        distMinimale <- dist
        coord <- c(i,j)
      }
    }
  }
  return(coord)
}

apprentissage <- function(hauteur,largeur,v,myArray,to,r0,n0,t,tailleVecteur,listeParam)
{
  t <<- t+1
  n <- calcul_taux_apprentissage(n0,t,to)
  r <- calcul_rayon(r0,t,to)
  winner <- recherche_distance_minimale(v,myArray,tailleVecteur)
  newArray <- init_matrice_vide(hauteur,largeur,tailleVecteur)
  for (i in 1:dim(myArray)[1])
  {
    for (j in 1:dim(myArray)[2])
    {
      coord <- c(i,j)
      for (k in 1:tailleVecteur)
      {
        wt <- myArray[i,j,k]
        newArray[i,j,k] <- wt + ((v[k]-wt)*n*exp(-(1/(2*r^2))*(sqrt((coord[2]-winner[2])^2+(coord[1]-winner[1])^2)^2)))
      }
    }
  }
  if((t%%100)==0)
  {
    title <- paste("check_plot_t=",t,".png",sep="")
    chemin<-paste("results/Kohonen/interm/img",listeParamInit,title,sep="/")
    png(file=chemin,width=1200,height=1500)
    schema <- c(c(0,0,0,0,1,0,0,0,0,0),c(2:101))
    layout(matrix(schema,11,10,byrow=TRUE))
    plot_vecteur(v,"V",0)
    afficher_plot(newArray,0)
    dev.off()
  }
  if((t%%10000)==0)
  {
    print(t)
    sauvegarde_matrice_kohonen(newArray,t,listeParam)
  }
  return(newArray)
}

#_____________________Determination pentapeptides______________________________#
remplissage_matrice_comptage <- function(myArray,myComptage,v)
{
  coord <- c()
  coord <- recherche_distance_minimale(v,myArray,length(v))
  myComptage[coord[1],coord[2]] <- myComptage[coord[1],coord[2]] + 1
  return(myComptage)
}

pourcentage <- function(x,nbVecteurTotal)
{
  x <- x/nbVecteurTotal*100
}

trier_matrice <- function(mat)
{
  mat3ColPourcentage <- matrix(0,nrow=100,ncol=3)
  colnames(mat3ColPourcentage)<-c("i","j","pourcentage")
  mat3ColPourcentage[,1] <- c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10))
  mat3ColPourcentage[,2] <- rep(c(1:10),10)
  mat3ColPourcentage[,3] <- as.vector(t(mat))
  
  mat3ColTriee <- mat3ColPourcentage[order(mat3ColPourcentage[,"pourcentage"],decreasing=TRUE),]
  return(mat3ColTriee)
}

creation_matrice_pentapeptide <- function(matPourcentage,matKohonen)
{
  alphabet<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T")
  matPentapeptides <- c()
  for(i in 1:length(alphabet))
  {
    matPentapeptides <- rbind(matPentapeptides,matKohonen[matPourcentage[i,"i"],matPourcentage[i,"j"],])
  }
  rownames(matPentapeptides) <- alphabet
  return(matPentapeptides)
}

#_____________________________Fichiers_________________________________________#
sauvegarde_matrice_kohonen <- function(myArray,t,listeParam)
{
  repertoire <- "results/Kohonen/interm/csv"
  nomDossier <- paste("sauvegarde-kohonen-t=",t,sep="")
  nomNouveauDossier <- paste(repertoire,listeParam,nomDossier,sep="/")
  dir.create(nomNouveauDossier)
  for(i in 1:dim(myArray)[3])
  {
    nomFichier <- paste(nomNouveauDossier,"/matrice-kohonen-W[",i,"].csv",sep="")
    write.table(myArray[,,i],file=nomFichier,sep=";",row.names=FALSE,col.names=FALSE)
  }
}

importation_matrice_kohonen <- function(cheminDossier)
{
  myArray <- init_matrice_vide(10,10,8)
  listeFichiers <- list.files(cheminDossier)
  for(i in 1:length(listeFichiers))
  {
    cheminFichier <- paste(cheminDossier,listeFichiers[i],sep="/")
    table <- read.csv(file=cheminFichier,sep=";",header=FALSE)
    table <- as.matrix(table)
    myArray[,,i]<-table
  }
  return(myArray)
}

#_______________________Fonctions Graphiques_____________________________#
plot_vecteur<-function(vecteur,title1,text)
{
  if(text==1)
  {
    plot(vecteur,ylim=c(-180,180),cex.main=2,type="l",xlab="",ylab="",xaxt='n',yaxt='n',main=title1)
    abline(h=0)
    text(position_texte(vecteur),labels=round(vecteur,2))
  }
  else
  {
    plot(vecteur,ylim=c(-180,180),cex.main=2,type="l",xlab="",ylab="",xaxt='n',yaxt='n',main=title1)
    abline(h=0)
  }
}

position_texte <- function(vecteur)
{
  positionTexte <- c()
  moyenne <- mean(vecteur)
  for(i in 1:length(vecteur))
  {
    if(vecteur[i] < moyenne)
    {
      positionTexte <- c(positionTexte,(vecteur[i]-15))
    }
    else
    {
      positionTexte <- c(positionTexte,(vecteur[i]+15))
    }
  }
  return(positionTexte)
}

afficher_plot <- function(myArray,text)
{
  for (i in 1:dim(myArray)[1])
  {
    for (j in 1:dim(myArray)[2])
    {
      title1 <- paste("w[",i,",",j,"]",sep="")
      plot_vecteur(myArray[i,j,],title1,text)  
    }
  }
}

# MAIN ##########################################################################

#_________________________Initialisations_______________________________________#
## Initialisation des parametres
# rayonInitial : doit être compris entre 1 et la taille de la matrice
# tauxApprInitial : doit être compris entre 0 et 1
rayonInitial <- 10
#taux d'apprentissage = 0.01 -> convergence trop rapide (mm pas en 1/5 de cycle)
tauxApprInitial <- 0.2
t <- 0
nbCycle <- 10
listeParamInit <- paste("Run-R=",rayonInitial,"-TxAppr=",tauxApprInitial,"-nbCycle=",nbCycle,sep="")

#Creation des dossiers de sauvegarde
#results/Kohonen/
cheminResultats <- "results/Kohonen/"
# Images intermédiaires
cheminIntermImg <- paste("results/Kohonen/interm/img",listeParamInit,sep="/")
dir.create(cheminIntermImg)
#results/Kohonen/csv/
cheminCsv <- paste("results/Kohonen/interm/csv",listeParamInit,sep="/")
dir.create(cheminCsv)

# Recuperation du tableau de vecteur
tableVecteur <- read.csv("results/vectors.csv",sep=";",header=FALSE)
tableVecteur <- as.matrix(tableVecteur)
tailleVect <- dim(tableVecteur)[2]
nbVecteurTotal <- dim(tableVecteur)[1]

# Skip le probleme des Nas
lesNa <- which(is.na(tableVecteur),arr.ind=TRUE)
if(!is.null(lesNa))
{
  for(i in 1:dim(lesNa)[1])
  {
    tableVecteur[lesNa[i,1],lesNa[i,2]] <- 0
  }
}

# Initialisation de la matrice de Kohonen
matKohonen <- initialisation_matrice(10,10,8)
png(file=paste(cheminResultats,"avant_apprentissage.png",sep="/"),width=1200,height=1200)
par(mfrow=c(10,10),pin=c(1,1))
afficher_plot(matKohonen,1)
dev.off()

#_______________________________Apprentissage___________________________________#
for (k in 1:nbCycle)
{
  # Creation d'une liste de nombre aleatoire qui correspond au tirage au sort des vecteurs
  listeVecteurAleatoire <- sample(1:nbVecteurTotal,nbVecteurTotal,replace=F)
  for (i in 1:nbVecteurTotal)
  {
    v <- tableVecteur[listeVecteurAleatoire[i],]
    testVecteur <- which(is.na(v),arr.ind=TRUE)
    if(length(testVecteur)==0)
    {
      matKohonen <- apprentissage(10,10,v,matKohonen,nbVecteurTotal,rayonInitial,tauxApprInitial,t,tailleVect,listeParamInit)
    }
  }
}

png(file=paste(cheminResultats,"apres_apprentissage.png",sep="/"),width=1200,height=1200)
par(mfrow=c(10,10),pin=c(1,1))
afficher_plot(matKohonen,1)
dev.off()

#_____________________Determination des 26 peptides_____________________________#
matriceComptage <- matrix(0,nrow=10,ncol=10)

for(i in 1:nbVecteurTotal)
{
  v <- tableVecteur[i,]
  matriceComptage<-remplissage_matrice_comptage(matKohonen,matriceComptage,v)
}

matriceComptagePourcentage <- apply(matriceComptage,c(1,2),pourcentage,nbVecteurTotal=nbVecteurTotal)
matricePourcentageTriee <- trier_matrice(matriceComptagePourcentage)
matricePentatpeptides <- creation_matrice_pentapeptide(matricePourcentageTriee,matKohonen)

# Graphiques avec la matrice de pentapeptides
png(file=paste(cheminResultats,"matrice-pentapeptides.png",sep="/"),width=1200,height=1200)
par(mfrow=c(5,4),pin=c(2,2))
for(i in 1:dim(matricePentatpeptides)[1])
{
  vecteur <- matricePentatpeptides[i,]
  plot(vecteur,ylim=c(-180,180),cex.main=2,type="l",xlab="",ylab="",xaxt='n',yaxt='n',
       main=rownames(matricePentatpeptides)[i],
       sub=paste(round(matricePourcentageTriee[i,3],2),"%",sep=" "))
  abline(h=0)
  text(position_texte(vecteur),labels=round(vecteur,2))
}
dev.off()
write.table(matricePentatpeptides,file=paste(cheminResultats,"matrice-pentapeptide.csv",sep="/"),sep=";",row.names=TRUE,col.names=FALSE)
