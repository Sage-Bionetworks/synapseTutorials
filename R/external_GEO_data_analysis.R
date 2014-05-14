###############################
# 
# originally written by Brian Bot @ Sage Bionetworks
###############################


# LOAD THE REQUIRED LIBRARIES
require(synapseClient)
require(affy)
require(simpleaffy)
require(limma)


#LOGIN TO SYNAPSE
synapseLogin()

#CREATE A PROJECT ON SYNAPSE
random_string = paste(sample(LETTERS,4),collapse="")
proj_name = paste('Demo_external_GEO_data_analysis',random_string,sep='-')
myProj <- Project(name=proj_name)
myProj <- synStore(myProj)



folderId = myProj$properties$id

## ALTERNATIVELY PROVIDE THE FOLDER ID WHERE YOU WANT TO STORE YOUR FILES
#folderId <- ""



## REFERENCE A FILE STORED IN GEO AND REGISTER WITH SYNAPSE
f <- File(path="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4757/suppl/GSE4757_RAW.tar", 
          parentId=folderId, synapseStore=FALSE)
f <- synStore(f)

f <- synGet(f@properties$id)

## UNTAR THE LOCAL COPY OF THE FILE
celDir <- file.path(tempdir(), "celfiles")
dir.create(celDir)
untar(getFileLocation(f), exdir=celDir, tar="/usr/bin/tar")

## READ IN THE CEL FILES
ab <- ReadAffy(celfile.path=celDir)
scanDates <- as.Date(substr(ab@protocolData@data$ScanDate, 1, 8), format="%m/%d/%y")
ab <- ab[, order(scanDates)]
scanDates <- scanDates[order(scanDates)]

## PRE-NORMALIZATION BOXPLOTS
boxPath <- file.path(tempdir(), "gse4757-boxplot-prenormalization.png")
png(boxPath, width=600, height=400)
boxplot(log2(exprs(ab)), range=0, col=factor(scanDates), xlab="CEL File", ylab="log2(expr)", main="Pre-Normalization colored by Scan Date")
dev.off()

boxFile <- synStore(File(path=boxPath, parentId=folderId), used=f)
onWeb(boxFile)

## SIMPLE RMA NORMALIZATION
eset <- rma(ab)
## PRE-NORMALIZATION BOXPLOTS
boxPathPost <- file.path(tempdir(), "gse4757-boxplot-postnormalization.png")
png(boxPathPost, width=600, height=400)
boxplot(exprs(eset), range=0, col=factor(scanDates), xlab="CEL File", ylab="log2(expr)", main="Post-Normalization (RMA) colored by Scan Date")
dev.off()

boxFilePost <- synStore(File(path=boxPathPost, parentId=folderId), used=f)
onWeb(boxFilePost)

