# Downloads the zip file (e. g. Agilent_Human1_cDNA.zip) defined by
# "which" from TIGR's RESOURCERER and then reads the data into R as a
# matrix
# which - a character string for the name of the zip file to be
#         downloaded from TIGR Resourcerer.
# AnnotationDbi is required
#
# Copyright 2004, J. Zhang. All rights reserved
#

getResourcerer <- function(which, organism,
    destDir = file.path(path.package("Resourcerer"), "temp"),
    baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
                clean = TRUE, exten = "zip"){

    fileName <- loadResourcerer(which, organism, destDir, baseUrl)
    resData <- as.matrix(read.table(fileName, skip = 1, sep = "\t",
                       header = FALSE, quote = "", comment.char = "",
                       colClasses = "character", strip.white = TRUE))
    resColName <- resData[1,]
    resData <- resData[2:nrow(resData), ]
    resData[resData == ""] <- NA
    colnames(resData) <- resColName
    if(clean){
        unlink(fileName)
    }

    return(resData)
}

# Downloads the zip file (e. g. Agilent_Human1_cDNA.zip) defined by
# "which" from TIGR's RESOURCERER and then returns the name of the
# expanded file.
# which - a character string for the name of the zip file to be
#         downloaded from TIGR.
#
# AnntationDbi is required
loadResourcerer <- function(which, organism, 
    destDir = file.path(path.package("Resourcerer"), "temp"),
    baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
                           exten = "zip"){
    switch(tolower(organism),
           human = url <- paste(baseUrl, "Human", which, sep = "/"),
           mouse = url <- paste(baseUrl, "Mouse", which, sep = "/"),
           rat = url <- paste(baseUrl, "Rat", which, sep = "/"),
           arabidopsis = url <- paste(baseUrl, "Arabidopsis", which,
             sep = "/"),
           barley = url <- paste(baseUrl, "Barley", which, sep = "/"),
           celegans = url <- paste(baseUrl, "C.elegans", which, sep = "/"),
           cattle = url <- paste(baseUrl, "Cattle", which, sep = "/"),
           grape = url <- paste(baseUrl, "Grape", which, sep = "/"),
           pine = url <- paste(baseUrl, "Pine", which, sep = "/"),
           potato = url <- paste(baseUrl, "Potato", which, sep = "/"),
           rice = url <- paste(baseUrl, "Rice", which, sep = "/"),
           sorghum = url <- paste(baseUrl, "Sorghum", which, sep = "/"),
           soybean = url <- paste(baseUrl, "soybean", which, sep = "/"),
           Xenopus = url <- paste(baseUrl, "Xenopus", which, sep = "/"),
           yeast = url <- paste(baseUrl, "Yeast", which, sep = "/"),
           zebrafish = url <- paste(baseUrl, "Zebrafish", which, sep = "/"),
           stop("Unknown organism name"))
    # destDir is not used untile the functions are included in a package
    temp <- loadFromUrl(url, destDir)
    unlink(temp)

    return(file.path(destDir, gsub(paste("\\.", exten, sep = ""), "", which)))
}

# Takes the name for an expanded TIGR cDNA file stored locally and
# returns a matrix mapping TIGE cDNA probe ids to either GenBank,
# LocusLink, or UniGene ids as defined by "which".
# tigrFile - a character string for an expanded TIGR Resourcerer
# annotation file (e. g. Agilent_Human1_cDNA).
# This function is intended to be used to get a base file ready for
# AnnotationDbi to use to build an annotation package
#

getProbe2ID <- function(tigrFile, baseMapType = c("gbNRef", "gb", "ug", "ll")){
    which <- match.arg(baseMapType)
    switch(which,
           gbNRef = map <- 3,
           gb = map <- 3,
           ll = map <- 5,
           ug = map <- 4)
    tigr <- as.matrix(read.table(tigrFile, skip = 2, sep = "\t",
                       header = FALSE, quote = "", comment.char = "",
                       colClasses = "character", strip.white = TRUE))
    tigr <- tigr[, c(1, map)]
    tigr[tigr == ""] <- "NA"
    colnames(tigr) <- c("Probe", which) 

    return(tigr)
}

# Build a BioC data package using probe ids and their mappings to GB,
# UG, or LL ids derived from Resourcerer
resourcerer2BioC <- function(which, organism = c("human", "mouse", "rat"),
       destDir =  file.path(path.package("Resourcerer"), "temp"),
       pkgName, pkgPath,
       baseMapType = c("gbNRef", "gb", "ug", "ll"), version = "1.1.0",
       baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
               check = FALSE, author = list(authors = "Anonymous",
               maintainer = "Anonymous <anonymous@email.com>"),
                             exten = "zip"){
    organism <- match.arg(organism)
    baseMapType = match.arg(baseMapType)

    if(missing(pkgName)){
        pkgName <- gsub(paste("\\.", exten, sep = ""), "", which)
    }
    if(missing(pkgPath)){
        pkgPath <- file.path(path.package("Resourcerer"), "temp")
    }

    baseFile <- file.path(path.package("Resourcerer"), "temp",
                          basename(tempfile()))
    unziped <- loadResourcerer(which, organism, destDir, baseUrl, exten)
    base <- getProbe2ID(unziped, baseMapType)
    write.table(base, baseFile, quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE)

    schema <- getSchema(tolower(organism))

    tryMe <- try(makeDBPackage(schema, fileName = baseFile,
                              affy = FALSE, prefix = pkgName, 
                  baseMapType = baseMapType, version = "1.0.0",
                  chipName = "Resourcerer annotation file",
                              outputDir = pkgPath))
    unlink(unziped)
    if(inherits(tryMe, "try-error")){   
        return("Unsuccessful")
    }
    return("Successful")
}

getSchema <- function(organism){
  if(tolower(organism) == "human"){
      if(!require(human.db0)){
        cat("human.db0 is not available. Install from BioC")
        source("http://www.bioconductor.org/biocLite.R")
        biocLite("human.db0")
      }
      return("HUMANCHIP_DB")
  }
  if(tolower(organism) == "mouse"){
      if(!require(mouse.db0)){
        cat("mouse.db0 is not available. Install from BioC")
        source("http://www.bioconductor.org/biocLite.R")
        biocLite("mouse.db0")
      }
      return("MOUSECHIP_DB")
  }
  if(tolower(organism) == "rat"){
     if(!require(rat.db0)){
        cat("mouse.db0 is not available. Install from BioC")
        source("http://www.bioconductor.org/biocLite.R")
        biocLite("mouse.db0")
      }
      return("RATCHIP_DB")
   }
}


# map2LL - a character string to a tab separated file with the first
# column for probe id and second column for LL id
# llRda - a character string for the name of an rda file for an
# environment with keys being probe ids and values for matching LL
# ids.

checkMapping <- function(pkgName, map2LL, llRda, outFile =
    file.path(path.package("Resourcerer"), "temp", "checkMapping.out")){

    output <- NULL
    load(llRda)
    bioC <- as.list(get(gsub("\\.rda", "", basename(llRda))))
    bioC <- cbind(names(bioC), unlist(bioC))
    tigr <- read.table(map2LL, header = FALSE, sep = "\t", as.is = TRUE,
                       strip.white = TRUE, colClasses = "character")
    combined <- as.matrix(merge(bioC, tigr, by.x = "V1", by.y = "V1",
                                all.x = TRUE))
    combined[is.na(combined)] <- "NA"
    diff <- combined[, 2] == combined[,3]
    output <- c(output, paste(pkgName, paste("Total =", nrow(combined)),
                paste("Match =", length(diff[diff])),
                paste("Mismatch =", length(diff[!diff])), sep = "\t"))
    write(output, outFile, append = TRUE)
}












