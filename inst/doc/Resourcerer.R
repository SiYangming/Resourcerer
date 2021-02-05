### R code from vignette source 'Resourcerer.Rnw'

###################################################
### code chunk number 1: Resourcerer.Rnw:40-42
###################################################
require("AnnotationDbi", quietly = TRUE)
require("Resourcerer", quietly = TRUE)


###################################################
### code chunk number 2: Resourcerer.Rnw:56-60
###################################################
agilent <- getResourcerer("Agilent_Human1_cDNA.zip", organism = "human",
destDir = file.path(path.package("Resourcerer"), "temp"), baseUrl =
            "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
            clean = TRUE )


###################################################
### code chunk number 3: Resourcerer.Rnw:69-70
###################################################
agilent[1:5, 1:4]


###################################################
### code chunk number 4: Resourcerer.Rnw:75-76
###################################################
as.vector(colnames(agilent))


###################################################
### code chunk number 5: Resourcerer.Rnw:94-109
###################################################

if(interactive()) {
    resourcerer2BioC("Agilent_Human1_cDNA.zip", organism = "human",
       destDir =  file.path(path.package("Resourcerer"), "temp"),
       pkgName = "AgilentHsa1",
       srcUrls = getSrcUrl("all", "Homo sapiens"),
       pkgPath =  file.path(path.package("Resourcerer"), "temp"),
       otherSrc = NULL, baseMapType = "gb",
       version = "1.1.0", fromWeb = TRUE,
       baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
               check = TRUE, author = list(authors = "Anonymous",
                              maintainer = "Anonymous <anonimous@email.com>"))
   }else{
       print("Code is executed only when invoked interactively")
   }


