\name{resourcerer2BioC}
\alias{resourcerer2BioC}
\title{A function that downloads an annotation file from TIGR
  Resourcerer and then creates a bioC annotation data package}
\description{
  TIGR Resourcerer maintains various annotation files for Affymetrix or
  cDNA chips. This function allows users to create a bioC annotation
  data package for the probes contained in the Resourcerer annotation
  file.
}
\usage{ 
resourcerer2BioC(
which, organism = c("human", "mouse", "rat"),
    destDir = file.path(path.package("Resourcerer"), "temp"),
    pkgName, pkgPath, baseMapType = c("gbNRef", "gb", "ug", "ll"),
    version = "1.1.0", baseUrl =
    "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer", check =
    FALSE, author = list(authors = "Anonymous", maintainer
     = "Anonymous <anonymous@email.com>"), exten = "zip")

}
\arguments{
  \item{which}{\code{which} a character string indicating which
    Resourcerer annotation file to be read in}
  \item{destDir}{\code{destDir} a character string for the path of a
    directory where the downloaded file will be stored. If missing, the
    temp directory will be the default}
  \item{baseUrl}{\code{baseUrl} a character string for the url of
    Resourcerer ftp site where directories containing annotation files
    for human, rat, mouse ... are stored}
  \item{baseMapType}{\code{baseMapType} a character string that is
    either "gb","ug", or "ll" to indicate whether the probe ids in
    baseName are mapped to GenBack accession numbers, UniGene ids, or
    LocusLink ids}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built (e. g. hgu95a, rgu34a)}
  \item{pkgPath}{\code{pkgPath} a character string for the full path of
    an existing directory where the built backage will be stored}
  \item{organism}{\code{organism} a character string for the name of the
    organism of concern (now can only be "human", "mouse", or "rat")}
  \item{version}{\code{version} a character string for the version number}
  \item{author}{\code{author} a list of character strings with an author
    element for the name of the author and maintainer element for the
    email address of the author}
  \item{check}{\code{check} a boolean incicating whether to check the
    mappings between probe ids and LocusLink ids obtained from Resourcer
    and bioC AnnotationDbi when baseMapType is 'll'}
  \item{exten}{\code{exten} a character string for the extension
    (e. g. zip) of the source data file to be processed}
}
\details{

  \code{baseUrl} is the root directory of TIGR ftp site for Resourcerer
  that contains subdirectories holding data for different organism.
}
\value{

  Function resourcerer2BioC returns invisible() if successfully executed.
}
\references{\url{http://pga.tigr.org/tigr-scripts/magic/r1.pl}}
\author{Jianhua Zhang}

\seealso{\code{\link{getResourcerer}}}
\examples{
  #############################################################
  ## The example takes a loooong time (about an hour) to run ##
  #############################################################
  if(interactive()){
    resourcerer2BioC("Agilent_Human1_cDNA.zip",  baseMapType = "gbNRef")
    unlink(file.path(path.package("Resourcerer"), "temp",
           "Agilent_Human1_cDNA.db"), TRUE)
  }    
}

\keyword{manip}




