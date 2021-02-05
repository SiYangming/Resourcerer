\name{getResourcerer}
\alias{getResourcerer}
\title{A function that downloads an annotation file from TIGR
  Resourcerer and then read the downloaded file in as a matrix}
\description{
  TIGR Resourcerer maintains various annotation files for Affymetrix or
  cDNA chips. This function allows users to read an annotation file into R
  as a matrix.
}
\usage{
getResourcerer(which, organism, destDir =
file.path(path.package("Resourcerer"), "temp"), baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer", clean = TRUE, exten = "zip")  
}

\arguments{
  \item{which}{\code{which} a character string indicating which
    annotation file to be read in. The annotation files are stored in
    subdirectories for various organisms under \code{baseUrl} (see
    below)} 
  \item{organism}{\code{organism} a character string for the name of the
    organism of interests}
  \item{destDir}{\code{destDir} a character string for the path of a
    directory where the downloaded file will be stored. If missing, the
    temp directory will be the default}
  \item{baseUrl}{\code{baseUrl} a character string for the url of
    Resourcerer ftp site where directories containing annotation files
    for human, rat, mouse ... are stored}
  \item{clean}{\code{clean} a boolean indicating whether the file
    downloaded from Resourcerer will be removed after the data contained
    have been read in}
  \item{exten}{\code{exten} a character string for the extension
    (e. g., zip) of the source data file to be processed}
}

\details{
  \code{baseUrl} is the root directory of TIGR ftp site for Resourcerer
  that contains subdirectories holding data for different organism.  
}
\value{
  Function getResourcerer returns a matrix derived from the source
  data. Column names of the returned matrix are taken directly from the
  source file provided by Resourcerer. Users are advised to visit the
  Resourcerer web site for more information about the source data files. 
}

\references{\url{http://pga.tigr.org/tigr-scripts/magic/r1.pl}}

\author{Jianhua Zhang}

\seealso{\code{"\link{resourcerer2BioC}"}}

\examples{
  resourcerer <- getResourcerer("Agilent_Human1_cDNA.zip", organism = "human",
       destDir = file.path(path.package("Resourcerer"), "temp"),
       baseUrl = "ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer",
       clean = TRUE, exten = "zip")
  resourcerer[1:3, ]
}

\keyword{manip}






