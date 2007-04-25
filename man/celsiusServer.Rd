\name{celsiusServer} %__CUT__%
\docType{class}
\alias{celsiusServer}
\alias{celsiusServer-class}
\title{ Object to hold values pertaining to a Celsius Server }
\description{ Holds values like the URL needed to reach the server or
  the default timings for accessing the features and samples for the
  server etc.}%- maybe also 'usage' for other objects documented here.
\section{Creating Objects}{
  \code{new('celsiusServer',}\cr
  \code{celsiusUrl = ...., #URL}\cr
  \code{platform = ...., #Platform for data retrieval }\cr
  \code{protocol = ...., #Normalization protocol for the data to be
    retrieved }\cr
  \code{stringency = ...., #Stringency setting for annotation retrieval }\cr
  \code{verbose = ...., #Verbosity preference }\cr    
  \code{transform = ...., #Transform the annotations vector into 0's and
    1's }\cr                
  \code{retry = ...., #Number of retry attempts with web server }\cr 
  \code{sampleTime  = ...., #estimated time to retrieve a sample}\cr
  \code{featureTime = ...., #estimated time to retrieve a feature}\cr
  \code{ )}}  
\section{Slots}{
  \describe{
    \item{\code{celsiusUrl}}{ Provide an approriate URL
      \code{"character"} to talk to a celsius server }
    \item{\code{platform}}{ Provide an approriate affymetrix platform
      \code{"character"} for data retrieval.  The preferred format is
      the affymetrix format, but you can also enter Bioconductor
      formatted platform names. Default is 'HG-U133A' }
    \item{\code{protocol}}{ Provide an approriate normalization protocol
      \code{"character"} for the data to be retrieved. Values can be:
      'rma','vsn','gcrma','plier' and 'mas5'.  Default value is 'rma'. }
    \item{\code{stringency}}{ Provide an approriate stringency setting
      \code{"numeric"}.  This is the quality of annotation desired.  A
      higher value results in more NA (un-annotated) values for
      annotation retrieval. Possible values include:
      600,500,400,300,200,100.  Default is 500. }
    \item{\code{verbose}}{ Indicate \code{"logical"} whether to print
      output to the screen whenever data has been retrieved from the
      DB. Default is FALSE. }    
    \item{\code{transform}}{ Indicate \code{"logical"} whether or not to
      format all annotation mask data into a list of only 0's and 1's
      (ie. to list NAs and and annotation conficts as 0's) Default is
      FALSE. }                
    \item{\code{retry}}{ Provide the number of retry attempts
      \code{"character"} for the client to attempt to get data from the
      web server. Default is 3. }                
    \item{\code{sampleTime}}{ Provide an expected wait time
      \code{"numeric"} to retrieve a sample from this platform (can
      determine this by using time and the getSampleData() in this
      library if you REALLY hate waiting) }
    \item{\code{featureTime}}{ Provide an expected wait time
      \code{"numeric"} to retrieve a feature from this platform (can
      determine this by using time and the getFeatureData() in this
      library if you REALLY hate waiting)}
  }
}
\author{ Marc Carlson }
\keyword{ classes }
