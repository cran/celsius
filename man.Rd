%$Id: man.Rd,v 1.45 2007/05/24 18:04:32 mcarlson Exp $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{celsius} %__CUT__%
\alias{celsius}
\alias{celsius-package}
\docType{package}
\title{
Retrieve Affymetrix microarray measurements and metadata from Celsius
}
\description{
Retrieve Affymetrix microarray measurements and metadata from Celsius
web services, see http://genome.ucla.edu/projects/celsius.
}
\details{
\tabular{ll}{
Package: \tab celsius\cr
Type: \tab Package\cr
Version: \tab 1.0.7\cr
Date: \tab 2007-05-21\cr
License: \tab Artistic License\cr
}
  This library is intended to provide a set of accessor methods to a
  celsius style database.  To use any of these functions, you must 1st
  create a Celsius Server object that points to an appropriate instance
  of a Celsius Database.  The default parameters point to the initial
  instance.

  Once this is created, there are a host of usable functions which can
  be called as needed Here is a brief overview:

  
  getPlatformData(): returns all affy platforms supported by the server

  listSamplesForPlatform(): lists all Sample IDs for a particular platform

  mapIdToCelsius(): this maps foreign Sample IDs to the Celsius Server instance

  getSampleData(): Returns data based on Sample IDs

  getFeatureData(): Returns data based on features (affy probeset IDs)

  getAssayData(): Returns an intersection of data for a set of both
  features and sample IDs
  
  getCvterms(): returns information (such as ID) about available
  controlled vocabulary terms

  getOntologyData(): Returns phenotypic data based on list of controlled
  voculary term IDs
  
  getPhenoData(): returns an intersection of data for controlled
  vocabulary terms and specific Sample IDs

  makeExprSet(): makes an expression set object from supplied IDs,
  features and/or controlled vocabulary terms

}
\author{
Allen Day and Marc Carlson
Maintainer: Marc Carlson <mcarlson@fhcrc.org>
}
\keyword{ data }

\seealso{
  There is also individual documentation for the major functions listed above.
}
%~~ simple examples of the most important functions ~~
\examples{

#first make an instance of a celsius Server object.
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu" );

#list all of the platforms supported by the server
platformTable = getPlatformData( celsius );

#list all the samples in the server that correspond to the
#HG-U133_Plus_2 platform...
celsius@platform = "HG-U133_Plus_2";
idList = listSamplesForPlatform( celsius );

#map some foreign IDs to the celsius servers internal sample IDs
celsius@platform = 'HG-U133A';
id = c( "AEX:E-MEXP-402", "GSE6210" );
mapTable = mapIdToCelsius( celsius, id = id );

#get some samples data based on sample IDs
celsius@verbose = TRUE;
data = getSampleData( celsius, id = c( "SN:1005595" ) );

#get some data based on features (affymetrix probesets).
data = getFeatureData( celsius, feature = c( "204412_s_at", "1316_at" ) );

#get some data that corresponds to both a set of Sample IDs AND some
#features (affymetrix probesets)
celsius@verbose = FALSE;
data = getAssayData( celsius,
  id = c( "SN:1018351", "SN:1043233" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at" ) );

#get some information about controlled vocabulary terms
#(ID is probably most important here)
lungTerms = getCvterms( celsius, q = "lung", db = "MA" );

#get data on all samples (on one platform) for a particular controlled
#vocabulary term
data = getOntologyData( celsius, cvterm=c( "MA:0000415" ) );

#get data on all samples (on one platform) for a controlled vocabulary
#term on specific samples
data = getPhenoData( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039", "SN:1004962" ),
  cvterm = c( "MA:0000415" ) );

#finally, we can make an expression set object straight from the
#database using IDs and controlled vocabulary terms
exprSet = makeExprSet( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039" ),
  cvterm = c( "MA:0000415" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at" ) );

#Or you can also make an expression set using a foreign ID:
celsius@platform = "Mouse430_2";
exprSet = makeExprSet(celsius,
  id = c( "GSE6210" ),
  cvterm = c( "MA:0000415" ),
  feature = c( "1415670_at", "1430537_at" ), mapids = TRUE );

#Finally, you can get raw annotations for an expression set by using the
#annotations parameter
celsius@platform = 'HG-U133A';
exprSet = makeExprSet( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039" ),
  cvterm = c( "MA:0000415" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at" ),
  annotations = TRUE );
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getOntologyData} %__CUT__%
\alias{getOntologyData}
\title{ Retrieve Annotations into R from Celsius webserver }
\description{
Retrieves annotation data for Affymetrix arrays in multi-graph (MGR) format
}
\usage{
getOntologyData(celsiusObject, cvterm=c());
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{cvterm}{ character vector of ontology term accessions }
}
\value{
  Return value is a table with rows labels as SN identifiers and column
  labels as probeset identifiers.
}
\author{ Jun Dong, Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  platform = "HG-U133A",
  stringency = 500 );

#get the data
data = getOntologyData( celsius,
  cvterm = c( "MA:0000168","MPATH:458" ) );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getFeatureData} %__CUT__%
\alias{getFeatureData}
\title{ Retrieve data (probeset-wise) into R from Celsius webserver }
\description{
Retrieves quantified Affymetrix data in multi-graph (MGR) format by
feature }
\usage{
getFeatureData(celsiusObject, feature = c() );
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{feature}{ character vector of feature IDs, database accession
    identifiers for probeset records used in Celsius }
}
\value{
  Return value is a table with rows labels as SN identifiers and column
  labels as probeset identifiers.
}
\author{ Jun Dong, Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  platform = "HG-U133A",
  protocol = "rma" );

#get the data
data = getFeatureData(celsius, feature = c( "204412_s_at", "1316_at") );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getSampleData} %__CUT__%
\alias{getSampleData}
\title{ Retrieve data (sample-wise) into R from Celsius webserver }
\description{
Retrieves quantified Affymetrix data in multi-graph (MGR) format by
array }
\usage{
getSampleData(celsiusObject, id = c() );
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ character vector of IDs, database accession identifiers for
    CEL records used in Celsius }
}
\value{
  Return value is a table with rows labels as SN identifiers and column
  labels as probeset identifiers. }
\author{ Jun Dong, Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  verbose = TRUE );

#get the data
data = getSampleData( celsius, id = c( "SN:1005595" ) );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getAssayData} %__CUT__%
\alias{getAssayData}
\title{ Get Celsius data using IDs, probesets OR the intersection of
  data that corresponds to BOTH IDs and probesets.}
\description{
  This function is to intelligently get the data from the DB so that the
  data is retrieved in the fastest possible way.  IF only IDs are
  given, then it will only return the full data on THOSE, if only an
  feature list is given, then it will return ALL the data for that
  platform on only those elements, and finally, if both a ID and an
  feature list are given then it will decide whether to retrieve on
  IDs and restrict by elements, or vice versa (just depending on what
  is faster).
}
\usage{
getAssayData( celsiusObject, id = NULL, feature = NULL );
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ character vector of IDs, database accession identifiers for
    CEL records used in Celsius }
  \item{feature}{ character vector of feature IDs, database accession
    identifiers for probeset records used in Celsius }
}
\value{
  Return value is a table with col labels as SN identifiers and row
  labels as probeset identifiers. }
\author{ Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu" );

#ways you can get the data
data = getAssayData( celsius,
  id = c( "SN:1018351" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at",
  "203400_s_at" ) );

data = getAssayData( celsius,
  id = c( "SN:1018351", "SN:1043233" ),
  feature = c( "201820_at" ) );

data = getAssayData( celsius,
  id = c( "SN:1018351", "SN:1043233" ) );

data = getAssayData( celsius,
  feature = c( "201820_at", "202222_s_at" ) );

data = getAssayData( celsius,
  id = c( "SN:1018351", "SN:1043233" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at",
  "203400_s_at" ) );

data = getAssayData( celsius,
  id = c( "SN:1018351", "SN:1043233", "SN:1018331", "SN:1041970",
  "SN:1008042" ),
  feature = c( "201820_at","202222_s_at" ));

data = getAssayData( celsius,
  id = c( "SN:1018351" ),
  feature = c( "201820_at" ) );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getPhenoData} %__CUT__%
\alias{getPhenoData}
\title{ Retrieve a matrix of annotations from controlled vocabulary
  terms }
\description{
  The following is a function to just retrieve a list of annotation
  terms ON a set list of IDs, or for all IDS (if none are provided).
  All Annots retrieved are specific to one platform.}
\usage{
getPhenoData(celsiusObject, id = NULL, cvterm = NULL);
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ Provide a vector of ids  }
  \item{cvterm}{ Provide a list of cvterms }
}
\value{
  Return value is a table with col labels as SN identifiers and row
  labels as annotation terms.
}
\author{ Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu" );

#ways you can get the data
data = getPhenoData( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039", "SN:1004962",
  "SN:1005722", "SN:1006573" ),
  cvterm = c( "MA:0000168", "MPATH:218" ) ); 

celsius@transform = TRUE;
data = getPhenoData( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039", "SN:1004962",
  "SN:1005722", "SN:1006573" ),
  cvterm = c( "MA:0000168", "MPATH:218" ) );

data = getPhenoData( celsius, cvterm = c( "MA:0000168", "MPATH:218" ) );

celsius@transform = FALSE;
data = getPhenoData( celsius,
  id = c( "SN:1000019", "SN:1000170", "SN:1005039", "SN:1004962",
  "SN:1005722", "SN:1006573" ),
  cvterm = c( "MPATH:218" ) ); 

data = getPhenoData( celsius,
  id = c( "SN:1000019" ),
  cvterm = c( "MA:0000168", "MPATH:218" ) ); 

data = getPhenoData( celsius,
  id = c( "SN:1000019" ), cvterm = c( "MPATH:218" ) );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{makeExprSet} %__CUT__%
\alias{makeExprSet}
\title{ Retrieve data and or phenotype annotations from Celsius and
  package into a BIOC expression set object }
\description{
  The following function retrieves a matrix of annotation data and also
  a matrix of expression data and then packages them into a BIOC
  expression object for the user.  The function requires a list of
  annotation terms along with either a list of IDs OR a list of elements
  (or both). If a list of elements is also provided, then the function
  will return the slice of data that corresponds to the intersection of
  those elements on the IDs listed.
}
\usage{
makeExprSet( celsiusObject, id = NULL, cvterm = NULL, feature = NULL,
  annotations = FALSE, mapids = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ Provide a list of SN ID's }
  \item{cvterm}{ Provide a list of cvterms }
  \item{feature}{ Provide a vector of elements } 
  \item{annotations}{ if TRUE, this includes supplementary annotations
    as terms }
  \item{mapids}{ if TRUE, this will map any foreign keys for the user
    and allow access to alternate data }
}
\value{
  Return value is a BIOC expression set object.
}
\author{ Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  platform = "HG-U133A" );

#get the data
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595"),
  cvterm = c("MPATH:458") );

#get the data
exprSet = makeExprSet( celsius, id = c( "SN:1005595" ) );

#get the data
exprSet = makeExprSet( celsius, id = c( "SN:1005595" ),
  feature = c("201820_at", "202222_s_at") );

#change the transform option to TRUE
celsius@transform = TRUE;

#get data (with transform now = TRUE)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ) );

#get the data (transform still = TRUE)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ) );

#change the transform option to FALSE
celsius@transform = FALSE;

#get the data (now with transform = FALSE)
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "222152_at", "201820_at", "202222_s_at", "202508_s_at",
  "203296_s_at", "203400_s_at" ) );

#get the data
exprSet = makeExprSet( celsius, feature = c( "201820_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  feature = c( "201820_at", "202222_s_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415" ),
  feature = c( "201820_at" ) );

#change the affy platform being sought
celsius@platform = "Mouse430_2";

#get the data (now for Mouse 430 2.0 arrays)
exprSet = makeExprSet( celsius,
  id = c( "GSE6210" ),
  cvterm = c( "MPATH:458", "MA:0000415" ),
  mapids = TRUE); 

#change the platform back to HG-U133A
celsius@platform = "HG-U133A";

#now get the data
#(mapids = TRUE, so now you can use foreign keys as sample ID's )
exprSet = makeExprSet( celsius,
  id = c( "26490", "26491", "26492" ),
  cvterm = c( "MPATH:458", "MA:0000415" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at",
  "203400_s_at" ),
  mapids = TRUE );

#get the data (annotations = TRUE, so you will get some extra data about
#each sample)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ),
  annotations = TRUE );

}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{mapIdToCelsius} %__CUT__%
\alias{mapIdToCelsius}
\title{ Expands vector of SN or external identifiers to one of SN
  identifiers }
\description{ Useful to convert, e.g. GSM or GSE IDs to coresponding SN
  identifier(s). }
\usage{
 mapTable = mapIdToCelsius( celsiusObject, id = c());
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ Provide a list of SN or external identifiers }
}
\value{
  Returns a two-column table.  Column 1 contains the "query" identifiers
  passed in as argument 1, Column 2 contains the "target" identifiers
  that correspond to the "query" identifiers.  Each row pair of
  identifiers is unique.
}
\author{ Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu" );

#define IDs you which to search on
id = c( "AEX:E-MEXP-402", "GSE6210" );

#map the IDs to the celsius database
mapTable = mapIdToCelsius( celsius, id = id );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getCvterms} %__CUT__%
\alias{getCvterms}
\title{ Retrieve ontology term (aka CV term) data structures (lists) by
  query phrase }
\description{ Useful for searching ontologies by ontology term name}
\usage{
  lungTerms = getCvterms(celsiusObject, q="lung", db="MA");
}
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with appropriate
    URLs }
  \item{q}{ Query phrase.  Meta-character "*" may be used as a wildcard }
  \item{db}{ DBspace of CV to be queried, e.g. "GO" for Gene Ontology,
    "MA" for Mouse Adult Anatomy Ontology }
}
\value{
  Returns a vector of lists.  Each list corresponds to an ontology term,
  with attributes:
  * accession - ontology term accession, e.g. "MA:0000415"
  * name - ontology term name, e.g. "lung"
  * db - db of ontology term, e.g. "MA"
  * ontology - name of ontology, e.g. "Mouse\_anatomy\_by\_time\_xproduct"
}
\author{ Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl="http://celsius.genomics.ctrl.ucla.edu" );

#get annotation term information
lungTerms = getCvterms( celsius, q = "lung", db = "MA" );
ovuTerms = getCvterms( celsius, q = "ovu*" );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{getPlatformData} %__CUT__%
\alias{getPlatformData}
\title{ lists all the platforms available on a Celsius server }
\description{ Useful to know which affy platforms are available.
  Presents a table with 2 cols: 1) platform name, and 2) number of cel
  files available}
\usage{
 platformTable = getPlatformData( celsiusObject );
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Providea celsius server object with approriate
    URLs }
}
\value{
  Returns a two-column table.  Column 1 contains the supported
  platforms, Column 2 contains the number of samples representing said platform.
}
\author{ Marc Carlson }
\examples{
#make a celsius instance
celsius = new("celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu" );

#get platform availability information
platformTable = getPlatformData( celsius );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{listSamplesForPlatform} %__CUT__%
\alias{listSamplesForPlatform}
\title{ Lists all the internal IDs for a given platform }
\description{ Useful for enumerating all the IDs present on a platform. }
\usage{
 idList = listSamplesForPlatform( celsiusObject);
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
}
\value{
  Returns a list of IDs that match up with a particular platform.
}
\author{ Marc Carlson }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  platform = "HG-U133_Plus_2" );

#get ids for the platform
idList = listSamplesForPlatform( celsius );
}
\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\name{.getWebData} %__CUT__%
%\alias{.getWebData}
%\title{ a private Function to retrieve data from the web server }
%\description{ a private Function }
%\usage{
% thing = .getWebData( celsiusObject, urlstring, useScan );
%}
%\arguments{
%  \item{celsiusObject}{ Provide a celsius server object with approriate
%    URLs }
%  \item{urlstring}{ Provide the other information to compose the rest of
%    the URL string }
%  \item{useScan}{ indicate whether scan needs to be run (only ever
%    planned for gathering ID data) }
%}
%\author{ Marc Carlson }
%%\examples{}
%\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\name{.getAsTable} %__CUT__%
%\alias{.getAsTable}
%\title{ a private Function to retrieve data from the web server using
%  read.table() }
%\description{ a private Function }
%\usage{
% thing = .getWebData( celsiusObject, urlstring);
%}
%\arguments{
%  \item{celsiusObject}{ Provide a celsius server object with approriate
%    URLs }
%  \item{urlstring}{ Provide the other information to compose the rest of
%    the URL string }
%}
%\author{ Marc Carlson }
%%\examples{}
%\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\name{.getWithScan} %__CUT__%
%\alias{.getWithScan}
%\title{ a private Function to retrieve data from the web server using
%  scan() this helper function is specialized for mapping IDs and is
%  not intended for general usage at this time }
%\description{ a private Function }
%\usage{
% thing = .getWithScan( celsiusObject, urlstring);
%}
%\arguments{
%  \item{celsiusObject}{ Provide a celsius server object with approriate
%    URLs }
%  \item{urlstring}{ Provide the other information to compose the rest of
%    the URL string }
%}
%\author{ Marc Carlson }
%%\examples{}
%\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\name{.getAllAnnots} %__CUT__%
%\alias{.getAllAnnots}
%\title{ a private Function to retrieve all annotation data from the web
%  server and map this to the requested IDs }
%\description{ a private Function }
%\usage{
% thing = .getAllAnnots( celsiusObject, id);
%}
%\arguments{
%  \item{celsiusObject}{ Provide a celsius server object with approriate
%    URLs }
%  \item{id}{ Provide a list of ids to get all annotation data for }
%}
%\author{ Marc Carlson }
%%\examples{}
%\keyword{ data }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$Log: man.Rd,v $
%Revision 1.45  2007/05/24 18:04:32  mcarlson
%Corrected build errors in documentation.
%
%Revision 1.44  2007/05/23 01:48:24  mcarlson
%edited out synopsis from the man.Rd file.
%
%Revision 1.43  2007/05/22 01:59:22  allenday
%syntax error
%
%Revision 1.42  2007/05/04 19:48:49  mcarlson
%updated the man.Rd file to format and make the help more consistent and...
%helpful.
%
%Revision 1.41  2007/05/04 18:12:22  mcarlson
%killed some bugs relating to annotations and possibly the refactor of
%makeExprSet().  The function is MUCH easier to maintain now.  And
%everything appears to work well.  Need guys to do testing.  Also, need
%to finish cleanup of the manual pages...
%
%Revision 1.40  2007/05/04 02:26:08  mcarlson
%Critical bug fix for annotations parameter, and some minor refactoring on
%the huge makeExprSet().  Which I still think needs to shrink a lot.
%
%Revision 1.39  2007/05/03 21:56:11  mcarlson
%Bug fix for the ID mapper to accomodate STRANGE IDs that begin with "0".
%Also some refactoring to reduce redundancy in the code.  makeExprSet()
%still needs to be refactored...
%
%Revision 1.38  2007/05/02 23:51:30  mcarlson
%Added code to retry when data is being grabbed (needs more testing).  And
%also refactored the settings out of the functions and into the celsius
%object.  THIS WILL MEAN THAT WE HAVE TO CHANGE THE ONLINE TUTORIAL.  So
%please exercise the necessart precautions.
%
%Revision 1.37  2007/04/28 02:35:50  mcarlson
%Improved the error handling of the function to factor out function that
%gets data from the web server.  Fixed a broken example.
%
%Revision 1.36  2007/04/28 01:22:38  mcarlson
%Added preliminary code to retry a get (if it comes back empty from the
%web server...).  But before I continue, I need to verify that the new func
%will be documented ok in the man file.
%
%Revision 1.35  2007/04/28 00:08:02  mcarlson
%Finally the id mapping function will have functional error handling.  A
%parameter has also been added to the highest level function so that ID
%mapping will be toggled off by default.  This should save time for people
%who are using default SNIDs.
%
%Revision 1.34  2007/04/26 01:45:29  allenday
%version bump!
%
%Revision 1.33  2007/04/26 01:26:38  mcarlson
%changed progress to verbose
%
%Revision 1.32  2007/04/26 01:20:17  mcarlson
%Added error handling code for the other two "get" functions and added a verbose
%option that now goes all the way to the top.
%
%Revision 1.31  2007/04/25 01:44:35  mcarlson
%Added Jun code for error handling for getSampleData().  Made minor edits to
%the documentation.  Still needs to have documentation improved...
%
%Revision 1.30  2007/04/21 02:02:45  mcarlson
%
%Added edit to the man.Rd file to try and stamp out error.
%
%Revision 1.29  2007/04/21 01:55:20  mcarlson
%Added param to makeExprSet() so that it can output all annot terms if requested.
%
%Revision 1.28  2007/04/20 23:26:05  mcarlson
%Added another clever example to the top set of functions.
%
%Revision 1.27  2007/04/20 23:09:39  mcarlson
%Fixed some small bugs. Changed examples to match Allens web tutorial.
%
%Revision 1.26  2007/04/20 01:26:19  allenday
%1.1 check passes
%
%Revision 1.25  2007/04/20 01:11:34  mcarlson
%
%Updated man.Rd
%
%Revision 1.24  2007/04/19 22:42:23  allenday
%docs and code work
%
%Revision 1.23  2007/04/19 02:22:31  allenday
%getCvterms function
%
%Revision 1.22  2007/04/19 01:20:00  mcarlson
%Added functions to get the available platforms and the IDs for a particular
%platform along with relevant documentation.
%
%Revision 1.21  2007/04/13 01:37:31  mcarlson
%Just putting more names in...
%
%Revision 1.20  2007/04/13 01:26:09  allenday
%doc fixes, make target for distributable file
%
%Revision 1.19  2007/04/13 00:46:05  mcarlson
%more edits to try to get the class to read in...
%
%Revision 1.18  2007/04/13 00:19:21  mcarlson
%Made changes to the code and the man.Rd file to try to get this thing to be
%compliant with documentaion standards.
%
%Revision 1.17  2007/04/12 23:05:45  mcarlson
%Cleaned up code, added prelim manual info for server object.  Fixes in to
%accomodate server URL changes and parameter name changes.
%
%Revision 1.16  2007/04/12 18:56:29  mcarlson
%Updated usage and synopsis for man.Rd to accomodate the new celsius object.
%
%Revision 1.15  2007/04/11 23:48:07  mcarlson
%Added code to make a small R object with get and set functions to allow
%access to the URL Strings contained within it.  Functions now use the
%contents of this object to get the data that they need.  This will allow
%separate instanced of the celsius DB to be instantiated and accessed later
%on.  Also, this object defaults to the "original" celsius DB.
%
%Revision 1.14  2007/03/30 00:31:25  mcarlson
%Dropped extraneous fucntion.  Refactored the 3 base functions for speed and
%efficiency.
%
%Revision 1.13  2007/03/28 20:59:39  mcarlson
%More name changes.
%
%Revision 1.12  2007/03/28 20:38:59  mcarlson
%corrected naming error.
%
%Revision 1.11  2007/03/28 19:25:59  mcarlson
%Changed names of variables for consistency and to try and increase clarity
%by making names more similar to the existing BIOC names.
%
%Revision 1.10  2007/03/28 02:43:43  mcarlson
%I think that most of the evil bugs are finally dead...
%All tests pass now.  But the code could still be tidied (renamed etc).
%
%Revision 1.9  2007/03/27 00:50:23  mcarlson
%Cleaned up .Rd file to lose a lot of warnings. Removed superfluous annot.manager()
%
%Revision 1.8  2007/03/24 02:02:52  mcarlson
%Now you can use IDs from alternate sources, and the function to make a
%BIOC object will just figure out what you mean, and then place the ID mapping
%into the pheno object along with the annotations for any CVTERMS that you
%request...  Next: this code needs to be cleaned up.  And good examples
%need to be made.
%
%Revision 1.7  2007/03/23 03:07:20  mcarlson
%Many many bug fixes.  Most notably, you no longer are required to have a
%cvterm to make a BIOC object, and a lot of coner cases where cvterm, ID,
%or feature == 1 no longer fail due to inately bad "R bahavoirs"...  Lets
%just say that R sucks in the assumptions department and that this code now
%works.
%
%Revision 1.6  2007/03/23 01:10:09  allenday
%id mapping code and docs
%
%Revision 1.5  2007/03/23 00:52:41  allenday
%added RCS Id and Log
%
