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
