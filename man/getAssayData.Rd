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
