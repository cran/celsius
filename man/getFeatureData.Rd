\name{getFeatureData} %__CUT__%
\alias{getFeatureData}
\title{ Retrieve data (probeset-wise) into R from Celsius webserver }
\description{
Retrieves quantified Affymetrix data in multi-graph (MGR) format by
feature }
\usage{
getFeatureData(celsiusObject, feature = c() );
}
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
