\name{getSampleData} %__CUT__%
\alias{getSampleData}
\title{ Retrieve data (sample-wise) into R from Celsius webserver }
\description{
Retrieves quantified Affymetrix data in multi-graph (MGR) format by
array }
\usage{
getSampleData(celsiusObject, id = c() );
}
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
