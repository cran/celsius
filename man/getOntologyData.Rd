\name{getOntologyData} %__CUT__%
\alias{getOntologyData}
\title{ Retrieve Annotations into R from Celsius webserver }
\description{
Retrieves annotation data for Affymetrix arrays in multi-graph (MGR) format
}
\usage{
getOntologyData(celsiusObject, cvterm=c());
}
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
