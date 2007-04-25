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
