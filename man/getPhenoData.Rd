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
