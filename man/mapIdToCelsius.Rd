\name{mapIdToCelsius} %__CUT__%
\alias{mapIdToCelsius}
\title{ Expands vector of SN or external identifiers to one of SN
  identifiers }
\description{ Useful to convert, e.g. GSM or GSE IDs to coresponding SN
  identifier(s). }
\usage{
 mapTable = mapIdToCelsius( celsiusObject, id = c());
}
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
