\name{getPlatformData} %__CUT__%
\alias{getPlatformData}
\title{ lists all the platforms available on a Celsius server }
\description{ Useful to know which affy platforms are available.
  Presents a table with 2 cols: 1) platform name, and 2) number of cel
  files available}
\usage{
 platformTable = getPlatformData( celsiusObject );
}
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
