\name{listSamplesForPlatform} %__CUT__%
\alias{listSamplesForPlatform}
\title{ Lists all the internal IDs for a given platform }
\description{ Useful for enumerating all the IDs present on a platform. }
\usage{
 idList = listSamplesForPlatform( celsiusObject);
}
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
