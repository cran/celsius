\name{makeExprSet} %__CUT__%
\alias{makeExprSet}
\title{ Retrieve data and or phenotype annotations from Celsius and
  package into a BIOC expression set object }
\description{
  The following function retrieves a matrix of annotation data and also
  a matrix of expression data and then packages them into a BIOC
  expression object for the user.  The function requires a list of
  annotation terms along with either a list of IDs OR a list of elements
  (or both). If a list of elements is also provided, then the function
  will return the slice of data that corresponds to the intersection of
  those elements on the IDs listed.
}
\usage{
makeExprSet( celsiusObject, id = NULL, cvterm = NULL, feature = NULL,
  annotations = FALSE, mapids = FALSE)
}
\arguments{
  \item{celsiusObject}{ Provide a celsius server object with approriate
    URLs }
  \item{id}{ Provide a list of SN ID's }
  \item{cvterm}{ Provide a list of cvterms }
  \item{feature}{ Provide a vector of elements } 
  \item{annotations}{ if TRUE, this includes supplementary annotations
    as terms }
  \item{mapids}{ if TRUE, this will map any foreign keys for the user
    and allow access to alternate data }
}
\value{
  Return value is a BIOC expression set object.
}
\author{ Marc Carlson and Allen Day }
\examples{
#make a celsius instance
celsius = new( "celsiusServer",
  celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
  platform = "HG-U133A" );

#get the data
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595"),
  cvterm = c("MPATH:458") );

#get the data
exprSet = makeExprSet( celsius, id = c( "SN:1005595" ) );

#get the data
exprSet = makeExprSet( celsius, id = c( "SN:1005595" ),
  feature = c("201820_at", "202222_s_at") );

#change the transform option to TRUE
celsius@transform = TRUE;

#get data (with transform now = TRUE)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ) );

#get the data (transform still = TRUE)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ) );

#change the transform option to FALSE
celsius@transform = FALSE;

#get the data (now with transform = FALSE)
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "222152_at", "201820_at", "202222_s_at", "202508_s_at",
  "203296_s_at", "203400_s_at" ) );

#get the data
exprSet = makeExprSet( celsius, feature = c( "201820_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  feature = c( "201820_at", "202222_s_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ) );

#get the data
exprSet = makeExprSet( celsius,
  cvterm = c( "MA:0000415" ),
  feature = c( "201820_at" ) );

#change the affy platform being sought
celsius@platform = "Mouse430_2";

#get the data (now for Mouse 430 2.0 arrays)
exprSet = makeExprSet( celsius,
  id = c( "GSE6210" ),
  cvterm = c( "MPATH:458", "MA:0000415" ),
  mapids = TRUE); 

#change the platform back to HG-U133A
celsius@platform = "HG-U133A";

#now get the data
#(mapids = TRUE, so now you can use foreign keys as sample ID's )
exprSet = makeExprSet( celsius,
  id = c( "26490", "26491", "26492" ),
  cvterm = c( "MPATH:458", "MA:0000415" ),
  feature = c( "201820_at", "202222_s_at", "202508_s_at", "203296_s_at",
  "203400_s_at" ),
  mapids = TRUE );

#get the data (annotations = TRUE, so you will get some extra data about
#each sample)
exprSet = makeExprSet( celsius,
  id = c( "SN:1005595", "SN:1000170" ),
  cvterm = c( "MA:0000415", "MPATH:458" ),
  feature = c( "201820_at", "202222_s_at" ),
  annotations = TRUE );

}
\keyword{ data }
