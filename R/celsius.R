#$Id: celsius.R,v 1.40 2007/05/04 18:12:23 mcarlson Exp $
#Copyright 2006 Allen Day, Marc Carlson.  Released under Artistic License
library("Biobase");
setClass("celsiusServer", representation(
	celsiusUrl = "character", 
	platform = "character",
	protocol = "character",
	stringency = "numeric",
	verbose = "logical",
	transform = "logical",
	retry = "numeric",	
	sampleTime = "numeric", 
	featureTime = "numeric"),
	prototype = list(
	celsiusUrl = "http://celsius.genomics.ctrl.ucla.edu",
	platform = "HG-U133A",
	protocol = "rma",
	stringency = 500,
	verbose = FALSE,
	transform = FALSE, 
	retry = 3, 
	sampleTime = 1.8, 
	featureTime = 2.6 ))
#
#####HOWTO: make get() and set()
#setGeneric("getcelsiusUrl", function(object) {
# standardGeneric("getcelsiusUrl")
#})
#setMethod("getcelsiusUrl", "celsiusServer", function(object) {
#  object@celsiusUrl
#})
#setGeneric("setcelsiusUrl", function(object, String) {
# standardGeneric("setcelsiusUrl")
#})
#setMethod("setcelsiusUrl", "celsiusServer", function(object, String) {
#  eval(eval(substitute(expression(object@celsiusUrl<<-String))))
#})

.getWebData = function(
  celsiusObject,
  urlString,
  useScan = FALSE
){
    if(celsiusObject@retry > 0){ counter = 1;}else{counter = 0;}
    thing = NULL;
    while( (counter<=celsiusObject@retry || counter==0) && is.null(thing) ){
      if(useScan == TRUE){thing = tryCatch(.getWithScan(celsiusObject,urlString), error=function(e) NULL);}
	  else{thing = tryCatch(.getAsTable(celsiusObject,urlString), error=function(e) NULL)};
	  if(is.null(thing)){cat("Tried and Failed" ,counter,"times to retrieve data...\n")}
	  if(counter!=0){counter = counter + 1;}
    }
    return(thing);      
}

.getAsTable = function(
  celsiusObject,
  urlString
){
    return(read.table(paste(celsiusObject@celsiusUrl, urlString, sep=""), sep="\t", comment.char="!"));
}

.getWithScan = function( 
  celsiusObject,
  urlString
){
    table = NULL;
	linescan = scan(paste(celsiusObject@celsiusUrl, urlString, sep=""),sep="\t",what="character");
    for ( i in seq( 1,length( linescan ),by=2 ) ) {
		table = rbind(table,linescan[i:(i+1)]);
    }
    return(table);
}

getOntologyData = function(
  celsiusObject,
  cvterm = c()
){
    cvterm = as.matrix(cvterm)[, 1];
    mat = NULL;
    j = 1;
    check = rep(TRUE, length(cvterm));
    sampleNames = NULL;
    sampleNamesDone = FALSE;
    for (s in cvterm) {
        cvData = .getWebData(celsiusObject, paste("/service/cel/annots?p=",celsiusObject@platform,";f=mask;m=",s,";l=",celsiusObject@stringency, sep=""));
        if (celsiusObject@verbose) {
            cat(j, " term ID ", s, " was retrieved for ", dim(cvData)[1], " samples", "\n");
        }
        if(is.null(cvData)){
            check[j]= FALSE;
        } else if (dim(cvData)[2] == 2 && dim(cvData)[1] > 1) {
            if(!sampleNamesDone){
                sampleNames=cvData[-1, 1];
                sampleNamesDone=TRUE;
            }
            cvData = as.numeric(as.matrix(cvData[-1, 2]))
            if(celsiusObject@transform == TRUE){
				cvData[cvData=="!"]=NA;
				cvData[is.na(cvData)]=0;
			}
        } else {
                cvData = NULL;
                check[j] = FALSE;
        }
        mat = cbind(mat, cvData);
        j = j + 1;
    }
    if(any(check)){
      colnames(mat) = cvterm[check];
      rownames(mat) = sampleNames;
    }
    badCVTerm = sum(!check);
    if (badCVTerm > 0) {
        if(badCVTerm == length(cvterm) ){
          cat("\nThere is no data retrieved from the database! Please check your list of cvterms!");
        } 
        else {
            cat("\n", badCVTerm, "Some CVTerms were discarded from your list as they were not in the database.");
            cat("\nAnd they are:\n", cvterm[!check], "\n");
        }
    }
    return(t(mat));
}

getFeatureData = function(
  celsiusObject,
  feature = c()
){
    feature = as.matrix(feature)[, 1];
    mat = NULL;
    j = 1;
    check = rep(TRUE, length(feature));
    sampleNames = NULL;
    sampleNamesDone = FALSE;
    for (s in feature) {
		featureData = .getWebData(celsiusObject, paste("/service/element/element?p=",celsiusObject@platform,";id=",s,";p=",celsiusObject@protocol, sep = ""));            
        if (celsiusObject@verbose) {
            cat(j, " feature ID ", s, " was retrieved for " , dim(featureData)[1], " samples", "\n");
        }
        if(is.null(featureData)){
             check[j]= FALSE;
        }else if (dim(featureData)[2] == 2 && dim(featureData)[1] > 1) {
            if(!sampleNamesDone){
                sampleNames=featureData[-1, 1];
                sampleNamesDone=TRUE;
            }
            featureData = as.numeric(as.matrix(featureData[-1, 2]))
        } else {
                featureData = NULL;
                check[j] = FALSE;
        }
        mat = cbind(mat, featureData);
        j = j + 1;
    }
    if(any(check)){
      colnames(mat) = feature[check];
      rownames(mat) = sampleNames;
    }
    badFeature = sum(!check);
    if (badFeature > 0 ){
       if(badFeature == length(feature) ){
         cat("There is no data retrieved from the database! Please check your list of ProbeSets!");
       } else {
         cat("\n", badFeature, "Probe sets were discarded from your list as they were not in the database for this platform.");
         cat("\nAnd they are:\n", feature[!check], "\n");
       }
    }
    return(t(mat));
}

getSampleData = function (
  celsiusObject,
  id = c()
){
   id = as.matrix(id)[, 1];
   mat = NULL;
   j = 1;
   check = rep(TRUE, length(id));
   probesetNames = NULL;
   probesetNamesDone = FALSE;
   for (s in id) {
	   egr = .getWebData(celsiusObject, paste("/service/cel/cel?id=",s, ";f=egr;p=", celsiusObject@protocol, sep = ""));   
       if (celsiusObject@verbose) {
           cat(j, " sample ID ", s, " was retrieved for ", dim(egr)[1], " features", "\n")
       }
       if(is.null(egr)){
           check[j]= FALSE;
       } else if (dim(egr)[2] == 2 && dim(egr)[1] > 1) {
           if (!probesetNamesDone) {
               probesetNames = egr[-1, 1];
               probesetNames = gsub( pattern = '.+?:', replacement = "", x = probesetNames);
               probesetNamesDone = TRUE;
           }
           egr = as.numeric(as.matrix(egr[-1, 2]));
       } else {
           egr = NULL;
           check[j] = FALSE;
       }
       mat = cbind(mat, egr);
       j = j + 1;
   }
   if(any(check)){
       colnames(mat) = id[check];
       rownames(mat) = probesetNames;
   }
   badID = sum(!check);
   if (badID > 0) {
       if (badID == length(id)) {
           cat("\nThere is no data retrieved from the database! Please check your list of IDs!\n\n");
       }
       else {
           cat("\n", badID, " IDs were discarded from your list as they were either not in the database or did not contain any probesets.", sep="");
           cat("\nAnd they are:\n", id[!check], "\n\n");
       }
   }
   return(mat)
}

getAssayData = function(
  celsiusObject,
  id = NULL,
  feature = NULL
){
  if( is.null(feature) && !is.null(id) ){
    exprMatrix = getSampleData(celsiusObject,id=id);
  } 
  if( is.null(id) && !is.null(feature) ){
    exprMatrix = getFeatureData(celsiusObject,feature=feature);
  }
  if( !is.null(id) && !is.null(feature)){
    numCols = length(id);
    numRows = length(feature);
    if(numCols *celsiusObject@sampleTime >= numRows * celsiusObject@featureTime) {
      exprMatrix = getFeatureData(celsiusObject,feature=feature);
      matchVec = match(colnames(exprMatrix), id);
      exprMatrix = as.matrix(exprMatrix[,!is.na(matchVec),drop=FALSE]);
    }
    else if(numRows * celsiusObject@featureTime > numCols *celsiusObject@sampleTime){
      exprMatrix = getSampleData(celsiusObject,id=id);
      matchVec = match(rownames(exprMatrix), feature);
      exprMatrix = as.matrix(exprMatrix[!is.na(matchVec),,drop=FALSE]);
    }
    rownames(exprMatrix) = feature;
    colnames(exprMatrix) = id;    
  }
  return(exprMatrix);
}

getPhenoData = function(
  celsiusObject,  
  id = NULL,
  cvterm = NULL
){
  annotMat = getOntologyData(celsiusObject,cvterm=cvterm);
  if( !is.null(id) ){
    matchVec = match(colnames(annotMat), id);
    annotMat = annotMat[,!is.na(matchVec),drop=FALSE];
  }
  return(annotMat)
}

.getAllAnnots = function(
  celsiusObject,
  id = NULL
){
	if(is.null(id)){id = as.character(listSamplesForPlatform(celsiusObject))}
    annotTable = read.delim(paste(celsiusObject@celsiusUrl,'/service/cel/annots?p=',celsiusObject@platform,';f=tab', sep=''));
    rownames(annotTable) = sub(':','.',annotTable[,1]);
    idSymbols = sub(':','.',id);
    matchVec = match(rownames(annotTable), idSymbols);
    annotTable = annotTable[!is.na(matchVec),];
    return(annotTable);
}

makeExprSet = function(
  celsiusObject,
  id = NULL,
  cvterm = NULL,
  feature = NULL,
  annotations = FALSE,
  mapids = FALSE
){
  annotTable = NULL;
  phenoData  = NULL;
  exprData = NULL;

  platformName = gsub(pattern='[^0-9A-Za-z]+', replacement='', x=celsiusObject@platform, perl=TRUE);
  platformName = sub(pattern="(\\w+)", replacement="\\L\\1", x=platformName, perl=TRUE);

  if(!is.null(id)){
    if(mapids == TRUE){idMap = mapIdToCelsius(celsiusObject,id);}else{idMap = cbind(id,id);}
    id = as.vector(idMap[,2]);
    idPheno = as.matrix(idMap[,1]);
    rownames(idPheno) = idMap[,2];
    colnames(idPheno) = "original_ID";
  }

  exprData = getAssayData(celsiusObject,id = id, feature = feature);
  rownames(exprData) = gsub( pattern = '.+?:', replacement = "", x = rownames(exprData));
  colnames(exprData) = sub(':','.',colnames(exprData));

  if(!is.null(cvterm)){
    phenoData = getPhenoData(celsiusObject,id = id, cvterm = cvterm)
    colnames(phenoData) = sub(':','.',colnames(phenoData)); 
    phenoData = t(phenoData);
    matchVec = match(rownames(phenoData), colnames(exprData));
    phenoData = phenoData[!is.na(matchVec),,drop=FALSE];
  }

  if(!is.null(id)){phenoData = cbind(phenoData,idPheno);}

  if(annotations == TRUE){
    annotTable = .getAllAnnots(celsiusObject,id);
    matchVec = match(rownames(annotTable), colnames(exprData));
    annotTable = annotTable[!is.na(matchVec),,drop=FALSE];
#    phenoData = cbind(phenoData,annotTable);
    if(is.null(id) && is.null(cvterm)){phenoData = annotTable}else{phenoData = cbind(phenoData,annotTable);}
  }

  if(is.null(id) && is.null(cvterm) && annotations == FALSE){  
      exprSet = new('ExpressionSet', exprs = exprData, annotation = platformName);
  }else{
    phenoData = data.frame(phenoData, row.names = colnames(exprData));
    phenoData = new("AnnotatedDataFrame", data = phenoData);
    exprSet = new('ExpressionSet', phenoData = phenoData, exprs = exprData, annotation = platformName);
  }
  return(exprSet);
}

mapIdToCelsius = function(
  celsiusObject,
  id = c()
){
  id = as.character(id);  #checkme: maybe not needed
  by = 50;
  buf = NULL;
  check = rep(TRUE, length(id));
  j = 1;
  for ( i in seq( 1,length( id ),by=by ) ) {
    if( i + by > length( id ) ) {
      t = .getWebData(celsiusObject, paste('/service/search/cel_accession?f=tab;q=',paste(id[i:length(id)],collapse='+'),sep=""),useScan=TRUE);
    }
    else {
	  t = .getWebData(celsiusObject, paste('/service/search/cel_accession?f=tab;q=',paste(id[i:((i+by)-1)],collapse='+'),sep=""),useScan=TRUE);
    }
    if(is.null(t)){
      check[j]= FALSE;
    }else {
      buf = rbind(buf,t);
    }
    j = j + 1;
  }
  badID = sum(!check);
   if (badID == length(id)) {
     cat("\nThere are no IDs in your list that map to the database!\n\n");
   }
  uniqueIDs =  as.vector(unique(id));
  colnames(buf) = c('query','target');
  uniqueRetIDs = as.vector(unique(buf[,1]));
  if(length(uniqueRetIDs) != length(uniqueIDs) ){
	matchVec = match(uniqueIDs, uniqueRetIDs);  
	cat("\nThe following IDs could not be found:", id[is.na(matchVec)],"\n\n");
  }
  return(buf);
}

getCvterms = function(
  celsiusObject,
  q = NULL,
  db = NULL
){
  z = read.table(paste(celsiusObject@celsiusUrl, '/service/search/ontology?q=',q,sep=""),sep="\t",header=FALSE,col.names=c('db','accession','name','ontology'));
  r = vector(length=1,mode="list");
  j = 1;
  for ( i in 1:dim(z)[1] ) {
    if ( is.null(db) || db == z$db[i] ) {
      k = as.list(x=new.env());
      k$db = as.character(z$db[i]);
      k$accession = as.character(sub('/',':',z$accession[i]));
      k$name = as.character(z$name[i]);
      k$ontology = as.character(z$ontology[i]);
      r[[j]] = k;
      j = j + 1;
    }
  }
  return(r);
}

getPlatformData = function(
  celsiusObject
){
  platTable = read.delim(paste(celsiusObject@celsiusUrl, '/service/platform/platforms?f=tab', sep=''),header=FALSE,col.names=c('platform','samples','field'));
  platTable = platTable[,1:2];
  return(platTable);
}

listSamplesForPlatform = function(
  celsiusObject 
){
  platTable = read.delim(paste(celsiusObject@celsiusUrl,'/service/cel/cels?p=',celsiusObject@platform,';f=tab', sep=''),header=FALSE);
  platTable = platTable[,1];
  return(platTable);
}

#$Log: celsius.R,v $
#Revision 1.40  2007/05/04 18:12:23  mcarlson
#killed some bugs relating to annotations and possibly the refactor of
#makeExprSet().  The function is MUCH easier to maintain now.  And
#everything appears to work well.  Need guys to do testing.  Also, need
#to finish cleanup of the manual pages...
#
#Revision 1.39  2007/05/04 06:23:14  mcarlson
#FINALLY! Refactored that bloated ugly makeExprsSet().
#
#Revision 1.38  2007/05/04 02:26:09  mcarlson
#Critical bug fix for annotations parameter, and some minor refactoring on
#the huge makeExprSet().  Which I still think needs to shrink a lot.
#
#Revision 1.37  2007/05/03 21:56:12  mcarlson
#Bug fix for the ID mapper to accomodate STRANGE IDs that begin with "0".
#Also some refactoring to reduce redundancy in the code.  makeExprSet()
#still needs to be refactored...
#
#Revision 1.36  2007/05/03 02:04:34  mcarlson
#Bug Fixes.  Also eliminated many lines of unecessary checking in middle tier
#functions.  makeExprsSet still needs a major refactoring.
#
#Revision 1.35  2007/05/02 23:51:31  mcarlson
#Added code to retry when data is being grabbed (needs more testing).  And
#also refactored the settings out of the functions and into the celsius
#object.  THIS WILL MEAN THAT WE HAVE TO CHANGE THE ONLINE TUTORIAL.  So
#please exercise the necessart precautions.
#
#Revision 1.34  2007/04/28 02:35:52  mcarlson
#Improved the error handling of the function to factor out function that
#gets data from the web server.  Fixed a broken example.
#
#Revision 1.33  2007/04/28 01:22:40  mcarlson
#Added preliminary code to retry a get (if it comes back empty from the
#web server...).  But before I continue, I need to verify that the new func
#will be documented ok in the man file.
#
#Revision 1.32  2007/04/28 00:15:55  jundong
#test it
#
#Revision 1.31  2007/04/28 00:08:02  mcarlson
#Finally the id mapping function will have functional error handling.  A
#parameter has also been added to the highest level function so that ID
#mapping will be toggled off by default.  This should save time for people
#who are using default SNIDs.
#
#Revision 1.30  2007/04/26 06:54:11  allenday
#id mapping paging
#
#Revision 1.29  2007/04/26 01:26:39  mcarlson
#changed progress to verbose
#
#Revision 1.28  2007/04/26 01:20:18  mcarlson
#Added error handling code for the other two "get" functions and added a verbose
#option that now goes all the way to the top.
#
#Revision 1.27  2007/04/25 01:44:35  mcarlson
#Added Jun code for error handling for getSampleData().  Made minor edits to
#the documentation.  Still needs to have documentation improved...
#
#Revision 1.26  2007/04/21 01:55:22  mcarlson
#Added param to makeExprSet() so that it can output all annot terms if requested.
#
#Revision 1.25  2007/04/20 23:09:42  mcarlson
#Fixed some small bugs. Changed examples to match Allens web tutorial.
#
#Revision 1.24  2007/04/20 01:11:34  mcarlson
#
#Updated man.Rd
#
#Revision 1.23  2007/04/19 02:22:32  allenday
#getCvterms function
#
#Revision 1.22  2007/04/19 01:20:01  mcarlson
#Added functions to get the available platforms and the IDs for a particular
#platform along with relevant documentation.
#
#Revision 1.21  2007/04/13 01:29:48  allenday
#license
#
#Revision 1.20  2007/04/13 00:19:21  mcarlson
#Made changes to the code and the man.Rd file to try to get this thing to be
#compliant with documentaion standards.
#
#Revision 1.19  2007/04/12 23:05:45  mcarlson
#Cleaned up code, added prelim manual info for server object.  Fixes in to
#accomodate server URL changes and parameter name changes.
#
#Revision 1.18  2007/04/11 23:48:08  mcarlson
#Added code to make a small R object with get and set functions to allow
#access to the URL Strings contained within it.  Functions now use the
#contents of this object to get the data that they need.  This will allow
#separate instanced of the celsius DB to be instantiated and accessed later
#on.  Also, this object defaults to the "original" celsius DB.
#
#Revision 1.17  2007/03/30 00:31:26  mcarlson
#Dropped extraneous fucntion.  Refactored the 3 base functions for speed and
#efficiency.
#
#Revision 1.16  2007/03/28 20:59:39  mcarlson
#More name changes.
#
#Revision 1.15  2007/03/28 20:39:00  mcarlson
#corrected naming error.
#
#Revision 1.14  2007/03/28 19:25:59  mcarlson
#Changed names of variables for consistency and to try and increase clarity
#by making names more similar to the existing BIOC names.
#
#Revision 1.13  2007/03/28 02:43:46  mcarlson
#I think that most of the evil bugs are finally dead...
#All tests pass now.  But the code could still be tidied (renamed etc).
#
#Revision 1.12  2007/03/27 01:23:18  mcarlson
#Cleaned out superfluous comments in code.
#
#Revision 1.11  2007/03/27 00:50:27  mcarlson
#Cleaned up .Rd file to lose a lot of warnings. Removed superfluous annot.manager()
#
#Revision 1.10  2007/03/24 02:02:52  mcarlson
#Now you can use IDs from alternate sources, and the function to make a
#BIOC object will just figure out what you mean, and then place the ID mapping
#into the pheno object along with the annotations for any CVTERMS that you
#request...  Next: this code needs to be cleaned up.  And good examples
#need to be made.
#
#Revision 1.9  2007/03/23 03:07:20  mcarlson
#Many many bug fixes.  Most notably, you no longer are required to have a
#cvterm to make a BIOC object, and a lot of coner cases where cvterm, ID,
#or feature == 1 no longer fail due to inately bad "R bahavoirs"...  Lets
#just say that R sucks in the assumptions department and that this code now
#works.
#
#Revision 1.8  2007/03/23 01:10:09  allenday
#id mapping code and docs
#
#Revision 1.7  2007/03/23 00:51:22  allenday
#added RCS Log
#
