/*============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 3.1, June 2016
* 
* Observation: This file is apparently leaking memory like crazy - SJ 030616
============================================
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <regex.h>
#include <RLexact.h>

using std::cerr;
using std::cout;
using std::ios;
using std::ifstream;
using std::endl;


/* functions defined in this file*/
double atod(char*);
long long regexperr (long long, const char*, bool);
long long multimatch(char*,long long,const char*,double**,long long);
void filereader(char*,char*,long long);
long long filesizer(char*);


/* global variables defined in RLexact.C */
extern char* infilename;
extern char* filedata;
extern long long filesize;

/*
int main () {

  int filesize = filesizer(infilename);
  filedata = (char*) malloc(filesize*sizeof(char));
  filereader(infilename,filedata,filesize);

  double **res = (double**) malloc(10*sizeof(double*));

  multimatch(filedata,filesize,"ostehaps",res,2);

  double **res2 = (double**) malloc(10*sizeof(double*));
  multimatch(filedata,filesize,"oaflygel",res2,3);

  for (int i =0;i<3;i++) {
    for (int j=0;j<5;j++) {
      cout << res2[i][j] << " ";
    }
    cout << endl;
  }

  double *res5 = (double*) malloc(MAXARRAYSIZE*sizeof(double));

  matchlines(filedata,"gruk",res5,true);

}
*/

long long filesizer (char* filename) {
  /* function for finding the size of filename. Silly! */
  ifstream datafile (filename, ios::in|ios::binary|ios::ate);
  if (datafile.fail()) {
    cerr << "Error while opening input file " << filename << endl;
    exit(-1);
  }
  datafile.seekg(0,ios::beg); // move to beginning of file
  long long begin = datafile.tellg(); // get position
  datafile.seekg(0,ios::end); // move to end of file
  long long end = datafile.tellg(); // get position
  long long filesize = end-begin; // calculate size
  datafile.close();
  return filesize; // and return
}

void filereader (char* filename, char* filedata, long long filesize) {
  /* function for reading binary file into char array */
  ifstream datafile (filename, ios::in|ios::binary|ios::ate);
  if (datafile.fail()) {
    cerr << "Error while opening input file " << filename << endl;
    exit(-1);
  }
  datafile.seekg(0,ios::beg);
  datafile.read(filedata,filesize); // here we must have filesize from caller

#ifdef FILEREAD_TEST
  // write the read content of the file to stdout
  cerr << "Read the following from datafile " << filename << endl;
  for (long long i =0;i<filesize;i++) {
    cerr << filedata[i];
  }
  cerr << endl;
#endif /* FILEREAD_TEST */
  datafile.close();
}

long long matchlines (char* input, const char* pattern, double* result, bool strict) {

  /* this function takes a pattern (usually a keyword) and searches for it in
     the input string. If found, the value following the keyword is translated
     into an array of doubles, returned in the result variable. If strict is
     set, the function will return error (-1) if the keyword is found, but no
     usable values could be filled into result - if strict _isn't_ set, it
     silently passes an empty array to the caller. The function returns the
     number of elements in the result array if succesful, or -1 if an error
     occured */

  // the two variables for holding the regular expressions for the keyword and the array-matching pattern
  regex_t* comppat;
  regex_t* comppat2;
  comppat = new regex_t;
  comppat2 = new regex_t;

  long long count=0; // the number of elements put into the result array

  regmatch_t resarray[2], resarray2[2]; // arrays used in the regexec expressions.

  const char *pattern2 = "\\(\\-\\?[0-9|\\.|\\-]\\+\\).*[0-9\\-\\.]*"; // the pattern that matches one double followed by any number of further entries

  // take the keyword, and form the appropriate regular expression from it, ie. keyword -> ^keyword\(.*\)
  long long l = 0;
  while (pattern[l]!=0) { l++;}
  char* realpat=(char*) malloc(sizeof(char)*(l+14));
  realpat[0]='^';
  for (long long m=0;m<l;m++) {
    realpat[m+1]=pattern[m];
  }
  realpat[l+1]='\\';
  realpat[l+2]='(';
  realpat[l+3]='.';
  realpat[l+4]='*';
  realpat[l+5]='\\';
  realpat[l+6]=')';
  realpat[l+7]=0;

#ifdef FILEREAD_TEST
  cerr << "Transforming pattern ";
  for (long long n=0; n< l;n++) {
    cout << pattern[n];
  }
  cout << " into ";
  for (long long n=0; n< l+7;n++) {
    cout << realpat[n];
  }
  cout << endl;
#endif /* FILEREAD_TEST */

  // compile the patterns into comppat (the keyword), and comppat2 (the array matcher)
  regcomp (comppat, realpat,REG_NEWLINE|REG_ICASE);
  regcomp (comppat2, pattern2,0);

#ifdef FILEREAD_TEST
  cerr << "Now trying to match following text" << endl << "---INPUT BEGINS---" << endl;
  long long n = 0;
  while(input[n]!=0) {
    cout << input[n++];
  }
  cerr << "---INPUT ENDS---" << endl << "with pattern " << realpat << endl;
#endif /* FILEREAD_TEST */

  // find a line containing the keyword
  long long found = regexec(comppat,input,2,resarray,0);

  long long size,size2; // size of the line and an entry in the array, respectively
  double number; // temporary holder of the ascii-to-double-translated entry of the array
  char* parentes = (char*) malloc (MAXARRAYSIZE*sizeof(char)); // the array correspondig to the supplied keyword
  char* charnumber = (char*) malloc (MAXARRAYSIZE*sizeof(char)); // un-translated entry of the array

  if (found) {
    // search for the keyword was unsuccesful
    regexperr(found,realpat,strict); // regexperr kills the program if the error is fatal.
    // otherwise, there just wasn't a relevant entry in the datafile, and that's OK
    return 0;
  } else {
    // we found the keyword

    // first, find the size of the array (in chars) and copy it to the parentes variable
    size = resarray[1].rm_eo - resarray[1].rm_so+1;
    free(parentes);
    parentes = (char*) malloc (sizeof(char)*(size+1)); // +1 since we include the endline char
    long long i;
    for (i =0; i<size;i++) {
      parentes[i]=input[i+resarray[1].rm_so];
    }
    parentes[i]=0;

#ifdef FILEREAD_TEST
    cerr << "Found ";
    long long n = 0;
    while(parentes[n]!=0) {
      cout << parentes[n++];
    }
    cout << endl;
#endif /* FILEREAD_TEST */

    // then, start the parsing of the parentes variable and build the result array
    long long found2=0;
    while (!found2) { // as there is still more numbers left to parse

#ifdef FILEREAD_TEST
      cout << "Trying to match \"" << parentes << "\" with \"" << pattern2 << "\"" << endl;
#endif /* FILEREAD_TEST */

      found2 = regexec(comppat2,parentes,2,resarray2,0); // find the first number in parentes
      if (found2 && (regexperr(found2,pattern2,false)==-1)) {
	// something went wrong (other than us simply finishing the line)
	return -1;
      }
      if (!found2) {
	// we have a candidate number

	size2 = resarray2[1].rm_eo-resarray2[1].rm_so; // the size (in chars of the number

	// copy out the number chars to charnumber variable
	long long j;
	for (j =0; j<size2;j++) {
	  charnumber[j] = parentes[resarray2[1].rm_so+j];
	}
	charnumber[j]=0;
	
#ifdef FILEREAD_TEST
	cout << "sending \"" << charnumber <<"\" to the ascii to double parser"<< endl;
#endif /* FILEREAD_TEST */

	// translate the number from ascii to double
	number = atod(charnumber);

#ifdef FILEREAD_TEST
	cout << "Received " << number << " from ascii to double parser" << endl;
#endif /* FILEREAD_TEST */

	result[	count++] = number; // put the number into the result array

	// chop the now parsed number of the beginning of the character array.
	long long k;
	for (k = 0; k< size-resarray2[1].rm_eo;k++) {
	  parentes[k]=parentes[k+resarray2[1].rm_eo+1];
	}
	parentes[k]=0;
      }
    }
    // if we met no fatal errors, count now holds the number of entries in the array
    regfree(comppat);
    regfree(comppat2);
    delete(comppat);
    delete(comppat2);
    free(realpat);
    free(parentes);
    free(charnumber);
    return count;
  }
}

long long matchlines (char* input, const char* pattern, long long* result, bool strict) {

  /* function overloading: This will call the general double array matcher,
     and translate the resulting array to integers */

  double *intermideate = (double*) malloc(MAXARRAYSIZE*sizeof(double));
  long long matches=matchlines(input,pattern,intermideate,strict); // does the actual matching
  if (matches < 0) {
    // something went wrong. Propagate error to caller
    return matches;
  }
  for (long long i=0;i<matches;i++) {
    result[i]=(long long)intermideate[i];
#ifdef FILEREAD_VERBOSE
    if (result[i]!=intermideate[i]) { // report all non-integer elements to user
      cout << "Warning! truncating  " << intermideate[i] << "to" << result[i] << " while scanning for " << pattern << endl;
    }
#endif /* FILEREAD_VERBOSE */
  }
  free(intermideate);
  return matches;
}

long long multimatch (char* input, long long length, const char* pattern, double** result, long long* sizes, long long count) {
  /* this function is specifically built to match a series of lines, all with the same
     pattern. Caller must specify the number of hits expected, to force stringent
     thinking. Returns actual number of lines matched (should always be count). */
  
  regex_t* compline = new regex_t;

  const char* linepat = "\\(.*\\)\n.*"; // this matches a line in the subexpression

  // we will alter the input, so we must make a seperate copy
  char* myInput = (char*) malloc(length*sizeof(char)); 
  for (long long i = 0;i<length;i++) {
    myInput[i]=input[i];
  }

  regmatch_t resarray[2];

  regcomp (compline, linepat,REG_NEWLINE); // compile the regular expression

  long long hits=0; // counts the number of lines matched

  while (length>0) { // while there is still input left
    double *tempres = (double*) malloc(MAXARRAYSIZE*sizeof(double));
    long long regerr = regexec(compline,myInput,2,resarray,0); // match next line

    long long linesize = resarray[1].rm_eo-resarray[1].rm_so+1; // length of line

    char *line = (char*) malloc((1+linesize)*sizeof(char));
    long long j =0;
#ifdef FILEREAD_TEST
    cout << "Now processing line: ";
#endif /* FILEREAD_TEST */
    for (j=0;j<linesize;j++) {
      line[j]=myInput[resarray[1].rm_so+j]; // make a variable with the current line
#ifdef FILEREAD_TEST
      cout << line[j];
#endif /* FILEREAD_TEST */
    }
#ifdef FILEREAD_TEST
      cout << endl;
#endif /* FILEREAD_TEST */

    line[j]=0; // make sure its terminated

#ifdef FILEREAD_TEST
    cout << "Beginning matching of " << line << " with pattern " << pattern << endl;
#endif /* FILEREAD_TEST */
    long long matches = matchlines(line,pattern,tempres,false); // match patern on line
#ifdef FILEREAD_TEST
    cout << "Ended matching of " << line << " with pattern " << pattern << endl;
#endif /* FILEREAD_TEST */
    if (matches>0) { // matches positive if succesful
      if ((hits)==count) { // we already matched hits lines - error!
	cerr << "too many matches for pattern " << pattern << " - " << hits << "count must be specified correctly" << endl;
	exit(-1);
      }
      sizes[hits]=matches; // now holds the size of this array
      result[hits++]=tempres; // result is passed by reference
    }

    // chop the first line of myInput
    for (long long k=0;k<length-(resarray[1].rm_eo+1);k++) {
      myInput[k]=myInput[resarray[1].rm_eo+1+k];
    }
    length-=resarray[1].rm_eo+1; // modify length acordingly
    free(line);

  }

  // check if we have too few matches after going through all of input
  if ((hits)<count) {
    cerr << "too few matches for pattern " << pattern << " - " << hits << "count must be specified correctly" << endl;
    exit(-1);
  }

  delete(compline);
  free(myInput);


  return hits;
}


long long multimatch (char* input, long long length, const char* pattern, long long** result, long long* sizes, long long count) {
  double **tempresults=(double**)malloc(count*sizeof(double*));
#ifdef TEST_MULTIMATCH
  cerr << "Entering multimatch (long long)" << endl;
#endif
  long long res = multimatch (input,length,pattern,tempresults,sizes,count);
#ifdef TEST_MULTIMATCH
  cerr << "Multimatch i regc.cpp " << res << endl;
#endif
  
  for (long long i=0;i<count;i++) {
    result[i]=(long long*)malloc(sizes[i]*sizeof(long long));
    for (long long j=0;j<sizes[i];j++) {
      result[i][j]=(long long)tempresults[i][j];
      if (result[i][j]!=tempresults[i][j]) {
	cout << "WARNING: truncating " << tempresults[i][j] << " to " << result[i][j] << endl;
      }
    }
  }

  free(tempresults);
#ifdef TEST_MULTIMATCH
  cerr << "Ending multimatch (long long)" << endl;
#endif
  return res;
}


long long regexperr (long long code, const char *pattern,bool strict) {
  /* this function will identify various errors from the regexp parsing, and kill the
     program if neccesary. Otherwise we return 0. */

  switch (code) { // code is originally a return value from regexec
  case REG_ESPACE:
    cerr << "regular expression matching of " << pattern << "resulted in out of memory error.\n";
    exit(-1);
    break;
  case REG_NOMATCH:
    if (strict) {
      cerr << "regular expression matching of " << pattern << " didn't match any line in input, and you specified strict.\n";
      exit(-1);
    }
    else {
      return 0;
    }
    break;
  default:
    // could be improved...
    cerr << "Unforseen regexp error. Dying...";
    exit(-1);
  }
}

double atod (char* myString) {
  /* this function tries as hard as possible to translate myString from ascii to double */

  char *errormessage= new char[MAXARRAYSIZE+40];
  double result = 0; // will hold the resulting double
  long long point = 0; // we will divide by this number, if there is a . in the string
  for (long long i=0; myString[i]!=0;i++) {
    if ((myString[i]>47) && (myString[i]<58)) { // this character is a digit. Take it into account
      point*=10;
      result*=10;
      result+=myString[i]-48;
    } else if (myString[i] == '-') {
      if (i) { // if "-" is anything but the first character, croak and die
	cerr << "Tried to translate \""  <<  myString <<  "\" to double, but it is not a number." << endl;
	exit(-1);
      }
    } else if (myString[i] == '.') {
      if (!point) { // if this is the first comma, start counting places after it
	point = 1;
      } else { // only one comma allowed. Croak and die
	cerr << "Tried to translate \""  <<  myString <<  "\" to double, but it is not a number." << endl;
	exit(-1);
      }
    } else { // Should never happen: Regexp must have been odd!
	cerr << "Tried to translate \""  <<  myString <<  "\" to double, but it is not a number." << endl;
	exit(-1);
    }
  }
  if (point) { // if we encoutered a comma, we must divide by point
    result/=point;
  }
  if (myString[0] == '-') { // if we start with -, the number is negative. 
    result*=-1;
  }
  // return our result
  delete[](errormessage);
  return result;
}

