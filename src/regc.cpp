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
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <regex.h>
#include <RLexact.h>
#include <regex>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include <cstring>

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

#ifdef TEST_FILEREAD
  // write the read content of the file to stdout
  cerr << "Read the following from datafile " << filename << endl;
  for (long long i =0;i<filesize;i++) {
    cerr << filedata[i];
  }
  cerr << endl;
#endif /* TEST_FILEREAD */
  datafile.close();
}





long long matchlines(char* input, const char* pattern, double* result, bool strict)
{
    if (!input || !pattern || !result) {
        return strict ? -1 : 0;
    }

    // Build POSIX regex: ^[[:space:]]*pattern[[:space:]]+(.*)
    // (We’ll manually strip trailing // comments after capture.)
    char regexbuf[512];
    // NOTE: If 'pattern' may contain regex metacharacters, you must escape them.
    // Given your "very strict input" requirement, we assume literal keyword.
    std::snprintf(regexbuf, sizeof(regexbuf),
                  "^[[:space:]]*%s[[:space:]]+(.*)", pattern);

    regex_t re;
    if (regcomp(&re, regexbuf, REG_EXTENDED | REG_ICASE) != 0) {
        return -1; // regex compile failed
    }

    regmatch_t matches[2];
    long long count = 0;
    bool patternFound = false;

    // Iterate input as raw char* safely
    const char* p = input;
    while (*p != '\0') {
        // Find end of current line (or end of buffer)
        const char* nl = std::strchr(p, '\n');
        size_t lineLen = nl ? static_cast<size_t>(nl - p) : std::strlen(p);

        // Create a temporary, null-terminated buffer for the line
        // (Avoids reading past memory when calling regexec.)
        char* linebuf = static_cast<char*>(std::malloc(lineLen + 1));
        if (!linebuf) {
            regfree(&re);
            return -1;
        }
        std::memcpy(linebuf, p, lineLen);
        linebuf[lineLen] = '\0';

        // Try to match the pattern on this line
        if (regexec(&re, linebuf, 2, matches, 0) == 0) {
            patternFound = true;

            // Extract captured group (values after the pattern)
            int start = matches[1].rm_so;
            int end   = matches[1].rm_eo;

            if (start >= 0 && end > start) {
                // Create a values buffer (null-terminated)
                size_t valsLen = static_cast<size_t>(end - start);
                char* values = static_cast<char*>(std::malloc(valsLen + 1));
                if (!values) {
                    std::free(linebuf);
                    regfree(&re);
                    return -1;
                }
                std::memcpy(values, linebuf + start, valsLen);
                values[valsLen] = '\0';

                // Remove any trailing // comment
                char* comment = std::strstr(values, "//");
                if (comment) {
                    *comment = '\0'; // truncate at comment
                }

                // Parse space-separated numbers with strtod
                const char* s = values;
                char* endptr = NULL;
                while (*s != '\0') {
                    // skip leading spaces
                    while (std::isspace(static_cast<unsigned char>(*s))) ++s;
                    if (*s == '\0') break;

                    double val = std::strtod(s, &endptr);
                    if (endptr == s) {
                        // token wasn't a valid number; stop parsing this line
                        break;
                    }
                    result[count++] = val;

                    s = endptr;
                    // skip spaces before next token
                    while (std::isspace(static_cast<unsigned char>(*s))) ++s;
                }

                std::free(values);
            }

            std::free(linebuf);
            break; // One pattern per line; stop after first match
        }

        std::free(linebuf);

        // Advance to next line (or end)
        if (nl) {
            p = nl + 1;
        } else {
            break;
        }
    }

    regfree(&re);

    if (patternFound && count == 0 && strict) {
        return -1; // keyword found but no usable values
    }
    // If not found, return     // If not found, return 0 (original behavior)
    cout << "\n Count is " << count << "\n";
    for (int i = 0; i < count; i++){
      cout << "\n Result " << i << " = " << result[i] << "\n";
    }
    return count;
}




long long matchlines_wrapper(char* input, const char* pattern, long long* result, bool strict) {

  /* function overloading: This will call the general double array matcher,
     and translate the resulting array to integers */

  double *intermediate = (double*) malloc(MAXARRAYSIZE*sizeof(double));
  
    if (!intermediate) {
        std::cerr << "matchlines(long long*): malloc failed\n";
        return -1;
    }

  long long matches=matchlines(input,pattern,intermediate,strict); // does the actual matching
  cout << "Matchlines 1 has been run! intermediate is " << intermediate[0] << "\n";
  cout << "Matches is: " << matches << "\n";
  if (matches < 0) {
    // something went wrong. Propagate error to caller
    return matches;
  }
  for (long long i=0;i<matches;i++) {
    cout << "If matches is ";
    result[i]=(long long)intermediate[i];
#ifdef FILEREAD_VERBOSE
    if (result[i]!=intermediate[i]) { // report all non-integer elements to user
      cout << "Warning! truncating  " << intermediate[i] << "to" << result[i] << " while scanning for " << pattern << endl;
    }
#endif /* FILEREAD_VERBOSE */
  }
  free(intermediate);
  return matches;
}

long long multimatch (char* input, long long length, const char* pattern, double** result, long long* sizes, long long count) {
  /* this function is specifically built to match a series of lines, all with the same
     pattern. Caller must specify the number of hits expected, to force stringent
     thinking. Returns actual number of lines matched (should always be count). */
  
  regex_t* compline = new regex_t;

  const char* linepat = "\\(.*\\)\n.*"; // this matches a line in the subexpression

  // we will alter the input, so we must make a seperate copy,
  char* myInput = (char*) malloc((length + 1)*sizeof(char)); 
  for (long long i = 0;i<length;i++) {
    myInput[i]=input[i];
  }
  myInput[length] = '\0'; // Add null terminator. ESSENTIAL

  regmatch_t resarray[2];

  regcomp (compline, linepat,REG_NEWLINE); // compile the regular expression

  long long hits=0; // counts the number of lines matched

  while (length>0) { // while there is still input left
    double *tempres = (double*) malloc(MAXARRAYSIZE*sizeof(double));
    long long regerr = regexec(compline,myInput,2,resarray,0); // match next line

    long long linesize = resarray[1].rm_eo-resarray[1].rm_so+1; // length of line

    char *line = (char*) malloc((1+linesize)*sizeof(char));
    long long j =0;
#ifdef TEST_FILEREAD
    cout << "Now processing line: ";
#endif /* TEST_FILEREAD */
    for (j=0;j<linesize;j++) {
      line[j]=myInput[resarray[1].rm_so+j]; // make a variable with the current line
#ifdef TEST_FILEREAD
      cout << line[j];
#endif /* TEST_FILEREAD */
    }
#ifdef TEST_FILEREAD
      cout << endl;
#endif /* TEST_FILEREAD */

    line[j]=0; // make sure its terminated

#ifdef TEST_FILEREAD
    cout << "Beginning matching of " << line << " with pattern " << pattern << endl;
#endif /* TEST_FILEREAD */
    long long matches = matchlines(line,pattern,tempres,false); // match patern on line
#ifdef TEST_FILEREAD
    cout << "Ended matching of " << line << " with pattern " << pattern << endl;
#endif /* TEST_FILEREAD */
    if (matches>0) { // matches positive if succesful
      if ((hits)==count) { // we already matched hits lines - error!
	cerr << "too many matches for pattern " << pattern << " - " << hits << " count must be specified correctly" << endl;
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
    cerr << "too few matches for pattern " << pattern << ". " << hits << " found, "<< count << " required. count must be specified correctly" << endl;
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

