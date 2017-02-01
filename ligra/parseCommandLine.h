// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef _PARSE_COMMAND_LINE
#define _PARSE_COMMAND_LINE

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <vector>
using namespace std;

struct commandLine {
  int argc;
  char** argv;
  string comLine;
  commandLine(int _c, char** _v, string _cl) 
    : argc(_c), argv(_v), comLine(_cl) {}

  commandLine(int _c, char** _v) 
    : argc(_c), argv(_v), comLine("bad arguments") {}

  void badArgument() {
    cout << "usage: " << argv[0] << " " << comLine << endl;
    abort();
  }

  // get an argument
  // i is indexed from the last argument = 0, second to last indexed 1, ..
  char* getArgument(int i) {
    if (argc < 2+i) badArgument();
    return argv[argc-1-i];
  }

  // looks for two filenames
  pair<char*,char*> IOFileNames() {
    if (argc < 3) badArgument();
    return pair<char*,char*>(argv[argc-2],argv[argc-1]);
  }

  pair<int,char*> sizeAndFileName() {
    if (argc < 3) badArgument();
    return pair<int,char*>(atoi(argv[argc-2]),(char*) argv[argc-1]);
  }

  bool getOption(string option) {
    for (int i = 1; i < argc; i++)
      if ((string) argv[i] == option) return true;
    return false;
  }

  char* getOptionValue(string option) {
    for (int i = 1; i < argc-1; i++)
      if ((string) argv[i] == option) return argv[i+1];
    return NULL;
  }

  void setOptionValue(string option, string value) {
    for (int i = 1; i < argc-1; i++) {
      if ((string) argv[i] == option)
        argv[i+1] = (char *)value.c_str();
    }
  }

  string getOptionValue(string option, string defaultValue) {
    for (int i = 1; i < argc-1; i++)
      if ((string) argv[i] == option) return (string) argv[i+1];
    return defaultValue;
  }

  int getOptionIntValue(string option, int defaultValue) {
    for (int i = 1; i < argc-1; i++)
      if ((string) argv[i] == option) {
	int r = atoi(argv[i+1]);
	return r;
      }
    return defaultValue;
  }

  vector<long> getOptionLongVector(string option) {
    char *r = NULL;
    for (int i = 1; i < argc-1; i++) {
      if ((string) argv[i] == option) {
        r = argv[i+1];
        break;
      }
    }
    if (r == NULL) return vector<long>();
    vector<long> arr;
    long val;
    int n;
    char* tok = r;
    n = sscanf(tok, "%ld", &val);
    if (!n) return vector<long>();
    while (n) {
      arr.push_back(val);
      tok = strchr(tok, ',' );
      if (tok == NULL)
        break;
      tok += 1;
      n = sscanf(tok, "%ld", &val);
    }
    return arr;
  }

  long getOptionLongValue(string option, long defaultValue) {
    for (int i = 1; i < argc-1; i++)
      if ((string) argv[i] == option) {
	long r = atol(argv[i+1]);
	return r;
      }
    return defaultValue;
  }

  double getOptionDoubleValue(string option, double defaultValue) {
    for (int i = 1; i < argc-1; i++)
      if ((string) argv[i] == option) {
	double val;
	if (sscanf(argv[i+1], "%lf",  &val) == EOF) {
	  badArgument();
	}
	return val;
      }
    return defaultValue;
  }

};
 
#endif // _PARSE_COMMAND_LINE
