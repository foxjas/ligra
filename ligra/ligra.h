// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <math.h>
#include <vector>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"

#include "papi.h"
#include <omp.h>
using namespace std;

typedef struct{
  int num;
  int *events;
  long long int *counters;
}papi_t;

//*****START FRAMEWORK*****

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options { DENSE, DENSE_FORWARD };

template <class vertex, class F>
bool* edgeMapDense(graph<vertex> GA, bool* vertexSubset, F &f, bool parallel = 0) {
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for (long i=0; i<numVertices; i++) {
    next[i] = 0;
    if (f.cond(i)) {
      G[i].decodeInNghBreakEarly(i, vertexSubset, f, next, parallel);
    }
  }}
  return next;
}

template <class vertex, class F>
bool* edgeMapDenseForward(graph<vertex> GA, bool* vertexSubset, F &f) {
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for(long i=0;i<numVertices;i++) next[i] = 0;}
  {parallel_for (long i=0; i<numVertices; i++){
    if (vertexSubset[i]) {
      G[i].decodeOutNgh(i, vertexSubset, f, next);
    }
  }}
  return next;
}

template <class vertex, class F>
pair<long,uintE*> edgeMapSparse(vertex* frontierVertices, uintE* indices, 
        uintT* degrees, uintT m, F &f, 
        long remDups=0, uintE* flags=NULL) {
  uintT* offsets = degrees;
  long outEdgeCount = sequence::plusScan(offsets, degrees, m);
  uintE* outEdges = newA(uintE,outEdgeCount);
  {parallel_for (long i = 0; i < m; i++) {
      uintT v = indices[i], o = offsets[i];
      vertex vert = frontierVertices[i]; 
      vert.decodeOutNghSparse(v, o, f, outEdges);
    }}
  uintE* nextIndices = newA(uintE, outEdgeCount);
  if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);
  // Filter out the empty slots (marked with -1)
  long nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonMaxF());
  free(outEdges);
  return pair<long,uintE*>(nextM, nextIndices);
}

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class vertex, class F>
vertexSubset edgeMap(graph<vertex> GA, vertexSubset &V, F f, intT threshold = -1, 
		 char option=DENSE, bool remDups=false) {
  long numVertices = GA.n, numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  long m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(vertex,m);
  {parallel_for (long i=0; i < m; i++){
    vertex v = G[V.s[i]];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  uintT outDegrees = sequence::plusReduce(degrees, m);
  if (outDegrees == 0) return vertexSubset(numVertices);
  if (m + outDegrees > threshold) { 
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R = (option == DENSE_FORWARD) ? 
      edgeMapDenseForward(GA,V.d,f) : 
      edgeMapDense(GA, V.d, f, option);
    vertexSubset v1 = vertexSubset(numVertices, R);
    //cout << "size (D) = " << v1.m << endl;
    return v1;
  } else { 
    pair<long,uintE*> R = 
      remDups ? 
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    //cout << "size (S) = " << R.first << endl;
    free(degrees);
    free(frontierVertices);
    return vertexSubset(numVertices, R.first, R.second);
  }
}

//*****VERTEX FUNCTIONS*****

//Note: this is the optimized version of vertexMap which does not
//perform a filter
template <class F>
void vertexMap(vertexSubset V, F add) {
  long n = V.numRows(), m = V.numNonzeros();
  if(V.isDense) {
    {parallel_for(long i=0;i<n;i++)
	if(V.d[i]) add(i);}
  } else {
    {parallel_for(long i=0;i<m;i++)
	add(V.s[i]);}
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
    }

    //cond function that always returns true
    inline bool cond_true (intT d) { return 1; }

    template<class vertex>
    void Compute(graph<vertex>&, commandLine);

    int parallel_main(int argc, char* argv[]) {
    commandLine P(argc,argv," [-s] <inFile>");
    char* iFile = P.getArgument(0);
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool binary = P.getOptionValue("-b");
    bool verbose = P.getOptionValue("-v");
    long rounds = P.getOptionLongValue("-rounds",1);
    printf("Rounds: %ld\n", rounds);
    double avgTime;
    double totalTime;
    double tDelta;
    vector<double> avgTimeBySrc;
    
    papi_t papi;
    const int num_events = 2;
    int events[num_events] = {PAPI_TLB_DM, PAPI_L2_TCM};
    //const int num_events = 4;
    //int events[num_events] = {PAPI_LD_INS, PAPI_SR_INS, PAPI_BR_INS, PAPI_TOT_INS}; 
    long long int counters[num_events] = {0};
    papi.num = num_events;
    papi.events = events;
    papi.counters = counters;
    //papi_struct_set(papi, num_events, events, counters);
    int retval;
    #ifdef PERF
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fprintf(stderr,"PAPI library init error!\\n");
        exit(1); 
    }
    for (int i=0; i<num_events; i++) {
       retval = PAPI_query_event(papi.events[i]);
       if (retval != PAPI_OK) {
           fprintf(stderr, "Counter listed at %d is not available; exitting\n", i);
           exit(1);
       }
    }
    #endif
    /*
    if (PAPI_thread_init((unsigned long (*)())omp_get_thread_num) != PAPI_OK){
        fprintf(stderr, "PAPI thread init failed\n");
        exit(1);
    }   
    */
  vector<long> srcList = P.getOptionLongVector("--R");
  if (!srcList.size()) {
    srcList.push_back(P.getOptionLongValue("-r", 0));
  }
  if (verbose) {
      printf("Sources:");
      for (vector<long>::const_iterator it = srcList.begin(); it != srcList.end(); it++) {
        printf(" %ld", *it);
      }
      printf("\n");
  }

  if (compressed) {
    if (symmetric) {
      graph<compressedSymmetricVertex> G =
        readCompressedGraph<compressedSymmetricVertex>(iFile,symmetric); //symmetric graph
      Compute(G,P);
      for(int r=0;r<rounds;r++) {
        startTime();
        Compute(G,P);
        nextTime("Running time");
      }
      G.del();
    } else {
      graph<compressedAsymmetricVertex> G =
        readCompressedGraph<compressedAsymmetricVertex>(iFile,symmetric); //asymmetric graph
      Compute(G,P);
      if(G.transposed) G.transpose();
      for(int r=0;r<rounds;r++) {
        startTime();
        Compute(G,P);
        nextTime("Running time");
        if(G.transposed) G.transpose();
      }
      G.del();
    }
  } else {
    if (symmetric) {
      graph<symmetricVertex> G =
        readGraph<symmetricVertex>(iFile,compressed,symmetric,binary); //symmetric graph
      Compute(G,P);
      int code;
      for (vector<long>::const_iterator it = srcList.begin(); it != srcList.end(); it++) {
        P.setOptionValue("-r", to_string(*it));
        tDelta = 0;
        for(int r=0;r<rounds;r++) {
          startTime();

          #ifdef PERF
          PAPI_start_counters(papi.events, papi.num);
          #endif
          Compute(G,P);
          //PAPI_strerror(code);

          #ifdef PERF
          code = PAPI_accum_counters(papi.counters, papi.num);
          #endif
          tDelta += stopT();
        }
        avgTimeBySrc.push_back(tDelta/rounds);
      }
      totalTime = totalTime();
      #ifdef PERF
      printf("PAPI_read_counters exit code: %d\n", code);
      for (int i=0; i<papi.num; i++) {
        printf("Counter %d=%lld\n", i, (papi.counters[i]/rounds));
      }
      #endif
      avgTime = totalTime/(rounds*srcList.size());
      std::cout << "Average time: " << avgTime << std::endl;
      if (verbose) {
        std::cout << "Average time by source:";
        for (vector<double>::const_iterator it = avgTimeBySrc.begin(); it != avgTimeBySrc.end(); it++) {
          std::cout << " " << *it;
        }
        std::cout << std::endl;
      }
      G.del();
    } else {
      graph<asymmetricVertex> G =
        readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary); //asymmetric graph
      Compute(G,P);
      if(G.transposed) G.transpose();
      int code;
      for (vector<long>::const_iterator it = srcList.begin(); it != srcList.end(); it++) {
        P.setOptionValue("-r", to_string(*it));
        tDelta = 0;
        for(int r=0;r<rounds;r++) {
          startTime();

          #ifdef PERF
          PAPI_start_counters(papi.events, papi.num);
          #endif
          Compute(G,P);
          //PAPI_strerror(code);

          #ifdef PERF
          code = PAPI_accum_counters(papi.counters, papi.num);
          #endif
          tDelta += stopT();
        }
        avgTimeBySrc.push_back(tDelta/rounds);
      }
      totalTime = totalTime();
      #ifdef PERF
      printf("PAPI_read_counters exit code: %d\n", code);
      for (int i=0; i<papi.num; i++) {
        printf("Counter %d=%lld\n", i, (papi.counters[i]/rounds));
      }
      #endif
      avgTime = totalTime/(rounds*srcList.size());
      std::cout << "Average time: " << avgTime << std::endl;
      if (verbose) {
        std::cout << "Average time by source:";
        for (vector<double>::const_iterator it = avgTimeBySrc.begin(); it != avgTimeBySrc.end(); it++) {
          std::cout << " " << *it;
        }
        std::cout << std::endl;
      }
      //print_stats(papi, time_ms);
      G.del();
    }
  }
}
#endif
