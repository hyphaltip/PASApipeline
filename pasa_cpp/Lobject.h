#ifndef __Lobject__
#define __Lobject__

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdint>
#include <algorithm>

extern bool DEBUG;

using namespace std;

class Lobject {
 public: 
  Lobject(int index, int num_alignments);
    
  int index;
  int num_alignments;
  
  vector<uint64_t> contained_bits;
  
  int LscoreF;
  int LscoreR;
  int combined_score;
  Lobject* toLptr;
  Lobject* fromLptr;
  string toString();

  void setTraceIndices (vector<int>);
  vector<int> getTraceIndices();

  void setContainedBit(int i);
  void setContainedIndices(vector<int>);
  int num_contained_indices;

  int num_unique_contained(Lobject& other);
 private:

  vector<int> traceIndices;
  
};

#endif
