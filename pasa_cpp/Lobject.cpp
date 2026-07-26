#include "Lobject.h"

Lobject::Lobject (int index, int num_alignments) {
  this->index = index;
  this->num_alignments = num_alignments;
  LscoreF = 0; 
  LscoreR = 0;
  combined_score = 0;
  toLptr = 0;
  fromLptr = 0;
  num_contained_indices = 0;
  int num_words = (num_alignments + 63) / 64;
  contained_bits.resize(num_words, 0);
}

void Lobject::setContainedBit(int i) {
  int word = i / 64;
  uint64_t bit = 1ULL << (i % 64);
  contained_bits[word] |= bit;
  num_contained_indices++;
}

void Lobject::setContainedIndices(vector<int> indices) {
  fill(contained_bits.begin(), contained_bits.end(), 0);
  num_contained_indices = 0;
  LscoreF = 0;
  LscoreR = 0;

  for (int i=0; i < (int)indices.size(); i++) {
    int idx = indices[i];
    int word = idx / 64;
    uint64_t bit = 1ULL << (idx % 64);
    contained_bits[word] |= bit;
    LscoreF++;
    LscoreR++;
    num_contained_indices++;
  }
}

int Lobject::num_unique_contained (Lobject& other) {
  int num = 0;
  int min_words = min((int)contained_bits.size(), (int)other.contained_bits.size());
  
  // Process 4 words at a time for instruction-level parallelism (POPCNT has 3-cycle latency)
  int i = 0;
  #pragma GCC ivdep
  for (; i + 4 <= min_words; i += 4) {
    num += __builtin_popcountll(contained_bits[i]   & ~other.contained_bits[i]);
    num += __builtin_popcountll(contained_bits[i+1] & ~other.contained_bits[i+1]);
    num += __builtin_popcountll(contained_bits[i+2] & ~other.contained_bits[i+2]);
    num += __builtin_popcountll(contained_bits[i+3] & ~other.contained_bits[i+3]);
  }
  #pragma GCC ivdep
  for (; i < min_words; i++) {
    uint64_t diff = contained_bits[i] & ~other.contained_bits[i];
    num += __builtin_popcountll(diff);
  }
  for (int i = min_words; i < (int)contained_bits.size(); i++) {
    num += __builtin_popcountll(contained_bits[i]);
  }
  return num;
}

string Lobject::toString () {
  ostringstream os;
  
  os << "Lobject index: [" << index << "] has LscoreF: " << LscoreF 
     << ", LscoreR: " << LscoreR 
     << ", combinedScore: " << combined_score 
     << endl << "contains the following alignment indices:" << endl;
  
  for (int i = 0; i < num_alignments; i++) {
    int word = i / 64;
    uint64_t bit = 1ULL << (i % 64);
    if (contained_bits[word] & bit) {
      os << "\tindex: " << index << " contains " << i << endl;
    }
  }
  return (os.str());
}

void Lobject::setTraceIndices(vector<int> traces) {
  traceIndices = traces;
}

vector<int> Lobject::getTraceIndices() {
  return (traceIndices);
}

