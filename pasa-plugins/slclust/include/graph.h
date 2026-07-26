#ifndef __GRAPH__
#define __GRAPH__

#include <vector>
#include <unordered_map>
#include <string>
#include "graphnode.h"
#include <iostream>
#include "common.h"

using namespace std;

class Graph {

 public:

  // all vertices
  vector<Graphnode*> allNodes;

  // map of vertex name to index in above vector
  unordered_map<string,int> nodeLookup;

  // destructor
  ~Graph(); // release memory allocated for graphnodes.


  // add a pair of linked vertices
  void addLinkedNodes(string a, string b);

  // helper function, get a graph node.
  Graphnode* getGraphnode (string s);

  // generate cluster output
  void printClusters ();


  // populates cluster list by finding connected components in graph.
  // iterative (explicit stack), not recursive: a chain of O(n) linked
  // singletons used to blow the call stack (SIGSEGV even with
  // `ulimit -s unlimited`, since that's set by the Perl wrapper and silently
  // swallowed if it fails -- see PerlLib/SingleLinkageClusterer.pm).
  void traverseGraph (Graphnode* start, vector<string>& cluster);


  /* Apply the Jaccard Coefficient */
  Graph* applyJaccardCoeff (float coeff);

  /* compare two nodes to get their link score (jaccard coeff) */
  float calclinkcoeff (Graphnode* a_node, Graphnode* b_node);

};


#endif
