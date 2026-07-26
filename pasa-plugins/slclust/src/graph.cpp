#include "graph.h"


// constructor
//currently using default


// destructor
// rlease memory allocated for graph nodes:
Graph::~Graph () {
  for (unsigned int i=0; i < allNodes.size(); i++) {
    delete allNodes[i];
  }
}
  

void Graph::addLinkedNodes (string a, string b) {
  if (VERBOSE >= debug)
    cout << "Adding nodes: " << a << ", " << b << endl;
  
  Graphnode* a_node = getGraphnode(a);
  Graphnode* b_node = getGraphnode(b);
  
  a_node->addLinkedNode(b_node);
  b_node->addLinkedNode(a_node);
  
}


Graphnode* Graph::getGraphnode (string s) {
  if (VERBOSE >= debug)
    cout << "Getting graphnode for : " << s << endl;
  
  unordered_map<string, int>::iterator p;
  // see if (s) already exists:
  p = nodeLookup.find(s);
  Graphnode* s_node;
  if (p != nodeLookup.end()) {
    if (VERBOSE >= debug) 
      cout << "Found it." << endl;
    
    int s_index = p->second;
    s_node = allNodes[s_index];
  
  } else {
    if (VERBOSE >= debug)
      cout << "Didn't find it, inserting node." << endl;
  
    s_node = new Graphnode(s);
    int s_index = allNodes.size();
    allNodes.push_back(s_node);
    if (VERBOSE >= debug)
      cout << "Added new node with " << s << " at position: " << s_index << endl;
    nodeLookup[s] = s_index;
  }
  
  return (s_node);
  
}


void Graph::printClusters () {
  vector<string> cluster;
  for (unsigned int i=0; i < allNodes.size(); i++) {
    Graphnode* g = allNodes[i];
    if (g->marked) {
      continue; // already in another cluster
    }
        
    if (!g->marked) {
      
      if (VERBOSE >= debug)
        cout << "Starting graph traversal from node[" << i << "]: " << g->toString() << endl;
      
      traverseGraph(g, cluster);
    }
    // cout << "Cluster : ";
    for (unsigned int j=0; j < cluster.size(); j++) {
      cout << cluster[j] << " ";
    }
    cout << endl;
    cluster.clear(); // init for next cluster
  }
}


void Graph::traverseGraph (Graphnode* start, vector<string>& cluster) {
  // Explicit-stack DFS, pushed in reverse per node so pop order reproduces
  // exactly what the old recursive version visited: it fully explored
  // linkedNodes[0]'s subtree before moving to linkedNodes[1], depth-first.
  vector<Graphnode*> stack;
  stack.push_back(start);

  while (! stack.empty()) {
    Graphnode* g = stack.back();
    stack.pop_back();

    if (g->marked) {
      continue;
    }
    cluster.push_back(g->getNodename());
    g->marked = true;

    vector<Graphnode*>& linkednodes = g->getLinkedNodes();
    for (int i = (int)linkednodes.size() - 1; i >= 0; i--) {
      if (! linkednodes[i]->marked) {
        stack.push_back(linkednodes[i]);
      }
    }
  }
}


 /* Apply the Jaccard Coefficient */
Graph* Graph::applyJaccardCoeff (float coeff) {
  
  Graph* gph = new Graph();
  
  int num_nodes = allNodes.size();
  
  
  for (unsigned int i=0; i < allNodes.size(); i++) {
    
    if (VERBOSE >= info) {
      cerr << "\rJaccard: examining node " << (i+1)  << " of " << num_nodes << "       "; 
    }
    
    Graphnode* g = allNodes[i];
    
    // examine each pair of linked nodes
    vector<Graphnode*>& linkednodes = g->getLinkedNodes();
    for (unsigned int j=0; j < linkednodes.size(); j++) {
      Graphnode* linkednode = linkednodes[j];
      if (calclinkcoeff(g, linkednode) >= coeff) {
        gph->addLinkedNodes(g->getNodename(), linkednode->getNodename());
      }
      
    }
    
  }
  
  return (gph);
  
}


float Graph::calclinkcoeff (Graphnode* a_node, Graphnode* b_node) {
  
  int num_a_links = a_node->numLinkedNodes();
  int num_b_links = b_node->numLinkedNodes();

  // determine number shared vertices. Both nodes' adjacency carries an O(1)
  // membership set (linkedNodesSet, via isLinkedTo), so there's no need to
  // build a fresh map<string,bool> of one side's neighbor names per edge
  // examined -- that map construction, and the string hashing/compares it
  // required, was the dominant cost of Jaccard filtering (measured: -j was
  // 25x slower than the unfiltered pass on the same input before this
  // change). Probe from the smaller adjacency list into the larger one's
  // set, same as the Rust implementation, to bound the work by
  // min(|A|,|B|) rather than |A|.
  Graphnode* smaller = (num_a_links <= num_b_links) ? a_node : b_node;
  Graphnode* larger  = (num_a_links <= num_b_links) ? b_node : a_node;
  vector<Graphnode*>& probe_links = smaller->getLinkedNodes();

  int num_common = 0;
  for (unsigned int i=0; i < probe_links.size(); i++) {
    if (larger->isLinkedTo(probe_links[i])) {
      num_common++;
    }
  }
  
  int total_num_vertices = num_a_links + num_b_links - num_common;
  int total_in_common = num_common + 2; // a and b aren't included in common count
  
  float linkcoeff = (float)total_in_common/total_num_vertices;
  
  if (VERBOSE >= debug) {
    cerr << "Link score between " << a_node->getNodename() << " and " 
         << b_node->getNodename() << " = " << linkcoeff << endl;
  }
  
  return (linkcoeff);
  
  
}
  







