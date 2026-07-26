#include "graphnode.h"


// constructor
Graphnode::Graphnode (string s) {
  this->nodeName = s;
  marked = false;
}



const string& Graphnode::getNodename () {
  return (this->nodeName);
}

void Graphnode::addLinkedNode (Graphnode* a) {
  // Graph::getGraphnode hands out exactly one Graphnode* per unique name, so
  // pointer identity is equivalent to the name comparison this replaced, and
  // is O(1) average instead of an O(degree) linear scan of string compares.
  if (linkedNodesSet.insert(a).second) {
    linkedNodes.push_back(a);
  }
}



int Graphnode::numLinkedNodes () {
  return (linkedNodes.size());
}


string Graphnode::toString () {
  string marked = (this->marked) ? "true" : "false";
  string ret = "(" + nodeName + ", marked: " + marked + ") linked to the following: ";
  for (unsigned int i=0; i < linkedNodes.size(); i++) {
    Graphnode* g = linkedNodes[i];
    string othermarked = (g->marked) ? "true" : "false";
    ret += "(" + g->getNodename() + ", marked: " + othermarked + ") ";
    if (i < linkedNodes.size() - 1) {
      ret += ", ";
    }
  }

  return (ret);
}



bool Graphnode::isLinkedTo (Graphnode* a) {
  return (linkedNodesSet.count(a) > 0);
}

vector<Graphnode*>& Graphnode::getLinkedNodes () {
  return (linkedNodes);
}
