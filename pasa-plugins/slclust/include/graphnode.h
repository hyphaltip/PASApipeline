#ifndef __GRAPHNODE__
#define __GRAPHNODE__

#include <vector>
#include <string>
#include <unordered_set>
#include "common.h"

using namespace std;

class Graphnode {
	public:

	string nodeName;
	vector<Graphnode*> linkedNodes;
	bool marked;

    // constructor
	Graphnode (string s);


    // get name of node
	const string& getNodename ();


    // add a linked node
	void addLinkedNode (Graphnode* a);


    // get number of linked nodes.
	int numLinkedNodes ();


    // describe node
    string toString ();


    // is `a` a linked node? O(1) average, backed by linkedNodesSet.
	bool isLinkedTo (Graphnode* a);

    // get the actual linked nodes list
	vector<Graphnode*>& getLinkedNodes ();

	private:

	// mirrors linkedNodes; exists only for O(1) duplicate-detection in
	// addLinkedNode and O(1) membership tests in isLinkedTo. linkedNodes
	// stays a vector so iteration order (and hence traversal/cluster output
	// order) is unaffected by this change.
	unordered_set<Graphnode*> linkedNodesSet;
};

#endif
