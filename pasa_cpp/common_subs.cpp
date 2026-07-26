#include "common_subs.h"
#include <iostream>
using namespace std;

void swap_ints(int& i, int& j) {
  int temp = i;
  i = j;
  j = temp;
}


int min (int a, int b) {
  if (a < b) {
    return (a);
  } else {
    return (b);
  }
}


int min (vector<int>& myList) {
  if (myList.size() == 0) {
    cerr << "Can't determine minimum of an empty vector." << endl;
    exit(1);
  }
  int Min = myList[0];
  for (int i = 1; i < myList.size(); i++) {
    if (myList[i] < Min) {
      Min = myList[i];
    }
  }
  return (Min);
}

int max (int a, int b) {
  if (a > b) {
    return (a);
  } else {
    return (b);
  }
}




int max (vector<int>& myList) {
  if (myList.size() == 0) {
    cerr << "Can't determine maximum of an empty vector." << endl;
    exit(1);
  }
  int Max = myList[0];
  for (int i = 1; i < myList.size(); i++) {
    if (myList[i] > Max) {
      Max = myList[i];
    }
  }
  return (Max);
}


bool overlap (struct coordset& a, struct coordset& b) {
  if (a.lend <= b.rend && a.rend >= b.lend) { //overlap
    return (true);
  } else {
    return (false);
  }
}




