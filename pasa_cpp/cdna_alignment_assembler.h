#ifndef __CDNA_ALIGNMENT_ASSEMBLER__
#define __CDNA_ALIGNMENT_ASSEMBLER__

#include "Lobject.h"
#include "cdna_alignment.h"
#include "alignment_segment.h"
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

extern bool DEBUG;

class CDNA_alignment_assembler {
 public:
  CDNA_alignment_assembler(vector<CDNA_alignment>& incomingAlignments);

  ~CDNA_alignment_assembler();  
  
  void assembleAlignments();
  
  vector<CDNA_alignment> get_assemblies ();  

  void set_fuzzlength(int);
  
  string toAlignIllustration(int lineLength);


 private:
  vector<CDNA_alignment>& alignments;
  vector<CDNA_alignment> assemblies;

  vector<vector<int> > assembly_containment_list;
  int fuzzlength;
  
  bool canMerge(CDNA_alignment&, CDNA_alignment&);
  CDNA_alignment mergeAlignments (CDNA_alignment& A, CDNA_alignment& B); 
  
  vector<Lobject> Lobjects;

  void do_full_Fscan();
  void do_full_Rscan();

  bool encapsulates(CDNA_alignment& A, CDNA_alignment& B);
  void determine_compatibilities_and_encapsulations();
  vector<int> forwardTrace (int);
  vector<int> backTrace(int); 
  CDNA_alignment create_assembly(vector<int>); 
  vector<int> get_top_scoring_alignment();
  vector<int> get_alignment_assembly_nucleating_at_alignment_index(int index);
  vector<int> unique_entries(vector<vector<int> >);
  void populateLobjects();
  
  vector<vector<int>> compatibilities;
  vector<vector<int>> encapsulations;

  /* scratch buffers reused by mergeAlignments so the per-call splice-coordinate
     sets cost no allocations; mergeAlignments is only ever called serially. */
  vector<int> merge_left_splicecoords;
  vector<int> merge_right_splicecoords;
  int num_alignments;
  
  Lobject* get_max_missing_Lobj(vector<Lobject*>&, vector<bool>&);
  
};

#endif
