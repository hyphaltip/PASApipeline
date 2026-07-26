#include "cdna_alignment_assembler.h"
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <climits>

#ifdef _OPENMP
#include <omp.h>
#endif

bool sort_CDNA_alignments (const CDNA_alignment& a, const CDNA_alignment& b) {
  return a.get_coords().lend < b.get_coords().lend;
}


CDNA_alignment_assembler::CDNA_alignment_assembler (vector<CDNA_alignment>& incomingAlignments) : alignments (incomingAlignments) {
  
  if (DEBUG) {
    cout << "-sorting alignments by lend position." << endl;
  }
  
  sort (alignments.begin(), alignments.end(), sort_CDNA_alignments);
  if (DEBUG) {
    cout << "-done sorting alignments." << endl;
    cout << "Here's the alignments orderd by lend pos." << endl;
    for (int i=0; i < alignments.size(); i++) {
      cout << i << " " << alignments[i].toString() << endl;
    }
  }
  
  fuzzlength = 20;
  
  num_alignments = incomingAlignments.size();
  
}


CDNA_alignment_assembler::~CDNA_alignment_assembler () {
}


bool sort_Lobjects_via_combined_Lscore_R_F (Lobject* a, Lobject* b) {
  
  if (a->combined_score < b->combined_score) {
    return(true);
  } else {
    return (false);
  }
}



void CDNA_alignment_assembler::assembleAlignments() {
  compatibilities.resize(num_alignments);
  encapsulations.resize(num_alignments);
  
  if (DEBUG) {
    cout << "Assembling alignments." << endl;
    cout << "Determining compatibilities and encapsulations." << endl;
  }

  determine_compatibilities_and_encapsulations();
  
  if (DEBUG) { cout << "Populating Lobjects." << endl; }
  
  populateLobjects();
  
  if (DEBUG) { cout << "Doing full Fscan." << endl; }
  
  do_full_Fscan();
  
  if (DEBUG) { cout << "Getting maximum scoring alignment assembly." << endl; }
  
  vector<int> topAssemblyIndices = get_top_scoring_alignment();
  CDNA_alignment assembly = create_assembly(topAssemblyIndices);
  
  ostringstream assemblyTitle;
  assemblyTitle << "assembly_" << assemblies.size();
  assembly.set_title(assemblyTitle.str());
  assemblyTitle.str("");
  assemblies.push_back(assembly);
  assembly_containment_list.push_back(topAssemblyIndices);
  
  if (topAssemblyIndices.size() == num_alignments) {
    return;
  }
  
  assemblies.clear();
  assembly_containment_list.clear();
  
  if (DEBUG) { cout << "Doing full Rscan." << endl; }
  do_full_Rscan();
  
  vector<bool> accountedFor(num_alignments, false);
  
  vector<Lobject*> untraversedLobjs;
  for (int i=0; i < num_alignments; i++) {
    if (!accountedFor[i]) {
      untraversedLobjs.push_back(&Lobjects[i]);
    }
  }
  
  for (int i = untraversedLobjs.size() - 1; i >= 0; i--) {
    Lobject& nextBestLobj = *(untraversedLobjs[i]);
    int nucleatingIndex = nextBestLobj.index;
    int combined_score = nextBestLobj.LscoreF + nextBestLobj.LscoreR - nextBestLobj.num_contained_indices;
  
    if (DEBUG) { 
      cout << "Combined Lscore for index: " << nucleatingIndex << " is " 
           << combined_score 
           <<" (LscoreF: " << nextBestLobj.LscoreF 
           << ", LscoreR: " << nextBestLobj.LscoreR 
           << ")" <<  endl; 
    }
    
    nextBestLobj.combined_score = combined_score;

    vector<int> alignmentIndices = get_alignment_assembly_nucleating_at_alignment_index(nucleatingIndex);
    nextBestLobj.setTraceIndices(alignmentIndices);
    
  }
  
  sort (untraversedLobjs.begin(), untraversedLobjs.end(), sort_Lobjects_via_combined_Lscore_R_F);
  
  vector< vector<Lobject*> > lobjsSameSizeVec;
  int max_index = untraversedLobjs.size() - 1;
  int curr_score = untraversedLobjs[ max_index ]->combined_score;
  vector<Lobject*> curr_bin;
  
  if (DEBUG) { cout << "\n\n\n**** Binning lobjs with same combined scores.\n";}
  
  for (int i = untraversedLobjs.size() - 1; i >= 0; i--) {
    Lobject* lobj = untraversedLobjs[i];
    if (DEBUG) { cout << "-adding lobj[" << i << "] to bin. " << endl; }
    int score = lobj->combined_score;
    
    if (DEBUG) { cout << "curr_score: " << curr_score << ", lobj score: " << score << endl; }
    
    if (score == curr_score) {
      if (DEBUG) { cout << "\tscores same, adding to current bin." <<endl; }
      curr_bin.push_back(lobj);
    }
    else {
      if (DEBUG) { cout << "\tscores diff, archiving current bin, creating new bin with this lobj." << endl; }
      lobjsSameSizeVec.push_back(curr_bin);
      curr_bin.clear();
      curr_bin.push_back(lobj);
      curr_score = score;
    }
  }
  
  if (curr_bin.size() > 0) {
    if (DEBUG) { cout << "-pushing last lobj bin onto the stack." << endl; }
    lobjsSameSizeVec.push_back(curr_bin);
    curr_bin.clear();
  }
  
  if (DEBUG) {
    cout << "\n-describing binned Lobjs:\n\n";
    cout << "There are " << lobjsSameSizeVec.size() << " bins of Lobjs.\n\n";
    for (int i = 0; i < lobjsSameSizeVec.size(); i++) {
      vector<Lobject*> lobj_bin = lobjsSameSizeVec[i];
      cout << "\tbin: " << i << " of size: " << lobj_bin.size() << endl; 
      for (int j=0; j < lobj_bin.size(); j++) {
        cout << "\t\tLobj: " << lobj_bin[j]->index << endl; }
    }
  }
  
  for (int i = 0; i < lobjsSameSizeVec.size(); i++) {
    vector<Lobject*> lobj_bin = lobjsSameSizeVec[i];
    
    if (DEBUG) { cout << "-**** Analyzing lobj_bin at bin pos: " << i << endl; }
        
    while (lobj_bin.size() > 0) {
      Lobject* max_missing_Lobj = get_max_missing_Lobj(lobj_bin, accountedFor);
      
      if (max_missing_Lobj == NULL) {
        if (DEBUG) { cout << "no more missing alignments in bin " << i << ". Trying next bin." << endl << endl; }
        break;
      }
      
      int nucleatingIndex = max_missing_Lobj->index;
      vector<int> alignmentIndices = max_missing_Lobj->getTraceIndices();
      bool hasUnconsumed = false;
      for (int j=0; j < alignmentIndices.size(); j++) {
        int index = alignmentIndices[j];
        if (! accountedFor[index]) {
          hasUnconsumed = true;
          break;
        }
      }
      if (hasUnconsumed) {
        CDNA_alignment newAssembly = create_assembly(alignmentIndices);
        assemblyTitle << "assembly_" << assemblies.size();
        newAssembly.set_title(assemblyTitle.str());
        assemblyTitle.str("");
        assemblies.push_back(newAssembly);
        assembly_containment_list.push_back(alignmentIndices);
        
        for (int j=0; j<alignmentIndices.size(); j++) {
          int index = alignmentIndices[j];
          accountedFor[index] = true;
        }
        
        bool all_accounted_for = true;
        for (int j=0; j < num_alignments; j++) {
          if (! accountedFor[j]) {
            all_accounted_for = false;
            break;
          }
        }
        if (all_accounted_for) {
          return;
        }
      } else {
        if (DEBUG) { cout << "Apparently, no unconsumed alignments in assembly nucleating at index: " << nucleatingIndex << endl; }
        break;
      }

      if (DEBUG) { 
        cout << "Now removing max Lobj from current bin\n";
        cout << "Before bin size: " << lobj_bin.size() << endl;
      }
      
      vector<Lobject*> new_lobj_bin;
      for (int i=0; i < lobj_bin.size(); i++) {
        Lobject* lobj = lobj_bin[i];
        if (lobj != max_missing_Lobj) {
          new_lobj_bin.push_back(lobj);
        }
      }
      
      if (DEBUG) { cout << "doing replacment.  new bin size: " << new_lobj_bin.size()  << ", old bin size: " << lobj_bin.size() << endl; }
      lobj_bin = new_lobj_bin;
      if (DEBUG) { cout << "after replacement, bin size: " << lobj_bin.size() << endl; }
            
    }

  }
  
  cerr << "Not all alignments were accounted for by alignment assemblies." << endl;
  exit(5);
}


vector<CDNA_alignment> CDNA_alignment_assembler::get_assemblies () {
  return (assemblies);
}

void CDNA_alignment_assembler:: set_fuzzlength(int length) {
  fuzzlength = length;
}

bool CDNA_alignment_assembler::canMerge(CDNA_alignment& a1, CDNA_alignment& a2) {
  if (! overlap(a1.get_coords(), a2.get_coords())) {
    if (DEBUG) { cout << "-can't merge: alignment coordsets don't overlap." << endl; }
    return (false);
  }
  
  if (a1.get_orientation() != a2.get_orientation()) {
    if (DEBUG) { cout << "-can't merge: diff orientations." << endl; }
    return (false);
  }
  
  vector<Alignment_segment>& a1_segments = a1.get_alignment_segments();
  vector<Alignment_segment>& a2_segments = a2.get_alignment_segments();
  
  int starting_a1 = -1;
  int starting_a2 = -1;
  for (int i=0; i < (int)a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i];
    for (int j=0; j < (int)a2_segments.size(); j++) {
      Alignment_segment& a2_seg = a2_segments[j];
      if (overlap(a1_seg.get_coords(), a2_seg.get_coords())) {
        starting_a1 = i;
        starting_a2 = j;
        break;
      }
    }
    if (starting_a1 != -1 && starting_a2 != -1) {
      break;
    }
  }
  
  if (starting_a1 == -1 || starting_a2 == -1) {
    if (DEBUG) { cout << "can't merge: couldn't align two segments of the alignments." << endl; }
    return(false);
  }
  
  if (! (starting_a1 == 0 || starting_a2 == 0)) {
    if (DEBUG) { cout << "can't merge: couldn't map first segment of either to the other." << endl; }
    return (false);
  }
  
  while (starting_a1 < (int)a1_segments.size() && starting_a2 < (int)a2_segments.size()) {
    Alignment_segment& a1_seg = a1_segments[starting_a1];
    Alignment_segment& a2_seg = a2_segments[starting_a2];
    struct coordset& a1_seg_coords = a1_seg.get_coords();
    struct coordset& a2_seg_coords = a2_seg.get_coords();
    int a1_lend = a1_seg_coords.lend;
    int a1_rend = a1_seg_coords.rend;
    int a2_lend = a2_seg_coords.lend;
    int a2_rend = a2_seg_coords.rend;
    
    if (overlap(a1_seg_coords, a2_seg_coords)) {
      
      if (a1_seg.get_left_splice_junction() || a2_seg.get_left_splice_junction()) {
        if (a1_seg.get_left_splice_junction() && a2_seg.get_left_splice_junction() && a1_lend != a2_lend) {
          if (DEBUG) { cout << "can't merge: diff left splice sites." << endl; }
          return (false);
        } else if (a1_seg.get_left_splice_junction() && (a2_lend + fuzzlength < a1_lend)) {
          if (DEBUG) { cout << "can't merge: left splice analysis, not within fuzz distance." << endl; }
          return (false);
        } else if (a2_seg.get_left_splice_junction() && (a1_lend + fuzzlength < a2_lend)) {
          if (DEBUG) { cout << "can't merge: left splice analysis, not within fuzz distance." << endl;}
          return (false);
        }
      }
      
      if (a1_seg.get_right_splice_junction() || a2_seg.get_right_splice_junction()) {
        
        if (a1_seg.get_right_splice_junction() && a2_seg.get_right_splice_junction() && a1_rend != a2_rend) {
          if (DEBUG) { cout << "can't merge: diff right splice sites." << endl; }
          return (false);
        } else if (a1_seg.get_right_splice_junction() && (a2_rend - fuzzlength > a1_rend)) {
          if (DEBUG) { cout << "can't merge: right splice analysis, not within fuzz distance." << endl; }
          return (false);
        } else if (a2_seg.get_right_splice_junction() && (a1_rend - fuzzlength > a2_rend)) {
          if (DEBUG) { cout << "can't merge: right splice analysis, not within fuzz distance." << endl; }
          return (false);
        }
      }
    } else {
      if (DEBUG) { cout << "can't merge: Two ordered segments do not overlap each other." << endl; }
      return (false);
    }
    starting_a1++;
    starting_a2++;
  }
  
  if (DEBUG) { cout << "-Merge possible. Passed all tests." << endl; }
  return (true);
}



CDNA_alignment CDNA_alignment_assembler::mergeAlignments(CDNA_alignment& A, CDNA_alignment& B) {
  unordered_set<int> leftsplicecoords;
  unordered_set<int> rightsplicecoords;
  
  char orientation = A.get_orientation();
  
  vector<Alignment_segment>& a1_segments = A.get_alignment_segments();
  for (int i=0; i < (int)a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i]; 
    struct coordset& a1_coordset = a1_seg.get_coords();
    int lend = a1_coordset.lend;
    int rend = a1_coordset.rend;
    if (a1_seg.get_left_splice_junction()) {
      leftsplicecoords.insert(lend);
    }
    if (a1_seg.get_right_splice_junction()) {
      rightsplicecoords.insert(rend);
    }
  }
  
  vector<Alignment_segment>& a2_segments = B.get_alignment_segments();
  for (int i=0; i < (int)a2_segments.size(); i++) {
    Alignment_segment& a2_seg = a2_segments[i]; 
    struct coordset& a2_coordset = a2_seg.get_coords();
    int lend = a2_coordset.lend;
    int rend = a2_coordset.rend;
    if (a2_seg.get_left_splice_junction()) {
      leftsplicecoords.insert(lend);
    }
    if (a2_seg.get_right_splice_junction()) {
      rightsplicecoords.insert(rend);
    }
  }
  
  vector<struct coordset> merged_coords;
  
  for (int i=0; i < (int)a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i];
    struct coordset& a1_coordset = a1_seg.get_coords();
    int a1_lend = a1_coordset.lend;
    int a1_rend = a1_coordset.rend;
    int merged_lend = -1;
    int merged_rend = -1;
    for (int j=0; j < (int)a2_segments.size(); j++) {
      Alignment_segment& a2_seg = a2_segments[j];
      struct coordset& a2_coordset = a2_seg.get_coords();
      int a2_lend = a2_coordset.lend;
      int a2_rend = a2_coordset.rend;
      if (overlap(a1_coordset, a2_coordset)) {
        
        if (leftsplicecoords.count(a1_lend)) {
          merged_lend = a1_lend;
        } else if (leftsplicecoords.count(a2_lend)) {
          merged_lend = a2_lend;
        } else {
          merged_lend = min(a1_lend, a2_lend);
        }
        
        if (rightsplicecoords.count(a1_rend)) {
          merged_rend = a1_rend;
        } else if (rightsplicecoords.count(a2_rend)) {
          merged_rend = a2_rend;
        } else {
          merged_rend = max(a1_rend, a2_rend);
        }
        break;
      }
    }
    struct coordset merged_coordset;
    if (merged_lend != -1 && merged_rend != -1) {
      merged_coordset.lend = merged_lend;
      merged_coordset.rend = merged_rend;
    } else {
      merged_coordset.lend = a1_lend;
      merged_coordset.rend = a1_rend;
    }
    merged_coords.push_back(merged_coordset);
  }
  
  for (int i=0; i < (int)a2_segments.size(); i++) {
    Alignment_segment& a2_seg = a2_segments[i];
    struct coordset& a2_coords = a2_seg.get_coords();
    int a2_lend = a2_coords.lend;
    int a2_rend = a2_coords.rend;
    bool overlapFlag = false;
    for (int j=0; j < (int)merged_coords.size(); j++) {
      struct coordset& m_coords = merged_coords[j];
      if (overlap(a2_coords, m_coords)) {
        overlapFlag = true;
        break;
      }
    }
    if (! overlapFlag) {
      struct coordset m_coords;
      m_coords.lend = a2_lend;
      m_coords.rend = a2_rend;
      merged_coords.push_back(m_coords);
    }
  }
  
  vector<Alignment_segment> new_seg_list;
  for (int i=0; i < (int)merged_coords.size(); i++) {
    struct coordset& coords = merged_coords[i];
    Alignment_segment new_seg (coords);
    new_seg_list.push_back(new_seg);
  }
  
  CDNA_alignment merged_alignment (new_seg_list, orientation);
  
  return (merged_alignment);
  
}


void CDNA_alignment_assembler::do_full_Fscan() {
  
  for (int i=1; i < num_alignments; i++) {
    Lobject& Lobj = Lobjects[i];
    int top_score = 0;
    int top_scoring_index = -1;
    
    for (int j : compatibilities[i]) {
      if (j >= i) continue;
      
      bool containment =  binary_search(encapsulations[i].begin(), encapsulations[i].end(), j)
                       || binary_search(encapsulations[j].begin(), encapsulations[j].end(), i);
      
      if (containment) continue;
      
      int curr_total_score = Lobjects[j].LscoreF + Lobj.num_unique_contained(Lobjects[j]);
      if (DEBUG) { cout << "FSCAN(" << i << "," << j << ") \tcurr total score: " << curr_total_score << endl;}
      
      if (curr_total_score > top_score) {
        top_scoring_index = j;
        top_score = curr_total_score;
      }
    }
    if (top_scoring_index > -1) {
      Lobject& topLobj = Lobjects[top_scoring_index];
      Lobj.fromLptr = &topLobj;
      Lobj.LscoreF = top_score;
      if (DEBUG) { cout << "-Assigning fromLptr from index " << i << " -> " << top_scoring_index << endl; }
    }
  }
}



void CDNA_alignment_assembler::do_full_Rscan() {
  for (int i= num_alignments - 2; i >= 0; i--) {
    Lobject& Lobj = Lobjects[i];
    int top_score = 0;
    int top_scoring_index = -1;
    
    for (int j : compatibilities[i]) {
      if (j <= i) continue;
      
      bool containment =  binary_search(encapsulations[i].begin(), encapsulations[i].end(), j)
                       || binary_search(encapsulations[j].begin(), encapsulations[j].end(), i);
      if (containment) continue;
      
      int curr_total_score = Lobjects[j].LscoreR + Lobj.num_unique_contained(Lobjects[j]);
      if (DEBUG) { cout << "RSCAN(" << i << "," << j << ") \tcurr total score: " << curr_total_score << endl;}
      
      if (curr_total_score > top_score) {
        top_scoring_index = j;
        top_score = curr_total_score;
      }
    }
    if (top_scoring_index > -1) {
      Lobject& topLobj = Lobjects[top_scoring_index];
      Lobj.toLptr = &topLobj;
      Lobj.LscoreR = top_score;
      if (DEBUG) { cout << "-Assigning toLptr from index " << i << " -> " << top_scoring_index << " with top score " << top_score << endl; }
    }
  }
}

bool CDNA_alignment_assembler::encapsulates (CDNA_alignment& A, CDNA_alignment& B) {
  struct coordset& Acoords = A.get_coords();
  int a_lend = Acoords.lend;
  int a_rend = Acoords.rend;
  
  struct coordset& Bcoords = B.get_coords();
  int b_lend = Bcoords.lend;
  int b_rend = Bcoords.rend;
  
  if (b_lend >= a_lend && b_rend <= a_rend) {
    return (true);
  } else {
    return (false);
  }
  
}

void CDNA_alignment_assembler::determine_compatibilities_and_encapsulations() {
  
  // Build interval tree: sorted array of (lend, index) pairs
  vector<pair<int32_t, int>> starts;
  starts.reserve(num_alignments);
  for (int i = 0; i < num_alignments; i++) {
    starts.emplace_back(alignments[i].get_coords().lend, i);
  }
  sort(starts.begin(), starts.end());
  
  compatibilities.resize(num_alignments);
  encapsulations.resize(num_alignments);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
  for (int i = 0; i < num_alignments; i++) {
    
    struct coordset& icoords = alignments[i].get_coords();
    int i_rend = icoords.rend;
    int i_lend = icoords.lend;
    
    // Binary search for first start > i_rend
    auto it = upper_bound(starts.begin(), starts.end(),
                          make_pair(i_rend, INT_MAX));
    int limit = it - starts.begin();
    
    for (int si = 0; si < limit; si++) {
      int j = starts[si].second;
      if (j <= i) continue;
      
      // Complete overlap check: rend must span back past i.lend
      if (alignments[j].get_coords().rend < i_lend) continue;
      
      if (DEBUG) { cout << "can merge " << i << " to " << j << " ?" << endl; }
      if (canMerge(alignments[i], alignments[j])) {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
        compatibilities[i].push_back(j);
        compatibilities[j].push_back(i);
        if (encapsulates(alignments[i],alignments[j])) {
          if (DEBUG) { cout << "alignment " << i << " encapsulates " << j << endl; }
          encapsulations[i].push_back(j);
        }
        if (encapsulates(alignments[j],alignments[i])) {
          if (DEBUG) { cout << "alignment " << j << " encapsulates " << i << endl; }
          encapsulations[j].push_back(i);
        }
        }
      }
    }
  }
  
  // Sort and deduplicate each row for binary-search-based lookups
  for (int i = 0; i < num_alignments; i++) {
    sort(compatibilities[i].begin(), compatibilities[i].end());
    compatibilities[i].erase(unique(compatibilities[i].begin(), compatibilities[i].end()), compatibilities[i].end());
    sort(encapsulations[i].begin(), encapsulations[i].end());
    encapsulations[i].erase(unique(encapsulations[i].begin(), encapsulations[i].end()), encapsulations[i].end());
  }
}


void CDNA_alignment_assembler::populateLobjects () {
  Lobjects.clear();
  Lobjects.reserve(num_alignments);
  for (int i=0; i < num_alignments; i++) {
    Lobject L (i, num_alignments);
    vector<int> contained;
    contained.push_back(i);
    for (int j : encapsulations[i]) {
      contained.push_back(j);
    }
    L.setContainedIndices(contained);
    Lobjects.push_back(L);
  }
}



vector<int> CDNA_alignment_assembler::forwardTrace (int startIndex) {
  if (DEBUG) { cout << "Beginning forwardTrace, starting at index: " << startIndex << endl; }
  vector<bool> tracker(num_alignments, false);
  Lobject* Lobj = & Lobjects[startIndex];
  while (Lobj != 0) {
    if (DEBUG) { cout << "trace index: " << Lobj->index << endl; }
    if (DEBUG) { cout << Lobj->toString(); }
    
    // Iterate only set bits using ctz
    for (int w = 0; w < (int)Lobj->contained_bits.size(); w++) {
      uint64_t word = Lobj->contained_bits[w];
      while (word) {
        int bit = __builtin_ctzll(word);
        int align_idx = w * 64 + bit;
        tracker[align_idx] = true;
        word &= word - 1;
      }
    }
    Lobj = Lobj->toLptr;
  }
  
  vector<int> unique;
  for (int i = 0; i < num_alignments; i++) {
    if (tracker[i]) {
      unique.push_back(i);
    }
  }
  return (unique);
}

vector<int> CDNA_alignment_assembler::backTrace(int startIndex) {
  if (DEBUG) { cout << "Beginning backTrace, starting at index: " << startIndex << endl; }
  vector<bool> tracker(num_alignments, false);
  Lobject* Lobj = & Lobjects[startIndex];
  while (Lobj != 0) {
    if (DEBUG) { cout << "trace index: " << Lobj->index << endl;
                 cout << Lobj->toString(); }
    
    for (int w = 0; w < (int)Lobj->contained_bits.size(); w++) {
      uint64_t word = Lobj->contained_bits[w];
      while (word) {
        int bit = __builtin_ctzll(word);
        int align_idx = w * 64 + bit;
        tracker[align_idx] = true;
        word &= word - 1;
      }
    }
    Lobj = Lobj->fromLptr;
  }
  
  vector<int> unique;
  for (int i = 0; i < num_alignments; i++) {
    if (tracker[i]) {
      unique.push_back(i);
    }
  }
  return (unique);
}

CDNA_alignment CDNA_alignment_assembler::create_assembly(vector<int> Alignment_index_listing) {
  if (Alignment_index_listing.empty()) {
    cerr << "empty list of indices, can't create assembly." << endl;
    exit(6);
  }
  
  sort(Alignment_index_listing.begin(), Alignment_index_listing.end());
  
  int alignment_index = Alignment_index_listing[0];
  
  CDNA_alignment assembly = alignments[alignment_index];
  for (int i = 1; i < (int)Alignment_index_listing.size(); i++) {
    alignment_index = Alignment_index_listing[i];
    CDNA_alignment& nextAlignment = alignments[alignment_index];
    CDNA_alignment newAssembly = mergeAlignments(assembly, nextAlignment);
    assembly = newAssembly;
  }
  
  return (assembly);
}


vector<int> CDNA_alignment_assembler::get_top_scoring_alignment() {
  int top_score = 0;
  int top_scoring_index = -1;
  for (int i=0; i < (int)Lobjects.size(); i++) {
    Lobject* L = & Lobjects[i];
    int Lscore = L->LscoreF;
    if (DEBUG) { cout << "LscoreF of alignment " << i << " is " << Lscore << endl; }
    if (Lscore > top_score) {
      top_scoring_index = i;
      top_score = Lscore;
    }
  }
  if (DEBUG) { cout << "Top score is alignment " << top_scoring_index << " with LscoreF of " << top_score << endl; }
  vector<int> indexListing = backTrace(top_scoring_index);
  return (indexListing);
}

vector<int> CDNA_alignment_assembler::get_alignment_assembly_nucleating_at_alignment_index(int index) {
  vector<vector<int> > vecvec;
  vecvec.push_back(backTrace(index));
  vecvec.push_back(forwardTrace(index));
  return (unique_entries(vecvec));
}

vector<int> CDNA_alignment_assembler::unique_entries(vector<vector<int> > vecvec) {
  vector<bool> uniqueMap(num_alignments, false);
  for (int j=0; j < (int)vecvec.size(); j++) {
    vector<int> myvec = vecvec[j];
    for (int i=0; i < (int)myvec.size(); i++) {
      int entry = myvec[i];
      uniqueMap[entry] = true;
    }
  }
  
  vector<int> uniqueEntries;
  for (int i = 0; i < num_alignments; i++) {
    if (uniqueMap[i]) {
      uniqueEntries.push_back(i);
    }
  }
  
  return(uniqueEntries);
}

string CDNA_alignment_assembler::toAlignIllustration (int lineLength) {
  
  int min_coord;
  int max_coord;
  vector<int> allCoords;
  for (int i=0; i < (int)alignments.size(); i++) {
    struct coordset& coords = alignments[i].get_coords();
    allCoords.push_back(coords.lend);
    allCoords.push_back(coords.rend);
  }
  sort (allCoords.begin(), allCoords.end());
  min_coord = allCoords[0];
  max_coord = allCoords[allCoords.size()-1];
  
  int rel_max = max_coord - min_coord;
  ostringstream alignment_text;
  ostringstream assembly_summary;
  alignment_text << "Individual Alignments: (" << num_alignments << ")" << endl;
  for (int i=0; i < (int)alignments.size(); i++) {
    alignment_text << alignments[i].toAlignIllustration(min_coord, rel_max, lineLength) << " index: [" << i << "]" << endl;
  }
  
  if (assemblies.size() != 0) {
    alignment_text << endl << "ASSEMBLIES: (" << assemblies.size() << ")" << endl;
    for (int i=0; i < (int)assemblies.size(); i++) {
      alignment_text << assemblies[i].toAlignIllustration(min_coord, rel_max, lineLength) << " score: (" << assembly_containment_list[i].size() << ") contains [";
      assembly_summary << "assembly: (" << i << ") contains alignments: [";
      vector<int> assemblyIndexList = assembly_containment_list[i];
      for (int j=0; j < (int)assemblyIndexList.size(); j++) {
        int alignmentIndex = assemblyIndexList[j];
        alignment_text << alignmentIndex;
        if (j != (int)assemblyIndexList.size() -1) {
          alignment_text << ",";
        }
        assembly_summary << alignments[alignmentIndex].get_title(); 
        if (j != (int)assemblyIndexList.size() -1) {
          assembly_summary << ",";
        }
      }
      alignment_text << "]" << endl;
      assembly_summary << "]" << " with structure [" << assemblies[i].toString() << "]" << " score: (" << assembly_containment_list[i].size() << ")" << endl;
    }
    alignment_text << assembly_summary.str();
  }
  return (alignment_text.str());
}


Lobject* CDNA_alignment_assembler::get_max_missing_Lobj (vector<Lobject*>& lobj_bin, vector<bool>& accountedFor) {
  
  int max_missing = 0;
  Lobject* lobj = NULL;
  
  if (DEBUG) { 
    cout << "\n\nget_max_missing_Lobj()" << endl
         << "sifting thru bin of size: " << lobj_bin.size() << endl;
  }
  
  for (int i=0; i < (int)lobj_bin.size(); i++) {
    Lobject* curr_lobj = lobj_bin[i];
    
    if (DEBUG) { cout << "\tanalyzing lobj[" << i << "] " << " at index: " << curr_lobj->index << endl; }
    
    vector<int> alignmentIndices = curr_lobj->getTraceIndices();
    
    int num_missing = 0;
    for (int j=0; j < (int)alignmentIndices.size(); j++) {
      int index = alignmentIndices[j];
      if (! accountedFor[index]) {
        num_missing++;
      }
    }
    if (num_missing > max_missing) {
      max_missing = num_missing;
      lobj = curr_lobj;
    }
  }

  if (DEBUG) { 
    if (max_missing > 0) {
      cout << "-reporting max missing as " << max_missing << " provided by Lobj index: " << lobj->index << endl;
    } else {
      cout << "-none missing in bin this round." << endl;
    }
  }
  
  return (lobj);
}
