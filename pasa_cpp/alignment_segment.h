
#include "common_structs.h"
#include "common_subs.h"
#include <string>

using namespace std;

#ifndef __ALIGNMENT_SEGMENT__
#define __ALIGNMENT_SEGMENT__

/* position of a segment within its alignment.  This was a std::string carried on
   every segment; alignments are copied often enough that the per-segment string
   dominated the profile. */
enum segment_type {
  SEGMENT_UNKNOWN,
  SEGMENT_SINGLE,
  SEGMENT_FIRST,
  SEGMENT_LAST,
  SEGMENT_INTERNAL
};

class Alignment_segment {
 public:

  Alignment_segment (int genomic_lend, int genomic_rend); //constructor
  Alignment_segment (struct coordset&); // copies coordset info over to local new coordset.

  struct coordset& get_coords();
  const struct coordset& get_coords() const;
  void set_coords(int lend, int rend);

  void set_type (segment_type);
  segment_type get_type ();

  void set_left_splice_junction(bool);
  bool get_left_splice_junction();
  
  void set_right_splice_junction(bool);
  bool get_right_splice_junction();

  string toString();
  segment_type type;
  
 private:
  void init();
  
  struct coordset coords; //genome coords of segment.
  bool has_left_splice_junction;
  bool has_right_splice_junction;
  
};

#endif

