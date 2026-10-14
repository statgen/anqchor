#ifndef __TSV_READER_H
#define __TSV_READER_H

#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <cstdint>
#include <map>
#include <set>

#define _FILE_OFFSET_BITS 64

extern "C" {
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
}
#include "qgen_error.h"

// a class to read tab-limited (tabixable) file using htsFile.h and kstring.h
class tsv_reader {
public:
  std::string filename;  // file name to read
  htsFile* hp;           // handles to htsFile 
  tbx_t* tbx;            // tabix index, if exists
  hts_itr_t* itr;        // iterator, if exists
  kstring_t str;         // kstring_t object to store the string
  int32_t lstr;          // length of the string
  int32_t nfields;       // number of tokenized fields (from kstring_t)
  int32_t* fields;       // indices of the starting point to each field 
  uint64_t nlines;        // total number of lines read
  int32_t delimiter;     // delimiter to tokenize
  
  bool open(const char* filename); // open a file, set up the file handle
  bool close();                    // close the file, returns false if fails
  int32_t read_line();             // read a line, returns the number of tokenzied fields (=nfields)
  const char* str_field_at(int32_t idx); // get a pointer to the string at index idx
  int32_t int_field_at(int32_t idx);     // get integer value at index idx
  uint64_t uint64_field_at(int32_t idx); // get unsigned long long (uint64_t) value at index idx
  int64_t int64_field_at(int32_t idx);   // get long long (int64_t) value at index idx  
  double double_field_at(int32_t idx);   // get double value at index idx
  int32_t store_to_vector(std::vector<std::string>& v); // store the tokenized values into string vectors
  bool jump_to(const char* reg);   // jump to a specific region using tabix
  bool jump_to(const char* chr, int32_t beg, int32_t end = INT_MAX); // jump to a specific region using tabix

  tsv_reader() : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0), delimiter(0) {
    str.l = str.m = 0; str.s = NULL;
  }

  tsv_reader(const char* filename) : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0), delimiter(0) {
    str.l = str.m = 0; str.s = NULL;    
    if ( !open(filename) )
      error("[E:%s:%s %s] Cannot open file %s for reading", __FILE__, __LINE__, __FUNCTION__, filename);    
  }  

  ~tsv_reader() {
    // if ( str.s ) free(str.s);
    // if ( fields ) free(fields);
    // if ( hp != NULL) close();
  }
};


class text_line_reader {
public:
  std::string filename;
  FILE* fp;
  char* buffer;
  int32_t max_line_length;
  int32_t last_line_length;
  int32_t count_lines;
  off_t cur_fp_offset;

  text_line_reader() : fp(NULL), buffer(NULL), max_line_length(1024), last_line_length(0), count_lines() {}

  bool open(const char* _filename, int32_t _max_line_length = 0);
  int32_t readline();
  int32_t close();
  ~text_line_reader();

  static int32_t load_to_set(const char* filename, std::set<std::string>& sset);
};

class dsv_hdr_reader {
protected:
  int32_t tokenize_fields();
  text_line_reader tlr;
  std::map<std::string,int32_t> col2idx;
  std::string delim;
  int32_t ret_hdr_line;
  std::vector<std::string> headers;
  std::vector<int32_t> fields;

  static bool is_substring_of_any(const char* str, std::vector<std::string>& queries) {
    for(int32_t i=0; i < (int32_t)queries.size(); ++i) {
      if ( strstr(str, queries[i].c_str()) != NULL )
        return true;
    }
    return false;
  }

public:
  int32_t nfields;

  bool open(const char* filename, int32_t _max_line_length = 65536, const char* _delim = " \t");
  bool read_hdr(const char* prefix_match, const char* prefix_skip = NULL, const char* ignore_chars = "#");
  bool read_hdr(std::vector<std::string>& prefixes_match, std::vector<std::string>& prefixes_skip, const char* ignore_chars = "#");
  bool close();
  inline bool has_header() { return col2idx.empty(); }
  int32_t read_line();

  inline const char* str_field_at(int32_t idx) {
    return &tlr.buffer[fields[idx]];
  }

  inline int32_t int_field_at(int32_t idx) {
    return atoi(&tlr.buffer[fields[idx]]);
  }

  inline double double_field_at(int32_t idx) {
    return atof(&tlr.buffer[fields[idx]]);
  }

  inline const char* str_field_colnames(const char* colname, const char* default_value = NULL) {
    std::map<std::string,int32_t>::iterator it = col2idx.find(colname);
    if ( it == col2idx.end() ) return default_value;
    else return &tlr.buffer[fields[it->second]];
  }

  inline int32_t int_field_colnames(const char* colname, int32_t default_value = 0) {
    std::map<std::string,int32_t>::iterator it = col2idx.find(colname);
    if ( it == col2idx.end() ) return default_value;
    else return atoi(&tlr.buffer[fields[it->second]]);
  }

  inline double double_field_colnames(const char* colname, double default_value = 0.0) {
    std::map<std::string,int32_t>::iterator it = col2idx.find(colname);
    if ( it == col2idx.end() ) return default_value;
    else return atof(&tlr.buffer[fields[it->second]]);
  }

  int32_t store_to_vector(std::vector<std::string>& v);

  dsv_hdr_reader(): nfields(0), delim(" \t"), ret_hdr_line(-1) {}
  dsv_hdr_reader(const char* filename) : nfields(0), delim(" \t") {
    if (!open(filename)) {
      error("Cannot open file %s for reading", filename);
    }
  }

  ~dsv_hdr_reader() { tlr.close(); }
};

#endif
