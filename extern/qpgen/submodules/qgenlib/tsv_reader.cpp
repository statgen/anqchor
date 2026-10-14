/* The MIT License

   Copyright (c) 2016-2022 Hyun Min Kang <hmkang@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "qgenlib/tsv_reader.h"
#define DSV_NOT_YET_PEEKED -9

bool tsv_reader::open(const char* filename) {
  this->filename = filename;
  
  htsFormat fmt;
  fmt.category = unknown_category;
  fmt.format = text_format;
  fmt.version.major = fmt.version.minor = -1;
  fmt.compression = no_compression;
  fmt.compression_level = -1;
  fmt.specific = NULL;  
  hp = hts_open_format(filename, "r", &fmt);
  if ( hp == NULL ) {
    return false;
  }
  return true;
}

bool tsv_reader::close() {
//  notice("foo %s", filename.c_str());
  if ( hp == NULL ) return false;
//  notice("bar");
  int32_t ret = hts_close(hp);
//  notice("baz %d", ret);
  hp = NULL;
  return ret == 0;
}

int32_t tsv_reader::read_line() {
  if ( itr == NULL ) {
    //if ( ( str.s != NULL ) && ( lstr > 0 ) ) free(str.s);
    lstr = hts_getline(hp, KS_SEP_LINE, &str);
  }
  else {
    lstr = tbx_itr_next(hp, tbx, itr, &str);
  }

  if ( lstr <= 0 ) {
    nfields = 0;
    fields = NULL;
    return 0; // lstr;
  }

  // need to free the previously allocated fields!!
  if ( fields != NULL ) { free(fields); fields = NULL; }  
  fields = ksplit(&str, delimiter, &nfields);

  //notice("lstr = %d, str = %s, delim = %d", lstr, str.s, delimiter);
    
  ++nlines;
  return nfields;
}

bool tsv_reader::jump_to(const char* chr, int32_t beg, int32_t end) {
  char buf[65536];
  snprintf(buf, 65536, "%s:%d-%d", chr, beg, end);
  return jump_to(buf);
}

bool tsv_reader::jump_to(const char* reg) {
  if ( tbx == NULL ) {
    tbx = tbx_index_load(filename.c_str());
    if ( !tbx ) error("[E:%s] Could not load .tbi/.csi index of %s\n", __PRETTY_FUNCTION__, filename.c_str());
  }

  if ( itr != NULL )
    tbx_itr_destroy(itr);

  itr = tbx_itr_querys(tbx, reg);
  
  if ( itr == NULL ) {
    notice("Failed jumping to %s, tbx = %x, itr = %x", reg, tbx, itr);
    return false;
  }
  else return true;
}

const char* tsv_reader::str_field_at(int32_t idx) {
  if ( idx >= nfields ) {
    error("[E:%s:%d %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);
    //return NULL;
  }
  return ( &str.s[fields[idx]] );
}

int32_t tsv_reader::int_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%d %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);    
  return ( atoi(&str.s[fields[idx]]) );
}

int64_t tsv_reader::int64_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%d %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);
  return ( strtoll(&str.s[fields[idx]], NULL, 10) );
}

uint64_t tsv_reader::uint64_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%d %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);
  return ( strtoull(&str.s[fields[idx]], NULL, 10) );
}

double tsv_reader::double_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%d %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);    
  return ( atof(&str.s[fields[idx]]) );
}

int32_t tsv_reader::store_to_vector(std::vector<std::string>& v) {
  v.resize(nfields);
  for(int32_t i=0; i < nfields; ++i) {
    v[i].assign(&str.s[fields[i]]);
  }
  return nfields;
}

// implementations for text_line_reader class

int32_t text_line_reader::load_to_set(const char* filename, std::set<std::string>& sset) {
  text_line_reader tlr;
  if ( tlr.open(filename) ) {
    while(tlr.readline()) {
      tlr.buffer[tlr.last_line_length-1] = '\0';
      sset.insert(tlr.buffer);
    }
    return tlr.count_lines;
  }
  else {
    return 0;
  }
}

bool text_line_reader::open(const char* _filename, int32_t _max_line_length) {
  filename.assign(_filename);
  if ( _max_line_length > 0 )
    max_line_length = _max_line_length;
  buffer = (char*)malloc(max_line_length);
  fp = fopen(filename.c_str(), "r");
  if ( fp == NULL ) {
    fprintf(stderr,"ERROR: Cannot open file %s for reading", filename.c_str());
    return false;
  }
  cur_fp_offset = ftello(fp);
  return true;
}

int32_t text_line_reader::readline() {
  const char* s = fgets(buffer, max_line_length, fp);
  off_t new_fp_offset = ftello(fp);
  while( new_fp_offset - cur_fp_offset == max_line_length -1 ) { // need to expand the current line length
    int32_t new_line_length = max_line_length * 2;
    buffer = (char*)realloc(buffer,new_line_length);
    s = fgets(&buffer[max_line_length-1], max_line_length+1, fp);
    max_line_length = new_line_length;
    new_fp_offset = ftello(fp);
  }
  last_line_length = new_fp_offset - cur_fp_offset;
  ++count_lines;
  cur_fp_offset = new_fp_offset;
  return s != NULL;
}

int32_t text_line_reader::close() {
  if ( fp != NULL ) {
    int32_t ret = fclose(fp);
    if ( ret == 0 )
      fp = NULL;
    return 0;
  }
  else
    return -1;
}

text_line_reader::~text_line_reader() {
  if ( fp != NULL )
    close();
  if ( buffer ) {
    free(buffer);
    buffer = NULL;
  }
}


// implementations for dsv_hdr_reader class

bool dsv_hdr_reader::open(const char* filename, int32_t _max_line_length, const char* _delim) {
  delim = _delim;
  if ( tlr.open(filename, _max_line_length) ) {
    ret_hdr_line = tlr.readline();
    return true;
  }
  else return true;
}

bool dsv_hdr_reader::close() {
  return tlr.close();
}

bool dsv_hdr_reader::read_hdr(const char* prefix_match, const char* prefix_skip, const char* ignore_chars) {
  std::vector<std::string> pmatch;
  pmatch.push_back(prefix_match);
  std::vector<std::string> pskip;
  pskip.push_back(prefix_skip);
  return read_hdr(pmatch, pskip, ignore_chars);
}

bool dsv_hdr_reader::read_hdr(std::vector<std::string>& prefixes_match, std::vector<std::string>& prefixes_skip, const char* ignore_chars) {
  while ( ret_hdr_line > 0 ) {
    if ( is_substring_of_any(tlr.buffer, prefixes_match) ) {
      headers.emplace_back(tlr.buffer);
      tokenize_fields();

      while ( (ignore_chars != NULL) && ( tlr.buffer[fields[0]] != '\0' ) && (strchr(ignore_chars,tlr.buffer[fields[0]]) != NULL) ) 
        ++fields[0];

      for(int32_t i=0; i < nfields; ++i) {
        col2idx[str_field_at(i)] = i;
      }
      
      ret_hdr_line = DSV_NOT_YET_PEEKED;
      return true;
    }
    else if ( is_substring_of_any(tlr.buffer, prefixes_skip ) ) {
      headers.emplace_back(tlr.buffer);
      ret_hdr_line = tlr.readline();
    }
    else { // no header line
      return false;
    }
  }
  return false;
}

int32_t dsv_hdr_reader::read_line() {
  if ( ret_hdr_line == DSV_NOT_YET_PEEKED ) { // 2nd line or forward. Have to read the line
    if ( tlr.readline() == 0 )
      return 0;
  }
  else if ( ret_hdr_line == 0 )
    return 0;
  else if ( ret_hdr_line > 0 ) 
    ret_hdr_line = DSV_NOT_YET_PEEKED;
  else {
    fprintf(stderr, "ERROR: ret_hdr_line == %d observed, which should never happen", ret_hdr_line);
    return 0;
  }
  return tokenize_fields();
}

int32_t dsv_hdr_reader::tokenize_fields() {
  tlr.buffer[tlr.last_line_length-1] = '\0';

  char* cursor = tlr.buffer;
  char* token = NULL;
  nfields = 0;
  while( ( token = strsep(&cursor, delim.c_str())) != NULL ) {
    if ( fields.size() <= nfields )
      fields.push_back(token-tlr.buffer);
    else
      fields[nfields] = token-tlr.buffer;
    ++nfields;
  }
  return nfields;
}

int32_t dsv_hdr_reader::store_to_vector(std::vector<std::string>& v) {
  v.resize(nfields);
  for(int32_t i=0; i < nfields; ++i) {
    v[i].assign(str_field_at(i));
  }
  return nfields;
}
    
  
    
