#ifndef __DATAFRAME_H
#define __DATAFRAME_H

#include <vector>
#include <cstdint>
#include <string>
#include <map>

class dataframe_t {
public:
  std::vector<std::string> colnames;
  std::map<std::string,uint32_t> col2idx;
  std::vector<std::vector<std::string> > columns;
  int32_t ncols;
  int32_t nrows;

  bool load(const char* tsvfile);

  inline const std::string& get_str_elem(int32_t row, int32_t col) {
    return columns[col][row];
  }
  
  inline int32_t get_int_elem(int32_t row, int32_t col) {
    return strtol(columns[col][row].c_str(), NULL, 10);
  }
  
  inline int64_t get_int64_elem(int32_t row, int32_t col) {
    return strtoll(columns[col][row].c_str(), NULL, 10);
  }
  
  inline uint64_t get_uint64_elem(int32_t row, int32_t col) {
    return strtoull(columns[col][row].c_str(), NULL, 10);
  }

  inline double get_double_elem(int32_t row, int32_t col) {
    return strtod(columns[col][row].c_str(), NULL);
  }

  inline std::vector<std::string>& get_column(int32_t col) {
    return columns[col];
  }
  
  inline std::vector<std::string>& get_column(const char* colname) {
    return columns[col2idx[colname]];
  }

  inline bool has_column(const char* colname) {
    return col2idx.find(colname) != col2idx.end();
  }

  inline int32_t get_colidx(const char* colname) {
    if ( has_column(colname) ) return col2idx[colname];
    else return -1;
  }

  inline const std::string& get_str_elem(int32_t row, const char* colname) {
    return get_column(colname).at(row);
  }  
  
  inline int32_t get_int_elem(int32_t row, const char* colname) {
    return strtol(get_column(colname).at(row).c_str(), NULL, 10);
  }

  inline int64_t get_int64_elem(int32_t row, const char* colname) {
    return strtoll(get_column(colname).at(row).c_str(), NULL, 10);
  }  

  inline uint64_t get_uint64_elem(int32_t row, const char* colname) {
    return strtoull(get_column(colname).at(row).c_str(), NULL, 10);
  }

  inline double get_double_elem(int32_t row, const char* colname) {
    return strtod(get_column(colname).at(row).c_str(), NULL);
  }

  int32_t add_empty_column(const char* colname);
  int32_t add_empty_row();
  void set_str_elem(const char *s, int32_t row, int32_t col);
  void set_str_elem(const char* s, int32_t row, const char* colname);

  dataframe_t() : nrows(0), ncols(0) {} // default constructor does not do anything
  dataframe_t(const char* tsvfile) { load(tsvfile); }
};

#endif
