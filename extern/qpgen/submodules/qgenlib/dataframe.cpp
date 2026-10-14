#include "qgenlib/dataframe.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"

bool dataframe_t::load(const char* tsvfile) {
  tsv_reader tr(tsvfile);

  // read header line;
  ncols = tr.read_line();
  for(int32_t i=0; i < ncols; ++i) {
    colnames.push_back(tr.str_field_at(i));
    if ( col2idx.find(colnames[i]) == col2idx.end() ) 
      col2idx[colnames[i]] = i;
    else
      error("Duplicated column name %s", colnames[i].c_str());      
  }
  nrows = 0;
  columns.resize(ncols);

  while( tr.read_line() ) {
    if ( tr.nfields != ncols )
      error("dataframe_t expects the same number of columns %d in each row, but observed %d", ncols, tr.nfields);
    for(int32_t i=0; i < ncols; ++i) {
      columns[i].push_back(tr.str_field_at(i));
    }
    ++nrows;
  }

  tr.close();
  return true;
}

int32_t dataframe_t::add_empty_column(const char* colname) {
  colnames.push_back(colname);
  columns.resize(ncols+1);
  columns[ncols].resize(nrows);
  if ( col2idx.find(colnames[ncols]) == col2idx.end() )
    col2idx[colnames[ncols]] = ncols;
  else
    error("Duplicated column name %s", colnames[ncols].c_str());
  ++ncols;
  return ncols-1;
}

int32_t dataframe_t::add_empty_row() {
  ++nrows;
  for(int32_t i=0; i < ncols; ++i)
    columns[i].resize(nrows);
  return nrows-1;
}

void dataframe_t::set_str_elem(const char*s, int32_t row, int32_t col) {
  columns[col][row].assign(s);
}

void dataframe_t::set_str_elem(const char*s, int32_t row, const char* colname) {
  get_column(colname).at(row).assign(s);
}


