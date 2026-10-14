#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"

int32_t main(int32_t argc, char** argv) {
  std::string inTSV;
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input TSV", NULL)
    LONG_STRING_PARAM("tsv",&inTSV, "InputTSV file. Plain, gzipped, or bgzipped formats are  allowed")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inTSV.empty() )
    error("[E:%s:%d %s] --tsv parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  //tsv_reader tr(inTSV.c_str());
  dsv_hdr_reader tr(inTSV.c_str());  
  for( int32_t i=0; tr.read_line() > 0; ++i ) {
    for( int32_t j=0; j < tr.nfields; ++j) {
      if (j > 0 ) printf("\t");
      printf("%lf", tr.double_field_at(j));
    }
    printf("\n");
  }
  return 0;
}
