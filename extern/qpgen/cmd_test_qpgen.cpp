#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qpgen.h"

int32_t test_qpgen(int32_t argc, char** argv) {
  std::string prefix;
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input file", NULL)
    LONG_STRING_PARAM("prefix", &prefix, "Input PLINK file prefix")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( prefix.empty() )
    error("[E:%s:%d %s] --prefix parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  PlinkReader pr;

  if ( pr.prep_prefix(prefix.c_str()) ) {
    notice("Successfully prepared the prefix %s", prefix.c_str());
  }
  else {
    error("Failed to prepare the prefix %s", prefix.c_str());
  }

  // print the sample size
  notice("Total sample size = %d", pr.samps.size());

  // read the first variant and genotype
  if ( pr.stream_genos() ) {
    notice("Successfully read the first variant");
    std::map<int32_t, uint32_t> cntmap;
    for(int32_t i=0; i < pr.int_buf.size(); ++i) {
      cntmap[pr.int_buf[i]]++;
    }
    for(std::map<int32_t, uint32_t>::iterator it = cntmap.begin(); it != cntmap.end(); ++it) {
      notice("Allele %d has count %u", it->first, it->second);
    }
  }
  else {
    error("Failed to read the first variant");
  }



  // PgenReader pgr;
  // int32_t nsamps = 500;
  // std::vector<int32_t> idx;
  // for(int32_t i=0; i < nsamps; ++i)
  //   idx.push_back(i+1);
  // pgr.Load(pgenf, nsamps, idx, 1);
  // //pgr.Load(pgenf, 0, idx, 0);

  // // extract a variant based on rsID
  // uint32_t variant_cnt = pgr.GetVariantCt();
  // notice("variant count is %u", variant_cnt);

  // uint32_t variant_idx = 0;

  // // get the allele count for this variant
  // uint32_t allele_ct = pgr.GetAlleleCt(variant_idx);

  // notice("allele count for variant %s is %d", rsid.c_str(), allele_ct);

  // // read the first ALT allele
  // std::vector<double> buf(1);
  // pgr.Read(buf.data(), 1, 0, variant_idx, 1);

  // for(int32_t i=0; i < buf.size(); ++i )
  //   printf("%lf\n", buf[i]);


  // //tsv_reader tr(inTSV.c_str());
  // dsv_hdr_reader tr(inTSV.c_str());  
  // for( int32_t i=0; tr.read_line() > 0; ++i ) {
  //   for( int32_t j=0; j < tr.nfields; ++j) {
  //     if (j > 0 ) printf("\t");
  //     printf("%lf", tr.double_field_at(j));
  //   }
  //   printf("\n");
  // }
  return 0;
}
