
#include <iostream>
#include <sstream>
#include "thirdparty/murmurhash/MurmurHash2.h"
#include "utility/gutil/jsonwriter.hpp"
#include "utility/gutil/gtimelib.hpp"
#include "core/gpe/gapi4/graphapi.hpp"
#include "core/gpe/gpelib4/enginedriver/engineservicerequest.hpp"
#include "thirdparty/boost/lexical_cast.hpp"


#ifndef _SORT_RPI_CPI_MATRIX_
#define _SORT_RPI_CPI_MATRIX_
namespace UDIMPL {
class sort_rpi_cpi_matrix {
public:

  int64_t key1;
  int64_t cpi;
  double value1;

  sort_rpi_cpi_matrix():
  key1(),
  cpi(),
  value1(){
  }

  sort_rpi_cpi_matrix(int64_t key1_, int64_t cpi_, double value1_){
    key1 = key1_;
    cpi = cpi_;
    value1 = value1_;
  }

 bool operator==(sort_rpi_cpi_matrix const &other) const {
    return key1 == other.key1
        && cpi == other.cpi
        && value1 == other.value1
    ;
  }

  friend std::size_t hash_value(const sort_rpi_cpi_matrix& other){
    std::stringstream os;
    os    << other.key1
    << other.cpi
    << other.value1
    ;
    std::string s = os.str();
    return MurmurHash64A(s.c_str(), s.size(), 0);
  }

 sort_rpi_cpi_matrix&  operator += (const sort_rpi_cpi_matrix& other){
    this->key1 += other.key1;
    this->cpi += other.cpi;
    this->value1 += other.value1;
    return *this;
  }

 bool operator<(sort_rpi_cpi_matrix const &other) const {
    if (key1 < other.key1) {
      return true;
    } else if (key1 > other.key1) {
      return false;
    }

    if (cpi < other.cpi) {
      return true;
    } else if (cpi > other.cpi) {
      return false;
    }

    if (value1 < other.value1) {
      return true;
    } else if (value1 > other.value1) {
      return false;
    }

    return false;
  }


  void json_printer (gutil::JSONWriter& writer, gpelib4::EngineServiceRequest& _request,
gapi4::UDFGraphAPI* graphAPI, bool verbose = false) const {

    writer.WriteStartObject();
    writer.WriteName("key1");
      writer.WriteInt(key1);
          writer.WriteName("cpi");
      writer.WriteInt(cpi);
          writer.WriteName("value1");
      writer.WriteFloat(value1);
      
    writer.WriteEndObject();
  }
  gutil::JSONWriter& json_write_name (
    gutil::JSONWriter& writer, gpelib4::EngineServiceRequest& _request,
    gapi4::UDFGraphAPI* graphAPI, bool verbose = false) const {
    
    std::string ss = boost::lexical_cast<std::string>(*this);
    return writer.WriteName(ss.c_str());
  }

  friend std::ostream& operator<<(std::ostream& os, const sort_rpi_cpi_matrix& m) {
    std::string tmp;
    os<<"[";
    os<<"key1 "<<m.key1<<"|";
    os<<"cpi "<<m.cpi<<"|";
    os<<"value1 "<<m.value1<<"]";
    return os ;
  }


  template <class ARCHIVE>
   void serialize(ARCHIVE& ar) {

    struct TempTuple {
      int64_t key1;
      int64_t cpi;
      double value1;
    };
    TempTuple tp;
    //initialize tp from this->data
    tp.key1 = key1;
    tp.cpi = cpi;
    tp.value1 = value1;
    //serialization for TempTuple
    ar (tp.key1,tp.cpi,tp.value1);
    //recover this->data from tp
    key1 = tp.key1;
    cpi = tp.cpi;
    value1 = tp.value1;

   }

}__attribute__((__packed__));
}//END namespace UDIMPL
#endif

