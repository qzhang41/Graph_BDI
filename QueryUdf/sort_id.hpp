
#include <iostream>
#include <sstream>
#include "thirdparty/murmurhash/MurmurHash2.h"
#include "utility/gutil/jsonwriter.hpp"
#include "utility/gutil/gtimelib.hpp"
#include "core/gpe/gapi4/graphapi.hpp"
#include "core/gpe/gpelib4/enginedriver/engineservicerequest.hpp"
#include "thirdparty/boost/lexical_cast.hpp"


#ifndef _SORT_ID_
#define _SORT_ID_
namespace UDIMPL {
class sort_id {
public:

  int64_t index1;
  double value1;

  sort_id():
  index1(),
  value1(){
  }

  sort_id(int64_t index1_, double value1_){
    index1 = index1_;
    value1 = value1_;
  }

 bool operator==(sort_id const &other) const {
    return index1 == other.index1
        && value1 == other.value1
    ;
  }

  friend std::size_t hash_value(const sort_id& other){
    std::stringstream os;
    os    << other.index1
    << other.value1
    ;
    std::string s = os.str();
    return MurmurHash64A(s.c_str(), s.size(), 0);
  }

 sort_id&  operator += (const sort_id& other){
    this->index1 += other.index1;
    this->value1 += other.value1;
    return *this;
  }

 bool operator<(sort_id const &other) const {
    if (index1 < other.index1) {
      return true;
    } else if (index1 > other.index1) {
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
    writer.WriteName("index1");
      writer.WriteInt(index1);
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

  friend std::ostream& operator<<(std::ostream& os, const sort_id& m) {
    std::string tmp;
    os<<"[";
    os<<"index1 "<<m.index1<<"|";
    os<<"value1 "<<m.value1<<"]";
    return os ;
  }


  template <class ARCHIVE>
   void serialize(ARCHIVE& ar) {

    struct TempTuple {
      int64_t index1;
      double value1;
    };
    TempTuple tp;
    //initialize tp from this->data
    tp.index1 = index1;
    tp.value1 = value1;
    //serialization for TempTuple
    ar (tp.index1,tp.value1);
    //recover this->data from tp
    index1 = tp.index1;
    value1 = tp.value1;

   }

}__attribute__((__packed__));
}//END namespace UDIMPL
#endif

