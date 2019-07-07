
#include <iostream>
#include <sstream>
#include "thirdparty/murmurhash/MurmurHash2.h"
#include "utility/gutil/jsonwriter.hpp"
#include "utility/gutil/gtimelib.hpp"
#include "core/gpe/gapi4/graphapi.hpp"
#include "core/gpe/gpelib4/enginedriver/engineservicerequest.hpp"
#include "thirdparty/boost/lexical_cast.hpp"


#ifndef _SORT_RPI_VERTEX_
#define _SORT_RPI_VERTEX_
namespace UDIMPL {
class sort_rpi_vertex {
public:

  int64_t key1;
  int64_t Lp;
  int64_t Up;
  int64_t rp;
  int64_t cpi;
  double row_scaling;
  double col_scaling;

  sort_rpi_vertex():
  key1(),
  Lp(),
  Up(),
  rp(),
  cpi(),
  row_scaling(),
  col_scaling(){
  }

  sort_rpi_vertex(int64_t key1_, int64_t Lp_, int64_t Up_, int64_t rp_, int64_t cpi_, double row_scaling_, double col_scaling_){
    key1 = key1_;
    Lp = Lp_;
    Up = Up_;
    rp = rp_;
    cpi = cpi_;
    row_scaling = row_scaling_;
    col_scaling = col_scaling_;
  }

 bool operator==(sort_rpi_vertex const &other) const {
    return key1 == other.key1
        && Lp == other.Lp
        && Up == other.Up
        && rp == other.rp
        && cpi == other.cpi
        && row_scaling == other.row_scaling
        && col_scaling == other.col_scaling
    ;
  }

  friend std::size_t hash_value(const sort_rpi_vertex& other){
    std::stringstream os;
    os    << other.key1
    << other.Lp
    << other.Up
    << other.rp
    << other.cpi
    << other.row_scaling
    << other.col_scaling
    ;
    std::string s = os.str();
    return MurmurHash64A(s.c_str(), s.size(), 0);
  }

 sort_rpi_vertex&  operator += (const sort_rpi_vertex& other){
    this->key1 += other.key1;
    this->Lp += other.Lp;
    this->Up += other.Up;
    this->rp += other.rp;
    this->cpi += other.cpi;
    this->row_scaling += other.row_scaling;
    this->col_scaling += other.col_scaling;
    return *this;
  }

 bool operator<(sort_rpi_vertex const &other) const {
    if (key1 < other.key1) {
      return true;
    } else if (key1 > other.key1) {
      return false;
    }

    if (Lp < other.Lp) {
      return true;
    } else if (Lp > other.Lp) {
      return false;
    }

    if (Up < other.Up) {
      return true;
    } else if (Up > other.Up) {
      return false;
    }

    if (rp < other.rp) {
      return true;
    } else if (rp > other.rp) {
      return false;
    }

    if (cpi < other.cpi) {
      return true;
    } else if (cpi > other.cpi) {
      return false;
    }

    if (row_scaling < other.row_scaling) {
      return true;
    } else if (row_scaling > other.row_scaling) {
      return false;
    }

    if (col_scaling < other.col_scaling) {
      return true;
    } else if (col_scaling > other.col_scaling) {
      return false;
    }

    return false;
  }


  void json_printer (gutil::JSONWriter& writer, gpelib4::EngineServiceRequest& _request,
gapi4::UDFGraphAPI* graphAPI, bool verbose = false) const {

    writer.WriteStartObject();
    writer.WriteName("key1");
      writer.WriteInt(key1);
          writer.WriteName("Lp");
      writer.WriteInt(Lp);
          writer.WriteName("Up");
      writer.WriteInt(Up);
          writer.WriteName("rp");
      writer.WriteInt(rp);
          writer.WriteName("cpi");
      writer.WriteInt(cpi);
          writer.WriteName("row_scaling");
      writer.WriteFloat(row_scaling);
          writer.WriteName("col_scaling");
      writer.WriteFloat(col_scaling);
      
    writer.WriteEndObject();
  }
  gutil::JSONWriter& json_write_name (
    gutil::JSONWriter& writer, gpelib4::EngineServiceRequest& _request,
    gapi4::UDFGraphAPI* graphAPI, bool verbose = false) const {
    
    std::string ss = boost::lexical_cast<std::string>(*this);
    return writer.WriteName(ss.c_str());
  }

  friend std::ostream& operator<<(std::ostream& os, const sort_rpi_vertex& m) {
    std::string tmp;
    os<<"[";
    os<<"key1 "<<m.key1<<"|";
    os<<"Lp "<<m.Lp<<"|";
    os<<"Up "<<m.Up<<"|";
    os<<"rp "<<m.rp<<"|";
    os<<"cpi "<<m.cpi<<"|";
    os<<"row_scaling "<<m.row_scaling<<"|";
    os<<"col_scaling "<<m.col_scaling<<"]";
    return os ;
  }


  template <class ARCHIVE>
   void serialize(ARCHIVE& ar) {

    struct TempTuple {
      int64_t key1;
      int64_t Lp;
      int64_t Up;
      int64_t rp;
      int64_t cpi;
      double row_scaling;
      double col_scaling;
    };
    TempTuple tp;
    //initialize tp from this->data
    tp.key1 = key1;
    tp.Lp = Lp;
    tp.Up = Up;
    tp.rp = rp;
    tp.cpi = cpi;
    tp.row_scaling = row_scaling;
    tp.col_scaling = col_scaling;
    //serialization for TempTuple
    ar (tp.key1,tp.Lp,tp.Up,tp.rp,tp.cpi,tp.row_scaling,tp.col_scaling);
    //recover this->data from tp
    key1 = tp.key1;
    Lp = tp.Lp;
    Up = tp.Up;
    rp = tp.rp;
    cpi = tp.cpi;
    row_scaling = tp.row_scaling;
    col_scaling = tp.col_scaling;

   }

}__attribute__((__packed__));
}//END namespace UDIMPL
#endif

