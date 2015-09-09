#ifndef FILE_GRIB_H
#define FILE_GRIB_H
#include <grib_api.h>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "File.h"
#include "../Variable.h"

//! Represents a Netcdf data file
class FileGrib : public File {
   public:
      // Can be a glob
      FileGrib(std::string iFilename, bool iReadOnly, const Options& iOptions);
      ~FileGrib();

      static std::string description();
      std::string name() const {return "grib";};
      //static bool isValid(std::string iFilename);
   protected:
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      void writeCore(std::vector<Variable::Type> iVariables);
      bool hasVariableCore(Variable::Type iVariable) const;

      void setLatLons(std::string iFilename);
      int getTime(std::string iFilename);
      //typedef std::pair<std::string, int> VarTime; // variable, time
      //std::map<VarTime, std::vector<std::string> > mFiles; // varTime, filenames 
      std::vector<std::string> mFiles;

      static std::string getVariableName(grib_handle* iH);
      static float getOffset(grib_handle* iH);

      typedef std::map<std::string,std::string> VariableMap; // gridppName, gribName
      VariableMap mVariableMap;
};
#endif
