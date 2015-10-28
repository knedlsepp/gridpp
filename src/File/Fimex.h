#ifndef FILE_FIMEX_H
#define FILE_FIMEX_H
#include <grib_api.h>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "File.h"
#include "../Variable.h"
#include <fimex/CDM.h>
#include <fimex/CDMReader.h>

class FileFimex : public File {
   public:
      FileFimex(std::string iFilename, bool iReadOnly, const Options& iOptions);
      ~FileFimex();

      static std::string description();
      std::string name() const {return "fimex";};
      std::string getVariableName(Variable::Type iVariable) const;
   protected:
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      void writeCore(std::vector<Variable::Type> iVariables);
      bool hasVariableCore(Variable::Type iVariable) const;
   private:
      std::vector<std::string> mFiles;
      MetNoFimex::CDM mCdm;
      boost::shared_ptr<MetNoFimex::CDMReader> mReader;
};
#endif
