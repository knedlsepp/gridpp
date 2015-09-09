#include "Grib.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"

FileGrib::FileGrib(std::string iFilename, bool iReadOnly, const Options& iOptions) :
      File(iFilename) {
   mFiles = Util::glob(getFilename());

   if(mFiles.size() == 0) {
      std::stringstream ss;
      ss << "Found no files matching " << getFilename();
      Util::error(ss.str());
   }

   // Retrieve variable map
   std::string mapFile;
   if(!iOptions.getValue("map", mapFile)) {
      Util::error("Grib files require 'map' option");
   }
   std::ifstream ifs(mapFile.c_str());
   if(ifs.good()) {
      while(ifs.good()) {
         std::string gridppName, gribName;
         ifs >> gridppName >> gribName;
         mVariableMap[gridppName] = gribName;
      }
   }
   else {
      Util::error("Could not open file '" + mapFile + "'");
   }

   // Locations
   setLatLons(mFiles[0]);

   // Times
   mNTime = 1;
   std::vector<double> times(1,0);
   setTimes(times);

   // Assume for now that each file is a member
   mNEns = mFiles.size();
}

FileGrib::~FileGrib() {
}

FieldPtr FileGrib::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string gridppName = Variable::getTypeName(iVariable);
   VariableMap::const_iterator it = mVariableMap.find(gridppName);
   if(it == mVariableMap.end()) {
      std::stringstream ss;
      ss << "File does not contain variable '" << gridppName << "'";
      Util::error(ss.str());
   }
   std::string gribName = it->second;

   FieldPtr field = getEmptyField();
   for(int m = 0; m < mFiles.size(); m++) {
      FILE* fid = fopen(mFiles[m].c_str(), "r");
      bool foundVariable = false;
      if(fid) {
         // Loop over all messages
         while(true) {
            int err = 0;
            grib_handle* h = grib_handle_new_from_file(0,fid,&err);
            if(h == NULL)
               break; // No more messages to process

            std::string currVariable = getVariableName(h);
            std::cout << "Found " << currVariable << std::endl;
            if(currVariable == gribName) {
               std::cout << "Assuming that " << currVariable << " is " << gridppName << std::endl;
               foundVariable = true;
               size_t N;
               GRIB_CHECK(grib_get_size(h,"values",&N),0);
               assert(N == mNLat * mNLon);

               double* arr = new double[N];
               GRIB_CHECK(grib_get_double_array(h,"values",arr,&N),0);

               int index = 0;
               for(int i = 0; i < mNLat; i++) {
                  for(int j = 0; j < mNLon; j++) {
                     (*field)(i,j,m) = arr[index];
                     index++;
                  }
               }
               break;
            }
            if(h) {
               grib_handle_delete(h);
            }
         }
      }
      else {
         std::stringstream ss;
         ss << "Could not open " << mFiles[m];
         Util::error(ss.str());
      }

      if(!foundVariable) {
         std::stringstream ss;
         ss << "Could not find variable " << Variable::getTypeName(iVariable);
         Util::error(ss.str());
      }
      fclose(fid);
   }
   return field;
}
void FileGrib::writeCore(std::vector<Variable::Type> iVariables) {
   return;
}
bool FileGrib::hasVariableCore(Variable::Type iVariable) const {
   return true;
}

/*
bool FileGrib::isValid(std::string iFilename) {
   std::vector<std::string> files = Util::glob(iFilename);
   for(int m = 0; m < files.size(); m++) {
      FILE* fid = fopen(files[m].c_str(), "r");
      if(!fid) {
         return false;
      }
      fclose(fid);
   }
   return true;
}
*/

void FileGrib::setLatLons(std::string iFilename) {
   // Get the locations from the first file
   FILE* fid = fopen(iFilename.c_str(),"r");
   if(!fid) {
      Util::error("Could not open " + iFilename);
   }
   else {
      int err = 0;
      grib_handle* h = grib_handle_new_from_file(0,fid,&err);
      if (h == NULL) {
         Util::error("Could not create handle for " + iFilename);
      }

      // Iterate over all locations
      long nLon, nLat;
      GRIB_CHECK(grib_get_long(h,"Ni",&nLon),0);
      mNLon = nLon;
      GRIB_CHECK(grib_get_long(h,"Nj",&nLat),0);
      mNLat = nLat;

      grib_iterator* iter=grib_iterator_new(h,0,&err);                                     
      if (err != GRIB_SUCCESS) GRIB_CHECK(err,0);                       
      double lat, lon;
      double value;

      mLats.resize(mNLat);
      mLons.resize(mNLat);
      mElevs.resize(mNLat);

      std::cout << mNLat << " " << mNLon << std::endl;

      for(int i = 0; i < mNLat; i++) {
         mLats[i].resize(mNLon);
         mLons[i].resize(mNLon);
         mElevs[i].resize(mNLon);
         for(int j = 0; j < mNLon; j++) {
            int status = grib_iterator_next(iter,&lat,&lon,&value);
            if(status == 0)
               Util::error("Bad status");
            assert(lat <= 90 && lat >= -90);
            mLats[i][j] = lat;
            mLons[i][j] = lon;
            // std::cout << lat << " " << lon << std::endl;
         }
      }
      grib_iterator_delete(iter);
      grib_handle_delete(h);
   }
   fclose(fid);
}

std::string FileGrib::getVariableName(grib_handle* iH) {
   std::stringstream ss;
   char name1[512];
   char name2[512];
   char name3[512];
   char name4[512];
   size_t len1 = 512;
   size_t len2 = 512;
   size_t len3 = 512;
   size_t len4 = 512;

   int err = grib_get_string(iH, "shortName", name3, &len3); 
   assert(err == 0);
   std::string shortName = std::string(name3);

   err = grib_get_string(iH, "typeOfLevel", name1, &len1); 
   assert(err == 0);
   std::string typeOfLevel = std::string(name1);

   err = grib_get_string(iH, "levtype", name2, &len2); 
   assert(err == 0);
   std::string levtype = std::string(name2);
   if(typeOfLevel == "nominalTop") {
      levtype = "ntop";
   }

   err = grib_get_string(iH, "levelist", name4, &len4); 
   std::string level = "undef";
   if(err != GRIB_NOT_FOUND) {
      level = std::string(name4);
   }

   ss << shortName << "_" << levtype << "_" << level;
   return ss.str();
}

float FileGrib::getOffset(grib_handle* iH) {
   double offset;
   int err = grib_get_double(iH, "stepRange", &offset); 
   assert(err == 0);
   return (float) offset;
}

std::string FileGrib::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=grib", "GRIB file") << std::endl;
   ss << Util::formatDescription("", "Multiple GRIB files can be specified using wildcards, but the filenames must be put in quotes. Each file is assumed to be a different ensemble member.") << std::endl;
   ss << Util::formatDescription("   map=required", "A file containing the translation of GRIB variable names (shortName_levType_levelist) to gridpp names (those avilable with -v). Must be a text file as follows:") << std::endl;
   ss << Util::formatDescription("", "T t_sfc_2") << std::endl;
   ss << Util::formatDescription("", "RH rh_sfc_2") << std::endl;
   return ss.str();
}
