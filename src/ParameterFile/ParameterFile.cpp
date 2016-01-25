#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFile::ParameterFile(const Options& iOptions, bool iIsNew) :
      Scheme(iOptions),
      mFilename(""),
      mIsNew(iIsNew) {
   iOptions.getValue("file", mFilename);
}

ParameterFile* ParameterFile::getScheme(std::string iName, const Options& iOptions, bool iIsNew) {
   ParameterFile* p;
   if(iName == "metnoKalman") {
      p = new ParameterFileMetnoKalman(iOptions, iIsNew);
   }
   else if(iName == "text") {
      p = new ParameterFileText(iOptions, iIsNew);
   }
   else if(iName == "netcdf") {
      p = new ParameterFileNetcdf(iOptions, iIsNew);
   }
   else {
      Util::error("Parameter file type '" + iName + "' not recognized");
   }
   return p;
}

Parameters ParameterFile::getParameters(int iTime) const {
   // Find the right location to use
   std::map<Location, std::vector<Parameters> >::const_iterator it;
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
       it = mParameters.begin();
   }

   // Find the right time to use
   const std::vector<Parameters>& timeParameters = it->second;
   if(timeParameters.size() == 1) {
      // One set of parameters for all times
      return timeParameters[0];
   }
   else if(timeParameters.size() > iTime) {
      return timeParameters[iTime];
   }
   else {
      std::stringstream ss;
      ss << "Parameter file '" << mFilename << "' does not have values for time " << iTime << ".";
      Util::error(ss.str());
   }
}
Parameters ParameterFile::getParameters(int iTime, const Location& iLocation) const {
   if(mParameters.size() == 0)
      return Parameters();
   // Find the right location to use
   Location loc = getNearestLocation(iTime, iLocation);

   // Find the right time to use
   std::map<Location, std::vector<Parameters> >::const_iterator it = mParameters.find(loc);
   const std::vector<Parameters>& timeParameters = it->second;
   if(timeParameters.size() == 1) {
      // One set of parameters for all times
      return timeParameters[0];
   }
   else if(timeParameters.size() > iTime) {
      return timeParameters[iTime];
   }
   else {
      std::stringstream ss;
      ss << "Parameter file '" << mFilename << "' does not have values for time " << iTime << ".";
      Util::error(ss.str());
   }
}

Location ParameterFile::getNearestLocation(int iTime, const Location& iLocation) const {
   Location loc(Util::MV,Util::MV,Util::MV);
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
      std::map<Location, std::vector<Parameters> >::const_iterator it = mParameters.begin();
      loc = it->first;
   }
   else {
      std::map<Location, std::vector<Parameters> >::const_iterator it;
      it = mParameters.find(iLocation);
      // Nearest neighbour
      // TODO
      if(it == mParameters.end()) {
         std::map<Location, std::vector<Parameters> >::const_iterator it;
         float minDist = Util::MV;
         for(it = mParameters.begin(); it != mParameters.end(); it++) {
            const Location& currLoc = it->first;
            float dist = iLocation.getDistance(currLoc);
            if(!Util::isValid(minDist) || (Util::isValid(dist) && dist < minDist)) {
               // Check that this location actually has parameters available for this time
               bool hasAtThisTime = it->second.size() > iTime && it->second[iTime].size() != 0;
               if(hasAtThisTime) {
                  minDist = dist;
                  loc = it->first;
               }
            }
         }
         if(!Util::isValid(minDist)) {
            // TODO: No nearby location found
            Util::error("No parameter locations found");
         }
         else {
         }
      }
      else {
         loc = it->first;
      }
   }
   return loc;
}

void ParameterFile::setParameters(Parameters iParameters, int iTime, const Location& iLocation) {
   std::map<Location, std::vector<Parameters> >::const_iterator it = mParameters.find(iLocation);
   if(mParameters[iLocation].size() <= iTime) {
      mParameters[iLocation].resize(iTime+1);
   }
   mParameters[iLocation][iTime] = iParameters;
}
void ParameterFile::setParameters(Parameters iParameters, int iTime) {
   std::map<Location, std::vector<Parameters> >::iterator it = mParameters.begin();
   it->second[iTime] = iParameters;
}

std::string ParameterFile::getFilename() const {
   return mFilename;
}

std::vector<Location> ParameterFile::getLocations() const {
   std::vector<Location> locations;
   std::map<Location, std::vector<Parameters> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      locations.push_back(it->first);
   }
   return locations;
}

std::vector<int> ParameterFile::getTimes() const {
   std::set<int> times;
   std::map<Location, std::vector<Parameters> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      for(int i = 0; i < it->second.size(); i++) {
         if(it->second[i].size() != 0)
            times.insert(i);
      }
   }
   return std::vector<int>(times.begin(), times.end());
}

bool ParameterFile::isLocationDependent() const {
   return mParameters.size() > 1;
}
// TODO
bool ParameterFile::isTimeDependent() const {
   return true;
}

int ParameterFile::getNumParameters() const {
   int size = Util::MV;

   std::vector<Location> locations;
   std::map<Location, std::vector<Parameters> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      for(int i = 0; i < it->second.size(); i++) {
         int currSize = it->second[i].size();
         if(currSize != 0) {
            if(Util::isValid(size) && currSize != size)
               return Util::MV;
            size = currSize;
         }
         else
            return Util::MV;
      }
   }
   return size;
}

std::string ParameterFile::getDescription(bool iSpatialOnly) {
   std::stringstream ss;
   if(iSpatialOnly)
      ss << "What file type is the parameters in? One of 'text' and 'metnoKalman'.";
   else
      ss << "What file type is the parameters in? One of 'text' and 'metnoKalman'.";
   return ss.str();
}

void ParameterFile::setFilename(std::string iFilename) {
   mFilename = iFilename;
}

std::string ParameterFile::getDescriptions() {
   std::stringstream ss;
   ss << ParameterFileText::description() << std::endl;
   ss << ParameterFileMetnoKalman::description() << std::endl;
   ss << ParameterFileNetcdf::description() << std::endl;
   return ss.str();
}

long ParameterFile::getCacheSize() const {
   std::map<Location, std::vector<Parameters>, Location::CmpIgnoreElevation >::const_iterator it;
   long total = 0;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      for(int i = 0; i < it->second.size(); i++) {
         total++;
      }
   }
   return total;
}
