#include "Fimex.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"
#include <fimex/CDM.h>
#include <fimex/CDMVariable.h>
#include <fimex/GribCDMReader.h>
#include <fimex/CDMFileReaderFactory.h>
#include <fimex/XMLInput.h>
#include <fimex/Data.h>
#include <fimex/CDMconstants.h>
#include <fimex/SliceBuilder.h>
#include <boost/smart_ptr/shared_array.hpp>

FileFimex::FileFimex(std::string iFilename, bool iReadOnly, const Options& iOptions) :
      File(iFilename), mDebug(false) {
   // Turn out debug information
   iOptions.getValue("debug", mDebug);

   // Read config file

   if(mDebug)
      std::cout << "Reading: " << getFilename() << std::endl;

   // Determine filetype
   mifi_filetype fileType = MetNoFimex::CDMFileReaderFactory::detectFileType(getFilename());

   // Create reader
   std::string xmlFilename = "";
   if(iOptions.getValue("config", xmlFilename)) {
      MetNoFimex::XMLInputFile configXml(xmlFilename);
      mReader = MetNoFimex::CDMFileReaderFactory::create(fileType,getFilename(), configXml);
   }
   else {
      mReader = MetNoFimex::CDMFileReaderFactory::create(fileType, getFilename());
   }
   mCdm = mReader->getCDM();

   if(mDebug) {
      MetNoFimex::CDM::DimVec dims = mCdm.getDimensions();
      assert(dims.size() > 0);
      for(int d = 0; d < dims.size(); d++) {
         std::cout << dims[d].getName() << " " << dims[d].getLength() << std::endl;
      }
   }

   MetNoFimex::CDM::VarVec variables = mCdm.getVariables();
   if(mDebug) {
      for(int i = 0; i < variables.size(); i++) {
         std::cout << variables[i].getName() << std::endl;
      }
   }

   // Get the dimension sizes and the time vector values of the first variable with lat/lons
   std::string latName, lonName, timeName;
   std::string yName, xName;
   for(int i = 0; i < variables.size(); i++) {
      // Find names of latitude and longitude variables
      bool found = mCdm.getLatitudeLongitude(variables[i].getName(), latName, lonName);
      if(found) {
         xName = mCdm.getHorizontalXAxis(variables[i].getName());
         yName = mCdm.getHorizontalYAxis(variables[i].getName());
         mNLon = mCdm.getDimension(xName).getLength();
         mNLat = mCdm.getDimension(yName).getLength();
         mNEns = 1;
         if(mCdm.hasDimension("ensemble_member")) {
            mNEns = mCdm.getDimension("ensemble_member").getLength();
         }
         timeName = mCdm.getTimeAxis(variables[i].getName());
         mNTime = mCdm.getDimension(timeName).getLength();
         MetNoFimex::DataPtr dataPtr = mReader->getData(timeName);
         boost::shared_array<double> times = dataPtr->asDouble();
         std::vector<double> timesVec(mNTime, 0);
         for(int t = 0; t < mNTime; t++ ){
            timesVec[t] = times[t];
            if(mDebug)
               std::cout << "Found time: " << times[t] << std::endl;
         }
         setTimes(timesVec);
         break;
      }
   }

   // Get the lat/lon values
   assert(latName != "");
   MetNoFimex::CDMVariable latVar = mCdm.getVariable(latName);
   MetNoFimex::CDMVariable lonVar = mCdm.getVariable(lonName);
   boost::shared_array<float> lats = mReader->getData(latName)->asFloat();
   boost::shared_array<float> lons = mReader->getData(lonName)->asFloat();

   // Initialize to missing
   mLats.resize(mNLat);
   mLons.resize(mNLat);
   for(int i = 0; i < mNLat; i++) {
      mLats[i].resize(mNLon, Util::MV);
      mLons[i].resize(mNLon, Util::MV);
   }

   // Assign lat/lon values in internal variables
   MetNoFimex::SliceBuilder slice(mCdm, latName);
   std::vector< std::string > dimNames = slice.getDimensionNames();
   int numDims = dimNames.size();
   if(numDims == 1) {
      for(int i = 0; i < mNLat; i++) {
         for(int j = 0; j < mNLon; j++) {
            mLats[i][j] = lats[i];
            mLons[i][j] = lons[j];
         }
      }
   }
   else if(numDims == 2) {
      for(int i = 0; i < mNLat; i++) {
         slice.setStartAndSize(yName, i, 1);
         boost::shared_array<double> latData = mReader->getDataSlice(latName, slice)->asDouble();
         boost::shared_array<double> lonData = mReader->getDataSlice(lonName, slice)->asDouble();
         for(int j = 0; j < mNLon; j++) {
            mLats[i][j] = latData[j];
            mLons[i][j] = lonData[j];
         }
      }
   }
   else {
      Util::error("Cannot handle case where latitude variable has more than 2 dimensions");
   }

   if(mDebug) {
      for(int i = 0; i < mNLat; i++) {
         for(int j = 0; j < mNLon; j++) {
            std::cout << mLats[i][j] << " " << mLons[i][j] << std::endl;
         }
      }
      std::cout << "Dimensions: " << mNLat << " " << mNLon << " " << mNTime << " " << mNEns << std::endl;
   }
}

FileFimex::~FileFimex() {
}

// Get a slice for the current timestep
FieldPtr FileFimex::getFieldCore(Variable::Type iVariable, int iTime) const {
   if(mDebug)
      std::cout << "Retrieving time " << iTime << std::endl;
   double s = Util::clock();

   std::string varName = getVariableName(iVariable);

   // Dimension names
   std::string ensembleName = "ensemble_member";
   std::string timeName = "time";
   std::string xName = mCdm.getHorizontalXAxis(varName);
   std::string yName = mCdm.getHorizontalYAxis(varName);

   // Slice the current lead time
   // Slicing strategy: We will slice the time and ensemble dimension, but not the x dimension as
   // this has a significant performance penalty (about 3x slower). We will therefore need to know
   // the order of the x and y dimensions
   MetNoFimex::SliceBuilder slice(mCdm, varName);
   slice.setStartAndSize("time", iTime, 1);

   // Determine if we have an ensemble dimension
   MetNoFimex::CDMVariable var = mCdm.getVariable(varName);
   std::vector<std::string> dims = var.getShape();
   bool hasEnsembleDim = false;
   for(int i = 0; i < dims.size(); i++) {
      if(dims[i] == ensembleName) {
         hasEnsembleDim = true;
         break;
      }
   }

   // Determine if the X or Y dimension varies fastests in the data array
   bool doesXVaryFastest = true;
   bool found = false;
   for(int i = 0; i < dims.size(); i++) {
      if(dims[i] == yName) {
         doesXVaryFastest = false;
         found = true;
         break;
      }
      else if(dims[i] == xName) {
         doesXVaryFastest = true;
         found = true;
         break;
      }
   }
   assert(found);

   // Retrieve values and assign into internal data structure
   FieldPtr field = getEmptyField();
   for(int e = 0; e < mNEns; e++) {
      if(hasEnsembleDim)
         slice.setStartAndSize(ensembleName, e, 1);
      boost::shared_array<float> values = mReader->getDataSlice(varName, slice)->asFloat();
      if(doesXVaryFastest) {
         for(int j = 0; j < mNLon; j++) {
            for(int i = 0; i < mNLat; i++) {
               int index = j * mNLat + i;
               (*field)(i,j,e) = values[index];
            }
         }
      }
      else {
         for(int i = 0; i < mNLat; i++) {
            for(int j = 0; j < mNLon; j++) {
               int index = i * mNLon + j;
               (*field)(i,j,e) = values[index];
            }
         }
      }
   }
   double e = Util::clock();
   if(mDebug)
      std::cout << "Time: " << e-s << " seconds" << std::endl;
   return field;
}
void FileFimex::writeCore(std::vector<Variable::Type> iVariables) {
   Util::error("Cannot write fimex file. Not yet implemented.");
   return;
}
bool FileFimex::hasVariableCore(Variable::Type iVariable) const {
   bool status = mCdm.hasVariable(getVariableName(iVariable));
   return status;
}

std::string FileFimex::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=fimex", "Any file understood by the fimex library") << std::endl;
   ss << Util::formatDescription("   config=", "XML configuration file used by FIMEX to decode the files") << std::endl;
   ss << Util::formatDescription("   debug=0", "If 1, turn on debugging information") << std::endl;
   return ss.str();
}

std::string FileFimex::getVariableName(Variable::Type iVariable) const {
   if(iVariable == Variable::PrecipAcc) {
      return "precipitation_amount_acc";
   }
   else if(iVariable == Variable::Cloud) {
      return "cloud_area_fraction";
   }
   else if(iVariable == Variable::T) {
      return "air_temperature_2m";
   }
   else if(iVariable == Variable::Precip) {
      return "precipitation_amount";
   }
   else if(iVariable == Variable::Pop) {
      return "precipitation_amount_prob_low";
   }
   else if(iVariable == Variable::Pop6h) {
      return "precipitation_amount_prob_low_6h";
   }
   else if(iVariable == Variable::PrecipLow) {
      return "precipitation_amount_low_estimate";
   }
   else if(iVariable == Variable::PrecipMiddle) {
      return "precipitation_amount_middle_estimate";
   }
   else if(iVariable == Variable::PrecipHigh) {
      return "precipitation_amount_high_estimate";
   }
   else if(iVariable == Variable::U) {
      return "x_wind_10m";
   }
   else if(iVariable == Variable::V) {
      return "y_wind_10m";
   }
   else if(iVariable == Variable::W) {
      // TODO: Correct name?
      return "windspeed_10m";
   }
   else if(iVariable == Variable::WD) {
      // TODO: Correct name?
      return "winddirection_10m";
   }
   else if(iVariable == Variable::RH) {
      return "relative_humidity_2m";
   }
   else if(iVariable == Variable::Phase) {
      // TODO: Correct name?
      return "phase";
   }
   else if(iVariable == Variable::P) {
      return "surface_air_pressure";
   }
   else if(iVariable == Variable::MSLP) {
      return "air_pressure_at_sea_level";
   }
   else if(iVariable == Variable::QNH) {
      // TODO: What name to use?
      return "qnh";
   }
   else if(iVariable == Variable::SwinAcc) {
      return "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time";
   }
   else if(iVariable == Variable::LwinAcc) {
      return "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time";
   }
   else if(iVariable == Variable::Fake) {
      return "fake";
   }
   else {
      // TODO:
      Util::error("Variable '" + Variable::getTypeName(iVariable) + "' not defined for FileArome");
   }
   return "";
}
