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
      File(iFilename) {

   std::string filenameFormat;
   iOptions.getValue("filenameFormat", filenameFormat)
   mFiles = Util::glob(getFilename());
   if(mFiles.size() == 0) {
      std::stringstream ss;
      ss << "Found no files matching " << getFilename();
      Util::error(ss.str());
   }

   std::vector<std::pair<std::string, std::string> > members;
   std::vector<std::string> files;
   std::vector<std::string> args;
   // args.push_back("/disk1/mbr000.grb");
   args.push_back("memberName:mbr000");
   args.push_back("memberName:mbr001");
   /*
   MetNoFimex::CDMFileReaderFactory::parseGribArgs(args, members, files);
   for(int i = 0; i < members.size(); i++) {
      std::cout << "Member: " << members[i].first << " " << members[i].second << std::endl;
   }
   for(int i = 0; i < files.size(); i++) {
      std::cout << "Files: " << files[i] << std::endl;
   }
   */

   std::string xmlFilename = "/usr/share/fimex/cdmGribReaderConfig.xml";
   iOptions.getValue("config", xmlFilename);
   MetNoFimex::XMLInputFile configXml(xmlFilename);
   // mReader = new MetNoFimex::GribCDMReader(files, configXml, members);
   mifi_filetype fileType = MetNoFimex::CDMFileReaderFactory::detectFileType("/disk1/mbr000.grb");

   std::cout << fileType << std::endl;
   mReader = MetNoFimex::CDMFileReaderFactory::create(fileType, "glob:/disk1/mbr00*.grb", configXml, args);
   mCdm = mReader->getCDM();

   MetNoFimex::CDM::DimVec dims = mCdm.getDimensions();
   assert(dims.size() > 0);
   for(int d = 0; d < dims.size(); d++) {
      std::cout << dims[d].getName() << " " << dims[d].getLength() << std::endl;
   }

   std::string xdim, ydim;

   // Find names of latitude and longitude
   std::string latitudeName, longitudeName, timeName;
   MetNoFimex::CDM::VarVec variables = mCdm.getVariables();
   for(int i = 0; i < variables.size(); i++) {
      std::cout << variables[i].getName() << std::endl;
   }
   for(int i = 0; i < variables.size(); i++) {
      bool status = mCdm.getLatitudeLongitude(variables[i].getName(), latitudeName, longitudeName);
      // std::cout << status << " " << variables[i].getName() << " " << latitudeName << std::endl;
      if(status) {
         mNLon = mCdm.getDimension(mCdm.getHorizontalXAxis(variables[i].getName())).getLength();
         mNLat = mCdm.getDimension(mCdm.getHorizontalYAxis(variables[i].getName())).getLength();
         mNEns = mCdm.getDimension("ensemble_member").getLength();
         timeName = mCdm.getTimeAxis(variables[i].getName());
         mNTime = mCdm.getDimension(timeName).getLength();
         MetNoFimex::CDMVariable timeVar = mCdm.getVariable(timeName);
         boost::shared_array<double> times = timeVar.getData()->asDouble();
         std::vector<double> timesVec(mNTime, 0);
         for(int t = 0; t < mNTime; t++ ){
            timesVec[t] = times[t];
            std::cout << "Found time: " << times[t] << std::endl;
         }
         setTimes(timesVec);
         break;
      }
   }
   assert(latitudeName != "");
   MetNoFimex::CDMVariable latVar = mCdm.getVariable(latitudeName);
   MetNoFimex::CDMVariable lonVar = mCdm.getVariable(longitudeName);
   boost::shared_array<float> lats = latVar.getData()->asFloat();
   boost::shared_array<float> lons = lonVar.getData()->asFloat();
   mLats.resize(mNLat);
   mLons.resize(mNLat);
   int index = 0;
   for(int i = 0; i < mNLat; i++) {
      mLats[i].resize(mNLon);
      mLons[i].resize(mNLon);
      for(int j = 0; j < mNLon; j++) {
         int index = i + j*mNLat;
         mLats[i][j] = lats[index];
         mLons[i][j] = lons[index];
      }
   }

   // Assume for now that each file is a member
   // mNEns = 1;//mFiles.size();

   /*
   for(int i = 0; i < mNLat; i++) {
      for(int j = 0; j < mNLon; j++) {
         std::cout << i << " " << j << " " << mLats[i][j] << " " << mLons[i][j] << std::endl;
      }
   }
   */

   std::cout << "Dimensions: " << mNLat << " " << mNLon << " " << mNTime << " " << mNEns << std::endl;
}

FileFimex::~FileFimex() {
   // delete mReader;
}

FieldPtr FileFimex::getFieldCore(Variable::Type iVariable, int iTime) const {

   std::string varName = getVariableName(iVariable);
   MetNoFimex::SliceBuilder slice(mCdm, varName);
   MetNoFimex::DataPtr data = mReader->getDataSlice(varName, slice);
   // MetNoFimex::CDF::DimVec dims = mReader->getDimensions(varName);
   MetNoFimex::CDMVariable var = mCdm.getVariable(varName);
   std::vector<std::string> dims = var.getShape();
   for(int i = 0; i < dims.size(); i++) {
      std::cout << "shape: " << dims[i] << std::endl;
   }
   FieldPtr field = getEmptyField();
   boost::shared_array<float> values = data->asFloat();
   for(int i = 0; i < mNLat; i++) {
      for(int j = 0; j < mNLon; j++) {
         for(int e = 0; e < mNEns; e++) {
            int index = i + j*mNLat + e*mNLat*mNLon;//i*mNTime + j*mNLat*mNTime;
            (*field)(i,j,e) = values[index];
         }
      }
   }
   return field;
}
void FileFimex::writeCore(std::vector<Variable::Type> iVariables) {
   return;
}
bool FileFimex::hasVariableCore(Variable::Type iVariable) const {
   bool status = mCdm.hasVariable(getVariableName(iVariable));
   return status;
}

std::string FileFimex::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=grib", "GRIB file") << std::endl;
   ss << Util::formatDescription("", "Multiple GRIB files can be specified using wildcards, but the filenames must be put in quotes. Each file is assumed to be a different ensemble member.") << std::endl;
   ss << Util::formatDescription("   map=required", "A file containing the translation of GRIB variable names (shortName_levType_levelist) to gridpp names (those avilable with -v). Must be a text file as follows:") << std::endl;
   ss << Util::formatDescription("", "T t_sfc_2") << std::endl;
   ss << Util::formatDescription("", "RH rh_sfc_2") << std::endl;
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
