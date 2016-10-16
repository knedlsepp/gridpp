#include "Calibrator.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../Options.h"
#include "../ParameterFile/ParameterFile.h"
#include "../File/File.h"

Calibrator::Calibrator(const Options& iOptions) : Scheme(iOptions) {

}
std::unique_ptr<Calibrator> Calibrator::getScheme(std::string iName, const Options& iOptions) {

   if(iName == "zaga") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'zaga' needs variable");
      }
      return std::make_unique<CalibratorZaga>(Variable::getType(variable), iOptions);
   }
   else if(iName == "cloud") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'cloud' needs variable");
      }
      return std::make_unique<CalibratorCloud>(Variable::getType(variable), iOptions);
   }
   else if(iName == "accumulate") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'accumulate' needs variable");
      }
      return std::make_unique<CalibratorAccumulate>(Variable::getType(variable), iOptions);
   }
   else if(iName == "gaussian") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'gaussian' needs variable");
      }
      return std::make_unique<CalibratorGaussian>(Variable::getType(variable), iOptions);
   }
   else if(iName == "neighbourhood") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'neighbourhood' needs variable");
      }
      return std::make_unique<CalibratorNeighbourhood>(Variable::getType(variable), iOptions);
   }
   else if(iName == "phase") {
      std::unique_ptr<CalibratorPhase> c = std::make_unique<CalibratorPhase>(iOptions);
      float minPrecip;
      if(iOptions.getValue("minPrecip", minPrecip)) {
         c->setMinPrecip(minPrecip);
      }
      bool useWetbulb;
      if(iOptions.getValue("useWetbulb", useWetbulb)) {
         c->setUseWetbulb(useWetbulb);
      }

      return std::move(c);
   }
   else if(iName == "windDirection") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'windDirection' needs variable");
      }
      return std::make_unique<CalibratorWindDirection>(Variable::getType(variable), iOptions);
   }
   else if(iName == "kriging") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'kriging' needs variable");
      }
      return std::make_unique<CalibratorKriging>(Variable::getType(variable), iOptions);
   }
   else if(iName == "diagnose") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'diagnose' needs variable");
      }
      return std::make_unique<CalibratorDiagnose>(Variable::getType(variable), iOptions);
   }
   else if(iName == "window") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'window' needs variable");
      }
      return std::make_unique<CalibratorWindow>(Variable::getType(variable), iOptions);
   }
   else if(iName == "qc") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'qc' needs variable");
      }
      return std::make_unique<CalibratorQc>(Variable::getType(variable), iOptions);
   }
   else if(iName == "qnh") {
      return std::make_unique<CalibratorQnh>(iOptions);
   }
   else if(iName == "qq") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'regression' needs variable");
      }
      return std::make_unique<CalibratorQq>(Variable::getType(variable), iOptions);
   }
   else if(iName == "regression") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'regression' needs variable");
      }
      return std::make_unique<CalibratorRegression>(Variable::getType(variable), iOptions);
   }
   else if(iName == "bct") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'bct' needs variable");
      }
      return std::make_unique<CalibratorBct>(Variable::getType(variable), iOptions);
   }
   else if(iName == "sort") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'sort' needs variable");
      }
      return std::make_unique<CalibratorSort>(Variable::getType(variable), iOptions);
   }
   else if(iName == "altitude") {
      return std::make_unique<CalibratorAltitude>(iOptions);
   }
   else {
      Util::error("Could not instantiate calibrator with name '" + iName + "'");
      return NULL;
   }
}
bool Calibrator::calibrate(File& iFile, const ParameterFile* iParameterFile) const {
   if(requiresParameterFile() && iParameterFile == NULL) {
      std::stringstream ss;
      ss << "Calibrator '" << name() << "' requires a parameter file";
      Util::error(ss.str());
   }
   return calibrateCore(iFile, iParameterFile);
}

void Calibrator::shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter) {
   if(iBefore.size() != iAfter.size()) {
      return;
   }
   if(iBefore.size() == 0)
      return;

   int N = iBefore.size();
   std::vector<std::pair<float,int> > pairs(N);
   for(int e = 0; e < N; e++) {
      if(!Util::isValid(iBefore[e]) || !Util::isValid(iAfter[e])) {
         return;
      }
      pairs[e].first = iBefore[e];
      pairs[e].second = e;
   }
   std::vector<float> afterCopy = iAfter;
   // Sort values so that the rank of a member is the same before and after calibration
   std::sort(pairs.begin(), pairs.end(), Util::sort_pair_first<float,int>());
   std::sort(afterCopy.begin(), afterCopy.end());
   for(int e = 0; e < N; e++) {
      int ei = pairs[e].second;
      float valueCal = afterCopy[e];
      iAfter[ei] = valueCal;
   }
}

Parameters Calibrator::train(const std::vector<ObsEns>& iData) const {
   Util::error("Cannot train method. Not yet implemented.");
   return Parameters();
}

std::string Calibrator::getDescriptions() {
   std::stringstream ss;
   ss << CalibratorAccumulate::description() << std::endl;
   ss << CalibratorAltitude::description() << std::endl;
   ss << CalibratorBct::description() << std::endl;
   ss << CalibratorCloud::description() << std::endl;
   ss << CalibratorDiagnose::description() << std::endl;
   ss << CalibratorGaussian::description() << std::endl;
   ss << CalibratorKriging::description() << std::endl;
   ss << CalibratorNeighbourhood::description() << std::endl;
   ss << CalibratorPhase::description() << std::endl;
   ss << CalibratorQc::description() << std::endl;
   ss << CalibratorQnh::description() << std::endl;
   ss << CalibratorQq::description() << std::endl;
   ss << CalibratorRegression::description() << std::endl;
   ss << CalibratorSort::description() << std::endl;
   ss << CalibratorWindDirection::description() << std::endl;
   ss << CalibratorZaga::description() << std::endl;
   return ss.str();
}
