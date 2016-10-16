#include <boost/scoped_ptr.hpp>
#include <cmath>

#include "Downscaler.h"
#include "../File/File.h"
#include "../KDTree.h"

std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > > Downscaler::mNeighbourCache;

Downscaler::Downscaler(Variable::Type iVariable, const Options& iOptions) : Scheme(iOptions),
      mVariable(iVariable) {
}

bool Downscaler::downscale(const File& iInput, File& iOutput) const {
   if(iInput.getNumTime() != iOutput.getNumTime())
      return false;

   downscaleCore(iInput, iOutput);
   return true;
}

std::unique_ptr<Downscaler> Downscaler::getScheme(std::string iName, Variable::Type iVariable, const Options& iOptions) {
   if(iName == "nearestNeighbour") {
      return std::make_unique<DownscalerNearestNeighbour>(iVariable, iOptions);
   }
   else if(iName == "gradient") {
      return std::make_unique<DownscalerGradient>(iVariable, iOptions);
   }
   else if(iName == "smart") {
      std::unique_ptr<DownscalerSmart> d = std::make_unique<DownscalerSmart>(iVariable, iOptions);
      float searchRadius = Util::MV;
      if(iOptions.getValue("searchRadius", searchRadius)) {
         d->setSearchRadius(searchRadius);
      }
      float numSmart = Util::MV;
      if(iOptions.getValue("numSmart", numSmart)) {
         d->setNumSmart(numSmart);
      }
      float minElevDiff = Util::MV;
      if(iOptions.getValue("minElevDiff", minElevDiff)) {
         d->setMinElevDiff(minElevDiff);
      }
      return std::move(d);
   }
   else if(iName == "bypass") {
      return std::make_unique<DownscalerBypass>(iVariable, iOptions);
   }
   else if(iName == "pressure") {
      return std::make_unique<DownscalerPressure>(iVariable, iOptions);
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iName + "'");
      return NULL;
   }
}

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);

   // Check if the grid is the same
   if(iFrom.getNumLat() == iTo.getNumLat() && iFrom.getNumLon() == iTo.getNumLon()) {
      if(ilats == olats && ilons == olons) {
         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
         Util::status("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            getNearestNeighbourBruteForce(iFrom, olons[i][j], olats[i][j], iI[i][j], iJ[i][j]);
         }
      }
   }
   addToCache(iFrom, iTo, iI, iJ);
}

void Downscaler::getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   // Check if the grid is the same
   if(iFrom.getNumLat() == iTo.getNumLat() && iFrom.getNumLon() == iTo.getNumLon()) {
      vec2 ilats = iFrom.getLats();
      vec2 ilons = iFrom.getLons();
      vec2 olats = iTo.getLats();
      vec2 olons = iTo.getLons();

      if(ilats == olats && ilons == olons) {
         int nLon = iTo.getNumLon();
         int nLat = iTo.getNumLat();

         iI.resize(nLat);
         iJ.resize(nLat);

         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
         Util::status("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   KDTree searchTree(iFrom.getLats(), iFrom.getLons());
   searchTree.getNearestNeighbour(iTo, iI, iJ);

   addToCache(iFrom, iTo, iI, iJ);
}

bool Downscaler::isCached(const File& iFrom, const File& iTo) {
   std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > >::const_iterator it = mNeighbourCache.find(iFrom.getUniqueTag());
   if(it == mNeighbourCache.end()) {
      return false;
   }
   std::map<Uuid, std::pair<vec2Int, vec2Int> >::const_iterator it2 = it->second.find(iTo.getUniqueTag());
   if(it2 == it->second.end()) {
      return false;
   }
   return true;
}

void Downscaler::addToCache(const File& iFrom, const File& iTo, vec2Int iI, vec2Int iJ) {
   std::pair<vec2Int, vec2Int> pair(iI, iJ);
   mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()] = pair;
}
bool Downscaler::getFromCache(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(!isCached(iFrom, iTo))
      return false;
   iI = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].first;
   iJ = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].second;
   return true;
}

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, float iLon, float iLat, int& iI, int &iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();

   iI = Util::MV;
   iJ = Util::MV;

   if(Util::isValid(iLat) && Util::isValid(iLon)) {
      float minDist = Util::MV;
      for(int i = 0; i < ilats.size(); i++) {
         for(int j = 0; j < ilats[0].size(); j++) {
            if(Util::isValid(ilats[i][j]) && Util::isValid(ilons[i][j])) {
               float currDist = Util::getDistance(ilats[i][j], ilons[i][j], iLat, iLon);
               if(!Util::isValid(minDist) || currDist < minDist) {
                  iI = i;
                  iJ = j;
                  minDist = currDist;
               }
            }
         }
      }
   }
}

std::string Downscaler::getDescriptions() {
   std::stringstream ss;
   ss << DownscalerNearestNeighbour::description();
   ss << DownscalerGradient::description();
   ss << DownscalerSmart::description();
   ss << DownscalerPressure::description();
   ss << DownscalerBypass::description();
   return ss.str();
}

void Downscaler::clearCache() {
   mNeighbourCache.clear();
}
