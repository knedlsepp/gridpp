#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscalerGradient : public ::testing::Test {
      public:
         void SetUp() {
            mFrom = std::make_unique<FileArome>("testing/files/10x10.nc");
            std::stringstream ss;
            ss << "nLat=1 nLon=4 nEns=1 nTime=" << mFrom->getNumTime();
            mTo = std::make_unique<FileFake>(Options(ss.str()));
            setLatLonElev(*mTo, (const float[]) {5}, (const float[]){2,2,12,20}, (const float[]){120, 1500, 600, -100});
         }
         void TearDown() {
         }
         vec2 makeVec2(int nLat, int nLon, const std::vector<float>& values) {
            vec2 grid;
            grid.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               grid[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  int index = i*nLon + j;
                  grid[i][j] = values[index];
               }
            }
            return grid;
         };
         // iElev must be a flat array with lon varying fastest
         void setLatLonElev(FileFake& iFile, const float iLat[], const float iLon[], const float iElev[]) {
            vec2 lat;
            vec2 lon;
            vec2 elev;
            int nLat = iFile.getNumLat(); 
            int nLon = iFile.getNumLon();
            lat.resize(nLat);
            lon.resize(nLat);
            elev.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               elev[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
                  elev[i][j] = iElev[j+i*nLon];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
            iFile.setElevs(elev);
         };
      protected:
         std::unique_ptr<FileArome> mFrom;
         std::unique_ptr<FileFake> mTo;
   };

   TEST_F(TestDownscalerGradient, 10x10) {
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) + gradient * (elev - elev(nn))
         EXPECT_FLOAT_EQ(301.31491, toT(0,0,0)); // 301 - 0.00797 * (120-160)
         EXPECT_FLOAT_EQ(290.34964, toT(0,1,0));
         EXPECT_FLOAT_EQ(301.29544, toT(0,2,0)); // 301 - 0.01068 * (600-346)
         EXPECT_FLOAT_EQ(308.77686, toT(0,3,0));
      }

      // Fix the gradient
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 constantGradient=0.01"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) + gradient * (elev - elev(nn))
         EXPECT_FLOAT_EQ(300.60367, toT(0,0,0));
         EXPECT_FLOAT_EQ(314.40369, toT(0,1,0));
         EXPECT_FLOAT_EQ(306.53052, toT(0,2,0));
         EXPECT_FLOAT_EQ(299.53052, toT(0,3,0));
      }
   }
   // Check that logTransform works
   TEST_F(TestDownscalerGradient, 10x10log) {
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 logTransform=1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) * exp(gradient * (elev - elev(nn)))
         EXPECT_FLOAT_EQ(301.30985, toT(0,0,0)); // 301 * exp(-2.59620e-5 * (120-159.63))
      }

      // Fix the gradient
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 logTransform=1 constantGradient=-0.01"));
         d.downscale(*mFrom, *mTo);
         const Field& toT  = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) * exp(gradient * (elev - elev(nn)))
         EXPECT_FLOAT_EQ(447.3916, toT(0,0,0));
      }
   }
   TEST_F(TestDownscalerGradient, 10x10minElevDiff) {
      {
         DownscalerGradient d(Variable::T ,Options("searchRadius=1 minElevDiff=1500"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) + gradient * (elev - elev(nn))
         EXPECT_FLOAT_EQ(301, toT(0,0,0));
         EXPECT_FLOAT_EQ(301, toT(0,1,0));
         EXPECT_FLOAT_EQ(304, toT(0,2,0));
         EXPECT_FLOAT_EQ(304, toT(0,3,0));
      }

      // Fix the gradient
      // This should not be affected by the minElevDiff
      {
         DownscalerGradient d(Variable::T ,Options("searchRadius=1 minElevDiff=1500 constantGradient=0.01"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // T = T(nn) + gradient * (elev - elev(nn))
         EXPECT_FLOAT_EQ(300.60367, toT(0,0,0));
         EXPECT_FLOAT_EQ(314.40369, toT(0,1,0));
         EXPECT_FLOAT_EQ(306.53052, toT(0,2,0));
         EXPECT_FLOAT_EQ(299.53052, toT(0,3,0));
      }
   }
   TEST_F(TestDownscalerGradient, 10x10negativeTemperatures) {
      setLatLonElev(*mTo, (const float[]) {5}, (const float[]){2,2,2,2}, (const float[]){100000, 10000, 0, 0});
      {
         // When the gradient goes outside the domain of the variable, use nearest neightbour
         DownscalerGradient d(Variable::T ,Options("searchRadius=1 minElevDiff=0"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // Gradient = -0.00797
         EXPECT_FLOAT_EQ(301,   toT(0,0,0)); // nearest neighbour
         EXPECT_FLOAT_EQ(222.80988, toT(0,1,0));
         EXPECT_FLOAT_EQ(302.2684, toT(0,2,0));
      }

      // Fix the gradient
      {
         DownscalerGradient d(Variable::T ,Options("searchRadius=1 minElevDiff=0 constantGradient=-0.1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         EXPECT_FLOAT_EQ(301, toT(0,0,0)); // nearest neighbour
         EXPECT_FLOAT_EQ(301, toT(0,1,0)); // nearest neighbour
         EXPECT_FLOAT_EQ(316.96323, toT(0,2,0));
      }
   }
   TEST_F(TestDownscalerGradient, minGradient) {
      DownscalerGradient d(Variable::T, Options("searchRadius=1 minGradient=-0.01"));
      bool status = d.downscale(*mFrom, *mTo);
      EXPECT_TRUE(status);
      const Field& toT   = *mTo->getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(4, toT.getNumLon());
      // Gradient within range:
      EXPECT_FLOAT_EQ(301.31491, toT(0,0,0)); // 301 - 0.00797 * (120-159.6324)
      // Gradient is too negative (-0.01068):
      EXPECT_FLOAT_EQ(301.4695,  toT(0,2,0)); // 304 - 0.01 * (600-346.9477)
   }
   TEST_F(TestDownscalerGradient, minmaxGradient) {
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 minGradient=-0.01 maxGradient=-0.008"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // Gradient too large (-0.00797):
         EXPECT_FLOAT_EQ(301.3171, toT(0,0,0)); // 301 - 0.008 * (120-159.6324)
         // Gradient is too negative (-0.01068):
         EXPECT_FLOAT_EQ(301.4695,  toT(0,2,0)); // 304 - 0.01 * (600-346.9477)
      }
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 minGradient=0.001 logTransform=1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // Gradient too low (-2.59620e-5)
         EXPECT_FLOAT_EQ(289.3039, toT(0,0,0)); // 301 * exp(0.001 * (120-159.6324))
      }
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=1 maxGradient=-0.01 logTransform=1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         ASSERT_EQ(1, toT.getNumLat());
         ASSERT_EQ(4, toT.getNumLon());
         // Gradient too high (-2.59620e-5)
         EXPECT_FLOAT_EQ(447.3916, toT(0,0,0)); // 301 * exp(-0.01 * (120-159.63))
      }
   }
   TEST_F(TestDownscalerGradient, defaultGradient) {
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=0 defaultGradient=0.1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         EXPECT_FLOAT_EQ(297.0368, toT(0,0,0)); // 301 + 0.1 * (120-159.6324)
         EXPECT_FLOAT_EQ(329.3052,  toT(0,2,0)); // 304 + 0.1 * (600-346.9477)
      }
      {
         DownscalerGradient d(Variable::T, Options("searchRadius=0 defaultGradient=-0.1"));
         bool status = d.downscale(*mFrom, *mTo);
         EXPECT_TRUE(status);
         const Field& toT   = *mTo->getField(Variable::T, 0);
         EXPECT_FLOAT_EQ(304.9632, toT(0,0,0)); // 301 - 0.1 * (120-159.6324)
         EXPECT_FLOAT_EQ(278.6948,  toT(0,2,0)); // 304 - 0.1 * (600-346.9477)
      }
   }
   TEST_F(TestDownscalerGradient, options) {
      DownscalerGradient d(Variable::T, Options("searchRadius=3 defaultGradient=0.1 constantGradient=0.3 minGradient=0 maxGradient=0.2 logTransform=1 minElevDiff=1500"));
      EXPECT_FLOAT_EQ(3, d.getSearchRadius());
      EXPECT_FLOAT_EQ(0.1, d.getDefaultGradient());
      EXPECT_FLOAT_EQ(0.3, d.getConstantGradient());
      EXPECT_FLOAT_EQ(0, d.getMinGradient());
      EXPECT_FLOAT_EQ(0.2, d.getMaxGradient());
      EXPECT_FLOAT_EQ(1, d.getLogTransform());
      EXPECT_FLOAT_EQ(1500, d.getMinElevDiff());
   }
   TEST_F(TestDownscalerGradient, missingValues) {

   }
   TEST_F(TestDownscalerGradient, description) {
      DownscalerGradient::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
