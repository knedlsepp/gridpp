#include "../File/Ec.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileEcTest : public ::testing::Test {
   };

   TEST_F(FileEcTest, projectionGrid) {
      FileEc file1("testing/files/validEc1.nc");
      vec2 lats = file1.getLats();
      vec2 lons = file1.getLons();
      vec2 elevs = file1.getElevs();
      ASSERT_EQ(3, lats.size());
      ASSERT_EQ(3, lons.size());
      ASSERT_EQ(3, elevs.size());
      ASSERT_EQ(2, file1.getNumEns());
      ASSERT_EQ(2, file1.getNumTime());
      EXPECT_FLOAT_EQ(1, lats[1][1]);
      EXPECT_FLOAT_EQ(2, lats[2][2]);
      EXPECT_FLOAT_EQ(1, lons[1][1]);
      EXPECT_FLOAT_EQ(2, lons[0][2]);
      EXPECT_FLOAT_EQ(-1.716614e-05, elevs[1][1]);
      EXPECT_FLOAT_EQ(11.73094, elevs[0][1]);
      EXPECT_FLOAT_EQ(3, lons.size());
      FieldPtr temp = file1.getField(Variable::T, 1);
      EXPECT_FLOAT_EQ(21, (*temp)(0,0,0));
      EXPECT_FLOAT_EQ(33, (*temp)(0,2,1));
   }

   TEST_F(FileEcTest, latlonGrid) {
      // The file is on a perfect lat/lon grid, so the lat/lon variables
      // only have one dimension.
      FileEc file1("testing/files/validEc2.nc");
      vec2 lats = file1.getLats();
      vec2 lons = file1.getLons();
      vec2 elevs = file1.getElevs();
      ASSERT_EQ(3, lats.size());
      ASSERT_EQ(3, lons.size());
      ASSERT_EQ(3, elevs.size());
      ASSERT_EQ(2, file1.getNumEns());
      ASSERT_EQ(2, file1.getNumTime());
      EXPECT_FLOAT_EQ(1, lats[1][1]);
      EXPECT_FLOAT_EQ(2, lats[2][2]);
      EXPECT_FLOAT_EQ(1, lons[1][1]);
      EXPECT_FLOAT_EQ(2, lons[0][2]);
      EXPECT_FLOAT_EQ(-1.716614e-05, elevs[1][1]);
      EXPECT_FLOAT_EQ(11.73094, elevs[0][1]);
      EXPECT_FLOAT_EQ(3, lons.size());
      FieldPtr temp = file1.getField(Variable::T, 1);
      EXPECT_FLOAT_EQ(21, (*temp)(0,0,0));
      EXPECT_FLOAT_EQ(33, (*temp)(0,2,1));
   }

   TEST_F(FileEcTest, order) {
      // In this test, the x and y dims are reversed for longitude and altitude
      FileEc file1("testing/files/validEc3.nc");
      vec2 lats = file1.getLats();
      vec2 lons = file1.getLons();
      vec2 elevs = file1.getElevs();
      ASSERT_EQ(3, lats.size());
      ASSERT_EQ(3, lons.size());
      ASSERT_EQ(3, elevs.size());
      ASSERT_EQ(2, file1.getNumEns());
      ASSERT_EQ(2, file1.getNumTime());
      EXPECT_FLOAT_EQ(1, lats[1][1]);
      EXPECT_FLOAT_EQ(2, lats[2][2]);
      EXPECT_FLOAT_EQ(1, lons[1][1]);
      EXPECT_FLOAT_EQ(2, lons[0][2]);
      EXPECT_FLOAT_EQ(-1.716614e-05, elevs[1][1]);
      EXPECT_FLOAT_EQ(11.73094, elevs[0][1]);
      EXPECT_FLOAT_EQ(3, lons.size());
      FieldPtr temp = file1.getField(Variable::T, 1);
      EXPECT_FLOAT_EQ(21, (*temp)(0,0,0));
      EXPECT_FLOAT_EQ(33, (*temp)(0,2,1));
   }

   TEST_F(FileEcTest, validFiles) {
      // Valid lat/lon ec file
      FileEc file1("testing/files/validEc1.nc");
      // Valid LQQT file
      FileEc file2("testing/files/validEc2.nc");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
