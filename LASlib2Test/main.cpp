#define CATCH_CONFIG_MAIN // This tells Catch to provide a main()

#include "include/catch.hpp"
#include "lasreader.hpp"
#include "Las2Rings.h"

using namespace std;


TEST_CASE( "Trying to read las", "[]" )
{
    string file_path = "../5points_14.las";
    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(file_path.c_str());


    SECTION("File exists") {
        bool result = std::ifstream(file_path.c_str()).good();
        REQUIRE( result );
    }

    SECTION("Can open with LASreader") {
        LASreader* lasreader = lasreadopener.open();
        REQUIRE( lasreader->read_point() );
        lasreader->close();
    }

    SECTION("There is data") {
        LASreader* lasreader = lasreadopener.open();
        lasreader->read_point();
        REQUIRE(lasreader->get_x() > 0);
        lasreader->close();
    }

    SECTION("LAS2Rings")
    {
        SECTION("Can create")
        {
            Las2Rings* las2rings = new Las2Rings();
            delete las2rings;
        }
    }
}
