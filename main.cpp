/* Project las2rings * */

/* *****************************************************
 * This program reads [TLS] LAS/LAZ files and classifies
 * returns supposedly pertaining to vertical surfaces
 * of elements in the laser scanned area
 *
 * Version: 0.0 (20160919)
 *     Basic structure + options and arguments
 *
 * Programmers:
 *     Luiz Estraviz, Caio Hamamura, Tiago do Conto
 * *****************************************************/

#include <getopt.h>    // getopt_long()
#include <stdlib.h>    // exit()
#include <iostream>    // std::cout
#include "lasreader.hpp"
#include "laswriter.hpp"

bool       is_square = false;
std::string inpf_las = "DummyInputFileName.txt";
std::string outf_las = "DummyOutputFileName.txt";


   /** \brief Function: PrintHelp()
    *
    * Prints help when argument -h or --help is used, or
    * when an undefined argument is used
    */

void PrintHelp()
{
    std::cout <<
            "--input <inpfname>  : input LAS/LAZ file name\n"
            "--output <outfname> : output LAS/LAZ file name\n"
            "--east <UTM num>    : UTM easting (float) of the scene central point\n"
            "--north <UTM num>   : UTM northing (float) of the scene central point\n"
            "--radius <val>      : sets the radius/diagonal value (float) for clipping\n"
            "--square            : sets search for rectangles, instead of rings (default)\n"
            "--bottomsize <val>  : sets the min radius/diagonal (float) of rings/rectangles\n"
            "--topsize <val>     : sets the max radius/diagonal (float) of rings/rectangles\n"
            "--help              : show this help\n";
    exit(1);
}

/* **************************************************
   Function: ProcessArgs(argc, argv[])
   las2rings takes options from the user at the
   shell. An option is typically of the form -h
   (short) or --help (long). Options can also have
   a following argument, for example -o someFile.txt
   or --open someFile.txt. All these types of options
   and arguments can be processed using the
   getopt_long function in #include <getopt.h>
   **************************************************/
void ProcessArgs(int argc, char** argv)
{
/* **************************************************
   define local variables and set default values
   **************************************************/
    int opt = -1;
    float radius           = -1.0,
          utm_east         = -1.0,
          utm_north        = -1.0,
          bottom_ring_diam = -1.0,
          top_ring_diam    = -1.0;

/* **************************************************
   short_opts is a string with shortcuts for the
   options. The colon in front of the shortcut
   identifies the options which demand an argument
   **************************************************/
    const char* const short_opts = "i:o:e:n:r:sb:t:h";

/* **************************************************
   long_opts is a struct for the options with
   correspondent shortcuts. The struct has the
   following fields:
   const char *name -> name of the option (string)
   int has_arg -> says whether the option takes an
                  argument. There are three
                  legitimate values: no_argument,
                  required_argument and
                  optional_argument.
   int *flag  ,  int val
   -> control how to report or act on the option
      when it occurs. If flag is a null pointer,
      then the val is a value which identifies this
      option. Otherwise, it should be the address
      of an int variable which is the flag for this
      option. The value in val is the value to store
      in the flag to indicate that the option was seen.
   **************************************************/
    const option long_opts[] = {
            {"input",      1, 0, 'i'},
            {"output",     1, 0, 'o'},
            {"east",       1, 0, 'e'},
            {"north",      1, 0, 'n'},
            {"radius",     1, 0, 'r'},
            {"square",     0, 0, 's'},
            {"bottomsize", 1, 0, 'b'},
            {"topsize",    1, 0, 't'},
            {"help",       0, 0, 'h'},
            {0,      0, 0, 0}
    };

/* **************************************************
   loop for setting the values defined by the user to
   the arguments of each option
   **************************************************/
    while (true)
    {
/* **************************************************
   Function getopt_long is used here. The last
   argument is set to zero because the indexing of
   long_opts is not required.
   **************************************************/
        opt = getopt_long(argc, argv, short_opts, long_opts, 0);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'i':
            inpf_las = std::string(optarg);
            std::cout << "Input file : " << inpf_las << std::endl;
            break;
        case 'o':
            outf_las = std::string(optarg);
            std::cout << "Output file: " << outf_las << std::endl;
            break;
        case 'e':
            utm_east  = atof(optarg);
            std::cout << "UTM easting : " << utm_east << std::endl;
            break;
        case 'n':
            utm_north = atof(optarg);
            std::cout << "UTM northing: " << utm_north << std::endl;
            break;
        case 'r':
            radius = atof(optarg);
            std::cout << "Clipping radius/diagonal: " << radius << std::endl;
            break;
        case 's':
            is_square = true;
            std::cout << "Search for rings is set off. Squared is on.\n";
            break;
        case 'b':
            bottom_ring_diam = atof(optarg);
            std::cout << "Min ring diameter (bottom size): " << bottom_ring_diam << " (or min diagonal if -s)" << std::endl;
            break;
        case 't':
            top_ring_diam = atof(optarg);
            std::cout << "Max ring diameter (top size): " << top_ring_diam << " (or max diagonal if -s)" << std::endl;
            break;

        case 'h':   // -h or --help
        case '?':   // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

/* **************************************************
 * Function: main function
 *
 * argc (argument count) is type int
 *    tells how many command-line arguments are there
 *    it is always at least 1, because the first
 *    string (argv[0]) is the program name
 *
 * argv (argument vector) is type char** or char* []
 *    contains the actual command-line arguments
 *    as an array of strings, the first of which is
 *    the program's name
 * ************************************************* */
#include "lasreader.hpp"
#include "laswriter.hpp"
#include "lasfilter.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <set>

typedef union I32I32I64 { I32 i32[2]; I64 i64; } I32I32I64;

using namespace std;
typedef set<I64> my_pixels;
typedef unordered_map<I64,U32> my_counts;
typedef unordered_map<I32,my_counts*> my_slices;

static void add_count(my_counts* counts, I32 pos_x, I32 pos_y)
{
    I32I32I64 i32i32i64;
    my_counts::iterator count_iterator;

    i32i32i64.i32[0] = pos_x;
    i32i32i64.i32[1] = pos_y;

    count_iterator = counts->find(i32i32i64.i64);
    if (count_iterator != counts->end())
    {
        (*count_iterator).second++;
    }
    else
    {
        counts->insert(my_counts::value_type(i32i32i64.i64, 1));
    }
}

int main(int argc, char **argv)
{
    F32 step_xy = 0.025f;
    F32 step_z = 1.0f;
    F32 c_density = 0.1f;
    F32 c_absolute = 4;
    F32 r_min = 0.05f;
    F32 r_max = 0.50f;
    F32 r_incr = step_xy;
    U32 min_hough_votes = 2;

    string file_path = "C:/LAStools/data/lcer.laz";

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(file_path.c_str());

 //   CHAR argv [] = { "-keep_z", "62.4", "65", "-filtered_transform", "-set_classification", "22" };

    CHAR* MY_other_argv[4];

    MY_other_argv[0] = (CHAR*)malloc(1000); strcpy(MY_other_argv[0], " ");
    MY_other_argv[1] = (CHAR*)malloc(1000); strcpy(MY_other_argv[1], "-keep_z");
    MY_other_argv[2] = (CHAR*)malloc(1000); strcpy(MY_other_argv[2], "1.0");
    MY_other_argv[3] = (CHAR*)malloc(1000); strcpy(MY_other_argv[3], "1.99");

    lasreadopener.parse(4, MY_other_argv);
    LASreader* lasreader = lasreadopener.open();

    cout << "Point reading" << endl;
    I32 p_count = 0;

    my_slices slices;
    my_slices::iterator slice_iterator;
    my_counts* slice_counts = 0;
    my_counts::iterator slice_count_iterator;

    // populate counters

    I32 min_z = I32_MAX;
    I32 max_z = I32_MIN;

    I32I32I64 i32i32i64;
    I32 pos_x;
    I32 pos_y;

    while (lasreader->read_point())
    {
        I32 pos_z = I32_FLOOR(lasreader->get_z()/step_z);

        if (pos_z < min_z) min_z = pos_z;
        if (pos_z > max_z) max_z = pos_z;

        // find the slice with correct pos_z

        slice_iterator = slices.find(pos_z);
        if (slice_iterator != slices.end())
        {
            slice_counts = (*slice_iterator).second;
        }
        else
        {
            slice_counts = new my_counts();
            slices.insert(my_slices::value_type(pos_z, slice_counts));
        }

        assert(slice_counts);

        // increment correct counter

        pos_x = I32_QUANTIZE(lasreader->get_x()/step_xy);
        pos_y = I32_QUANTIZE(lasreader->get_y()/step_xy);

        add_count(slice_counts, pos_x, pos_y);

        p_count++;
    }

    fprintf(stderr, "total points %d [%d,%d]\n", p_count, min_z, max_z);

    lasreader->close();

    #define DEBUG TRUE

    if (DEBUG)
    {
        // produce debug output for current slice

        if (slice_counts && (slice_counts->size()))
        {
            LASwriteOpener laswriteopener;
            laswriteopener.set_file_name(file_path.c_str());
            laswriteopener.set_appendix("_counters");
            laswriteopener.set_format(LAS_TOOLS_FORMAT_LAZ);

            LASheader lasheader;
            lasheader.point_data_format = 0;
            lasheader.point_data_record_length = 20;
            LASpoint laspoint;
            laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, &lasheader);
            LASwriter* laswriter = laswriteopener.open(&lasheader);

            slice_count_iterator = slice_counts->begin();

            while (slice_count_iterator != slice_counts->end())
            {
                i32i32i64.i64 = (*slice_count_iterator).first;

                pos_x = i32i32i64.i32[0];
                pos_y = i32i32i64.i32[1];

                laspoint.set_x(step_xy*pos_x+0.5f*step_xy);
                laspoint.set_y(step_xy*pos_y+0.5f*step_xy);
                laspoint.set_z((*slice_count_iterator).second);
                laswriter->write_point(&laspoint);
                laswriter->update_inventory(&laspoint);

                slice_count_iterator++;
            }
            laswriter->update_header(&lasheader, TRUE);
            laswriter->close(TRUE);
            delete laswriter;
        }
    }

    // process all slices

    I32 total = 0;

    I32 pos_z;

    for (pos_z = min_z; pos_z <= max_z; pos_z++)
    {
        slice_iterator = slices.find(pos_z);
        if (slice_iterator != slices.end())
        {
            U32 c;

            fprintf(stderr, "doing slice %g\n", step_z * pos_z);
            fprintf(stderr, "===============\n");

            slice_counts = (*slice_iterator).second;

            // find highest count per slice

            U32 c_max = 0;

            slice_count_iterator = slice_counts->begin();

            while (slice_count_iterator != slice_counts->end())
            {
                c = (*slice_count_iterator).second;
                if (c > c_max) c_max = c;
                slice_count_iterator++;
            }

            // those are the pixels we center circles on

            U32 c_min = c_absolute;

            if (c_min == 0) c_min = I32_CEIL(c_density*c_max);

            // store all pixels that succeed in the Hough transform

            my_pixels pixels;

            // loop over all radii and perform the Hough transform

            F32 r;
            for (r = r_min; r <= r_max; r += r_incr)
            {
                // construct circle with radius r

                I32 circle_x;
                I32 circle_x_num = I32_CEIL(r_min / step_xy);
                I32* circle_y = new I32[circle_x_num];
                for (circle_x = 0; circle_x < circle_x_num; circle_x++)
                {
                    circle_y[circle_x] = I32_CEIL((sin(acos(circle_x*step_xy/r)))*step_xy);
                }

                // empty grid

                my_counts hough_counts;

                // hough transform for all pixels with c >= c_min

                slice_count_iterator = slice_counts->begin();

                while (slice_count_iterator != slice_counts->end())
                {
                    if ((*slice_count_iterator).second >= c_min)
                    {
                        i32i32i64.i64 = (*slice_count_iterator).first;

                        pos_x = i32i32i64.i32[0];
                        pos_y = i32i32i64.i32[1];

                        for (circle_x = 0; circle_x < circle_x_num; circle_x++)
                        {
                            add_count(&hough_counts, pos_x + circle_x, pos_y + circle_y[circle_x]);
                            add_count(&hough_counts, pos_x + circle_x, pos_y - circle_y[circle_x]);
                        }
                        for (circle_x = 1; circle_x < circle_x_num; circle_x++)
                        {
                            add_count(&hough_counts, pos_x - circle_x, pos_y + circle_y[circle_x]);
                            add_count(&hough_counts, pos_x - circle_x, pos_y - circle_y[circle_x]);
                        }
                    }
                    slice_count_iterator++;
                }

                // delete circle with radius r

                delete [] circle_y;
                circle_y = 0;

                // find and store pixels that succeed in the Hough transform

                my_counts::iterator hough_count_iterator;

                hough_count_iterator = hough_counts.begin();
                while (hough_count_iterator != hough_counts.end())
                {
                    if ((*hough_count_iterator).second >= min_hough_votes)
                    {
                        I64 key = (*hough_count_iterator).first;
                        pixels.insert(key);
                    }
                    hough_count_iterator++;
                }
            }

            fprintf(stderr, "%d pixels passed min_hough_votes of %d\n", pixels.size(), min_hough_votes);

            if (pixels.size())
            {
                LASwriteOpener laswriteopener;
                laswriteopener.make_numbered_file_name(file_path.c_str(), 5);
                laswriteopener.set_format(LAS_TOOLS_FORMAT_LAZ);

                LASheader lasheader;
                lasheader.point_data_format = 0;
                lasheader.point_data_record_length = 20;
                LASpoint laspoint;
                laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, &lasheader);
                LASwriter* laswriter = laswriteopener.open(&lasheader);

                laspoint.set_z(step_z*pos_z+0.5f*step_z);

                my_pixels::iterator pixels_iterator;
                pixels_iterator = pixels.begin();
                while (pixels_iterator != pixels.end())
                {
                   i32i32i64.i64 = (*pixels_iterator);

                    pos_x = i32i32i64.i32[0];
                    pos_y = i32i32i64.i32[1];

                    laspoint.set_x(step_xy*pos_x+0.5f*step_xy);
                    laspoint.set_y(step_xy*pos_y+0.5f*step_xy);
                    laswriter->write_point(&laspoint);
                    laswriter->update_inventory(&laspoint);

                    pixels_iterator++;
                }
                laswriter->update_header(&lasheader, TRUE);
                laswriter->close(TRUE);
                delete laswriter;
            }
        }
        else
        {
            fprintf(stderr, "no point in slice %g\n", step_z * pos_z);
        }
    }

    fprintf(stderr, "total count %d\n", total);

    ProcessArgs(argc, argv);
    return 0;
}
