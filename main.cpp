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
 *     Luiz Estraviz, Caio Hakamada, Tiago do Conto
 * *****************************************************/

#include <getopt.h>    // getopt_long()
#include <stdlib.h>    // exit()
#include <iostream>    // std::cout

bool       is_square = false;
std::string inpf_las = "DummyInputFileName.txt";
std::string outf_las = "DummyOutputFileName.txt";

/* **************************************************
   Function: PrintHelp()
   Prints help when argument -h or --help is used, or
   when an undefined argument is used
   **************************************************/
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
int main(int argc, char **argv)
{
    ProcessArgs(argc, argv);
    return 0;
}
