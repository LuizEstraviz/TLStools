#include "Las2Rings.h"

Las2Rings::Las2Rings()
{
    //ctor
}

Las2Rings::~Las2Rings()
{
    //dtor
}

string Las2Rings::getHelp()
{
    return "--input <inpfname>  : input LAS/LAZ file name\n"
            "--output <outfname> : output LAS/LAZ file name\n"
            "--east <UTM num>    : UTM easting (float) of the scene central point\n"
            "--north <UTM num>   : UTM northing (float) of the scene central point\n"
            "--radius <val>      : sets the radius/diagonal value (float) for clipping\n"
            "--square            : sets search for rectangles, instead of rings (default)\n"
            "--bottomsize <val>  : sets the min radius/diagonal (float) of rings/rectangles\n"
            "--topsize <val>     : sets the max radius/diagonal (float) of rings/rectangles\n"
            "--help              : show this help\n";
}
