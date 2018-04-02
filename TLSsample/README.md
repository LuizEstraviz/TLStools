# TLSsample

A command line tool for resampling of dense Terrestrial Laser Scanning (TLS) point clouds.

## Getting Started

To run the executable just run `TLSsample -h` from shell and the following help page will be loaded:

```
# /*** TLSsample ***/
# /*** Command line arguments ***/

# -i or --input         : input file path
# -o or --odir          : output directory (defaults to the input file directory)
# -v or --voxel         : (default method) voxel side length (default = 0.05 m)
# -r or --random        : randomly sample a proportion (0-1] of the input cloud.
                         *this method is applied only if the argument is explicitly given. Otherwise, voxel sampling is performed.
# --odix                : output suffix for the sampled point cloud (default = sample)
# -? or -h or --help    : print this help

``` 

## Adopted methods

Two sampling methods are provided: random and systematic (or voxel-wise). The former performs a random selection according to a provided proportion of remaining points desired. The voxel-wise method (default) samples one random point per voxel of given size - this method outputs a point cloud of standardized density, i.e. every region of the point cloud will have approximately the same amount of pts/m³.

### Further details and recommendations on the command line usage

The command arguments follow the Unix style, specified by short `(-char)` and/or long `(--string)` parameters, followed by their respective values, when required.

* `-i (--input)`

the input path of the point cloud to be "thinned" in `las`, `laz` or `txt` format.

* `o (--odir)`

directory to save the point cloud output.

* `-v (--voxel)`

voxel side length for systematic resampling.

* `-r (--random)` 

proportion, between 0 and 1, to keep when random sampling the point cloud.

* `--odix`

suffix to append to the name of the point cloud output.