
# MCL

![Visualisation of MCL](/img/fa_header_800_200.png)

Markov CLustering or the Markov CLuster algorithm, MCL is a method
for clustering weighted or simple networks, a.k.a. graphs.  It is accompanied
in this source code by other network-related
programs, one of which is RCL (restricted contingency linkage) for fast
multi-resolution consensus clustering (see below).
**If you use this software, please cite**

van Dongen, Stijn, Graph clustering via a discrete uncoupling process, Siam
Journal on Matrix Analysis and Applications 30-1, p121-141, 2008.

[See https://doi.org/10.1137/040608635](https://doi.org/10.1137/040608635)

This MCL implementation is fast, threaded, and uses sparse matrices. It runs on a single
machine and can use multiple CPUs. On capable hardware it can cluster
networks with millions of nodes and billions of edges within hours.

## Installation
Installing MCL software requires a compilable source tree. The code in this
repository requires processing by autotools to produce such a tree.
Hence, to use MCL software this repository is not the right source unless
you are interested in development. For installing the current MCL release
from [micans.org/mcl](https://micans.org/mcl)
use [the script install-this-mcl.sh](install-this-mcl.sh). On Linux
and MacOS (if you have development tools installed on MacOS) the following
lines pasted in a terminal (or saved to file and sourced) will install MCL.

```
mkdir installmcl
cd installmcl
wget https://raw.githubusercontent.com/micans/mcl/main/install-this-mcl.sh -o install-this-mcl
chmod u+x install-this-mcl.sh
./install-this-mcl.sh
mcl --version        # test install
```

By default programs are installed in `$HOME/local/bin` and the multi-resolution
consensus clustering program `rcl` is enabled (see below). Edit the file `install-this-mcl.sh`
before executing it if you want to make changes.


## Applications and bioinformatics
MCL has been used a lot in the field of bioinformatics, starting with the TribeMCL
method published by Enright, van Dongen and Ouzounis.
A lot (too much) of information and documentation is available
at [micans.org/mcl](https://micans.org/mcl) .

## RCL, fast multi-resolution consensus clustering
RCL, (f)or Restricted Contingency Linkage, is a fast and *parameterless* method for integrating
multiple flat clusterings at different levels of resolutions. There is
no requirement on these clusterings; they can be made by any method or combination of methods,
for example by Leiden with different resolution parameters, or by MCL with
different inflation values.
The implementation provided here
in [the RCL directory](rcl) just requires the input to be in the native mcl
matrix format.
For Seurat results this is facilitated by the scripts
[rcl/srt2tab.sh](rcl/srt2tab.sh) (this establishes a mapping from labels to
indexes) and [rcl/srt2cls.sh](rcl/srt2cls.sh) (this translates a Seurat
`<LABEL><CLSID>` file to mcl matrix format).

[This preprint quite extensively describes RCL](https://www.biorxiv.org/content/10.1101/2022.10.09.511493v1),
including application of RCL to a large-scale single-cell kidney data set of 27k cells.


## Status and plans
The program MCL has been very stable or nearly unchanging for well over 15
years now. The last speed optimisations happened in 2010. I
aim to do some development in its sibling programs, including improving those
that implement (currently somewhat inelegant) mini-formats such as `mcx alter`
and `mcxsubs`.

RCL (see above) was recently added, so some focus will be to improve/extend its
implementation, as well as support and documentation.

A second are of development, tied to the first, will be
low-level fast loading, filtering and subsetting of large networks and matrices.
I have found occasional use for this in past projects, and under certain conditions
an mcl-edge recipe, if possible, will be a few times faster than a recipe using
Python panda or R. Challenge of this type are sometimes a reason
to expand mcl-edge's capabilities.

Another reason for new releases is that new compilers and
compiler settings have unearthed two or three blemishes in the code base that
needed fixing.
This code needs the C library in the github repo
[micans/cimfomfa](http://github.com/micans/cimfomfa),
hence the build procedure has been changed somewhat and needs more steps (see instructions above).

If you have questions, suggestions, or problems please open an issue or discussion.

