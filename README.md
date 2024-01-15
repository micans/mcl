
[Installation and MCL versions](#installation-and-mcl-versions)  
[Applications and bioinformatics](#applications-and-bioinformatics)  
[Quick pointers](#quick-pointers)  
[RCL, fast multi-resolution consensus clustering](#rcl-fast-multi-resolution-consensus-clustering)  
[Status and plans](#status-and-plans)  

# MCL

![Visualisation of MCL](/img/fa_header_800_200.png)

Markov CLustering or the Markov CLuster algorithm, MCL is a method
for clustering weighted or simple networks, a.k.a. graphs.  It is accompanied
in this source code by other network-related
programs, one of which is RCL (restricted contingency linkage) for fast
multi-resolution consensus clustering (see below).
**If you use this software, please cite**

van Dongen, Stijn, Graph clustering via a discrete uncoupling process, Siam
Journal on Matrix Analysis and Applications 30-1, p121-141, 2008. [https://doi.org/10.1137/040608635](https://doi.org/10.1137/040608635) .

The algorithm was conceived in 1998 and [first published in a technical report in 1998](https://ir.cwi.nl/pub/4604).
A PhD thesis and three more technical reports [followed in 2000](https://micans.org/mcl/index.html?sec_thesisetc).
The paper above is the result of a long-winded review process that started in
2000 and lay dormant for a long time, for reasons not entirely untypical within the realms of scientific publishing.
A lot more (too much) information and documentation is available at [micans.org/mcl](https://micans.org/mcl) .

This MCL implementation is fast, threaded, and uses sparse matrices. It runs on a single
machine and can use multiple CPUs. On capable hardware it can cluster
networks with millions of nodes and billions of edges within hours.

## Installation and MCL versions

Releases 14-137 and 22-282 of MCL are **available
on Bioconda and many flavours of Linux and BSD**, including **Debian**,
**Ubuntu** and **OpenBSD**.  Release 14-137 is a fine version; this MCL
implementation has not noticeably changed over the past decade, so for using
the clustering program `mcl` it does not matter which of these versions you have.
See [Status and Plans](#status-and-plans) below for more detail and why
the software is still being developed nonetheless.
Many thanks to Joost van baal, Kusalananda, Andreas Tille and other
maintainers who package(d) MCL for Debian and other Linux and BSD releases.
Example package installation commands:

```
conda install bioconda::mcl

apt-get install mcl           # Debian, Ubuntu
```

Installing MCL software without a package manager requires a compilable source tree. The code in this
repository requires processing by [autotools](https://en.wikipedia.org/wiki/GNU_Autotools) to produce such a tree.
Hence, to use MCL software this repository is not the right source unless you are interested in development.
This code additionally needs the C library in the github repository
[micans/cimfomfa](http://github.com/micans/cimfomfa) (previously cimfomfa was included within the mcl distribution).
The build procedure has been changed accordingly and coordinates installation of both cimfomfa and mcl.

For installing the current MCL release from [micans.org/mcl](https://micans.org/mcl)
use [the script install-this-mcl.sh](install-this-mcl.sh).
On Linux and MacOS (if you have development tools installed on MacOS) the following
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

MCL's build environment was created by Joost van Baal - many thanks Joost!

The current release is 22-282, which is without open issues that relate to mcl.
You'll need this release if you want to experiment with RCL (consensus clustering
integrating results for different granularities / inflation values / resolution
values).  The release also [fixes some mcxarray issues](ChangeLog).


## Applications and bioinformatics
MCL has been used a lot in the field of bioinformatics, starting with the TribeMCL
method published by Enright, van Dongen and Ouzounis.
For bioinformatic applications, please cite additionally

Enright A.J., Van Dongen S., Ouzounis C.A.
An efficient algorithm for large-scale detection of protein families,
Nucleic Acids Research 30(7):1575-1584 (2002).
[https://pubmed.ncbi.nlm.nih.gov/11917018/](https://pubmed.ncbi.nlm.nih.gov/11917018/)


## Quick pointers
The quickest way to try out MCL is to provide it with a file that has three tab-separated columns,
where each line is of the form `LABEL1<tab>LABEL2<tab>VALUE`. Such a line represents an edge
from `LABEl1` to `LABEL2` with weight `VALUE`. If the file is called `MYFILE` you can run MCL
like this:
```
mcl MYFILE --abc -I 2.0
```
Output will be in the file `out.MYFILE.I20`, where each line is cluster written as a list of labels.
It is recommended to try a few different inflation values (the `-I` parameter), e.g.
```
mcl MYFILE --abc -I 1.4
mcl MYFILE --abc -I 2.0
mcl MYFILE --abc -I 3.0
mcl MYFILE --abc -I 5.0
```

How you construct the network is important.
[Some recipes can be found here](https://micans.org/mcl/man/clmprotocols.html), and
[some Frequently Answered Questions here](https://micans.org/mcl/man/mclfaq.html) (the latter is a bit over the top).
For large data the 'abc' format just described becomes very slow to load.
[Use these instructions on the recipe page](https://micans.org/mcl/man/clmprotocols.html#large)
to convert 'abc' format to a binary format that is orders of magnitude faster to load.


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

If you have questions, suggestions, or problems please open an issue or discussion.

