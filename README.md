# mcl
{M}arkov {CL}uster algorithm, a tool for clustering undirected weighted or simple networks, a.k.a. graphs

Currently only development mcl releases have been made available from this repository.
Stable releases will follow in the not too distant future.
This code needs the C library in the github repo
[micans/cimfomfa](http://github.com/micans/cimfomfa),
hence the build procedure has been changed somewhat and needs more steps.

- [This script](build-mcl-21-257.sh) will pull http://micans.org/mcl/dev/mcl-21-257.tar.gz
  and http://micans.org/dev/cimfomfa-21-257.tar.gz and then (attempt to) build an mcl release.
- [installmcl.sh](installmcl.sh) is similar with a little more control and needs local tar archives,
  for example the two listed above.

