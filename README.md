# mcl
{M}arkov {CL}uster algorithm, a tool for clustering undirected weighted or simple networks, a.k.a. graphs

Currently only development mcl releases have been made available from this repository.
Stable releases will follow in the not too distant future.
This code needs the C library in the github repo
[micans/cimfomfa](http://github.com/micans/cimfomfa),
hence the build procedure has been changed somewhat and needs more steps.

- Use the script [buildrelease](buildrelease) to build an mcl release. Note; to build mcl both an mcl release and a cimfomfa release are required.
  The cimfomfa repository also has a buildrelease script.
- Look at [installmcl.sh](installmcl.sh) for an elementary build example. 
- These [mcl](http://micans.org/mcl/dev/mcl-21-257.tar.gz) and [cimfomfa](http://micans.org/dev/cimfomfa-21-257.tar.gz) releases can be used to
  test the installmcl.sh script and build an mcl release.

