# mcl

{M}arkov {CL}uster algorithm, a method and software for clustering undirected
weighted or simple networks, a.k.a. graphs

Markov CLustering or the Markov CLuster algorithm, mcl is a method and program
for clustering weighted or simple networks, a.k.a. graphs.  It is accompanied
in this source code by other network construction and (trait) analysis
programs.

## Applications and bioinformatics
MCL has been used a lot in the field of bioinformatics, starting with the TribeMCL
method published by Enright, van Dongen and Ouzounis.
A lot (too much) of information and documentation is available
at [micans.org](http://micans.org/mcl) .

## Mcl-edge
This is the name I use sometimes to refer to the collection of sibling
programs to mcl for network and matrix loading, filtering and subsetting,
including computation of network traits such as centrality and clustering coefficient.

## Status and plans
The program mcl has been very stable or nearly unchanging for well over 15
years now, with the last speed optimisations occurring about a decade ago. I
aim to do some development in its sibling programs, including improving those
that implement (currently somewhat inelegant) mini-formats such as `mcx alter`
and `mcxsubs`.

A second are of development, tied to the first, will be
low-level fast loading, filtering and subsetting of large networks and matrices.
I have found occasional use for this in past projects, and under certain conditions
an mcl-edge recipe, if possible, will be a few times faster than a recipe using
Python panda or R. Challenge of this type are sometimes a reason
to expand mcl-edge's capabilities.

Another reason for new releases is that new compilers and
compiler settings have unearthed two or three blemishes in the code base that
needed fixing.

I aim to make a new release, the first since mcl-14-137 (released at about day 137
in the year 2014), in autumn 2021. A development release (mcl-21-257) is available,
please have a look at the dev branch.

If you have suggestions please get in touch or open an issue.

