## Process this file with automake to produce Makefile.in
## $Id: Makefile.am,v 1.16 2007-05-18 11:02:25 flux Exp $

TROFF =  groff
COL   =  col
ZOEM  =  zoem

## tidy is a HTML syntax checker and reformatter, available from
## http://www.w3.org/People/Raggett/tidy/ , by Dave Raggett <dsr@w3.org>
TIDY = tidy

docdir = $(datadir)/doc/$(PACKAGE)

# mclfamily.html mclfamily.ps

doc_DATA_dep = mcl.html mcl.ps mcl.txt \
 mclcm.html mclcm.ps mclcm.txt \
 mclfaq.html mclfaq.ps mclfaq.txt \
 mclfamily.html mclfamily.ps \
 mcx.html mcx.ps \
 mcxi.html mcxi.ps \
 clm.html clm.ps \
 mcxconvert.html mcxconvert.ps \
 mcxarray.html mcxarray.ps \
 mcxdump.html mcxdump.ps \
 mcxdiameter.html mcxdiameter.ps \
 mcxerdos.html mcxerdos.ps \
 mcxalter.html mcxalter.ps \
 mcxclcf.html mcxclcf.ps \
 mcxctty.html mcxctty.ps \
 mcxrand.html mcxrand.ps \
 mcxload.html mcxload.ps \
 mcxsubs.html mcxsubs.ps \
 mcxmap.html mcxmap.ps \
 mcxquery.html mcxquery.ps \
 mcxassemble.html mcxassemble.ps \
 index.html mclindex.ps \
 clmimac.html clmimac.ps \
 clmoptics.html clmoptics.ps \
 clmclose.html clmclose.ps \
 clmresidue.html clmresidue.ps \
 clmdist.html clmdist.ps \
 clmvol.html clmvol.ps \
 clminfo.html clminfo.ps \
 clminfo2.html clminfo2.ps \
 clmmate.html clmmate.ps \
 clmmeet.html clmmeet.ps \
 clmformat.html clmformat.ps \
 clmorder.html clmorder.ps \
 mcxio.html mcxio.ps \
 clmprotocols.html clmprotocols.ps \
 clmprotocols2.html clmprotocols2.ps \
 tingea.log.html tingea.log.ps

doc_DATA = $(doc_DATA_dep) mclmcl.txt

man_MANS = mcl.1 mclcm.1 clm.1 mclindex.7 mclfamily.7 mclfaq.7 mcx.1 mcxi.1 mcxio.5 clmprotocols2.5 clmprotocols.5 \
   clmresidue.1 clmclose.1 clmoptics.1 clmimac.1 clmmate.1 clmmeet.1 clmformat.1 \
   clmvol.1 clmdist.1 clminfo2.1 clminfo.1 clmorder.1 \
   mcxconvert.1 mcxassemble.1 mcxmap.1 mcxsubs.1 mcxload.1 mcxdump.1 mcxarray.1 mcxrand.1 \
	mcxalter.1 mcxerdos.1 mcxdiameter.1 mcxctty.1 mcxclcf.1 mcxquery.1 \
	tingea.log.7



MAINTAINERCLEANFILES = $(man_MANS) $(doc_DATA_dep)

# ALIENFILES =  mcxdeblast.ps mclblastline.html mclblastline.ps
## This is quite hackish (though not deep).
## In the alien files, these are copied to this directory.
## This is needed so that hyperlinks work correctly when this
## documentation is accessed directly from the tarball, rather than
## after an installation.
## So this introduces the dependency that alien documentation must
## be processed and created first.

EXTRA_DIST = mclmcl.txt mcx.zmm\
 stamp.day stamp.month stamp.stamp stamp.year\
 clmformat.azm clmvol.azm clmdist.azm clmmate.azm clminfo2.azm clminfo.azm clmmeet.azm\
 clm.azm clmorder.azm mcl.azm mclcm.azm\
 mclfamily.azm mclfaq.azm mcxi.azm mcx.azm mcxquery.azm \
 mcxdump.azm mcxconvert.azm clmprotocols2.azm clmprotocols.azm mcxio.azm clmresidue.azm clmclose.azm clmoptics.azm clmimac.azm mcxassemble.azm\
 mcxalter.azm mcxerdos.azm mcxdiameter.azm mcxctty.azm mcxclcf.azm mcxrand.azm\
 mcxload.azm mcxmap.azm mcxarray.azm mcxsubs.azm index.in.azm index.html distindex.html\
 tingea.log.azm\
 $(MAINTAINERCLEANFILES)

SUFFIXES = .1 .5 .7 .azm .html .ps .txt

MYPATH = $(srcdir)

mclindex.7: $(MYPATH)/index.in.azm
	$(ZOEM) -i $(MYPATH)/index.in -d roff -o mclindex.7 -s install=1
	$(TROFF) -man mclindex.7 > mclindex.ps

index.html: $(MYPATH)/index.in.azm
	$(ZOEM) -i $(MYPATH)/index.in -o index.html -d html -s install=1
	$(ZOEM) -i $(MYPATH)/index.in -o distindex.html -d html -s dist=1

.azm.1:
	$(ZOEM) -d roff -i $< -o $@
	$(ZOEM) -d roff -i $< -o $@

.azm.5:
	$(ZOEM) -d roff -i $< -o $@
	$(ZOEM) -d roff -i $< -o $@

.azm.7:
	$(ZOEM) -d roff -i $< -o $@
	$(ZOEM) -d roff -i $< -o $@

.1.ps:
	$(TROFF) -man $< > $@

.5.ps:
	$(TROFF) -man $< > $@

.7.ps:
	$(TROFF) -man $< > $@

.1.txt:
	$(TROFF) -t -e -mandoc -Tascii $< | $(COL) -bx > $@

.5.txt:
	$(TROFF) -t -e -mandoc -Tascii $< | $(COL) -bx > $@

.7.txt:
	$(TROFF) -t -e -mandoc -Tascii $< | $(COL) -bx > $@

.azm.html:
	$(ZOEM) -d html -i $< -o $@
	$(ZOEM) -d html -i $< -o $@
	- $(TIDY) -e $@

SHARE1 = "--exclude=version,h,help"
SHARE2 = "--exclude=nop,version,set,h,help,debug,test,progress,write-binary"


# clm omitted
# mcx omitted
# mcxi omitted

# ../scripts/docme --silent --exclude=version,h,help ../src/shmcl/mcl mcl.azm

.PHONY: docme-test
docme-test:
	$(top_srcdir)/scripts/docme --doc-exceed=21 --silent $(SHARE1) ../src/shmcl/mcl $(top_srcdir)/doc/mcl.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmcl/mclcm        $(top_srcdir)/doc/mclcm.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxdump       $(top_srcdir)/doc/mcxdump.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxload       $(top_srcdir)/doc/mcxload.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxarray      $(top_srcdir)/doc/mcxarray.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxrand       $(top_srcdir)/doc/mcxrand.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shcl/clmformat     $(top_srcdir)/doc/clmformat.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxsubs       $(top_srcdir)/doc/mcxsubs.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxmap        $(top_srcdir)/doc/mcxmap.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE1) ../src/shmx/mcxassemble   $(top_srcdir)/doc/mcxassemble.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm close"   $(top_srcdir)/doc/clmclose.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm dist"    $(top_srcdir)/doc/clmdist.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm order"   $(top_srcdir)/doc/clmorder.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm mate"    $(top_srcdir)/doc/clmmate.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm meet"    $(top_srcdir)/doc/clmmeet.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm info"    $(top_srcdir)/doc/clminfo.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm info2"   $(top_srcdir)/doc/clminfo2.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm vol"     $(top_srcdir)/doc/clmvol.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm residue" $(top_srcdir)/doc/clmresidue.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm imac"    $(top_srcdir)/doc/clmimac.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shcl/clm optics"  $(top_srcdir)/doc/clmoptics.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx convert" $(top_srcdir)/doc/mcxconvert.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx diameter"  $(top_srcdir)/doc/mcxdiameter.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx erdos"   $(top_srcdir)/doc/mcxerdos.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx alter"   $(top_srcdir)/doc/mcxalter.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx clcf"    $(top_srcdir)/doc/mcxclcf.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx ctty"    $(top_srcdir)/doc/mcxctty.azm
	$(top_srcdir)/scripts/docme --silent $(SHARE2) "../src/shmx/mcx query"   $(top_srcdir)/doc/mcxquery.azm


mcl: mcl.txt            mcl.1            mcl.html            mcl.ps
mclcm: mclcm.txt      mclcm.1          mclcm.html          mclcm.ps

clm:                    clm.1            clm.html            clm.ps
clminfo2:          clminfo2.1       clminfo2.html       clminfo2.ps
clminfo:            clminfo.1        clminfo.html        clminfo.ps
clmmeet:            clmmeet.1        clmmeet.html        clmmeet.ps
clmformat:        clmformat.1      clmformat.html      clmformat.ps
clmdist:            clmdist.1        clmdist.html        clmdist.ps
clmvol:              clmvol.1         clmvol.html         clmvol.ps
clmmate:            clmmate.1        clmmate.html        clmmate.ps
clmresidue:      clmresidue.1     clmresidue.html     clmresidue.ps
clmimac:            clmimac.1        clmimac.html        clmimac.ps
clmoptics:        clmoptics.1      clmoptics.html      clmoptics.ps
clmclose:          clmclose.1       clmclose.html       clmclose.ps
clmorder:          clmorder.1       clmorder.html       clmorder.ps
                                                      
mcx:                    mcx.1            mcx.html            mcx.ps
mcxi:                  mcxi.1           mcxi.html           mcxi.ps
mcxassemble:    mcxassemble.1    mcxassemble.html    mcxassemble.ps
mcxmap:              mcxmap.1         mcxmap.html         mcxmap.ps
mcxsubs:            mcxsubs.1        mcxsubs.html        mcxsubs.ps
mcxarray:          mcxarray.1       mcxarray.html       mcxarray.ps
mcxdump:            mcxdump.1        mcxdump.html        mcxdump.ps
mcxload:            mcxload.1        mcxload.html        mcxload.ps
mcxconvert:      mcxconvert.1     mcxconvert.html     mcxconvert.ps

mcxdiameter:    mcxdiameter.1    mcxdiameter.html    mcxdiameter.ps
mcxerdos:          mcxerdos.1       mcxerdos.html       mcxerdos.ps
mcxalter:          mcxalter.1       mcxalter.html       mcxalter.ps
mcxclcf:            mcxclcf.1        mcxclcf.html        mcxclcf.ps
mcxctty:            mcxctty.1        mcxctty.html        mcxctty.ps
mcxrand:            mcxrand.1        mcxrand.html        mcxrand.ps
mcxquery:           mcxquery.1      mcxquery.html       mcxquery.ps
                                                     
mcxio:                mcxio.5          mcxio.html          mcxio.ps
clmprotocols:  clmprotocols.5   clmprotocols.html   clmprotocols.ps
clmprotocols2: clmprotocols2.5 clmprotocols2.html  clmprotocols2.ps
tingea.log:		  tingea.log.7		 tingea.log.html		tingea.log.ps

mclfaq:  mclfaq.txt  mclfaq.7         mclfaq.html         mclfaq.ps	minimcl
mclfamily:        mclfamily.7      mclfamily.html      mclfamily.ps

index:             mclindex.7          index.html distindex.html

minimcl:
	cp ../scripts/minimcl .

