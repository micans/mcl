/*   (C) Copyright 2001, 2002, 2003, 2004, 2005 Stijn van Dongen
 *   (C) Copyright 2006, 2007, 2008, 2009, 2010 Stijn van Dongen
 *
 * This file is part of MCL.  You can redistribute and/or modify MCL under the
 * terms of the GNU General Public License; either version 3 of the License or
 * (at your option) any later version.  You should have received a copy of the
 * GPL along with MCL, in the file COPYING.
*/

#include <string.h>
#include <stdio.h>

#include "clm.h"
#include "report.h"
#include "clmclose.h"

#include "tingea/io.h"
#include "tingea/err.h"
#include "tingea/types.h"
#include "tingea/alloc.h"
#include "tingea/opt.h"
#include "tingea/minmax.h"
#include "tingea/compile.h"

#include "impala/edge.h"
#include "impala/matrix.h"
#include "impala/vector.h"
#include "impala/io.h"
#include "impala/app.h"
#include "impala/iface.h"

#include "mcl/interpret.h"
#include "mcl/transform.h"

#include "clew/scan.h"
#include "clew/clm.h"


#include "impala/matrix.h"
#include "impala/io.h"
#include "impala/iface.h"
#include "impala/compose.h"
#include "impala/ivp.h"
#include "impala/app.h"
#include "impala/stream.h"

#include "clew/clm.h"

#include "tingea/types.h"
#include "tingea/err.h"
#include "tingea/opt.h"

static const char* me = "clmclose";

enum
{  MY_OPT_IMX = CLM_DISP_UNUSED
,  MY_OPT_ABC
,  MY_OPT_DOMAIN
,  MY_OPT_OUTPUT
,  MY_OPT_READASIS
,  MY_OPT_WRITECC
,  MY_OPT_WRITECOUNT
,  MY_OPT_WRITESIZES
,  MY_OPT_WRITESIZECOUNTS
,  MY_OPT_LEVELS
,  MY_OPT_LEVELS_NORM
,  MY_OPT_SL
,  MY_OPT_SLLIST
,  MY_OPT_SL_RCL_CUTOFF
,  MY_OPT_WRITEGRAPH
,  MY_OPT_WRITEGRAPHC
,  MY_OPT_CCBOUND
,  MY_OPT_TABIN
,  MY_OPT_MXOUT
,  MY_OPT_TABOUT
,  MY_OPT_TABXOUT
,  MY_OPT_MAPOUT
,  MY_OPT_TF
,  MY_OPT_DEBUG
,  MY_OPT_CAN
}  ;


mcxOptAnchor closeOptions[] =
{  {  "-o"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_OUTPUT
   ,  "<fname>"
   ,  "output file name"
   }
,  {  "-imx"
   ,  MCX_OPT_HASARG | MCX_OPT_REQUIRED
   ,  MY_OPT_IMX
   ,  "<fname>"
   ,  "input matrix file, presumably dumped mcl iterand or dag"
   }
,  {  "-abc"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_ABC
   ,  "<fname>"
   ,  "specify input using label pairs"
   }
,  {  "--is-undirected"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_READASIS
   ,  NULL
   ,  "use if graph is known to be symmetric (slightly faster)"
   }
,  {  "--write-cc"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITECC
   ,  NULL
   ,  "output cluster/connected-component file"
   }
,  {  "--write-size-counts"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITESIZECOUNTS
   ,  NULL
   ,  "output compmressed list of component sizes"
   }
,  {  "--write-sizes"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITESIZES
   ,  NULL
   ,  "output list of component sizes"
   }
,  {  "-levels"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_LEVELS
   ,  "low/step/high[/prefix]"
   ,  "write cluster size distribution for each (edge weight cut-off) level\n"
      "                if prefix is specified, write each to file"
   }
,  {  "-levels-norm"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_LEVELS_NORM
   ,  "<num>"
   ,  "divide each level defined by -levels by <num> to define cutoff"
   }
,  {  "--sl"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_SL
   ,  NULL
   ,  "output single linkage tree encoded as list of joins"
   }
,  {  "-sl-rcl-cutoff"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_SL_RCL_CUTOFF
   ,  "<num>"
   ,  "A value inbetween 0-1000 (suggested:100) at which to stop joining"
   }
,  {  "--write-count"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITECOUNT
   ,  NULL
   ,  "output number of components"
   }
,  {  "--write-block"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITEGRAPH
   ,  NULL
   ,  "write graph restricted to -dom argument"
   }
,  {  "--write-blockc"
   ,  MCX_OPT_DEFAULT
   ,  MY_OPT_WRITEGRAPHC
   ,  NULL
   ,  "write the complement of graph restricted to -dom argument"
   }
,  {  "-cc-bound"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_CCBOUND
   ,  "<num>"
   ,  "select selected components of size at least <num>"
   }
,  {  "-tab"
   ,  MCX_OPT_HASARG | MCX_OPT_HIDDEN
   ,  MY_OPT_TABIN
   ,  "<fname>"
   ,  "read tab file"
   }
,  {  "-write-sl-list"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_SLLIST
   ,  "<fname>"
   ,  "write list of join order with weights"
   }
,  {  "-write-tab"
   ,  MCX_OPT_HASARG | MCX_OPT_HIDDEN
   ,  MY_OPT_TABOUT
   ,  "<fname>"
   ,  "write tab file of selected domain"
   }
,  {  "-write-tabx"
   ,  MCX_OPT_HASARG | MCX_OPT_HIDDEN
   ,  MY_OPT_TABXOUT
   ,  "<fname>"
   ,  "write tab file of deselected domain"
   }
,  {  "-write-matrix"
   ,  MCX_OPT_HASARG | MCX_OPT_HIDDEN
   ,  MY_OPT_MXOUT
   ,  "<fname>"
   ,  "write matrix of selected domain"
   }
,  {  "-write-map"
   ,  MCX_OPT_HASARG | MCX_OPT_HIDDEN
   ,  MY_OPT_MAPOUT
   ,  "<fname>"
   ,  "write mapping"
   }
,  {  "-dom"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_DOMAIN
   ,  "<fname>"
   ,  "input domain file"
   }
,  {  "-tf"
   ,  MCX_OPT_HASARG
   ,  MY_OPT_TF
   ,  "<tf-spec>"
   ,  "first apply tf-spec to matrix"
   }
,  {  "--debug"
   ,  MCX_OPT_DEFAULT | MCX_OPT_HIDDEN
   ,  MY_OPT_DEBUG
   ,  NULL
   ,  "set debug"
   }
,  {  "--canonical"
   ,  MCX_OPT_DEFAULT | MCX_OPT_HIDDEN
   ,  MY_OPT_CAN
   ,  NULL
   ,  "make result matrix canonical"
   }
,  {  NULL ,  0 ,  0 ,  NULL, NULL}
}  ;


static mcxIO*  xfout    =  (void*) -1;
static mcxIO*  xfmx     =  (void*) -1;
static mcxIO*  xfabc    =  (void*) -1;
static mcxIO*  xftabin  =  (void*) -1;
static mcxIO*  xftabout =  (void*) -1;
static mcxIO*  xftabxout=  (void*) -1;
static mcxIO*  xfmapout =  (void*) -1;
static mcxIO*  xfmxout  =  (void*) -1;
static mcxIO*  xfdom    =  (void*) -1;
static mcxTing* tfting  =  (void*) -1;
static dim     ccbound_num  =  -1;
static mcxbool canonical=  -1;
static mcxbool debug_g    =  -1;
static mcxbool make_symmetric=  -1;
static mcxmode write_mode = -1;

static ofs     hi_g     =  -1;
static ofs     lo_g     =  -1;
static ofs     st_g     =  -1;
static const char* levels_pfx = NULL;

static double  sgl_rcl_thr_g = 0.0;
static double  norm_g   =  0.0;
static mcxbool sgl_g    =  FALSE;      /* once there was a reason for the -1 initialisations,
                                        * but TBH I forgot.
                                       */
const char* fn_nodelist =  "nodes.list";


static mcxstatus closeInit
(  void
)
   {  xfout    =  mcxIOnew("-", "w")
   ;  write_mode = MY_OPT_WRITESIZES
   ;  xfmapout =  NULL
   ;  xfmxout  =  NULL
   ;  xftabout =  NULL
   ;  xftabxout=  NULL
   ;  xftabin  =  NULL
   ;  xfmx     =  mcxIOnew("-", "r")
   ;  xfabc    =  NULL
   ;  xfdom    =  NULL
   ;  tfting   =  NULL
   ;  ccbound_num  =  0
   ;  canonical=  FALSE
   ;  debug_g  =  FALSE
   ;  make_symmetric =  TRUE
   ;  hi_g     =  0;
   ;  lo_g     =  0;
   ;  st_g     =  1;
   ;  return STATUS_OK
;  }


static mcxstatus closeArgHandle
(  int optid
,  const char* val
)
   {  switch(optid)
      {  case MY_OPT_IMX
      :  mcxIOnewName(xfmx, val)
      ;  break
      ;

         case MY_OPT_READASIS
      :  make_symmetric = FALSE
      ;  break
      ;

         case MY_OPT_ABC
      :  xfabc = mcxIOnew(val, "r")
      ;  break
      ;

         case MY_OPT_OUTPUT
      :  mcxIOnewName(xfout, val)
      ;  break
      ;

         case MY_OPT_WRITEGRAPHC
      :  write_mode = MY_OPT_WRITEGRAPHC
      ;  break
      ;

         case MY_OPT_WRITEGRAPH
      :  write_mode = MY_OPT_WRITEGRAPH
      ;  break
      ;

         case MY_OPT_WRITECC
      :  write_mode = MY_OPT_WRITECC
      ;  break
      ;

         case MY_OPT_WRITECOUNT
      :  write_mode = MY_OPT_WRITECOUNT
      ;  break
      ;

         case MY_OPT_WRITESIZECOUNTS
      :  write_mode = MY_OPT_WRITESIZECOUNTS
      ;  break
      ;

         case MY_OPT_LEVELS_NORM
      :  norm_g = atof(val)
      ;  break
      ;

         case MY_OPT_LEVELS
      :  {  unsigned long l, s, h
         ;  static char cbuf[50] = { 0 }
         ;  if
            (  4 != sscanf(val, "%lu/%lu/%lu/%49s", &l, &s, &h, cbuf)
            && 3 != sscanf(val, "%lu/%lu/%lu", &l, &s, &h)
            )
            mcxDie(1, me, "cannot parse -levels low/step/high or low/step/high/FILEPREFIX")
         ;  lo_g = l
         ;  hi_g = h
         ;  st_g = s
         ;  if (cbuf[0])
            levels_pfx = cbuf
      ;  }
         break
      ;

         case MY_OPT_SL
      :  sgl_g = TRUE
      ;  break
      ;

         case MY_OPT_SL_RCL_CUTOFF
      :  sgl_rcl_thr_g = atof(val)
      ;  break
      ;

         case MY_OPT_SLLIST
      :  fn_nodelist = val
      ;  break
      ;

         case MY_OPT_WRITESIZES
      :  write_mode = MY_OPT_WRITESIZES
      ;  break
      ;

         case MY_OPT_DOMAIN
      :  xfdom= mcxIOnew(val, "r")
      ;  break
      ;

         case MY_OPT_MAPOUT
      :  xfmapout = mcxIOnew(val, "w")
      ;  break
      ;

         case MY_OPT_MXOUT
      :  xfmxout = mcxIOnew(val, "w")
      ;  break
      ;

         case MY_OPT_CCBOUND
      :  ccbound_num = atoi(val)
      ;  break
      ;

         case MY_OPT_TABXOUT
      :  xftabxout = mcxIOnew(val, "w")
      ;  break
      ;

         case MY_OPT_TABOUT
      :  xftabout = mcxIOnew(val, "w")
      ;  break
      ;

         case MY_OPT_TABIN
      :  xftabin = mcxIOnew(val, "r")
      ;  break
      ;

         case MY_OPT_DEBUG
      :  debug_g = TRUE
      ;  break
      ;

         case MY_OPT_CAN
      :  canonical = TRUE
      ;  break
      ;

         case MY_OPT_TF
      :  tfting = mcxTingNew(val)
      ;  break
      ;

         default
      :  return STATUS_FAIL
   ;  }
      return STATUS_OK
;  }


static double mclv_check_ccbound
(  const mclv* vec
,  void* data
)
   {  dim bound = *((dim*) data)
   ;  return vec->n_ivps >= bound ? 1.0 : 0.0
;  }


static int edge_val_cmp
(  const void* x
,  const void* y
)
   {  const mcle* e = x
   ;  const mcle* f = y
   ;  return e->val < f->val ? 1 : e->val > f->val ? -1 : 0
;  }


                        /* the role of cid / cluster id is to identify set membership
                         * of nodes. It is re-used throughout linking; when linking two
                         * sets the largest set gets to keep its ID, the smaller set
                         * is updated to the same ID.
                         * The ID is used when inspecting an edge to see if its endpoints
                         * are in different clusters.
                        */
struct slnode
{  mcxTing* name        /* Name that's written to the join-order file      */
;  struct slnode* next  /* so that we can iterate through a set of leaf nodes to update their current cluster ID */
;  struct slnode* last  /* so that we can quickl merge two linked lists    */
;  dim      lid         /* leaf ID, not strictly necessary; equal to offset in array */
;  dim      cid         /* current cluster ID; starts out identical to leaf id */
;  dim      size        /* current count of all leaf nodes below this node */
;  dim      lss         /* current largest sub split below this node       */
;  dim      nsg         /* number of singletons joining a bigger cluster   */
;
}  ;

void* node_init(void* v)
{  struct slnode* node = v
;  node->name = mcxTingNew("")
;  node->next = NULL
;  node->last = NULL
;  node->lid  = 0
;  node->cid  = 0
;  node->size = 1
;  node->lss  = 0
;  node->nsg  = 0
;  return node
;
}


static mcxstatus closeMain
(  int          argc_unused      cpl__unused
,  const char*  argv_unused[]    cpl__unused
)
   {  mclx *dom =  NULL, *cc = NULL, *ccbound = NULL
   ;  mclx *mx  =  NULL
   ;  mclx *map  =  NULL
   ;  const mclv *ccbound_rows = NULL
   ;  const mclv *ccbound_cols = NULL
   ;  mclTab* tab = NULL
   ;  dim N_start = 0
   ;  dim N_bound = 0

   ;  mclxIOstreamer streamer = { 0 }

   ;  if ((xftabout || xftabxout) && !xftabin)
      mcxDie(1, me, "-write-tab currently requires -tab or -abc")

   ;  if (xftabin)
      tab = streamer.tab_sym_in = mclTabRead(xftabin, NULL, EXIT_ON_FAIL)

   ;  mcxIOopen(xfout, EXIT_ON_FAIL)

   ;  if (xfabc)
      mx
      =  mclxIOstreamIn
         (  xfabc
         ,     MCLXIO_STREAM_ABC
            |  MCLXIO_STREAM_MIRROR
            |  MCLXIO_STREAM_SYMMETRIC
            |  (tab ? MCLXIO_STREAM_GTAB_RESTRICT : 0)
         ,  NULL
         ,  mclpMergeMax
         ,  &streamer
         ,  EXIT_ON_FAIL
         )
   ;  else
      mx = mclxReadx(xfmx, EXIT_ON_FAIL, MCLX_REQUIRE_GRAPH)

   ;  dom =    xfdom
            ?  mclxRead(xfdom, EXIT_ON_FAIL)
            :  NULL

   ;  if (write_mode == MY_OPT_WRITEGRAPH && !dom)
      mcxDie(1, me, "--write-graph requires -dom option")
   ;  else if (write_mode == MY_OPT_WRITEGRAPHC && !dom)
      mcxDie(1, me, "--write-graphc requires -dom option")
   ;  else if (dom && !MCLD_EQUAL(dom->dom_rows, mx->dom_cols))
      mcxDie(1, me, "domains not equal")

   ;  N_start = N_ROWS(mx)

   ;  if (xftabout || xftabxout)
      {  if (streamer.tab_sym_out)
         tab = streamer.tab_sym_out
      ;  else if (streamer.tab_sym_in)
         tab = streamer.tab_sym_in
      ;  else
         mcxDie(1, me, "no tab read, no tab created")
   ;  }

      if (tfting)
      {  mclgTF* tfar = mclgTFparse(NULL, tfting)
      ;  if (!tfar)
         mcxDie(1, me, "errors in tf-spec")
      ;  mclgTFexec(mx, tfar)
   ;  }

      /* clm close has three main modes.
       * - provide granularity info for different levels and optionally write clusterings
       * - output single-linkage join-order and join-values
       * - other modes: + different ways of outputting granularity info
       *                + block and block complemement networks
       *                + tab related things
       *   The latter two are perhaps fairly exploratory and may be considered for purge
      */

      if (hi_g)
      {  int i
      ;  mcxbool dedup = write_mode == MY_OPT_WRITESIZECOUNTS ? TRUE : FALSE

      ;  for (i=lo_g; i<= hi_g; i+=st_g)
         {  double cutoff = norm_g > 0.0 ? i / norm_g : 1.0 * i
         ;  dim prevsize = 0
         ;  dim n_same   = 1, j

         ;  mclxUnary(mx, fltxGQ, &cutoff)
         ;  mclx* mycc = clmComponents(mx, dom)

         ;  if (levels_pfx)
            {  mcxTing* name = mcxTingPrint(NULL, "%s.L%d", levels_pfx, (int) i)
            ;  mcxIO* xflevel = mcxIOnew(name->str, "w")
            ;  mcxIOopen(xflevel, EXIT_ON_FAIL)
            ;  mclxaWrite(mycc, xflevel, MCLXIO_VALUE_NONE, RETURN_ON_FAIL)
            ;  mcxIOclose(xflevel)
            ;  mcxIOfree(&xflevel)
            ;  mcxTingFree(&name)
         ;  }

            fprintf(xfout->fp, "%2d:", i)

         ;  if (dedup)
            {  for (j=0;j<N_COLS(mycc);j++)
               {  dim thissize = mycc->cols[j].n_ivps
               ;  if (thissize == prevsize)
                  n_same++
               ;  else
                  {  if (n_same > 1)
                     fprintf(xfout->fp, "(%d)", (int) n_same)
                  ;  n_same = 1
                  ;  fprintf(xfout->fp, " %lu", (ulong) thissize)
               ;  }
                  prevsize = thissize
            ;  }
               if (n_same > 1)
               fprintf(xfout->fp, "(%d)", (int) n_same)
         ;  }
            else
            for (j=0;j<N_COLS(mycc);j++)
            fprintf(xfout->fp, " %lu", (ulong) mycc->cols[j].n_ivps)

         ;  fputc('\n', xfout->fp)
         ;  mclxFree(&mycc)
      ;  }
         return STATUS_OK
   ;  }

                  /* We require a canonical domain.
                   * Make this a function.
                   * There is E log(E) factor due to edge sorting - Not an issue I think.
                   * The tree merge operation is done using linked lists; only the nodes
                   * in the smaller of the two children branches needs updating, so overall
                   * that is probably O(N).
                  */
      else if (sgl_g)
      {  dim L, U, D
      ;  mclxNrofEntriesLUD(mx, &L, &U, &D)
      ;  mcxTell(me, "Input matrix entries: lower=%d, upper=%d, diagonal=%d", (int) L, (int) U, (int) D)

      ;  dim i, e=0, E, N = MCX_MAX(U,L), n_linked = 0
      ;  mcle* edges       =  mcxAlloc(sizeof edges[0] * N, EXIT_ON_FAIL)
      ;  double sumszsq    =  N_COLS(mx)
      ;  mcxIO* xflist     =  mcxIOnew(fn_nodelist, "w")
      ;  mcxTing* upname   =  mcxTingNew("")
      ;  struct slnode *NODE  =  mcxNAlloc(N_COLS(mx), sizeof NODE[0], node_init, EXIT_ON_FAIL)
      ;  /* struct slleaf *LEAF  =  mcxNAlloc(N_COLS(mx), sizeof NODE[0], leaf_init, EXIT_ON_FAIL) */
      ;  int n_singleton = 0

      ;  if (!mclxDomCanonical(mx))
         mcxDie(1, me, "I need canonical domains in link mode")

      ;  mcxIOopen(xflist, EXIT_ON_FAIL)

      ;  for (i=0;i<N_COLS(mx);i++)
         {  NODE[i].lid = i
         ;  NODE[i].cid = i
         ;  NODE[i].size = 1
         ;  NODE[i].last = NODE+i
         ;  mcxTingPrint(NODE[i].name, "leaf_%d", (int) i)
      ;  }

         for (i=0;i<N_COLS(mx);i++)
         {  mclv* v = mx->cols+i
         ;  dim j
         ;  for (j=0;j<v->n_ivps;j++)
            {  mclp* d = v->ivps+j
            ;  if (d->idx > i)
               {  mcle* edge = edges+(e++)
               ;  edge->src = i
               ;  edge->dst = d->idx
               ;  edge->val = d->val
            ;  }
            }
         }
         E = e
      ;  mcxTell(me, "have %d edges ..", (int) E)
      ;  qsort(edges, E, sizeof edges[0], edge_val_cmp)
      ;  mcxTell(me, "sorted")
      ;  e = 0
      ;  n_linked = 1

      ;  fprintf
         (  xfout->fp
         ,  "link\tval\tNID\tANN\tBOB\txcsz\tycsz\txycsz\tnedge\tctr\tlss\tnsg\n"
         )
      ;  while (e<E)
         {  pnum s = edges[e].src      /* edge source node              */
         ;  pnum d = edges[e].dst      /* edge destination node         */
         ;  pval v = edges[e].val
         ;  pnum si = NODE[s].cid      /* source (cluster) index        */
         ;  pnum di = NODE[d].cid      /* destination (cluster) index   */

         ;  pnum ni = NODE[si].size >= NODE[di].size ? si : di          /* New Index (re-used)*/
         ;  pnum ui = ni == si ? di : si                                /* this one needs Updating */
         ;  char sbuf[50]
         ;  char dbuf[50]
         ;  dim j
         ;  e++

         ;  if (v < sgl_rcl_thr_g)
            {  mcxTell(me, "Threshold reached at link %d", (int) n_linked)
            ;  break
         ;  }

            if (si == di)                 /* already linked / same cluster */
            continue

         ;  snprintf(sbuf, 50, "%d", (int) s)
         ;  snprintf(dbuf, 50, "%d", (int) d)

         ;  {  dim sz1 = NODE[si].size
            ;  dim sz2 = NODE[di].size
            ;  dim lss_sub = MCX_MAX(NODE[si].lss, NODE[di].lss)
            ;  dim sgl_sub = NODE[si].nsg + NODE[di].nsg
            ;  dim sz_sum  = sz1 + sz2

            ;  sumszsq +=
                  (sz1 + sz2) * 1.0 * (sz1 + sz2)
               -  sz1 * 1.0 * sz1
               -  sz2 * 1.0 * sz2

            ;  if (sz1 == 1) fprintf(xflist->fp, "%s\t%.2f\n", tab ? mclTabGet(tab, s, NULL) : sbuf, v)
            ;  if (sz2 == 1) fprintf(xflist->fp, "%s\t%.2f\n", tab ? mclTabGet(tab, d, NULL) : dbuf, v)

            ;  mcxTingPrint(upname, "L%d_%d", (int) n_linked, (int) (sz_sum))

            ;  NODE[ni].lss = MCX_MAX(lss_sub, MCX_MIN(sz1, sz2))    /* overwrites si or di */
            ;  NODE[ni].nsg = sgl_sub + ((sz1 == 1) ^ (sz2 == 1))    /* overwrites si or di */

            ;  fprintf
               (  xfout->fp, "%d\t%.2f\t" "%s\t%s\t%s\t" "%d\t%d\t%d\t" "%.2f\t%.0f\t%lu\t%lu\n"
               ,  (int) n_linked, (double) v
               ,  upname->str, NODE[si].name->str, NODE[di].name->str
               ,  (int) sz1, (int) sz2, (int) sz_sum

               ,  (double) (e * 100.0 / E)
               ,  (0.5 + sumszsq / N_COLS(mx))
               ,  (long unsigned) NODE[ni].lss
               ,  (long unsigned) NODE[ni].nsg
               )

            ;  NODE[ni].size = sz1 + sz2
            ;  mcxTingWrite(NODE[ni].name, upname->str)
            ;
               {  struct slnode* node_ui = NODE+ui    /* all these dudes need setting to ni */
               ;  NODE[ni].last->next = node_ui       /* link them to the ni node as well */
               ;  NODE[ni].last       = node_ui->last /* the last ni node thus is the last ui node */

               ;  while (node_ui)
                  {  node_ui->cid = ni
                  ;  node_ui = node_ui->next          /* After this nothing points to ui anymore */
               ;  }
               }
               if (++n_linked == N_COLS(mx))
               break
         ;  }
         }

         if (E)
         mcxTell(me, "Finished linking at %.1f of edges", e * 100.0 / E)

      ;  for (i=0;i<N_COLS(mx);i++)
         {  struct slnode* node = NODE+i

                        /* Detect/write singletons: if a node is linked and has cid == lid
                         * then it will have a next; if it does not have a next and is linked
                         * then it will have cid != lid
                        */
         ;  if (node->cid == node->lid && ! node->next)
            {  char ibuf[50]
            ;  snprintf(ibuf, 50, "%d", (int) i)
            ;  fprintf(xflist->fp, "%s\t0.0\n", tab ? mclTabGet(tab, i, NULL) : ibuf)
            ;  fprintf
               (  xfout->fp, "%d\t%.2f\t" "sgl_%d\t%s\t%s\t" "%d\t%d\t%d\t" "%.2f\t%.0f\t%lu\t%lu\n"
               ,  (int) n_linked++
               ,  1000.0

               ,  (int) i
               ,  NODE[i].name->str
               ,  NODE[i].name->str

               ,  (int) 1
               ,  (int) 1
               ,  (int) 1

               ,  (double) (e * 100.0 / E)
               ,  (0.5 + sumszsq / N_COLS(mx))
               ,  (long unsigned) 0
               ,  (long unsigned) 0
               )
            ;  n_singleton++
         ;  }
         }
         if (n_singleton)
         mcxTell(me, "%d singletons in data", (int) n_singleton)
      ;  mcxIOclose(xflist)
      ;  mcxIOclose(xfout)
      ;  return STATUS_OK
   ;  }

      cc = make_symmetric ? clmComponents(mx, dom) : clmUGraphComponents(mx, dom)

                              /*
                               * thin out domain based on cc
                              */
   ;  if (ccbound_num)
      {  ccbound_cols = mclxColSelect(cc, mclv_check_ccbound, &ccbound_num)
      ;  ccbound_rows = mclgUnionv(cc, ccbound_cols, NULL, SCRATCH_READY, NULL)
   ;  }
      else
         ccbound_cols = cc->dom_cols
      ,  ccbound_rows = cc->dom_rows

   ;  N_bound = ccbound_rows->n_ivps

   ;  if
      (  canonical
      && (! (  map
            =  mclxMakeMap
               (  mclvClone(ccbound_rows)
               ,  mclvCanonical(NULL, ccbound_rows->n_ivps, 1.0)
               )
            )
         )
      )
      mcxDie(1, me, "cannot make a map")

   ;  if (N_bound < N_start)
      ccbound = mclxSub(cc, ccbound_cols, ccbound_rows)
   ;  else
      ccbound = cc

   ;  if (xfmxout)
      {  if (N_bound < N_start)      /* thin out matrix */
         {  mclx* sub = mclxSub(mx, ccbound_rows, ccbound_rows)
         ;  mclxFree(&mx)
         ;  mx = sub
      ;  }

         if (map)
         {  if (mclxMapRows(mx, map))
            mcxDie(1, me, "cannot map rows")

         ;  if (mclxMapCols(mx, map))
            mcxDie(1, me, "cannot map cols")
      ;  }
         mclxWrite(mx, xfmxout, MCLXIO_VALUE_GETENV, RETURN_ON_FAIL)
   ;  }

      if (xftabxout)
      {  mclv* deselect = mcldMinus(tab->domain, ccbound->dom_rows, NULL)
      ;  if (canonical)
         mcxErr(me, "--canonical and writing tab not yet implemented. beerware.")
      ;  else
         mclTabWrite(tab, xftabxout, deselect, RETURN_ON_FAIL)
      ;  mclvFree(&deselect)
   ;  }


      if (xftabout)
      {  mclTab* tabsel = mclTabSelect(tab, ccbound->dom_rows), *tabout
      ;  if (map)
            tabout = mclTabMap(tabsel, map)
         ,  mclTabFree(&tabsel)
      ;  else
         tabout = tabsel
      ;  if (!tabout)
         mcxDie(1, me, "no tab, baton")
      ;  mclTabWrite(tabout, xftabout, NULL, RETURN_ON_FAIL)

      ;  mclTabFree(&tabout)
   ;  }

      if (map)
      {  if (mclxMapRows(ccbound, map))
         mcxDie(1, me, "cannot map rows")

      ;  if (mclxMapCols(ccbound, NULL))
         mcxDie(1, me, "cannot map cols")
   ;  }

      if (write_mode == MY_OPT_WRITEGRAPH)
      {  mclx* bl = mclxBlockUnion(mx, cc)
      ;  mclxWrite(bl, xfout, MCLXIO_VALUE_GETENV, RETURN_ON_FAIL)
      ;  mclxFree(&mx)
      ;  mx = bl
   ;  }
      else if (write_mode == MY_OPT_WRITEGRAPHC)
      {  mclx* bl = mclxBlocksC(mx, cc)
      ;  mclxWrite(bl, xfout, MCLXIO_VALUE_GETENV, RETURN_ON_FAIL)
      ;  mclxFree(&mx)
      ;  mx = bl
   ;  }
      if (write_mode == MY_OPT_WRITECC)
      {  if (streamer.tab_sym_out)
         {  mclxIOdumper dumper
         ;  mclxIOdumpSet(&dumper, MCLX_DUMP_LINES | MCLX_DUMP_NOLEAD, NULL, NULL, NULL)
         ;  mclxIOdump
            (  ccbound
            ,  xfout
            ,  &dumper
            ,  NULL
            ,  streamer.tab_sym_out
            ,  MCLXIO_VALUE_NONE
            ,  EXIT_ON_FAIL
            )
      ;  }
         else
         mclxaWrite(ccbound, xfout, MCLXIO_VALUE_NONE, RETURN_ON_FAIL)
   ;  }

      else if (write_mode == MY_OPT_WRITECOUNT)
      fprintf(xfout->fp, "%lu\n", (ulong) N_COLS(ccbound))

   ;  else if (write_mode == MY_OPT_WRITESIZES || write_mode == MY_OPT_WRITESIZECOUNTS)
      {  dim j 
      ;  mcxbool dedup = write_mode == MY_OPT_WRITESIZECOUNTS ? TRUE : FALSE
         
      ;  if (dedup)
         {  dim prevsize = 0
         ;  dim n_same   = 1
         ;  for (j=0;j<N_COLS(ccbound);j++)
            {  dim thissize = ccbound->cols[j].n_ivps
            ;  if (thissize == prevsize)
               n_same++
            ;  else
               {  if (n_same > 1)
                  fprintf(xfout->fp, "(%d)", (int) n_same)
               ;  n_same = 1
               ;  fprintf(xfout->fp, "%s%lu", j ? " " : "", (ulong) thissize)
            ;  }
               prevsize = thissize
         ;  }
            if (n_same > 1)
            fprintf(xfout->fp, "(%d)", (int) n_same)
      ;  }
         else
         for (j=0;j<N_COLS(ccbound);j++)
         fprintf(xfout->fp, "%s%lu", j ? " " : "", (ulong) ccbound->cols[j].n_ivps)
      ;  fputc('\n', xfout->fp)
   ;  }

      if (xfmapout && map)
      mclxaWrite(map, xfmapout, MCLXIO_VALUE_NONE, RETURN_ON_FAIL)

   ;  mcxIOfree(&xfmx)
   ;  mcxIOfree(&xfabc)
   ;  mcxIOfree(&xfout)

   ;  mclTabFree(&(streamer.tab_sym_out))

   ;  mclxFree(&mx)
   ;  mclxFree(&cc)
   ;  mclxFree(&dom)
   ;  return STATUS_OK
;  }


mcxDispHook* mcxDispHookClose
(  void
)
   {  static mcxDispHook closeEntry
   =  {  "close"
      ,  "close [options] -imx <mx file>"
      ,  closeOptions
      ,  sizeof(closeOptions)/sizeof(mcxOptAnchor) - 1
      ,  closeArgHandle
      ,  closeInit
      ,  closeMain
      ,  0
      ,  0
      ,  MCX_DISP_MANUAL
      }
   ;  return &closeEntry
;  }


