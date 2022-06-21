/*   (C) Copyright 2006, 2007, 2008, 2009, 2010, 2011 Stijn van Dongen
 *   (C) Copyright 2012, 2013 Stijn van Dongen
 *
 * This file is part of MCL.  You can redistribute and/or modify MCL under the
 * terms of the GNU General Public License; either version 3 of the License or
 * (at your option) any later version.  You should have received a copy of the
 * GPL along with MCL, in the file COPYING.
*/


#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>
#include <time.h>

#include "impala/io.h"
#include "impala/matrix.h"
#include "impala/tab.h"
#include "impala/stream.h"
#include "impala/ivp.h"
#include "impala/compose.h"
#include "impala/vector.h"
#include "impala/edge.h"
#include "impala/app.h"
#include "impala/iface.h"
#include "impala/pval.h"

#include "tingea/types.h"
#include "tingea/ding.h"
#include "tingea/ting.h"
#include "tingea/io.h"
#include "tingea/err.h"
#include "tingea/equate.h"
#include "tingea/minmax.h"
#include "tingea/rand.h"
#include "tingea/array.h"
#include "tingea/opt.h"

#include "mcl/proc.h"
#include "mcl/procinit.h"
#include "mcl/alg.h"

#include "clew/clm.h"
#include "clew/cat.h"

const char* usagelines[];

const char* me = "mcxtest";

void usage
(  const char**
)  ;

static mclv* reduce_v
(  const mclv* u
)
   {  mclv* v = mclvClone(u)
   ;  dim n = v->n_ivps
   ;  double s = mclvSum(v)
   ;  double sq = mclvPowSum(v, 2.0)
   ;  if (s)
      mclvSelectGqBar(v, 0.25 * sq / s)
;fprintf(stderr, "from %d to %d entries\n", (int) n, (int) v->n_ivps)
   ;  return v
;  }


static dim update_meet_zip
(  mclVector*  v1
,  mclp* ivp1max
,  const mclVector*  v2
)
   {  mclp  *ivp1 = v1->ivps
         ,  *ivp2 = v2->ivps
         ,  *ivp2max = ivp2 + v2->n_ivps

   ;  dim n_zeroed = 0
   ;  while (ivp1 < ivp1max && ivp2 < ivp2max)
      {  if (ivp1->idx < ivp2->idx)
         ivp1++
      ;  else if (ivp1->idx > ivp2->idx)
         ivp2++
      ;  else
            ivp1++->val = 0
         ,  ivp2++
         ,  n_zeroed++
   ;  }
      return n_zeroed
;  }


static dim update_todo
(  mclv* todo
,  mclp* end
,  const mclv* strike
)
   {  dim i, n_zero = 0
   ;  mclp* p = strike->ivps, *pmax = p + strike->n_ivps
   ;  while (p < pmax)
      {  if (todo->ivps[p->idx].val)
         n_zero++
      ;  todo->ivps[p++->idx].val = 0.0
   ;  }
      return n_zero
;  }


static double mclv_inner
(  const mclv* a
,  const mclv* b
,  dim N
)
   {  dim j
   ;  double ip = 0.0

   ;  if (a->n_ivps < N || b->n_ivps < N)
      return mclvIn(a, b)

   ;  for (j=0;j<a->n_ivps;j++)
      ip += a->ivps[j].val * b->ivps[j].val
   ;  return ip
;  }


static double pearson
(  const mclv* v1
,  const mclv* v2
,  const mclv* sums
,  const mclv* Nssq
,  double N
)
   {  ofs    c       =  v1->vid
   ;  ofs    d       =  v2->vid 
   ;  double s1      =  sums->ivps[c].val
   ;  double s2      =  sums->ivps[d].val
   ;  double Nsq1    =  Nssq->ivps[c].val
   ;  double Nsq2    =  Nssq->ivps[d].val

   ;  double nom     =  sqrt(Nsq1 - s1 * s1) *  sqrt(Nsq2 - s2 * s2)
   ;  return nom ? (N* mclv_inner(v1, v2, N) - s1*s2) / nom : 0.0
;  }


static void add_edge
(  mclx* m
,  ofs c
,  ofs d
,  double v
)
   {  mclp* p1 = mclxInsertIvp(m, c, d)
   ;  mclp* p2 = mclxInsertIvp(m, d, c)
   ;  if (p1)
      p1->val = v
   ;  if (p2)
      p2->val = v
;  }



static double fltxCos
(  pval     v
,  void*    not_used
)
   {  return cos(v)
;  }


int main
(  int                  argc
,  const char*          argv[]
)  
   {  mclx* tbl, *nw, *buffer_nw, *buffer_long
   ;  mclv* Nssqs, *sums, *scratch_nw, *scratch_long, *todo, *order
   ;  mcxIO* xfdata
   ;  dim i, c, d
   ;  dim N
   ;  mcxIO* xfout = mcxIOnew("tst2", "w")
   ;  mclv themclv
   ;  mclp themclp
   ;  dim n_pearson = 0

   ;  double THRESHOLD = acos(0.7)

   ;  if (argc != 3) mcxDie(1, me, "need <matrix> <pearson threshold>")

   ;  xfdata = mcxIOnew(argv[1], "r")
   ;  tbl = mclxRead(xfdata, EXIT_ON_FAIL)

   ;  for (i=0;i<N_COLS(tbl);i++)
      mclvCanonicalEmbed(tbl->cols+i, tbl->cols+i, N_ROWS(tbl), 0.0)

   ;  THRESHOLD = acos(atof(argv[2]))

   ;  N = N_ROWS(tbl)

   ;  nw    = mclxAllocZero(mclvClone(tbl->dom_cols), mclvClone(tbl->dom_cols))

   ;  Nssqs    =  mclvCopy(NULL, tbl->dom_cols)
   ;  sums     =  mclvCopy(NULL, tbl->dom_cols)

   ;  for (c=0;c<N_COLS(tbl);c++)
      {  Nssqs->ivps[c].val = N * mclvPowSum(tbl->cols+c, 2.0)
      ;  sums->ivps[c].val = mclvSum(tbl->cols+c)
   ;  }

      for (i=0; i<N_COLS(nw); i++)
      {  dim new = i, j

      ;  for (j=0;j<i;j++)
         {  double p, arc

         ;  dim k = i - j - 1

         ;  p = pearson(tbl->cols+k, tbl->cols+new, sums, Nssqs, N * 1.0)
         ;  n_pearson++

         ;  if (p > 1.0)
            p = 1.0
         ;  arc = acos(p)
                                                /* fixme; k and new both need not be computed */
         ;  if (arc <= THRESHOLD)
            {  dim x

            ;  mclv* todo = nw->cols+new
            ;  mclgUnionv(nw, nw->cols+k, NULL, SCRATCH_READY, todo)

            ;  mclvInsertIdx(todo, new, 1.0)             /* will be translated later (need to store zeroes) */
            ;  mclvInsertIdx(todo, k, 1.0 + arc)         /* will be translated later */
            ;  mclvInsertIdx(nw->cols+k, new, arc)

            ;  for (x=0; x<todo->n_ivps; x++)
               {  double q, arc2
               ;  dim id = todo->ivps[x].idx
               ;  if (id == k || id == new)
                  continue
               ;  q = pearson(tbl->cols+id, tbl->cols+new, sums, Nssqs, N * 1.0)
               ;  n_pearson++
               ;  if (q > 1.0)
                  q = 1.0
               ;  arc2 = acos(q)
               ;  todo->ivps[x].val = arc2 <= THRESHOLD ? 1.0 + arc2 : 0.0    /* fixme tricky */
               ;  if (arc2 <= THRESHOLD)
                  mclvInsertIdx(nw->cols+id, new, arc2)
            ;  }
               mclvUnary(todo, fltxPositive, NULL)
            ;  mclvAffine(todo, 1.0, 0.0)
            ;  break
         ;  }
         }
         if (j == i) /* nothing hit */
         mclvInsertIdx(nw->cols+new, new, 0.0)
   ;  }

      mclxUnary(nw, fltxCos, NULL)
   ;  mclxWrite(nw, xfout, 8, RETURN_ON_FAIL)
   ;  mcxTell(me, "reduction %.3f", n_pearson * 2.0 / (1.0 * N_COLS(nw) * (N_COLS(nw)-1)))
   ;  return 0
;  }


const char* usagelines[] =
{  NULL
}  ;



