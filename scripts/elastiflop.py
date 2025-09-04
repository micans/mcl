
# Instigated and inspired by https://github.com/micans/mcl/discussions/34, opened by Zhitian-wu
#
#  Zhitia-wu found the following interesting class of states that have the remarkable
#  property of converging to a stable flip-flop state for a large set of starting conditions, with the equally
#  remarkable property of adaptation that when changing inflation those states converge to a different stable orbit.
#  Note that these states are not stable under perturbations of matrix entries.
#
#  Start with the following matrix, where x is initially quite a bit smaller than 1/2n (or zero).
#
#        <---     n + 1 copies    -->                          
#                                                              
#        1/n   1/n       ...      1/n       x         ^
#        1/n   1/n       ...      1/n       x         |        
#         ..    ..       ...       ..      ..         n copies 
#         ..    ..       ...       ..      ..         |        
#        1/n   1/n       ...      1/n       x         v        
#          -     -       ...        -       y      [ here y = ½ (1-nx) ]
#          -     -       ...        -       y
#
#  Then compute mcl process on this matrix. Note that it is defined by n and x (y is determined by x or vice versa)
#
#  For expansion and inflation, only the positions with x and y values need to be computed, all elements 1/n do not change.
#  Hence we only need a single expression for updating x and a single expression for updating y.
#  Computationally it is slightly more convenient to keep track of y rather than express it always as a function of x.
#  For some insight into expansion it is however useful to express the new x' as a function of x:
#
#  Expansion:
#    x' = x + y/n + x * y        = x + y(1/n + x) =  x + ½ (1 - nx)(1/n + x)
#    y' = y * y
#
#  So
#
#    x' =  x + ½ (1 - nx)(1/n + x)
#
#       =  x + ½ (1/n +x -x -n x^2)
#
#       =  x + 1/2n - ½ n x^2
#
#       =  1/2n +  (x - ½ n x^2)       where (x - ½ n x^2) is << 1/2n  if x << 1/2n
#
#       =  1/2n + x + O(x^2) 
#
#  Thus for very small x we get that x' is approxiately x away from 1/2n.
#  Then for high inflation values x'' (the inflated value of the expansion value x') will be very small again,
#  but crucially x' is about 1/2n and y' is about 1/4, so x'' is determined by n and the inflation value.
#  Whatever x'' ends up as, the next expansion value will be about x'' away again from 1/2n.
#
#  These counteracting forces settle into a stable orbit dependably and from various initial conditions.
#  Inflation can be changed, then expansion/inflation converge to again cancel each other in a different stable orbit.
#
#  Observed: Define k such that 1/k = x' (the result from expansion). Per above k gets close to 2n
#
#  - If we consider the difference value q = 1/(n-k'/2) then the ratio q(I+delta) / q(I) between the stable values for k'
#    as we increase inflation by delta seems to approach (N/2)**(delta).
#
#  Arguments to this program:
#
#     n       as above
#     k       x will be initialised as 1/k, unless k=0 in which case x is set to 0
#     infl    initial inflation value
#     lim     number of expansion/inflation cycles to compute
#     delta   inflation will change every 5 cycles. It will initially increment upwards by delta,
#              then slightly over half-way it will decrement by delta (every 5 cycles).
# 
# Example invocations:
# python3 elastiflop.py 31 215 2.0 50 1.0
#       achieves moving trampoline with inflation 2, higher N and k
# 
# python3 elastiflop.py 9 135 4.0 50 1.0
#       inflation 4, lower N and k
# 
# python3 elastiflop.py 5 35 8.0 50 1.0
#       inflation 8, lower N and k


import sys
import math

N, kstart, infl, lim, delta = 9, 130, 6, 100, 1.0

try:
    N         = int(sys.argv[1])
    kstart    = int(sys.argv[2])
    infl      = float(sys.argv[3])
    lim       = int(sys.argv[4])
    delta     = float(sys.argv[5])

except IndexError:
    print("Provide <N> <kstart> <start-inflation> <lim> - proceeding with pre-sets")

x = 1/kstart if kstart > 0 else 0.0
y = 0.5 - (N * x / 2.0)
ite = 0

sum = x * N + y * 2

print("At start x=%f y=%f sum=%f k=%d" % (x, y, sum, kstart))

qprev = 1.0
up = True

while ite < lim:
   if ite % 5 == 0:
      print("Inflation set to %.2f" % infl)
   else:
      print("---")

   # Compute expansion, new values for x and y
   x2 = x + y / N + x * y
   y2 = y * y
   x2b = x + 0.5 * (1.0/N - N*x*x)              # x computed as a function of x only.
   sum2 = x2 * N + y2 * 2                       # display as check

   try:
      contraction = math.log(y / x) / math.log(y2 / x2)
      kinv = 1/x2
      q = 1/(N-kinv/2)
      qratio = q / qprev if up else qprev / q
   except (ValueError, ZeroDivisionError):
      contraction, kinv, q, qratio = 0, 0, 0, 0
   print("Exp ite=%d x=%.16f y=%.16f sum=%f contraction=%f k'=%.10f kinv2Ndiff=%.4f ratio=%.10f" % (ite, x2, y2, sum2, contraction, kinv, q, qratio))
   # test = x2 - x2b                            # An old test. It checked out.
   # print("Test=%.12f" % test)

   # Inflate the expanded values ...
   x3 = x2 ** infl
   y3 = y2 ** infl

   sum = x3 * N + y3 * 2         # use for normalisation

   x4 = x3 / sum
   y4 = y3 / sum

   # Compute inflation from result; this is obviously the same as infl,
   # this is to show the computation is completely analogous to contraction above;
   # the output shows that contraction/expansion works as inverse to inflation.
   try:
      inflation = math.log(y4 / x4) / math.log(y2 / x2) 
      kinv = 1/x4
   except ValueError:
      inflation = 0
      kinv = 0

   sum4 = x4 * N + y4 * 2         # display as check

   print("Inf ite=%d x=%.16f y=%.16f sum=%f   inflation=%f k'=%.10f" % (ite, x4, y4, sum4, inflation, kinv))

   ite += 1
   if ite % 5 == 0:
      if 2 * ite - 5 > lim:
         infl -= delta
         up = False
      else:
         infl += delta
      print("")
      qprev = q

   x,y = x4, y4

