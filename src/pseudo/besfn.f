      function besfn(x,l)
      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

      dimension f(nang0)
      xs = 1.0
      if (x.lt.0.0) then
         xs = -1.0
      endif

c Note: if x<0 this function now returns il(abs(x)), whereas is
c          x>0 this function now returns jl(x)
c My convention is that il(x) = i**-l * jl(ix)
      if (x.eq.0.0d0) then
         besfn = 1.0d0
         if (l.gt.0) besfn = 0.0d0
         if (l.lt.0) then
            write (7,*)'besfn(0.0,',l,') is infinite.'
            stop
         endif
         return
      endif

      if (abs(x).lt.0.5d0) then
         prod = 1.0d0
         sum = 0.0d0
         do 3 j = 1,l
            prod = prod * dfloat(j+j+1)
 3       continue 
         prod = abs(x)**l/prod
         xmult = -1.0d0 * xs * x * x/2.0d0
         sum = sum + prod
         do 2 j = 1,6
            prod = prod * xmult/dfloat(j)/dfloat(j+j+l+l+1)
            sum = sum + prod
 2       continue 
         besfn = sum
      else
         f(1) = cos(x)/x
         f(2) = sin(x)/x
         if (xs.lt.0.0) then
            f(1) = cosh(x)/abs(x)
            f(2) = sinh(x)/x
         endif
         do 1 i = 3,l+2
            f(i) = (dfloat(i + i - 5)/abs(x) * f(i-1) - f(i-2)) * xs
 1       continue
         besfn = f(l+2)
      endif 

      return
      end
      
      
      function besder(x,l)
      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

      xs = 1.0
      if (x.lt.0.0) then
         xs = -1.0
      endif
c Note: if x<0 this function now returns il'(abs(x)), whereas if
c          x>0 this function now returns jl'(x)
c Within this comment, ' means d/dx
c My convention is that il(x) = i**-l * jl(ix)
      if (x.eq.0.0d0) then
         if (l.gt.1) besder = 0.0d0
         if (l.eq.0) besder = 0.0d0
         if (l.eq.1) besder = 1.0d0/3.0d0
         if (l.lt.0) then
            write (7,*) 'besder(0.0,',l,') is infinite.'
            stop
         endif
         return
      endif
      xl = dfloat(l)
      temp = xl * besfn(x,l-1) - xs * (xl + 1.0d0) * besfn(x,l+1)
      besder = temp/(xl + xl + 1.0d0)
      return
      end
      
      
      function besder2(x,l)
      implicit double precision(a-h,o-z)
      
#include "PARAMOPT"

      xs = 1.0
      if (x.lt.0.0) then
         xs = -1.0
      endif
c Note: if x<0 this function now returns il'(abs(x)), whereas if
c          x>0 this function now returns jl'(x)
c Within this comment, ' means d/dx
c My convention is that il(x) = i**-l * jl(ix)
c To get this recursion relation, use Laplacian in spherical coords.
      xl = dfloat(l)
      temp = (xl * (xl + 1.0d0) /x/x - xs) * besfn(x,l)
      besder2 = temp - 2.0d0/abs(x) * besder(x,l)
      return
      end


