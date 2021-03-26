      function pratt(d1,d2,d3,d4)
      implicit real*8 (a-h,o-z)

c pratt method applied to accelerate convergence of sc potential
c see f. herman and s. skillman, atomic structure calculations.
c the value of pratt must lie between 0 and 0.5.

      x1 = d1+d4
      x2 = d2+d3
      if(abs(x1-x2).lt.0.0001) goto 7
      a = (d4-d2)/(x1-x2)
      if(a.lt.0.) goto8
      if(a.lt.0.5) goto9
    7 a = 0.5
      goto 9
    8 a = 0.
    9 pratt = a
      return
      end
