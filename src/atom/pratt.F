c
c Copyright (c) 1998-2012 The OPIUM Group
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c
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
