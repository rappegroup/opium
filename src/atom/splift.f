      subroutine splift(x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)
      implicit double precision(a-h,o-z)
c
c
c     Revision 1.1  89/10/26  19:53:45  sverre
c     Initial revision
c     
c     sandia mathematical program library
c     applied mathematics division 2613
c     sandia laboratories
c     albuquerque, new mexico  87185
c     control data 6600/7600  version 7.2  may 1978
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    issued by sandia laboratories                     *
c  *                   a prime contractor to the                       *
c  *                united states department of energy                 *
c  * * * * * * * * * * * * * * * notice  * * * * * * * * * * * * * * * *
c  * this report was prepared as an account of work sponsored by the   *
c  * united states government.  neither the united states nor the      *
c  * united states department of energy nor any of their employees,    *
c  * nor any of their contractors, subcontractors, or their employees  *
c  * makes any warranty, express or implied, or assumes any legal      *
c  * liability or responsibility for the accuracy, completeness or     *
c  * usefulness of any information, apparatus, product or process      *
c  * disclosed, or represents that its use would not infringe          *
c  * owned rights.                                                     *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * the primary document for the library of which this routine is     *
c  * part is sand77-1441.                                              *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     written by rondall e. jones
c
c     abstract
c         splift fits an interpolating cubic spline to the n data points
c         given in x and y and returns the first and second derivatives
c         in yp and ypp.  the resulting spline (defined by x, y, and
c         ypp) and its first and second derivatives may then be
c         evaluated using splint.  the spline may be integrated using
c         spliq.  for a smoothing spline fit see subroutine smoo.
c
c     description of arguments
c         the user must dimension all arrays appearing in the call list,
c         e.g.   x(n), y(n), yp(n), ypp(n), w(3n)
c
c       --input--
c
c         nout - unit for error messages
c         x    - array of abscissas of data (in increasing order)
c         y    - array of ordinates of data
c         n    - the number of data points.  the arrays x, y, yp, and
c                ypp must be dimensioned at least n.  (n .ge. 4)
c         isx  - must be zero on the initial call to splift.
c                if a spline is to be fitted to a second set of data
c                that has the same set of abscissas as a previous set,
c                and if the contents of w have not been changed since
c                that previous fit was computed, then isx may be
c                set to one for faster execution.
c         a1,b1,an,bn - specify the end conditions for the spline which
c                are expressed as constraints on the second derivative
c                of the spline at the end points (see ypp).
c                the end condition constraints are
c                        ypp(1) = a1*ypp(2) + b1
c                and
c                        ypp(n) = an*ypp(n-1) + bn
c                where
c                        abs(a1).lt. 1.0  and  abs(an).lt. 1.0.
c
c                the smoothest spline (i.e., least integral of square
c                of second derivative) is obtained by a1=b1=an=bn=0.
c                in this case there is an inflection at x(1) and x(n).
c                if the data is to be extrapolated (say, by using splint
c                to evaluate the spline outside the range x(1) to x(n)),
c                then taking a1=an=0.5 and b1=bn=0 may yield better
c                results.  in this case there is an inflection
c                at x(1) - (x(2)-x(1)) and at x(n) + (x(n)-x(n-1)).
c                in the more general case of a1=an=a  and b1=bn=0,
c                there is an inflection at x(1) - (x(2)-x(1))*a/(1.0-a)
c                and at x(n) + (x(n)-x(n-1))*a/(1.0-a).
c
c                a spline that has a given first derivative yp1 at x(1)
c                and ypn at y(n) may be defined by using the
c                following conditions.
c
c                a1=-0.5
c
c                b1= 3.0*((y(2)-y(1))/(x(2)-x(1))-yp1)/(x(2)-x(1))
c
c                an=-0.5
c
c                bn=-3.0*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)/(x(n)-x(n-1))
c
c       --output--
c
c         yp   - array of first derivatives of spline (at the x(i))
c         ypp  - array of second derivatives of spline (at the x(i))
c         ierr - a status code
c              --normal code
c                 1 means that the requested spline was computed.
c              --abnormal codes
c                 2 means that n, the number of points, was .lt. 4.
c                 3 means the abscissas were not strictly increasing.
c
c       --work--
c
c         w    - array of working storage dimensioned at least 3n.
      dimension x(n),y(n),yp(n),ypp(n),w(n,3)
c
      if (n.lt.4) go to 200
      nm1  = n-1
      nm2  = n-2
      if (isx.gt.0) go to 40
      do 5 i=2,n
      if (x(i)-x(i-1)) 300,300,5
    5 continue
c
c     define the tridiagonal matrix
c
      w(1,3) = x(2)-x(1)
      do 10 i=2,nm1
      w(i,2) = w(i-1,3)
      w(i,3) = x(i+1)-x(i)
   10 w(i,1) = 2.D0*(w(i,2)+w(i,3))
      w(1,1) = 4.D0
      w(1,3) =-4.D0*a1
      w(n,1) = 4.D0
      w(n,2) =-4.D0*an
c
c     l u decomposition
c
      do 30 i=2,n
      w(i-1,3) = w(i-1,3)/w(i-1,1)
   30 w(i,1)   = w(i,1) - w(i,2)*w(i-1,3)
c
c     define *constant* vector
c
   40 ypp(1) = 4.D0*b1
      dold   = (y(2)-y(1))/w(2,2)
      do 50 i=2,nm2
      dnew   = (y(i+1) - y(i))/w(i+1,2)
      ypp(i) = 6.D0*(dnew - dold)
      yp(i)  = dold
   50 dold   = dnew
      dnew   = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(nm1) = 6.D0*(dnew - dold)
      ypp(n) = 4.D0*bn
      yp(nm1)= dold
      yp(n)  = dnew
c
c     forward substitution
c
      ypp(1) = ypp(1)/w(1,1)
      do 60 i=2,n
   60 ypp(i) = (ypp(i) - w(i,2)*ypp(i-1))/w(i,1)
c
c     backward substitution
c
      do 70 j=1,nm1
      i = n-j
   70 ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
c
c     compute first derivatives
c
      yp(1)  = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2.D0*ypp(1)
     1         + ypp(2))/6.D0
      do 80 i=2,nm1
   80 yp(i)  = yp(i) + w(i,2)*(ypp(i-1) + 2.D0*ypp(i))/6.D0
      yp(n)  = yp(n) + (x(n)-x(nm1))*(ypp(nm1) + 2.D0*ypp(n))/6.D0
c
      ierr = 1
      return
  200 ierr = 2
      write(nout,210)
  210 format(' in splift, there were less than 4 data values.')
      return
  300 ierr = 3
      write(nout,310)
  310 format(' in splift,',
     1' the abscissas were not strictly increasing.')
      return
      end
