      double precision function fsimpson(xa,xb,c)
      implicit none
      real*8 xa,xb,c
c
c     integrate dsqrt(c^2-x^2)/x from xa to xb
c
      integer*4 i,n
      real*8 x,dx,y,f1,f2
c
      integer*4 nmax
      real*8 eps
      data nmax/2048/
      data eps/1.0d-03/
c
      n=1
      dx=xb-xa
      y=0.5d0*(dsqrt(c*c-xa*xa)/xa+dsqrt(c*c-xb*xb)/xb)
      f1=0.d0
      f2=y*dx
10    n=2*n
      dx=0.5d0*dx
      do i=1,n-1,2
        x=xa+dble(i)*dx
        y=y+dsqrt(c*c-x*x)/x
      enddo
      f2=y*dx
      if(dabs(f2-f1).gt.eps*f2.and.n.lt.nmax)then
        f1=f2
        goto 10
      endif
      fsimpson=f2
c
      return
      end