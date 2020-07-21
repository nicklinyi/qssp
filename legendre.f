      subroutine legendre(raddis,plm,ldeg,ldegmax)
      implicit none
      integer*4 ldeg,ldegmax
      real*8 raddis
      real*8 plm(0:ldegmax,0:2)
c
c     calculate Plm(l,m,x)/(1-x^2)^(m/2)
c     where Plm are the associated Legendre polynomials
c
      integer*4 l,m
      real*8 x
c
      x=dcos(raddis)
      plm(0,0)=1.d0
      plm(0,1)=0.d0
      plm(1,1)=1.d0
      plm(0,2)=0.d0
      plm(1,2)=0.d0
      plm(2,2)=3.d0
      do m=0,2
        plm(m+1,m)=dble(2*m+1)*x*plm(m,m)
        do l=m+2,ldeg
          plm(l,m)=(dble(2*l-1)*x*plm(l-1,m)
     &             -dble(l+m-1)*plm(l-2,m))/dble(l-m)
        enddo
      enddo
      return
      end
