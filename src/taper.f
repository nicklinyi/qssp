      subroutine taper(ldeg,ldegmax,tap)
      implicit none
      integer*4 ldeg(4),ldegmax
      real*8 tap(0:ldegmax)
c
      integer*4 l
      real*8 fac
c
      real*8 PI
      data PI/3.14159265358979d0/
c
      do l=0,ldeg(1)-1
        tap(l)=0.d0
      enddo
      fac=0.5d0*PI/dble(1+ldeg(2)-ldeg(1))
      do l=ldeg(1),ldeg(2)-1
        tap(l)=dsin(fac*dble(l-ldeg(1)))**2
      enddo
      do l=ldeg(2),ldeg(3)
        tap(l)=1.d0
      enddo
      fac=0.5d0*PI/dble(1+ldeg(4)-ldeg(3))
      do l=ldeg(3)+1,ldeg(4)
        tap(l)=dcos(fac*dble(l-ldeg(3)))**2
      enddo
      do l=ldeg(4)+1,ldegmax
        tap(l)=0.d0
      enddo
c
      return
      end
