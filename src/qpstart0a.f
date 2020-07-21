      subroutine qpstart0a(ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c     for n = 0
c
      integer*4 ly,jh
      complex*16 cy(3)
c
      integer i
      real*8 peps
      complex*16 c2mu,cxi,cxp2,ca,cb,ya,yb
c
      real*8 eps
      complex*16 c2
      data eps/1.0d-08/
      data c2/(2.d0,0.d0)/
c
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
      cxp2=(kp(ly)*crrup(ly))**2
c
      if(jh.eq.1.and.ly.ge.lys)then
        cy(1)=zh12p(0,1,ly)
        cy(2)=cxi*cxp2-c2*c2mu*zh12p(0,1,ly)
        cy(3)=-cga(ly)
        peps=dabs(dimag(zh12p(0,1,ly)))/cdabs(zh12p(0,1,ly))
        if(peps.lt.eps)then
          ca=dcmplx((peps/eps)**2,0.d0)
          cb=(1.d0,0.d0)-ca
          ya=cy(1)
          yb=mas3x3up(1,1,ly)
          do i=1,3
            cy(i)=ca*cy(i)*yb+cb*mas3x3up(i,1,ly)*ya
          enddo
        endif
      else if(jh.eq.2.and.ly.le.lys)then
        cy(1)=zh12p(0,2,ly)
        cy(2)=cxi*cxp2-c2*c2mu*zh12p(0,2,ly)
        cy(3)=-cga(ly)
      else
        stop ' Error in qpstart0a: bad jh value!'
      endif
      return
      end