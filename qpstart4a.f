      subroutine qpstart4a(ldeg,ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly,jh
      complex*16 cy(4,2)
c
      integer*4 i
      real*8 peps
      complex*16 cldeg,cxp2,c2lp1,c2lp3,ca,cb,ya,yb
c
      real*8 eps
      complex*16 c0,c1,c2,c3
      data eps/1.0d-08/
      data c0,c1,c2,c3/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      cxp2=(kp(ly)*crrup(ly))**2
      c2lp1=c2*cldeg+c1
      c2lp3=c2*cldeg+c3
c
      if(jh.eq.1.and.ly.ge.lys)then
        cy(1,1)=cldeg-zh12p(ldeg,1,ly)
        cy(2,1)=-cro(ly)*comi2*crrup(ly)**2
        cy(3,1)=cga(ly)
        cy(4,1)=cga(ly)*(cldeg+c1)
c
        peps=dabs(dimag(zh12p(ldeg,1,ly)))/cdabs(zh12p(ldeg,1,ly))
        if(peps.lt.eps)then
          ca=dcmplx((peps/eps)**2,0.d0)
          cb=(1.d0,0.d0)-ca
          ya=cy(1,1)
          yb=mas4x4up(1,1,ly)
          do i=1,4
            cy(i,1)=ca*cy(i,1)*yb+cb*mas4x4up(i,1,ly)*ya
          enddo
        endif
        cy(1,2)=c0
        cy(2,2)=c0
        cy(3,2)=c1
        cy(4,2)=c2lp1
      else if(jh.eq.2.and.ly.le.lys)then
        cy(1,1)=cldeg-zh12p(ldeg,2,ly)
        cy(2,1)=-cro(ly)*comi2*crrup(ly)**2
        cy(3,1)=cga(ly)
        cy(4,1)=cga(ly)*(cldeg+c1)
c
        cy(1,2)=c0
        cy(2,2)=c0
        cy(3,2)=c1
        cy(4,2)=c0
      else
        stop ' Error in qpstart4a: bad jh value!'
      endif
c
c     Note:
c       y1 <- y1 (normal displacement)
c       y2 <- y2 (normal stress)       !!! Note definition different to the case with self-gravitation !!!
c       y3 <- y5 (potential)
c       y4 <- y6 (gravity)
c
      return
      end