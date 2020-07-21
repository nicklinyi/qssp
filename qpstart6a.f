      subroutine qpstart6a(ldeg,ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly,jh
      complex*16 cy(6,3)
c
      integer*4 i,j
      real*8 peps,seps
      complex*16 cldeg,cxi,c2mu,cxp2,cxs2
      complex*16 cllp1,cllm1,c2lp1,c2lm1,c2lp3
      complex*16 ca,cb,ya,yb
c
      real*8 eps
      complex*16 c0,c1,c2,c3
      data eps/1.0d-08/
      data c0,c1,c2,c3/(0.d0,0.d0),(1.d0,0.d0),
     &                 (2.d0,0.d0),(3.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
      cllm1=cldeg*(cldeg-c1)
      cllp1=cldeg*(cldeg+c1)
      c2lp1=c2*cldeg+c1
      c2lm1=c2*cldeg-c1
      c2lp3=c2*cldeg+c3
      cxp2=(kp(ly)*crrup(ly))**2
      cxs2=(ks(ly)*crrup(ly))**2
c
      if(jh.eq.1.and.ly.ge.lys)then
        cy(1,1)=cldeg-zh12p(ldeg,1,ly)
        cy(2,1)=-cxi*cxp2+c2mu*(cllm1+c2*zh12p(ldeg,1,ly))
        cy(3,1)=c1
        cy(4,1)=c2mu*(cldeg-c1-zh12p(ldeg,1,ly))
        cy(5,1)=cga(ly)
        cy(6,1)=cga(ly)*(cldeg+c1)
c
        peps=dabs(dimag(zh12p(ldeg,1,ly)))/cdabs(zh12p(ldeg,1,ly))
        if(peps.lt.eps)then
          ca=dcmplx((peps/eps)**2,0.d0)
          cb=(1.d0,0.d0)-ca
          ya=cdsqrt(cy(1,1)*dconjg(cy(1,1))+cy(3,1)*dconjg(cy(3,1)))
          yb=cdsqrt(mas6x6up(1,1,ly)*dconjg(mas6x6up(1,1,ly))
     &             +mas6x6up(3,1,ly)*dconjg(mas6x6up(3,1,ly)))
          do i=1,6
            cy(i,1)=ca*cy(i,1)*yb+cb*mas6x6up(i,1,ly)*ya
          enddo
        endif
c
        cy(1,2)=-cllp1
        cy(2,2)=c2mu*cllp1*(c1-cldeg+zh12sv(ldeg,1,ly))
        cy(3,2)=-(cldeg+c1)+zh12sv(ldeg,1,ly)
        cy(4,2)=cmu(ly)*(cxs2-c2*(cldeg**2-c1)-c2*zh12sv(ldeg,1,ly))
        cy(5,2)=c0
        cy(6,2)=cga(ly)*cllp1
c
        seps=dabs(dimag(zh12sv(ldeg,1,ly)))/cdabs(zh12sv(ldeg,1,ly))
        if(seps.lt.eps)then
          ca=dcmplx((seps/eps)**2,0.d0)
          cb=(1.d0,0.d0)-ca
          ya=cdsqrt(cy(1,2)*dconjg(cy(1,2))+cy(3,2)*dconjg(cy(3,2)))
          yb=cdsqrt(mas6x6up(1,3,ly)*dconjg(mas6x6up(1,3,ly))
     &             +mas6x6up(3,3,ly)*dconjg(mas6x6up(3,3,ly)))
          do i=1,6
            cy(i,2)=ca*cy(i,2)*yb+cb*mas6x6up(i,3,ly)*ya
          enddo
        endif
        cy(1,3)=c0
        cy(2,3)=c0
        cy(3,3)=c0
        cy(4,3)=c0
        cy(5,3)=c1
        cy(6,3)=c2lp1
      else if(jh.eq.2.and.ly.le.lys)then
        cy(1,1)=cldeg-zh12p(ldeg,2,ly)
        cy(2,1)=-cxi*cxp2+c2mu*(cllm1+c2*zh12p(ldeg,2,ly))
        cy(3,1)=c1
        cy(4,1)=c2mu*(cldeg-c1-zh12p(ldeg,2,ly))
        cy(5,1)=cga(ly)
        cy(6,1)=cga(ly)*(cldeg+c1)
c
        cy(1,2)=-cllp1
        cy(2,2)=c2mu*cllp1*(c1-cldeg+zh12sv(ldeg,2,ly))
        cy(3,2)=-(cldeg+c1)+zh12sv(ldeg,2,ly)
        cy(4,2)=cmu(ly)*(cxs2-c2*(cldeg**2-c1)-c2*zh12sv(ldeg,2,ly))
        cy(5,2)=c0
        cy(6,2)=cga(ly)*cllp1
c
        cy(1,3)=c0
        cy(2,3)=c0
        cy(3,3)=c0
        cy(4,3)=c0
        cy(5,3)=c1
        cy(6,3)=c0
      else
        stop ' Error in qpstart6a: bad jh value!'
      endif
c
      return
      end