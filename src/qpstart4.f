      subroutine qpstart4(ldeg,ly,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly
      complex*16 cy(4,2)
c
      integer*4 j
      complex*16 cldeg,cgl,ca2,cb2,cgrdr
      complex*16 kpg,ksg,cfp,cfs,chp,chs,spa,spb,cxup,cx2
      complex*16 cllm1,cllp1,c2lp1,c2lp3,c22lp3
      complex*16 ph0,ph1,ps0,ps1
      complex*16 spbphj,spbpsj
c
      complex*16 ci,c1,c2,c3,c4
      data ci,c1,c2,c3,c4/(0.d0,1.d0),
     &     (1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      cgrdr=cgrup(ly)/crrup(ly)
      c2lp1=c2*cldeg+c1
c
      if(dreal(comi).le.0.d0)then
c
c       imcompressible fluid in the quasi-static case
c
        cy(1,1)=cldeg
        cy(2,1)=c1
        cy(3,1)=(0.d0,0.d0)
        cy(4,1)=(0.d0,0.d0)
      else
        cllm1=cldeg*(cldeg-c1)
        cllp1=cldeg*(cldeg+c1)
        c2lp3=c2*cldeg+c3
        c22lp3=c2*(c2*cldeg+c3)
c
        ca2=(comi2+cga(ly)+cgrdr)/cvp(ly)**2
        cb2=comi2/cvs(ly)**2
        cgl=c4*cllp1*(cgrdr/(cvp(ly)*cvs(ly)))**2
c
        kpg=cdsqrt(comi2+cga(ly)+cgrdr-cllp1*cgrdr**2/comi2)/cvp(ly)
        cfp=-comi2/cgrdr
        chp=cfp-cldeg-c1
c
        if(cdabs(kpg*crrup(ly))**2.lt.0.5d0*dble(ldeg))then
c
          cxup=kpg*crrup(ly)
          ph0=spbphj(ldeg,cxup)
          ph1=spbphj(ldeg+1,cxup)
          ps0=spbpsj(ldeg,cxup)
          ps1=spbpsj(ldeg+1,cxup)
c
          cy(1,1)=-(cldeg*chp*ps0/c2+cfp*ph1)/c2lp3
          cy(2,1)=-(chp*ps0/c2-ph1)/c2lp3
          cy(3,1)=(cfp*cvp(ly)**2-(cldeg+c1)*cvs(ly)**2)/crrup(ly)**2
     &             -cga(ly)*cfp*ps0/c22lp3
          cy(4,1)=c2lp1*cy(3,1)+cga(ly)*cldeg*chp*ps0/c22lp3
        else
          call spbjh(ldeg,kpg,rrup(ly),rrlw(ly),
     &               zjupg(0),zjlwg(0),wjg(0),
     &               zhupg(0),zhlwg(0),whg(0))
c
          spa=c1-zjupg(ldeg)
          spb=-ci*zjupg(ldeg)
          cx2=(kpg*crrup(ly))**2
c
          cy(1,1)=cldeg*chp*spa-cfp*spb
          cy(2,1)=chp*spa+spb
          cy(3,1)=cga(ly)*cfp*spa
          cy(4,1)=cga(ly)*(c2lp1*cfp-cldeg*chp)*spa
        endif
      endif
c
      cy(1,2)=cldeg
      cy(2,2)=c1
      cy(3,2)=cldeg*cgrdr-comi2
      cy(4,2)=c2lp1*(cldeg*cgrdr-comi2)-cga(ly)*cldeg
c
c     Note:
c       y1 <- y1 (normal displacement)
c       y2 <- y3 (horizontal displacement)
c       y3 <- y5 (potential)
c       y4 <- y6 (gravity)
c
      return
      end