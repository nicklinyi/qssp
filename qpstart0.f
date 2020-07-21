      subroutine qpstart0(ly,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c     for n = 0
c
      integer*4 ly
      complex*16 cy(3)
c
      complex*16 cxi,c2mu,cgrdr,kpg,cxup,spa,spb
      complex*16 spbphj
c
      complex*16 ci,c1,c2,c3
      data ci,c1,c2,c3/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
      cgrdr=cgrup(ly)/crrup(ly)
c
      kpg=cdsqrt(comi2+cga(ly)+cgrdr)/cvp(ly)
      cxup=kpg*rrup(ly)
c
      if(cdabs(cxup).gt.0.1d0)then
        call spbjh(0,kpg,rrup(ly),rrlw(ly),
     &             zjupg(0),zjlwg(0),wjg(0),
     &             zhupg(0),zhlwg(0),whg(0))
c
        spa=c1-zjupg(0)
        spb=-ci*zjupg(0)
      else
        spa=spbphj(0,cxup)
        spb=spbphj(1,cxup)*cxup**2/c3
      endif
c
      cy(1)=spb
      cy(2)=cxi*cxup**2*spa+c2*c2mu*spb
      cy(3)=-cga(ly)*spa
      return
      end