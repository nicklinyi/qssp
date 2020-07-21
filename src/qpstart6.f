      subroutine qpstart6(ldeg,ly,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly
      complex*16 cy(6,3)
c
      integer*4 i
      complex*16 cldeg,cgl,ca2,cb2,cxi,c2mu,cgrdr
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
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
      cllm1=cldeg*(cldeg-c1)
      cllp1=cldeg*(cldeg+c1)
      c2lp1=c2*cldeg+c1
      c2lp3=c2*cldeg+c3
      c22lp3=c2*(c2*cldeg+c3)
      cgrdr=cgrup(ly)/crrup(ly)
c
      ca2=(comi2+cga(ly)+cgrdr)/cvp(ly)**2
      cb2=comi2/cvs(ly)**2
      cgl=c4*cllp1*(cgrdr/(cvp(ly)*cvs(ly)))**2
c
      kpg=cdsqrt((ca2+cb2-cdsqrt((cb2-ca2)**2+cgl))/c2)
      cfp=((kpg*cvs(ly))**2-comi2)/cgrdr
      chp=cfp-cldeg-c1
c
      ksg=cdsqrt((ca2+cb2+cdsqrt((cb2-ca2)**2+cgl))/c2)
      cfs=((ksg*cvs(ly))**2-comi2)/cgrdr
      chs=cfs-cldeg-c1
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
        cy(2,1)=-cxi*cfp*ph0+cmu(ly)*(-cllm1*chp*ps0
     &         +c2*(c2*cfp+cllp1)*ph1)/c2lp3
        cy(3,1)=-(chp*ps0/c2-ph1)/c2lp3
        cy(4,1)=cmu(ly)*(ph0-((cldeg-c1)*chp*ps0
     &         +c2*(cfp+c1)*ph1)/c2lp3)
        cy(5,1)=(cfp*cvp(ly)**2-(cldeg+c1)*cvs(ly)**2)/crrup(ly)**2
     &         -cga(ly)*cfp*ps0/c22lp3
        cy(6,1)=c2lp1*cy(5,1)+cga(ly)*cldeg*chp*ps0/c22lp3
      else
        call spbjh(ldeg,kpg,rrup(ly),rrlw(ly),
     &             zjupg(0),zjlwg(0),wjg(0),
     &             zhupg(0),zhlwg(0),whg(0))
c
        spa=c1-zjupg(ldeg)
        spb=-ci*zjupg(ldeg)
        cx2=(kpg*crrup(ly))**2
c
        cy(1,1)=cldeg*chp*spa-cfp*spb
        cy(2,1)=(-cxi*cfp*cx2+c2mu*cllm1*chp)*spa
     &         +c2mu*(c2*cfp+cllp1)*spb
        cy(3,1)=chp*spa+spb
        cy(4,1)=cmu(ly)*(cx2+c2*(cldeg-c1)*chp)*spa
     &         -c2mu*(cfp+c1)*spb
        cy(5,1)=cga(ly)*cfp*spa
        cy(6,1)=cga(ly)*(c2lp1*cfp-cldeg*chp)*spa
      endif
c
      if(cdabs(ksg*crrup(ly))**2.lt.0.5d0*dble(ldeg))then
c
        cxup=ksg*crrup(ly)
        ph0=spbphj(ldeg,cxup)
        ph1=spbphj(ldeg+1,cxup)
        ps0=spbpsj(ldeg,cxup)
        ps1=spbpsj(ldeg+1,cxup)
c
        cy(1,2)=-(cldeg*chs*ps0/c2+cfs*ph1)/c2lp3
        cy(2,2)=-cxi*cfs*ph0+cmu(ly)*(-cllm1*chs*ps0
     &         +c2*(c2*cfs+cllp1)*ph1)/c2lp3
        cy(3,2)=-(chs*ps0/c2-ph1)/c2lp3
        cy(4,2)=cmu(ly)*(ph0-((cldeg-c1)*chs*ps0
     &         +c2*(cfs+c1)*ph1)/c2lp3)
        cy(5,2)=(cfs*cvp(ly)**2-(cldeg+c1)*cvs(ly)**2)/crrup(ly)**2
     &         -cga(ly)*cfs*ps0/c22lp3
        cy(6,2)=c2lp1*cy(5,2)+cga(ly)*cldeg*chs*ps0/c22lp3
      else
c
        call spbjh(ldeg,ksg,rrup(ly),rrlw(ly),
     &             zjupg(0),zjlwg(0),wjg(0),
     &             zhupg(0),zhlwg(0),whg(0))
c
        spa=c1-zjupg(ldeg)
        spb=-ci*zjupg(ldeg)
        cx2=(ksg*crrup(ly))**2
c
        cy(1,2)=cldeg*chs*spa-cfs*spb
        cy(2,2)=(-cxi*cfs*cx2+c2mu*cllm1*chs)*spa
     &       +c2mu*(c2*cfs+cllp1)*spb
        cy(3,2)=chs*spa+spb
        cy(4,2)=cmu(ly)*(cx2+c2*(cldeg-c1)*chs)*spa
     &       -c2mu*(cfs+c1)*spb
        cy(5,2)=cga(ly)*cfs*spa
        cy(6,2)=cga(ly)*(c2lp1*cfs-cldeg*chs)*spa
      endif
c
      cy(1,3)=cldeg
      cy(2,3)=c2mu*cllm1
      cy(3,3)=c1
      cy(4,3)=c2mu*(cldeg-c1)
      cy(5,3)=cldeg*cgrdr-comi2
      cy(6,3)=c2lp1*(cldeg*cgrdr-comi2)-cga(ly)*cldeg
      return
      end