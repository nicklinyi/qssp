      subroutine qpdifmats(ly,ldeg,rr,mat)
      use qpalloc
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      complex*16 mat(6,6)
c
c     6x6 coefficient matrix for spheroidal mode l > 0 in solid media
c
      real*8 dro,rorr,ro1,mass
      complex*16 cldeg,clp1,cll1,cup,clw
      complex*16 crr,drr,crorr,clarr,cmurr,cxirr,cgarr,cgrrr
c
      complex*16 c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
      crr=dcmplx(rr,0.d0)
      cup=(crr-rrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      clarr=cup*claup(ly)+clw*clalw(ly)
      cmurr=cup*cmuup(ly)+clw*cmulw(ly)
      cxirr=clarr+c2*cmurr
c
      crorr=cup*croup(ly)+clw*crolw(ly)
      rorr=dreal(crorr)
      cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
      if(rr.le.rrlw(ly))then
        mass=0.d0
        crorr=crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
      else
        crorr=cup*croup(ly)+clw*crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
        dro=(rorr-rolw(ly))/(rr-rrlw(ly))
        ro1=rorr-dro*rrlw(ly) 
        mass=PI*(rr-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rr**2+rr*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rr**3+rr**2*rrlw(ly)
     &        +rr*rrlw(ly)**2+rrlw(ly)**3))
      endif
      cgrrr=cgrlw(ly)*(crrlw(ly)/crr)**2
     &     +dcmplx(BIGG*mass/rr**2,0.d0)
c
      mat(1,1)=(c1-c2*clarr/cxirr)/crr
      mat(1,2)=c1/cxirr/crr
      mat(1,3)=cll1*clarr/cxirr/crr
      mat(1,4)=(0.d0,0.d0)
      mat(1,5)=(0.d0,0.d0)
      mat(1,6)=(0.d0,0.d0)
c
      mat(2,1)=c4*(cmurr*(c1+c2*clarr/cxirr)/crr-crorr*cgrrr)
     &        -crorr*comi2*crr
      mat(2,2)=c2*clarr/cxirr/crr
      mat(2,3)=cll1*(crorr*cgrrr
     &        -c2*cmurr*(c1+c2*clarr/cxirr)/crr)
      mat(2,4)=cll1/crr
      mat(2,5)=crorr*clp1*crr
      mat(2,6)=-crorr*crr
c
      mat(3,1)=-c1/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=c2/crr
      mat(3,4)=c1/cmurr/crr
      mat(3,5)=(0.d0,0.d0)
      mat(3,6)=(0.d0,0.d0)
c
      mat(4,1)=crorr*cgrrr-c2*cmurr*(c1+c2*clarr/cxirr)/crr
      mat(4,2)=-clarr/cxirr/crr
      mat(4,3)=c2*cmurr*(c2*cll1*(c1-cmurr/cxirr)-c1)/crr
     &        -crorr*comi2*crr
      mat(4,4)=-c1/crr
      mat(4,5)=-crorr*crr
      mat(4,6)=(0.d0,0.d0)
c
      mat(5,1)=cgarr/crr
      mat(5,2)=(0.d0,0.d0)
      mat(5,3)=(0.d0,0.d0)
      mat(5,4)=(0.d0,0.d0)
      mat(5,5)=-clp1/crr
      mat(5,6)=c1/crr
c
      mat(6,1)=cgarr*clp1/crr
      mat(6,2)=(0.d0,0.d0)
      mat(6,3)=-cgarr*cll1/crr
      mat(6,4)=(0.d0,0.d0)
      mat(6,5)=(0.d0,0.d0)
      mat(6,6)=cldeg/crr
c
      return
      end