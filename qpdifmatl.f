      subroutine qpdifmatl(ly,ldeg,rr,mat)
      use qpalloc
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      complex*16 mat(6,6)
c
c     4x4 coefficient matrix for spheroidal mode l > 0 in liquid media
c
      real*8 mass,rorr,beta,dro,ro1
      complex*16 cldeg,clp1,cll1
      complex*16 cup,clw,crr,crorr,cgrrr,cvprr
      complex*16 cgarr,gamma,cn2rr
      logical*2 bv
c
      complex*16 c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
c
      crr=dcmplx(rr,0.d0)
      cup=(crr-crrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      cvprr=cup*cvpup(ly)+clw*cvplw(ly)
      beta=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
c
      if(rr.le.rrlw(ly))then
        mass=0.d0
        crorr=crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
      else if(ly.ge.lyos)then
        crorr=cup*croup(ly)+clw*crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        dro=(rorr-rolw(ly))/(rr-rrlw(ly))
        ro1=rorr-dro*rrlw(ly) 
        mass=PI*(rr-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rr**2+rr*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rr**3+rr**2*rrlw(ly)
     &        +rr*rrlw(ly)**2+rrlw(ly)**3))
      else
        rorr=rolw(ly)*dexp(dlog(roup(ly)/rolw(ly))
     &                  *(rr-rrlw(ly))/(rrup(ly)-rrlw(ly)))
        crorr=dcmplx(rorr,0.d0)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        mass=2.d0*PI2/beta**3
     &      *(rorr*(beta*rr*(beta*rr-2.d0)+2.d0)
     &       -rolw(ly)*(beta*rrlw(ly)*(beta*rrlw(ly)-2.d0)+2.d0))
      endif
      cgrrr=cgrlw(ly)*(crrlw(ly)/crr)**2
     &     +dcmplx(BIGG*mass/rr**2,0.d0)
c
c     stabilization at the static limit
c
      if(ly.lt.lyos)then
        bv=fbvatm.gt.0.d0
        if(cdabs(comi)/PI2.lt.fbvatm)then
          gamma=dcmplx((dsin(0.25d0*cdabs(comi)/fbvatm))**3,0.d0)
          cvprr=cvprr/cdsqrt(gamma
     &         +(gamma-c1)*dcmplx(beta,0.d0)*cvprr**2/cgrrr)
        endif
      else if(ly.ge.lyos.and.ly.lt.lyob)then
        bv=fbvocean.gt.0.d0
        if(cdabs(comi)/PI2.lt.fbvocean)then
          gamma=dcmplx((dsin(0.25d0*cdabs(comi)/fbvocean))**3,0.d0)
          cvprr=cvprr/cdsqrt(gamma
     &         +(gamma-c1)*dcmplx(beta,0.d0)*cvprr**2/cgrrr)
        endif
      else if(ly.ge.lycm.and.ly.lt.lycc)then
        bv=fbvcore.gt.0.d0
        if(cdabs(comi)/PI2.lt.fbvcore)then
          gamma=dcmplx((dsin(0.25d0*cdabs(comi)/fbvcore))**3,0.d0)
          cvprr=cvprr/cdsqrt(gamma
     &         +(gamma-c1)*dcmplx(beta,0.d0)*cvprr**2/cgrrr)
        endif
      else
        stop 'Error in qpdifmatl: unknown fiquid layer!'
      endif
c
c     Brunt-Väisälä frequency
c
      cn2rr=-cgrrr*(dcmplx(beta,0.d0)+cgrrr/cvprr**2)
c
      mat(1,1)=cgrrr/cvprr**2-c1/crr
      mat(1,2)=cll1/crr-comi2*crr/cvprr**2
      mat(1,3)=-crr/cvprr**2
      mat(1,4)=(0.d0,0.d0)
c
      if(bv.and.cdabs(comi).gt.0.d0)then
        mat(2,1)=(c1-cn2rr/comi2)/crr
        mat(2,2)=cn2rr/cgrrr
        mat(2,3)=cn2rr/(comi2*cgrrr)
        mat(2,4)=(0.d0,0.d0)
      else
        mat(2,1)=c1/crr
        mat(2,2)=(0.d0,0.d0)
        mat(2,3)=(0.d0,0.d0)
        mat(2,4)=(0.d0,0.d0)
      endif
c
      mat(3,1)=cgarr/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=-clp1/crr
      mat(3,4)=c1/crr
c
      mat(4,1)=cgarr*clp1/crr
      mat(4,2)=-cgarr*cll1/crr
      mat(4,3)=(0.d0,0.d0)
      mat(4,4)=cldeg/crr
      return
      end