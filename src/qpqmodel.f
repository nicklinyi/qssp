      subroutine qpqmodel(f)
      use qpalloc
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      real*8 f
c
      integer*4 ly
      real*8 fac,pup,plw,sup,slw
c
      if(.not.dispersion)then
        fac=0.d0
      else if(f.le.FSBLW)then
        fac=dlog(FSBLW/FSBREF)/PI
      else if(f.ge.FSBUP)then
        fac=dlog(FSBUP/FSBREF)/PI
      else
        fac=dlog(f/FSBREF)/PI
      endif
c
      if(1.d0+fac/qsmin.le.0.d0)then
        stop ' Error in qpqmodel: too small Qs value!'
      endif
      do ly=1,ly0
        if(ly.ge.lyob.and.ly.lt.lycm.or.ly.ge.lycc)then
          sup=1.d0+fac/qsup(ly)
          slw=1.d0+fac/qslw(ly)
          cvsup(ly)=dcmplx(vsup(ly),0.d0)
     &             *dcmplx(sup,0.5d0/qsup(ly))
          cvslw(ly)=dcmplx(vslw(ly),0.d0)
     &             *dcmplx(slw,0.5d0/qslw(ly))
          cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        else
          cvsup(ly)=(0.d0,0.d0)
          cvslw(ly)=(0.d0,0.d0)
          cvs(ly)=(0.d0,0.d0)
        endif
        pup=1.d0+fac/qpup(ly)
        plw=1.d0+fac/qplw(ly)
        cvpup(ly)=dcmplx(vpup(ly),0.d0)
     &           *dcmplx(pup,0.5d0/qpup(ly))
        cvplw(ly)=dcmplx(vplw(ly),0.d0)
     &           *dcmplx(plw,0.5d0/qplw(ly))
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
        cmuup(ly)=croup(ly)*cvsup(ly)**2
        claup(ly)=croup(ly)*cvpup(ly)**2-(2.d0,0.d0)*cmuup(ly)
        cmulw(ly)=crolw(ly)*cvslw(ly)**2
        clalw(ly)=crolw(ly)*cvplw(ly)**2-(2.d0,0.d0)*cmulw(ly)
        cmu(ly)=cro(ly)*cvs(ly)**2
        cla(ly)=cro(ly)*cvp(ly)**2-(2.d0,0.d0)*cmu(ly)
      enddo
      return
      end
