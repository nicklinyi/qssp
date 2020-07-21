      subroutine qpshkern(f,ldeg,ysh)
      use qpalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ysh(2,2): sh solution vector (complex)
c
      integer*4 ldeg
      real*8 f
      complex*16 ysh(2,2)
c
      integer*4 i,j,ly,istp,lyup,lylw,lwup
      real*8 xmin
      complex*16 cldeg,cll1,cruplw,cxt
c
      complex*16 c1,c2
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      lyup=lyupt(ldeg)
      lylw=min0(lylwb,lylwt(ldeg))
c
      if(ldeg.le.1.or.lyr.lt.lyup.or.lyr.gt.lylw.or.lylw.lt.lylwa)return
c
      if(.not.setzh12t)then
        if(lys.ge.lycc)then
          do ly=lycc,ly0
            cxt=kt(ly)*crrup(ly)
            call spbh12(ldegpsv(ly),ldegmax,cxt,zh12sh(0,1,ly))
          enddo
        else
          do ly=lyob,min0(lycm-1,ly0)
            cxt=kt(ly)*crrup(ly)
            call spbh12(ldegpsv(ly),ldegmax,cxt,zh12sh(0,1,ly))
          enddo
        endif
        setzh12t=.true.
      endif
c
      xmin=dsqrt(2.d0*dble(2*ldeg+1))
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      cll1=cldeg*(cldeg+c1)
c
      do ly=lyup,lylw
        if(vsup(ly).gt.0.d0)then
          if(rrlw(ly).gt.0.d0)then
            cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
          else
            cruplw=(0.d0,0.d0)
          endif
          if(xt(ly).le.xmin)then
            ksmallt(ly)=.true.
            wj(ldeg,ly,3)=cldeg*cruplw
            wh(ldeg,ly,3)=-(cldeg+c1)*cruplw
          else if(ksmallt(ly))then
            ksmallt(ly)=.false.
            call spbjh(ldegsh(ly),kt(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,3),zjlw(0,ly,3),wj(0,ly,3),
     &                 zhup(0,ly,3),zhlw(0,ly,3),wh(0,ly,3))
          endif
        endif
      enddo
c
      do ly=lyup,lylw
        if(vsup(ly).gt.0.d0)then
          if(ly.lt.lys)then
            lwup=0
            cpt(1,ly)=-wj(ldeg,ly,3)
            cpt(2,ly)=-wh(ldeg,ly,3)
          else
            lwup=1
            cpt(1,ly)=wj(ldeg,ly,3)
            cpt(2,ly)=wh(ldeg,ly,3)
          endif
          call qptmat(ldeg,ly,lylw,lwup)
        endif
      enddo
c
      call qptprop(ldeg,ysh,lyup,lylw)
c
      return
      end

