      subroutine qppsvkerng(f,ldeg,ypsv)
      use qpalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ypsv(6,4): psv solution vector (complex) with gravity effect
c
      integer*4 ldeg
      real*8 f
      complex*16 ypsv(6,4)
c
      integer*4 i,ia,na,ly,istp,lyup,lylw,lwup
      integer*4 lya(4)
      real*8 xmin
      complex*16 cldeg,cruplw,cxp,cxs
c
      complex*16 c1,c2
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      lyup=min0(lyupp(ldeg),lyups(ldeg))
      lylw=min0(lylwb,max0(lylwp(ldeg),lylws(ldeg)))
c
      if(cdabs(comi).le.0.d0)then
        if(vsup(lys).le.0.d0)return
        if(lys.ge.lyob.and.lys.lt.lycm)then
          lyup=max0(lyup,lyob)
          lylw=min0(lylw,lycm)
        else if(lys.ge.lycc)then
          lyup=max0(lyup,lycc)
          lylw=min0(lylw,ly0)
        endif
      endif
c
      if(lyr.lt.lyup.or.lyr.gt.lylw.or.lylw.lt.lylwa)return
c
      cldeg=dcmplx(dble(ldeg),0.d0)
c
      xmin=dsqrt(2.d0*dble(2*ldeg+1))
c
      if(.not.setzh12a)then
        do ly=1,ly0
          cxp=kp(ly)*crrup(ly)
          call spbh12(ldegpsv(ly),ldegmax,cxp,zh12p(0,1,ly))
          if(vsup(ly).gt.0.d0)then
            cxs=ks(ly)*crrup(ly)
            call spbh12(ldegpsv(ly),ldegmax,cxs,zh12sv(0,1,ly))
          endif
        enddo
        setzh12a=.true.
      endif
c
      na=0
      if(.not.freesurf.or.lyup.gt.1)then
        na=na+1
        lya(na)=lyup
      endif
      if(lylwa.gt.0)then
        na=na+1
        lya(na)=lylwa
      endif
      if(lylwb.le.ly0)then
        na=na+1
        lya(na)=lylwb
      endif
      if(lylw.ne.lylwa.and.lylw.ne.lylwb)then
        na=na+1
        lya(na)=lylw
      endif
c
      lwup=-1
c
      do ia=1,na
        ly=lya(ia)
        if(rrlw(ly).gt.0.d0)then
          cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
        else
          cruplw=(0.d0,0.d0)
        endif
        if(xp(ly).le.xmin)then
          ksmallp(ly)=.true.
          wj(ldeg,ly,1)=cldeg*cruplw
          wh(ldeg,ly,1)=-(cldeg+c1)*cruplw
        else if(ksmallp(ly))then
          ksmallp(ly)=.false.
          call spbjh(ldegpsv(ly),kp(ly),rrup(ly),rrlw(ly),
     &               zjup(0,ly,1),zjlw(0,ly,1),wj(0,ly,1),
     &               zhup(0,ly,1),zhlw(0,ly,1),wh(0,ly,1))
        endif
        if(vsup(ly).gt.0.d0)then
          if(xs(ly).le.xmin)then
            ksmalls(ly)=.true.
            wj(ldeg,ly,2)=(cldeg+c2)*cruplw
            wh(ldeg,ly,2)=-(cldeg-c1)*cruplw
          else if(ksmalls(ly))then
            ksmalls(ly)=.false.
            call spbjh(ldegpsv(ly),ks(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,2),zjlw(0,ly,2),wj(0,ly,2),
     &                 zhup(0,ly,2),zhlw(0,ly,2),wh(0,ly,2))
          endif
        endif
        if(ldeg.eq.0)then
          call qpsmat0(ly,lylw,lwup)
        else if(vsup(ly).gt.0.d0)then
          call qpsmat(ldeg,ly,lylw,lwup)
        else
          call qpsmatc(ldeg,ly,lylw,lwup)
        endif
      enddo
c
      if(ldeg.eq.0)then
        call qpspropg0(ypsv,lyup,lylw)
      else
        call qpspropg(ypsv,ldeg,lyup,lylw)
      endif
c
      return
      end