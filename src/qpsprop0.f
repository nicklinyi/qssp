      subroutine qpsprop0(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of spheroidal response for degree 0
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,j,istp,ly,lyswap,key
      complex*16 y0(3),c(3),yup(3),ylw(3)
      complex*16 wave(2),coef(3,3),b(3,2)
c
      if(lylwa.gt.lylw)return
c
c===============================================================================
c
c     propagation from surface to source
c
      if(freesurf.and.lyup.eq.1)then
        yup(1)=(1.d0,0.d0)
        yup(2)=(0.d0,0.d0)
        yup(3)=(0.d0,0.d0)
      else
        call qpstart0a(lyup,2,yup)
      endif
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,3)
c
      do ly=lyup,lys-1
        call caxcb(mas3x3inv(1,1,ly),yup,3,3,1,c)
        wave(1)=cdexp( cps(1,ly))
        wave(2)=cdexp(-cps(2,ly))
        if(ly.ge.lyr)then
          do i=1,3
            y0(i)=y0(i)*wave(2)
          enddo
        endif
c
	  c(1)=c(1)*wave(2)*wave(1)
c       c(2)=c(2)
        c(3)=c(3)*wave(2)
c
        call caxcb(mas3x3lw(1,1,ly),c,3,3,1,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,3)
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c     yup(3)=yup(3)
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lylwb)then
        lyswap=lylwb
        call qpstart0a(lyswap,1,ylw)
      else
        do i=1,3
          ylw(i)=mas3x3up(i,1,lylw)
        enddo
      endif
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylw-1,lys,-1
        call caxcb(mas3x3inv(1,1,ly),ylw,3,3,1,c)
        wave(1)=cdexp(-cps(1,ly))
        wave(2)=cdexp( cps(2,ly))
        if(ly.lt.lyr)then
          do i=1,3
            y0(i)=y0(i)*wave(1)
          enddo
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
        c(3)=c(3)*wave(1)
c
        call caxcb(mas3x3up(1,1,ly),c,3,3,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      y0(1)=y0(1)/crrup(lyr)
      y0(2)=y0(2)/crrup(lyr)**2
c     y0(3)=y0(3)
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(3,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      b(3,2)=(0.d0,0.d0)
      do j=1,3
        do i=1,3
          coef(i,j)=(0.d0,0.d0)
        enddo
      enddo
      do i=1,3
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
        coef(i,3)=(0.d0,0.d0)
      enddo
c
c     a constant will be added to potential for region below the source
c
      coef(3,3)=-(1.d0,0.d0)
c
      key=0
      call cdsvd500(coef,b,3,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsprop0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
          ypsv(3,istp)=ypsv(3,istp)+b(3,istp)
        enddo
      endif
c
      if(lylwa.le.0)then
        do istp=1,2
          ypsv(5,istp)=ypsv(3,istp)
          ypsv(6,istp)=ypsv(3,istp)/crrup(lyr)
          ypsv(3,istp)=(0.d0,0.d0)
        enddo
        return
      endif
c
c
c===============================================================================
c
c     propagation from bottom to source
c
      lyswap=lylwa
      call qpstart0a(lyswap,1,ylw)
c
      if(lylwa.eq.lyr.and.lylwa.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylwa-1,lys,-1
        call caxcb(mas3x3inv(1,1,ly),ylw,3,3,1,c)
        wave(1)=cdexp(-cps(1,ly))
        wave(2)=cdexp( cps(2,ly))
        if(ly.lt.lyr)then
          do i=1,3
            y0(i)=y0(i)*wave(1)
          enddo
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
        c(3)=c(3)*wave(1)
c
        call caxcb(mas3x3up(1,1,ly),c,3,3,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      if(lyr.gt.lys)then
        y0(1)=y0(1)/crrup(lyr)
        y0(2)=y0(2)/crrup(lyr)**2
c       y0(3)=y0(3)
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(3,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      b(3,2)=(0.d0,0.d0)
      do i=1,3
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
        coef(i,3)=(0.d0,0.d0)
      enddo
      coef(3,3)=-(1.d0,0.d0)
      key=0
      call cdsvd500(coef,b,3,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsprop0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=ypsv(i,istp)-b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=ypsv(i,istp)-b(2,istp)*y0(i)
          enddo
          ypsv(3,istp)=ypsv(3,istp)-b(3,istp)
        enddo
      endif
c
      do istp=1,2
        ypsv(5,istp)=ypsv(3,istp)
        ypsv(6,istp)=ypsv(3,istp)/crrup(lyr)
        ypsv(3,istp)=(0.d0,0.d0)
      enddo
c
      return
      end
