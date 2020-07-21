      subroutine qptprop(ldeg,ysh,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of toroidal response
c     ysh(2,2): solution vector (complex)
c
      integer*4 ldeg,lyup,lylw
      complex*16 ysh(2,2)
c
c     work space
c
      integer*4 i,istp,ly,lyswap,key
      complex*16 cldeg
      complex*16 y0(2),yup(2),ylw(2),wave(2),c(2)
      complex*16 coef(2,2),b(2,2)
c
      if(lylwa.gt.lylw)return
c
      cldeg=dcmplx(dble(ldeg),0.d0)
c
c===============================================================================
c
c     propagation from surface to source
c
      if(freesurf.and.lyup.eq.1.or.
     &   lyup.eq.lyob.and.lyob.gt.1.or.lyup.eq.lycc)then
        yup(1)=(1.d0,0.d0)
        yup(2)=(0.d0,0.d0)
      else
        call qpstart2t(ldeg,lyup,2,yup)
      endif
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,2)
c
      do ly=lyup,lys-1
        call caxcb(mat2x2inv(1,1,ly),yup,2,2,1,c)
        wave(1)=cdexp( cpt(1,ly))
        wave(2)=cdexp(-cpt(2,ly))
        if(ly.ge.lyr)then
          do i=1,2
            y0(i)=y0(i)*wave(2)
          enddo
        endif
c
	  c(1)=c(1)*wave(1)*wave(2)
c       c(2)=c(2)
c
        call caxcb(mat2x2lw(1,1,ly),c,2,2,1,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,2)
      enddo
c     yup(1)=yup(1)
      yup(2)=yup(2)/crrup(lys)
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lycm)then
        ylw(1)=(1.d0,0.d0)
        ylw(2)=(0.d0,0.d0)
      else if(lylw.eq.lylwb)then
        call qpstart2t(ldeg,lylw,1,ylw)
      else
        ylw(1)=mat2x2up(1,1,lylw)
        ylw(2)=mat2x2up(2,1,lylw)
      endif
      if(lylw.eq.lyr.and.lylw.gt.lys)then
        call cmemcpy(ylw,y0,2)
      endif
c
      do ly=lylw-1,lys,-1
        call caxcb(mat2x2inv(1,1,ly),ylw,2,2,1,c)
        wave(1)=cdexp(-cpt(1,ly))
        wave(2)=cdexp( cpt(2,ly))
c
        if(ly.lt.lyr)then
          do i=1,2
            y0(i)=y0(i)*wave(1)
          enddo
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
c
        call caxcb(mat2x2up(1,1,ly),c,2,2,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,2)
      enddo
c     ylw(1)=ylw(1)
      ylw(2)=ylw(2)/crrup(lys)
c
c     y0(1)=y0(1)
      y0(2)=y0(2)/crrup(lyr)
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qptprop: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(2,istp)*y0(i)
          enddo
        enddo
      endif
c
      if(lylwa.le.0)return
c
c===============================================================================
c
c     propagation from bottom to source
c
      lyswap=lylwa
      call qpstart2t(ldeg,lyswap,1,ylw)
      if(lylwa.eq.lyr.and.lylwa.gt.lys)then
        call cmemcpy(ylw,y0,2)
      endif
c
      do ly=lylwa-1,lys,-1
        call caxcb(mat2x2inv(1,1,ly),ylw,2,2,1,c)
        wave(1)=cdexp(-cpt(1,ly))
        wave(2)=cdexp( cpt(2,ly))
c
        if(ly.lt.lyr)then
          do i=1,2
            y0(i)=y0(i)*wave(1)
          enddo
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
c
        call caxcb(mat2x2up(1,1,ly),c,2,2,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,2)
      enddo
c     ylw(1)=ylw(1)
      ylw(2)=ylw(2)/crrup(lys)
c
      if(lyr.gt.lys)then
c       y0(1)=y0(1)
        y0(2)=y0(2)/crrup(lyr)
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qptprop: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ysh(i,istp)=ysh(i,istp)-b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,2
            ysh(i,istp)=ysh(i,istp)-b(2,istp)*y0(i)
          enddo
        enddo
      endif
      return
      end