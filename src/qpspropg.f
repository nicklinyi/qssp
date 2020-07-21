      subroutine qpspropg(ypsv,ldeg,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of speroidal response (with gravity)
c     ypsv(6,4): solution vector (complex)
c
      integer*4 ldeg,lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,j,j0,istp,ly,lyswap,ily,nly,key
      real*8 y4max,rr1,rr2,dlnr,h,f
      complex*16 cdet,alf,bet,cyabs,cyswap,ca,cb
      complex*16 y0(6,3),c(2)
      complex*16 yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      complex*16 coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
      external qpdifmatl,qpdifmats
c
      if(lylwa.gt.lylw)return
c
      f=dreal(comi)/PI2
c
c===============================================================================
c
c     propagation from surface to atmosphere/ocean bottom
c
      if(lyob.gt.lyup)then
        do j=1,2
          do i=1,4
            yupc(i,j)=(0.d0,0.d0)
          enddo
        enddo
c
        if(freesurf.and.lyup.eq.1)then
          yupc(1,1)=comi2*crrup(lyup)/cgrup(lyup)
          yupc(2,1)=(1.d0,0.d0)
c
          yupc(1,2)=(1.d0,0.d0)
          yupc(3,2)=cgrup(lyup)/crrup(lyup)
c
          if(ldeg.eq.1)then
            yupc(4,2)=(3.d0,0.d0)*yupc(3,2)
          endif
        else
          call qpstart4g(ldeg,lyup,2,yupc)
        endif
c
        if(lyr.eq.lyup)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=yupc(1,j)
            y0(2,j)=croup(lyup)*crrup(lyup)**2
     &             *(yupc(1,j)*cgrup(lyup)/crrup(lyup)
     &             -comi2*yupc(2,j)-yupc(3,j))
            y0(3,j)=yupc(2,j)
            y0(5,j)=yupc(3,j)
            y0(6,j)=yupc(4,j)
          enddo
        endif
      endif
c
      do ly=lyup,min0(lyob,lys)-1
        if(lyup.lt.lyos.and.lyob.gt.lyos.and.ly.eq.lyos)then
c
c         interface ocean-atmosphere
c
          ca=yupc(1,1)*cgrup(ly)/crrup(ly)-yupc(3,1)
          cb=yupc(1,2)*cgrup(ly)/crrup(ly)-yupc(3,2)
c
          if(cdabs(ca).gt.cdabs(cb))then
            do i=1,4
              yupc(i,2)=yupc(i,1)*cb-yupc(i,2)*ca
            enddo
            yupc(2,2)=yupc(2,2)*crolw(ly-1)/croup(ly)
c
            yupc(1,1)=yupc(1,1)*comi2
            yupc(2,1)=(yupc(2,1)*crolw(ly-1)*comi2
     &               +ca*(croup(ly)-crolw(ly-1)))/croup(ly)
            yupc(3,1)=yupc(3,1)*comi2
            yupc(4,1)=yupc(4,1)*comi2
            if(ly.ge.lyr)then
              do i=1,6
                y0(i,2)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,1)=y0(i,1)*comi2
              enddo
            endif
          else
            do i=1,4
              yupc(i,1)=yupc(i,1)*cb-yupc(i,2)*ca
            enddo
            yupc(2,1)=yupc(2,1)*crolw(ly-1)/croup(ly)
c
            yupc(1,2)=yupc(1,2)*comi2
            yupc(2,2)=(yupc(2,2)*crolw(ly-1)*comi2
     &               +cb*(croup(ly)-crolw(ly-1)))/croup(ly)
            yupc(3,2)=yupc(3,2)*comi2
            yupc(4,2)=yupc(4,2)*comi2
            if(ly.ge.lyr)then
              do i=1,6
                y0(i,1)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,2)=y0(i,2)*comi2
              enddo
            endif
          endif
        endif
c
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,1)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,1)=yupc(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+yupc(i,2)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            yupc(i,2)=yupc(i,2)-alf*yupc(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,2)*dconjg(yupc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,2)=yupc(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(yupc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr-1)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=yupc(1,j)
            y0(2,j)=crolw(ly)*crrlw(ly)**2
     &             *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &             -comi2*yupc(2,j)-yupc(3,j))
            y0(3,j)=yupc(2,j)
            y0(5,j)=yupc(3,j)
            y0(6,j)=yupc(4,j)
          enddo
        endif
      enddo
c
      if(lys.ge.lyob)then
        if(lyup.lt.lyob)then
          ly=lyob-1
          do j=1,2
            yup(1,j)=yupc(1,j)
            yup(2,j)=crolw(ly)*crrlw(ly)**2
     &              *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &              -comi2*yupc(2,j)-yupc(3,j))
            yup(3,j)=(0.d0,0.d0)
            yup(4,j)=(0.d0,0.d0)
            yup(5,j)=yupc(3,j)
            yup(6,j)=yupc(4,j)
          enddo
          yup(1,3)=(0.d0,0.d0)
          yup(2,3)=(0.d0,0.d0)
          yup(3,3)=(1.d0,0.d0)
          yup(4,3)=(0.d0,0.d0)
          yup(5,3)=(0.d0,0.d0)
          yup(6,3)=(0.d0,0.d0)
          if(lyr.eq.lyob)call cmemcpy(yup,y0,18)
        else if(lyup.lt.lycm)then
          if(freesurf.and.lyup.eq.1)then
            do j=1,3
              do i=1,6
                yup(i,j)=(0.d0,0.d0)
              enddo
            enddo
            yup(1,1)=(1.d0,0.d0)
            yup(3,2)=(1.d0,0.d0)
            yup(5,3)=(1.d0,0.d0)
            if(ldeg.eq.1)then
              yup(6,3)=(3.d0,0.d0)*yup(5,3)
            endif
          else
            call qpstart6g(ldeg,lyup,2,yup)
          endif
          if(lyr.eq.lyup)call cmemcpy(yup,y0,18)
        endif
      endif
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source
c
      do ly=max0(lyup,lyob),min0(lycm,lys)-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,1)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,1)=yup(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,2)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,2)=yup(i,2)-alf*yup(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,2)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,2)=yup(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,3)*dconjg(yup(i,1))/cypnorm(i,ly)**2
            bet=bet+yup(i,3)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,3)=yup(i,3)-alf*yup(i,1)-bet*yup(i,2)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,3)*dconjg(yup(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,3)=yup(i,3)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(yup,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,18)
      enddo
c===============================================================================
c
c     propagation from core-mantle boundary to source or solid core surface
c
      if(lys.ge.lycm)then
        if(lyup.lt.lycm)then
          y4max=0.d0
          do j=1,3
            if(y4max.lt.cdabs(yup(4,j)))then
              y4max=cdabs(yup(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=yup(i,j0)
            yup(i,j0)=yup(i,3)
            yup(i,3)=cyswap
          enddo
          do j=1,2
            yupc(1,j)=yup(4,3)*yup(1,j)-yup(4,j)*yup(1,3)
            yupc(2,j)=yup(4,3)*yup(2,j)-yup(4,j)*yup(2,3)
            yupc(3,j)=yup(4,3)*yup(5,j)-yup(4,j)*yup(5,3)
            yupc(4,j)=yup(4,3)*yup(6,j)-yup(4,j)*yup(6,3)
          enddo
          if(lycm.ge.lyr)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=yup(4,3)*y0(i,j)-yup(4,j)*y0(i,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
c
c         y2 = Ut
c
          do j=1,2
            c(j)=-yupc(2,j)/croup(lycm)/crrup(lycm)**2
     &          +yupc(1,j)*cgrup(lycm)/crrup(lycm)-yupc(3,j)
          enddo
          do i=1,4
            yupc(i,1)=c(2)*yupc(i,1)-c(1)*yupc(i,2)
          enddo
          yupc(2,1)=(0.d0,0.d0)
          do i=1,4
            yupc(i,2)=comi2*yupc(i,2)
          enddo
          yupc(2,2)=c(2)
          if(lycm.ge.lyr)then
            do i=1,6
              y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
              y0(i,2)=comi2*y0(i,2)
            enddo
          endif
        else if(lyup.lt.lycc)then
          call qpstart4g(ldeg,lyup,2,yupc)
          if(lyr.eq.lyup)then
            do j=1,2
              y0(1,j)=yupc(1,j)
              y0(2,j)=croup(lyup)*crrup(lyup)**2
     &                *(yupc(1,j)*cgrup(lyup)/crrup(lyup)
     &                -comi2*yupc(2,j)-yupc(3,j))
              y0(3,j)=yupc(2,j)
              y0(5,j)=yupc(3,j)
              y0(6,j)=yupc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        endif
      endif
c
      do ly=lycm,min0(lycc,lys)-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,1)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,1)=yupc(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+yupc(i,2)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            yupc(i,2)=yupc(i,2)-alf*yupc(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,2)*dconjg(yupc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,2)=yupc(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(yupc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr-1)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=yupc(1,j)
            y0(2,j)=crolw(ly)*crrlw(ly)**2
     &             *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &             -comi2*yupc(2,j)-yupc(3,j))
            y0(3,j)=yupc(2,j)
            y0(5,j)=yupc(3,j)
            y0(6,j)=yupc(4,j)
          enddo
        endif
      enddo
c
c===============================================================================
c
c     propagation from solid core surface to source
c
      if(lys.ge.lycc)then
        if(lyup.lt.lycc)then
          ly=lycc-1
          do j=1,2
            yup(1,j)=yupc(1,j)
            yup(2,j)=crolw(ly)*crrlw(ly)**2
     &            *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &            -comi2*yupc(2,j)-yupc(3,j))
            yup(3,j)=(0.d0,0.d0)
            yup(4,j)=(0.d0,0.d0)
            yup(5,j)=yupc(3,j)
            yup(6,j)=yupc(4,j)
          enddo
          yup(1,3)=(0.d0,0.d0)
          yup(2,3)=(0.d0,0.d0)
          yup(3,3)=(1.d0,0.d0)
          yup(4,3)=(0.d0,0.d0)
          yup(5,3)=(0.d0,0.d0)
          yup(6,3)=(0.d0,0.d0)
          if(lyr.eq.lycc)call cmemcpy(yup,y0,18)
        else if(lyup.lt.ly0)then
          call qpstart6g(ldeg,lyup,2,yup)
          if(lyr.eq.lyup)call cmemcpy(yup,y0,18)
        endif
      endif
c
      do ly=lycc,lys-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,1)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,1)=yup(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,2)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,2)=yup(i,2)-alf*yup(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,2)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,2)=yup(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,3)*dconjg(yup(i,1))/cypnorm(i,ly)**2
            bet=bet+yup(i,3)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,3)=yup(i,3)-alf*yup(i,1)-bet*yup(i,2)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,3)*dconjg(yup(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,3)=yup(i,3)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(yup,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,18)
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        do j=1,2
          yup(1,j)=yupc(1,j)
          yup(2,j)=croup(lys)*crrup(lys)**2
     &          *(yupc(1,j)*cgrup(lys)/crrup(lys)
     &          -comi2*yupc(2,j)-yupc(3,j))
          yup(3,j)=(0.d0,0.d0)
          yup(4,j)=(0.d0,0.d0)
          yup(5,j)=yupc(3,j)
          yup(6,j)=yupc(4,j)
        enddo
        do i=1,6
          yup(i,3)=(0.d0,0.d0)
        enddo
      endif
      do j=1,3
        yup(1,j)=yup(1,j)/crrup(lys)
        yup(2,j)=yup(2,j)/crrup(lys)**2
        yup(3,j)=yup(3,j)/crrup(lys)
        yup(4,j)=yup(4,j)/crrup(lys)**2
c        yup(5,j)=yup(5,j)
        yup(6,j)=yup(6,j)/crrup(lys)
      enddo
c
c===============================================================================
c
c     propagation within inner core
c
      if(lylw.ge.lycc)then
c
c       lowest layer is within inner core
c
        if(lylw.eq.lylwb)then
          call qpstart6g(ldeg,lylw,1,ylw)
        else
          call qpstart6g(ldeg,lylw,0,ylw)
        endif
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          call cmemcpy(ylw,y0,18)
        endif
      endif
c
      do ly=lylw-1,max0(lycc,lys),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lys.lt.lycc)then
        if(lylw.ge.lycc)then
c
c         interface conditions: solid to liquid
c
          y4max=0.d0
          do j=1,3
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(4,3)*ylw(1,j)-ylw(4,j)*ylw(1,3)
            ylwc(2,j)=ylw(4,3)*ylw(2,j)-ylw(4,j)*ylw(2,3)
            ylwc(3,j)=ylw(4,3)*ylw(5,j)-ylw(4,j)*ylw(5,3)
            ylwc(4,j)=ylw(4,3)*ylw(6,j)-ylw(4,j)*ylw(6,3)
          enddo
          if(lycc.le.lyr.and.lycc.gt.lys)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=ylw(4,3)*y0(i,j)-ylw(4,j)*y0(i,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
c
c         y2 = Ut
c
          do j=1,2
            c(j)=-ylwc(2,j)/crolw(lycc-1)/crrlw(lycc-1)**2
     &          +ylwc(1,j)*cgrlw(lycc-1)/crrlw(lycc-1)-ylwc(3,j)
          enddo
          do i=1,4
            ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
          enddo
          ylwc(2,1)=(0.d0,0.d0)
          do i=1,4
            ylwc(i,2)=comi2*ylwc(i,2)
          enddo
          ylwc(2,2)=c(2)
          if(lycc.le.lyr.and.lycc.gt.lys)then
            do i=1,6
              y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
              y0(i,2)=comi2*y0(i,2)
            enddo
          endif
        else if(lylw.ge.lycm)then
c
c         lowest layer is within the liquid core
c
          if(lylw.eq.lylwb)then
            call qpstart4g(ldeg,lylw,1,ylwc)
          else
            call qpstart4g(ldeg,lylw,0,ylwc)
          endif
c
          if(lylw.eq.lyr.and.lylw.gt.lys)then
            do j=1,3
              do i=1,6
                y0(i,j)=(0.d0,0.d0)
              enddo
            enddo
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=croup(lylw)*crrup(lylw)**2
     &               *(ylwc(1,j)*cgrup(lylw)/crrup(lylw)
     &               -comi2*ylwc(2,j)-ylwc(3,j))
              y0(3,j)=ylwc(2,j)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
          endif
        endif
      endif
c
      do ly=min0(lylw-1,lycc-1),max0(lycm,lys),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,1)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+ylwc(i,2)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,2)*dconjg(ylwc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(ylwc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=croup(ly)*crrup(ly)**2
     &             *(ylwc(1,j)*cgrup(ly)/crrup(ly)
     &             -comi2*ylwc(2,j)-ylwc(3,j))
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lys.lt.lycm)then
        if(lylw.ge.lycm)then
c
c         interface conditions: liquid to solid
c
          do j=1,2
            ylw(1,j)=ylwc(1,j)
            ylw(2,j)=croup(lycm)*crrup(lycm)**2
     &              *(ylwc(1,j)*cgrup(lycm)/crrup(lycm)
     &              -comi2*ylwc(2,j)-ylwc(3,j))
            ylw(3,j)=(0.d0,0.d0)
            ylw(4,j)=(0.d0,0.d0)
            ylw(5,j)=ylwc(3,j)
            ylw(6,j)=ylwc(4,j)
          enddo
          ylw(1,3)=(0.d0,0.d0)
          ylw(2,3)=(0.d0,0.d0)
          ylw(3,3)=(1.d0,0.d0)
          ylw(4,3)=(0.d0,0.d0)
          ylw(5,3)=(0.d0,0.d0)
          ylw(6,3)=(0.d0,0.d0)
c
          if(lycm.eq.lyr.and.lycm.gt.lys)then
            call cmemcpy(ylw,y0,18)
          endif
        else if(lylw.ge.lyob)then
c
c         lowest layer is within the mantle
c
          if(lylw.eq.lylwb)then
            call qpstart6g(ldeg,lylw,1,ylw)
          else
            call qpstart6g(ldeg,lylw,0,ylw)
          endif
          if(lylw.eq.lyr.and.lylw.gt.lys)then
            call cmemcpy(ylw,y0,18)
          endif
        endif
      endif
c
      do ly=min0(lylw-1,lycm-1),max0(lys,lyob),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lys.lt.lyob)then
        if(lylw.ge.lyob)then
c
c         interface conditions: solid to liquid
c
          y4max=0.d0
          do j=1,3
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(4,3)*ylw(1,j)-ylw(4,j)*ylw(1,3)
            ylwc(2,j)=ylw(4,3)*ylw(2,j)-ylw(4,j)*ylw(2,3)
            ylwc(3,j)=ylw(4,3)*ylw(5,j)-ylw(4,j)*ylw(5,3)
            ylwc(4,j)=ylw(4,3)*ylw(6,j)-ylw(4,j)*ylw(6,3)
          enddo
          if(lyob.le.lyr.and.lyob.gt.lys)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=ylw(4,3)*y0(i,j)-ylw(4,j)*y0(i,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
c
c         y2 = Ut
c
          do j=1,2
            c(j)=-ylwc(2,j)/crolw(lyob-1)/crrlw(lyob-1)**2
     &          +ylwc(1,j)*cgrlw(lyob-1)/crrlw(lyob-1)-ylwc(3,j)
          enddo
          do i=1,4
            ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
          enddo
          ylwc(2,1)=(0.d0,0.d0)
          do i=1,4
            ylwc(i,2)=comi2*ylwc(i,2)
          enddo
          ylwc(2,2)=c(2)
          if(lyob.le.lyr.and.lyob.gt.lys)then
            do i=1,6
              y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
              y0(i,2)=comi2*y0(i,2)
            enddo
          endif
        else if(lylw.lt.lyob)then
c
c         lowest layer is within the atmosphere
c
          if(lylw.eq.lylwb)then
            call qpstart4g(ldeg,lylw,1,ylwc)
          else
            call qpstart4g(ldeg,lylw,0,ylwc)
          endif
c
          if(lylw.eq.lyr.and.lylw.gt.lys)then
            do j=1,3
              do i=1,6
                y0(i,j)=(0.d0,0.d0)
              enddo
            enddo
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=croup(lylw)*crrup(lylw)**2
     &               *(ylwc(1,j)*cgrup(lylw)/crrup(lylw)
     &               -comi2*ylwc(2,j)-ylwc(3,j))
              y0(3,j)=ylwc(2,j)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
          endif
        endif
      endif
c
      do ly=min0(lylw-1,lyob-1),lys,-1
        if(lyob.gt.lyos.and.ly.eq.lyos-1)then
c
c         interface ocean-atmosphere
c
          ca=ylwc(1,1)*cgrlw(ly)/crrlw(ly)-ylwc(3,1)
          cb=ylwc(1,2)*cgrlw(ly)/crrlw(ly)-ylwc(3,2)
c
          if(cdabs(ca).gt.cdabs(cb))then
            do i=1,4
              ylwc(i,2)=ylwc(i,1)*cb-ylwc(i,2)*ca
            enddo
            ylwc(2,2)=ylwc(2,2)*croup(ly+1)/crolw(ly)
c
            ylwc(1,1)=ylwc(1,1)*comi2
            ylwc(2,1)=(ylwc(2,1)*croup(ly+1)*comi2
     &               +ca*(crolw(ly)-croup(ly+1)))/crolw(ly)
            ylwc(3,1)=ylwc(3,1)*comi2
            ylwc(4,1)=ylwc(4,1)*comi2
            if(ly.lt.lyr)then
              do i=1,6
                y0(i,2)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,1)=y0(i,1)*comi2
              enddo
            endif
          else
            do i=1,4
              ylwc(i,1)=ylwc(i,1)*cb-ylwc(i,2)*ca
            enddo
            ylwc(2,1)=ylwc(2,1)*croup(ly+1)/crolw(ly)
c
            ylwc(1,2)=ylwc(1,2)*comi2
            ylwc(2,2)=(ylwc(2,2)*croup(ly+1)*comi2
     &               +cb*(crolw(ly)-croup(ly+1)))/crolw(ly)
            ylwc(3,2)=ylwc(3,2)*comi2
            ylwc(4,2)=ylwc(4,2)*comi2
            if(ly.lt.lyr)then
              do i=1,6
                y0(i,1)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,2)=y0(i,2)*comi2
              enddo
            endif
          endif
        endif
c
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,1)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+ylwc(i,2)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,2)*dconjg(ylwc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(ylwc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=croup(ly)*crrup(ly)**2
     &             *(ylwc(1,j)*cgrup(ly)/crrup(ly)
     &             -comi2*ylwc(2,j)-ylwc(3,j))
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=croup(lys)*crrup(lys)**2
     &            *(ylwc(1,j)*cgrup(lys)/crrup(lys)
     &            -comi2*ylwc(2,j)-ylwc(3,j))
          ylw(3,j)=(0.d0,0.d0)
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        do i=1,6
          ylw(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
c
      do j=1,3
        y0(1,j)=y0(1,j)/crrup(lyr)
        y0(2,j)=y0(2,j)/crrup(lyr)**2
        y0(3,j)=y0(3,j)/crrup(lyr)
        y0(4,j)=y0(4,j)/crrup(lyr)**2
c        y0(5,j)=y0(5,j)
        y0(6,j)=y0(6,j)/crrup(lyr)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      if(vsup(lys).le.0.d0)then
        do istp=1,2
          do i=1,4
            b4(i,istp)=(0.d0,0.d0)
          enddo
          b4(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,2
          do i=1,2
            coef4(i,j)=yup(i,j)
            coef4(i,j+2)=-ylw(i,j)
          enddo
          do i=3,4
            coef4(i,j)=yup(i+2,j)
            coef4(i,j+2)=-ylw(i+2,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef4,b4,4,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpspropg: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.le.lys)then
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
        do istp=3,4
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            b6(i,istp)=(0.d0,0.d0)
          enddo
          b6(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup(i,j)
            coef6(i,j+3)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpspropg: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.le.lys)then
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j+3,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
      if(lylwa.le.0)return
c
c===============================================================================
c
c     propagation within inner core
c
      if(lylwa.ge.lycc)then
c
c       lowest layer is within inner core
c
        lyswap=lylwa
        call qpstart6g(ldeg,lyswap,1,ylw)
c
        if(lylwa.eq.lyr.and.lylwa.gt.lys)then
          call cmemcpy(ylw,y0,18)
        endif
      endif
c
      do ly=lylwa-1,max0(lycc,lys),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lys.lt.lycc)then
        if(lylwa.ge.lycc)then
c
c         interface conditions: solid to liquid
c
          y4max=0.d0
          do j=1,3
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(4,3)*ylw(1,j)-ylw(4,j)*ylw(1,3)
            ylwc(2,j)=ylw(4,3)*ylw(2,j)-ylw(4,j)*ylw(2,3)
            ylwc(3,j)=ylw(4,3)*ylw(5,j)-ylw(4,j)*ylw(5,3)
            ylwc(4,j)=ylw(4,3)*ylw(6,j)-ylw(4,j)*ylw(6,3)
          enddo
          if(lycc.le.lyr.and.lycc.gt.lys)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=ylw(4,3)*y0(i,j)-ylw(4,j)*y0(i,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
c
c         y2 = Ut
c
          do j=1,2
            c(j)=-ylwc(2,j)/crolw(lycc-1)/crrlw(lycc-1)**2
     &          +ylwc(1,j)*cgrlw(lycc-1)/crrlw(lycc-1)-ylwc(3,j)
          enddo
          do i=1,4
            ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
          enddo
          ylwc(2,1)=(0.d0,0.d0)
          do i=1,4
            ylwc(i,2)=comi2*ylwc(i,2)
          enddo
          ylwc(2,2)=c(2)
          if(lycc.le.lyr.and.lycc.gt.lys)then
            do i=1,6
              y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
              y0(i,2)=comi2*y0(i,2)
            enddo
          endif
        else if(lylwa.ge.lycm)then
c
c         lowest layer is within the liquid core
c
          lyswap=lylwa
          call qpstart4g(ldeg,lyswap,1,ylwc)
c
          if(lylwa.eq.lyr.and.lylwa.gt.lys)then
            do j=1,3
              do i=1,6
                y0(i,j)=(0.d0,0.d0)
              enddo
            enddo
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=croup(lylwa)*crrup(lylwa)**2
     &               *(ylwc(1,j)*cgrup(lylwa)/crrup(lylwa)
     &               -comi2*ylwc(2,j)-ylwc(3,j))
              y0(3,j)=ylwc(2,j)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
          endif
        endif
      endif
c
      do ly=min0(lylwa-1,lycc-1),max0(lycm,lys),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,1)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+ylwc(i,2)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,2)*dconjg(ylwc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(ylwc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=croup(ly)*crrup(ly)**2
     &             *(ylwc(1,j)*cgrup(ly)/crrup(ly)
     &             -comi2*ylwc(2,j)-ylwc(3,j))
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lys.lt.lycm)then
        if(lylwa.ge.lycm)then
c
c         interface conditions: liquid to solid
c
          do j=1,2
            ylw(1,j)=ylwc(1,j)
            ylw(2,j)=croup(lycm)*crrup(lycm)**2
     &              *(ylwc(1,j)*cgrup(lycm)/crrup(lycm)
     &              -comi2*ylwc(2,j)-ylwc(3,j))
            ylw(3,j)=(0.d0,0.d0)
            ylw(4,j)=(0.d0,0.d0)
            ylw(5,j)=ylwc(3,j)
            ylw(6,j)=ylwc(4,j)
          enddo
          ylw(1,3)=(0.d0,0.d0)
          ylw(2,3)=(0.d0,0.d0)
          ylw(3,3)=(1.d0,0.d0)
          ylw(4,3)=(0.d0,0.d0)
          ylw(5,3)=(0.d0,0.d0)
          ylw(6,3)=(0.d0,0.d0)
c
          if(lycm.eq.lyr.and.lycm.gt.lys)then
            call cmemcpy(ylw,y0,18)
          endif
        else if(lylwa.ge.lyob)then
c
c         lowest layer is within the mantle
c
          lyswap=lylwa
          call qpstart6g(ldeg,lyswap,1,ylw)
c
          if(lylwa.eq.lyr.and.lylwa.gt.lys)then
            call cmemcpy(ylw,y0,18)
          endif
        endif
      endif
c
      do ly=min0(lylwa-1,lycm-1),max0(lys,lyob),-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vsup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lys.lt.lyob)then
        if(lylwa.ge.lyob)then
c
c         interface conditions: solid to liquid
c
          y4max=0.d0
          do j=1,3
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(4,3)*ylw(1,j)-ylw(4,j)*ylw(1,3)
            ylwc(2,j)=ylw(4,3)*ylw(2,j)-ylw(4,j)*ylw(2,3)
            ylwc(3,j)=ylw(4,3)*ylw(5,j)-ylw(4,j)*ylw(5,3)
            ylwc(4,j)=ylw(4,3)*ylw(6,j)-ylw(4,j)*ylw(6,3)
          enddo
          if(lyob.le.lyr.and.lyob.gt.lys)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=ylw(4,3)*y0(i,j)-ylw(4,j)*y0(i,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
c
c         y2 = Ut
c
          do j=1,2
            c(j)=-ylwc(2,j)/crolw(lyob-1)/crrlw(lyob-1)**2
     &          +ylwc(1,j)*cgrlw(lyob-1)/crrlw(lyob-1)-ylwc(3,j)
          enddo
          do i=1,4
            ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
          enddo
          ylwc(2,1)=(0.d0,0.d0)
          do i=1,4
            ylwc(i,2)=comi2*ylwc(i,2)
          enddo
          ylwc(2,2)=c(2)
          if(lyob.le.lyr.and.lyob.gt.lys)then
            do i=1,6
              y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
              y0(i,2)=comi2*y0(i,2)
            enddo
          endif
        else if(lylwa.lt.lyob)then
c
c         lowest layer is within the atmosphere
c
          lyswap=lylwa
          call qpstart4g(ldeg,lyswap,1,ylwc)
c
          if(lylwa.eq.lyr.and.lylwa.gt.lys)then
            do j=1,3
              do i=1,6
                y0(i,j)=(0.d0,0.d0)
              enddo
            enddo
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=croup(lylwa)*crrup(lylwa)**2
     &               *(ylwc(1,j)*cgrup(lylwa)/crrup(lylwa)
     &               -comi2*ylwc(2,j)-ylwc(3,j))
              y0(3,j)=ylwc(2,j)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
          endif
        endif
      endif
c
      do ly=min0(lylwa-1,lyob-1),lys,-1
        if(lyob.gt.lyos.and.ly.eq.lyos-1)then
c
c         interface ocean-atmosphere
c
          ca=ylwc(1,1)*cgrlw(ly)/crrlw(ly)-ylwc(3,1)
          cb=ylwc(1,2)*cgrlw(ly)/crrlw(ly)-ylwc(3,2)
c
          if(cdabs(ca).gt.cdabs(cb))then
            do i=1,4
              ylwc(i,2)=ylwc(i,1)*cb-ylwc(i,2)*ca
            enddo
            ylwc(2,2)=ylwc(2,2)*croup(ly+1)/crolw(ly)
c
            ylwc(1,1)=ylwc(1,1)*comi2
            ylwc(2,1)=(ylwc(2,1)*croup(ly+1)*comi2
     &               +ca*(crolw(ly)-croup(ly+1)))/crolw(ly)
            ylwc(3,1)=ylwc(3,1)*comi2
            ylwc(4,1)=ylwc(4,1)*comi2
            if(ly.lt.lyr)then
              do i=1,6
                y0(i,2)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,1)=y0(i,1)*comi2
              enddo
            endif
          else
            do i=1,4
              ylwc(i,1)=ylwc(i,1)*cb-ylwc(i,2)*ca
            enddo
            ylwc(2,1)=ylwc(2,1)*croup(ly+1)/crolw(ly)
c
            ylwc(1,2)=ylwc(1,2)*comi2
            ylwc(2,2)=(ylwc(2,2)*croup(ly+1)*comi2
     &               +cb*(crolw(ly)-croup(ly+1)))/crolw(ly)
            ylwc(3,2)=ylwc(3,2)*comi2
            ylwc(4,2)=ylwc(4,2)*comi2
            if(ly.lt.lyr)then
              do i=1,6
                y0(i,1)=y0(i,1)*cb-y0(i,2)*ca
                y0(i,2)=y0(i,2)*comi2
              enddo
            endif
          endif
        endif
c
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*dmax1(f/vpup(ly),dble(ldeg)/rrlw(ly)))/5
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,1)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+ylwc(i,2)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,2)*dconjg(ylwc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(ylwc,4,2,ly,ldeg,qpdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=croup(ly)*crrup(ly)**2
     &             *(ylwc(1,j)*cgrup(ly)/crrup(ly)
     &             -comi2*ylwc(2,j)-ylwc(3,j))
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=croup(lys)*crrup(lys)**2
     &            *(ylwc(1,j)*cgrup(lys)/crrup(lys)
     &            -comi2*ylwc(2,j)-ylwc(3,j))
          ylw(3,j)=(0.d0,0.d0)
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        do i=1,6
          ylw(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
c
      if(lyr.gt.lys)then
        do j=1,3
          y0(1,j)=y0(1,j)/crrup(lyr)
          y0(2,j)=y0(2,j)/crrup(lyr)**2
          y0(3,j)=y0(3,j)/crrup(lyr)
          y0(4,j)=y0(4,j)/crrup(lyr)**2
c          y0(5,j)=y0(5,j)
          y0(6,j)=y0(6,j)/crrup(lyr)
        enddo
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      if(vsup(lys).le.0.d0)then
        do istp=1,2
          do i=1,4
            b4(i,istp)=(0.d0,0.d0)
          enddo
          b4(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,2
          do i=1,2
            coef4(i,j)=yup(i,j)
            coef4(i,j+2)=-ylw(i,j)
          enddo
          do i=3,4
            coef4(i,j)=yup(i+2,j)
            coef4(i,j+2)=-ylw(i+2,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef4,b4,4,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpspropg: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.le.lys)then
          do istp=1,2
            do i=1,6
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)-b4(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,2
            do i=1,6
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)-b4(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
        do istp=3,4
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            b6(i,istp)=(0.d0,0.d0)
          enddo
          b6(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup(i,j)
            coef6(i,j+3)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpspropg: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.le.lys)then
          do istp=1,4
            do i=1,6
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)-b6(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,6
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)-b6(j+3,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
      return
      end