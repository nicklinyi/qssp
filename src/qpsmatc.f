      subroutine qpsmatc(ldeg,ly,lylw,lwup)
      use qpalloc
      implicit none
c
c     calculate 4x4 spheroidal layer matrix for a solid shell
c
      integer*4 ldeg,ly,lylw,lwup
c
      integer*4 i,j,key
      complex*16 cldeg,cllm1,cllp1,c2lp1,c2lm1,c2lp3,cxp2,acc
      complex*16 ph0jup,ph1jup,ph0yup,ph1yup
      complex*16 ph0jlw,ph1jlw,ph0ylw,ph1ylw
      complex*16 spb(2,2),mas(4,4)
      complex*16 spbphj,spbphy
c
      complex*16 ci,c1,c2,c3
      data ci,c1,c2,c3/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      cllm1=cldeg*(cldeg-c1)
      cllp1=cldeg*(cldeg+c1)
      c2lp1=c2*cldeg+c1
      c2lm1=c2*cldeg-c1
      c2lp3=c2*cldeg+c3
      acc=cro(ly)*comi2
c
c     for upper radius
c
      cxp2=(kp(ly)*crrup(ly))**2
      if(ksmallp(ly))then
        ph0jup=spbphj(ldeg,kp(ly)*crrup(ly))
        ph1jup=spbphj(ldeg+1,kp(ly)*crrup(ly))
c
        ph0yup=spbphy(ldeg,kp(ly)*crrup(ly))
        ph1yup=spbphy(ldeg+1,kp(ly)*crrup(ly))
c
        mas4x4up(1,1,ly)=cldeg*ph0jup-cxp2*ph1jup/c2lp3
        mas4x4up(2,1,ly)=-acc*crrup(ly)**2*ph0jup
        mas4x4up(3,1,ly)=cga(ly)*ph0jup
        mas4x4up(4,1,ly)=cga(ly)*(cldeg+c1)*ph0jup
c
        mas4x4up(1,2,ly)=cldeg*ph0yup-ph1yup*c2lp1
        mas4x4up(2,2,ly)=-acc*crrup(ly)**2*ph0yup
        mas4x4up(3,2,ly)=cga(ly)*ph0yup
        mas4x4up(4,2,ly)=cga(ly)*(cldeg+c1)*ph0yup
      else
        spb(1,1)=c1-zjup(ldeg,ly,1)
        spb(2,1)=-ci*zjup(ldeg,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhup(ldeg,ly,1)
        do j=1,2
          mas4x4up(1,j,ly)=cldeg*spb(1,j)-spb(2,j)
          mas4x4up(2,j,ly)=-acc*crrup(ly)**2*spb(1,j)
          mas4x4up(3,j,ly)=cga(ly)*spb(1,j)
          mas4x4up(4,j,ly)=cga(ly)*(cldeg+c1)*spb(1,j)
        enddo
      endif
c
      do j=3,4
        do i=1,4
          mas4x4up(i,j,ly)=(0.d0,0.d0)
        enddo
      enddo
      mas4x4up(3,3,ly)=c1
      mas4x4up(4,3,ly)=c2lp1
      mas4x4up(3,4,ly)=c1
c
      if(ly.eq.lylw.or.lwup.lt.0)return
c
c     for lower radius
c
      cxp2=(kp(ly)*crrlw(ly))**2
      if(ksmallp(ly))then
        ph0jlw=spbphj(ldeg,kp(ly)*crrlw(ly))
        ph1jlw=spbphj(ldeg+1,kp(ly)*crrlw(ly))
c
        ph0ylw=spbphy(ldeg,kp(ly)*crrlw(ly))
        ph1ylw=spbphy(ldeg+1,kp(ly)*crrlw(ly))
c
        mas4x4lw(1,1,ly)=cldeg*ph0jlw-cxp2*ph1jlw/c2lp3
        mas4x4lw(2,1,ly)=-acc*crrlw(ly)**2*ph0jlw
        mas4x4lw(3,1,ly)=cga(ly)*ph0jlw
        mas4x4lw(4,1,ly)=cga(ly)*(cldeg+c1)*ph0jlw
c
        mas4x4lw(1,2,ly)=cldeg*ph0ylw-ph1ylw*c2lp1
        mas4x4lw(2,2,ly)=-acc*crrlw(ly)**2*ph0ylw
        mas4x4lw(3,2,ly)=cga(ly)*ph0ylw
        mas4x4lw(4,2,ly)=cga(ly)*(cldeg+c1)*ph0ylw
      else
        spb(1,1)=c1-zjlw(ldeg,ly,1)
        spb(2,1)=-ci*zjlw(ldeg,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhlw(ldeg,ly,1)
        do j=1,2
          mas4x4lw(1,j,ly)=cldeg*spb(1,j)-spb(2,j)
          mas4x4lw(2,j,ly)=-acc*crrlw(ly)**2*spb(1,j)
          mas4x4lw(3,j,ly)=cga(ly)*spb(1,j)
          mas4x4lw(4,j,ly)=cga(ly)*(cldeg+c1)*spb(1,j)
        enddo
      endif
c
      do j=3,4
        do i=1,4
          mas4x4lw(i,j,ly)=(0.d0,0.d0)
        enddo
      enddo
      mas4x4lw(3,3,ly)=c1
      mas4x4lw(4,3,ly)=c2lp1
      mas4x4lw(3,4,ly)=c1
c
      if(lwup.eq.0)then
c
c       calculate inverse matrix at upper radius
c
        do j=1,4
          do i=1,4
            mas(i,j)=mas4x4up(i,j,ly)
            mas4x4inv(i,j,ly)=(0.d0,0.d0)
          enddo
          mas4x4inv(j,j,ly)=c1
        enddo
      else
c
c       calculate inverse matrix at lower radius
c
        do j=1,4
          do i=1,4
            mas(i,j)=mas4x4lw(i,j,ly)
            mas4x4inv(i,j,ly)=(0.d0,0.d0)
          enddo
          mas4x4inv(j,j,ly)=c1
        enddo
      endif
      key=0
      call cdsvd500(mas,mas4x4inv(1,1,ly),4,4,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsmatc: anormal exit from cdsvd500!'
        return
      endif
c
      return
      end