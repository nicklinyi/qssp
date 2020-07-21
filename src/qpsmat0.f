      subroutine qpsmat0(ly,lylw,lwup)
      use qpalloc
      implicit none
c
c     calculate 3x3 spheroidal layer matrix for a solid shell in case of degree l = 0
c
      integer*4 ly,lylw,lwup
c
      integer*4 j
      complex*16 c2mu,cxi,cx,cx2,cdet
      complex*16 ph0jup,ph1jup,ph0yup,ph1yup,ps0jup,ps0yup
      complex*16 ph0jlw,ph1jlw,ph0ylw,ph1ylw,ps0jlw,ps0ylw
      complex*16 spb(2,2)
      complex*16 spbphj,spbphy,spbpsj
c
      complex*16 ci,c1,c2,c3,c6
      data ci,c1,c2,c3,c6/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),
     &                    (3.d0,0.d0),(6.d0,0.d0)/
c
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
c
      cx=kp(ly)*crrup(ly)
      cx2=cx*cx
      if(ksmallp(ly))then
        ph0jup=spbphj(0,cx)
        ph1jup=spbphj(1,cx)
c
        ph0yup=spbphy(0,cx)
        ph1yup=spbphy(1,cx)
c
        mas3x3up(1,1,ly)=cx2*ph1jup/c3
        mas3x3up(2,1,ly)=cx2*(cxi*ph0jup-c2*c2mu*ph1jup/c3)
        mas3x3up(3,1,ly)=-cga(ly)*ph0jup
        mas3x3up(1,2,ly)=ph1yup
        mas3x3up(2,2,ly)=cxi*cx2*ph0yup-c2*c2mu*ph1yup
        mas3x3up(3,2,ly)=-cga(ly)*ph0yup
      else
        spb(1,1)=c1-zjup(0,ly,1)
        spb(2,1)=-ci*zjup(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhup(0,ly,1)
c
        do j=1,2
          mas3x3up(1,j,ly)=spb(2,j)
          mas3x3up(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
          mas3x3up(3,j,ly)=-cga(ly)*spb(1,j)
        enddo
      endif
      mas3x3up(1,3,ly)=(0.d0,0.d0)
      mas3x3up(2,3,ly)=(0.d0,0.d0)
      mas3x3up(3,3,ly)=c1
c
      if(ly.eq.lylw.or.lwup.lt.0)return
c
      cx=kp(ly)*crrlw(ly)
      cx2=cx*cx
      if(ksmallp(ly))then
        ph0jlw=spbphj(0,cx)
        ph1jlw=spbphj(1,cx)
c
        ph0ylw=spbphy(0,cx)
        ph1ylw=spbphy(1,cx)
c
        mas3x3lw(1,1,ly)=cx2*ph1jlw/c3
        mas3x3lw(2,1,ly)=cx2*(cxi*ph0jlw-c2*c2mu*ph1jlw/c3)
        mas3x3lw(3,1,ly)=-cga(ly)*ph0jlw
        mas3x3lw(1,2,ly)=ph1ylw
        mas3x3lw(2,2,ly)=cxi*cx2*ph0ylw-c2*c2mu*ph1ylw
        mas3x3lw(3,2,ly)=-cga(ly)*ph0ylw
      else
        spb(1,1)=c1-zjlw(0,ly,1)
        spb(2,1)=-ci*zjlw(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhlw(0,ly,1)
c
        do j=1,2
          mas3x3lw(1,j,ly)=spb(2,j)
          mas3x3lw(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
          mas3x3lw(3,j,ly)=-cga(ly)*spb(1,j)
        enddo
      endif
      mas3x3lw(1,3,ly)=(0.d0,0.d0)
      mas3x3lw(2,3,ly)=(0.d0,0.d0)
      mas3x3lw(3,3,ly)=c1
c
      if(lwup.eq.0)then
c
c       calculate inverse masrix at upper radius
c
        cdet=mas3x3up(1,1,ly)*mas3x3up(2,2,ly)
     &      -mas3x3up(1,2,ly)*mas3x3up(2,1,ly)
        mas3x3inv(1,1,ly)=mas3x3up(2,2,ly)/cdet
        mas3x3inv(1,2,ly)=-mas3x3up(1,2,ly)/cdet
        mas3x3inv(1,3,ly)=(0.d0,0.d0)
        mas3x3inv(2,1,ly)=-mas3x3up(2,1,ly)/cdet
        mas3x3inv(2,2,ly)=mas3x3up(1,1,ly)/cdet
        mas3x3inv(2,3,ly)=(0.d0,0.d0)
        mas3x3inv(3,1,ly)=(mas3x3up(2,1,ly)*mas3x3up(3,2,ly)
     &                   -mas3x3up(3,1,ly)*mas3x3up(2,2,ly))/cdet
        mas3x3inv(3,2,ly)=(mas3x3up(1,2,ly)*mas3x3up(3,1,ly)
     &                   -mas3x3up(1,1,ly)*mas3x3up(3,2,ly))/cdet
        mas3x3inv(3,3,ly)=c1
      else
c
c       calculate inverse masrix at lower radius
c
        cdet=mas3x3lw(1,1,ly)*mas3x3lw(2,2,ly)
     &      -mas3x3lw(1,2,ly)*mas3x3lw(2,1,ly)
        mas3x3inv(1,1,ly)=mas3x3lw(2,2,ly)/cdet
        mas3x3inv(1,2,ly)=-mas3x3lw(1,2,ly)/cdet
        mas3x3inv(1,3,ly)=(0.d0,0.d0)
        mas3x3inv(2,1,ly)=-mas3x3lw(2,1,ly)/cdet
        mas3x3inv(2,2,ly)=mas3x3lw(1,1,ly)/cdet
        mas3x3inv(2,3,ly)=(0.d0,0.d0)
        mas3x3inv(3,1,ly)=(mas3x3lw(2,1,ly)*mas3x3lw(3,2,ly)
     &                   -mas3x3lw(3,1,ly)*mas3x3lw(2,2,ly))/cdet
        mas3x3inv(3,2,ly)=(mas3x3lw(1,2,ly)*mas3x3lw(3,1,ly)
     &                   -mas3x3lw(1,1,ly)*mas3x3lw(3,2,ly))/cdet
        mas3x3inv(3,3,ly)=c1
      endif
c
      return
      end