      subroutine qpsmat(ldeg,ly,lylw,lwup)
      use qpalloc
      implicit none
c
c     calculate 6x6 spheroidal layer matrix for a solid shell
c
      integer*4 ldeg,ly,lylw,lwup
c
      integer*4 i,j,key
      complex*16 cldeg,cxi,c2mu,cllm1,cllp1,c2lp1,c2lm1,c2lp3
      complex*16 cs2dp2,cxp2,cxs2
      complex*16 dphjup,dphjlw,dphyup,dphylw
      complex*16 ph0jup(2),ph1jup(2),ph0yup(2),ph1yup(2)
      complex*16 ph0jlw(2),ph1jlw(2),ph0ylw(2),ph1ylw(2)
      complex*16 ps0jup(2),ps0jlw(2),ps0yup(2),ps0ylw(2)
      complex*16 spb(2,4),mas(6,6)
      complex*16 spbphj,spbphy,spbpsj,spbpsy,spbdphj,spbdphy
c
      complex*16 ci,c1,c2,c3,c4
      data ci,c1,c2,c3,c4/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),
     &                    (3.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
      cllm1=cldeg*(cldeg-c1)
      cllp1=cldeg*(cldeg+c1)
      c2lp1=c2*cldeg+c1
      c2lm1=c2*cldeg-c1
      c2lp3=c2*cldeg+c3
c
c     for upper radius
c
      cxp2=(kp(ly)*crrup(ly))**2
      cxs2=(ks(ly)*crrup(ly))**2
      if(ksmallp(ly).or.ksmalls(ly))then
        ph0jup(1)=spbphj(ldeg,kp(ly)*crrup(ly))
        ph1jup(1)=spbphj(ldeg+1,kp(ly)*crrup(ly))
        ps0jup(1)=spbpsj(ldeg,kp(ly)*crrup(ly))
c
        ph0yup(1)=spbphy(ldeg,kp(ly)*crrup(ly))
        ph1yup(1)=spbphy(ldeg-1,kp(ly)*crrup(ly))
        ps0yup(1)=spbpsy(ldeg,kp(ly)*crrup(ly))
      endif
      if(ksmallp(ly))then
        mas6x6up(1,1,ly)=cldeg*ph0jup(1)-cxp2*ph1jup(1)/c2lp3
        mas6x6up(2,1,ly)=(-cxi*cxp2+c2mu*cllm1)*ph0jup(1)
     &                   +c2*c2mu*cxp2*ph1jup(1)/c2lp3
        mas6x6up(3,1,ly)=ph0jup(1)
        mas6x6up(4,1,ly)=c2mu*((cldeg-c1)*ph0jup(1)
     &                   -cxp2*ph1jup(1)/c2lp3)
        mas6x6up(5,1,ly)=cga(ly)*ph0jup(1)
        mas6x6up(6,1,ly)=cga(ly)*(cldeg+c1)*ph0jup(1)
c
        mas6x6up(1,2,ly)=-(cldeg+c1)*ph0yup(1)
     &                   +cxp2*ph1yup(1)/c2lm1
        mas6x6up(2,2,ly)=(-cxi*cxp2+c2mu*(cldeg+c1)*(cldeg+c2))
     &                  *ph0yup(1)-c2*c2mu*cxp2*ph1yup(1)/c2lm1
        mas6x6up(3,2,ly)=ph0yup(1)
        mas6x6up(4,2,ly)=c2mu*(-(cldeg+c2)*ph0yup(1)
     &                         +cxp2*ph1yup(1)/c2lm1)
        mas6x6up(5,2,ly)=cga(ly)*ph0yup(1)
        mas6x6up(6,2,ly)=cga(ly)*(cldeg+c1)*ph0yup(1)
      else
        spb(1,1)=c1-zjup(ldeg,ly,1)
        spb(2,1)=-ci*zjup(ldeg,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhup(ldeg,ly,1)
        do j=1,2
          mas6x6up(1,j,ly)=cldeg*spb(1,j)-spb(2,j)
          mas6x6up(2,j,ly)=(-cxi*cxp2+c2mu*cllm1)*spb(1,j)
     &                  +c2*c2mu*spb(2,j)
          mas6x6up(3,j,ly)=spb(1,j)
          mas6x6up(4,j,ly)=c2mu*((cldeg-c1)*spb(1,j)-spb(2,j))
          mas6x6up(5,j,ly)=cga(ly)*spb(1,j)
          mas6x6up(6,j,ly)=cga(ly)*(cldeg+c1)*spb(1,j)
        enddo
      endif
      if(ksmalls(ly))then
        ph0jup(2)=spbphj(ldeg,ks(ly)*crrup(ly))
        ph1jup(2)=spbphj(ldeg+1,ks(ly)*crrup(ly))
        ps0jup(2)=spbpsj(ldeg,ks(ly)*crrup(ly))
c
        ph0yup(2)=spbphy(ldeg,ks(ly)*crrup(ly))
        ph1yup(2)=spbphy(ldeg-1,ks(ly)*crrup(ly))
        ps0yup(2)=spbpsy(ldeg,ks(ly)*crrup(ly))
c
        cs2dp2=(cvs(ly)/cvp(ly))**2
        dphjup=spbdphj(ldeg,kp(ly)*crrup(ly),ks(ly)*crrup(ly))
     &            *(c1-(cvs(ly)/cvp(ly))**2)
        dphyup=spbdphy(ldeg,kp(ly)*crrup(ly),ks(ly)*crrup(ly))
     &        *(c1-(cvs(ly)/cvp(ly))**2)
c
        mas6x6up(1,3,ly)=cldeg*dphjup-c2*cs2dp2*ph1jup(1)
        mas6x6up(2,3,ly)=c2mu*(cllm1*dphjup-c2lp3*ph0jup(1)
     &                    +c4*cs2dp2*ph1jup(1)+c2*cldeg*ph1jup(2))
        mas6x6up(3,3,ly)=dphjup+c2*ph1jup(2)/(cldeg+c1)
        mas6x6up(4,3,ly)=c2mu*((cldeg-c1)*dphjup
     &                        +c2lp3*ph0jup(2)/(cldeg+c1)
     &                -c2*(cs2dp2*ph1jup(1)+ph1jup(2)/(cldeg+c1)))
        mas6x6up(5,3,ly)=-cga(ly)*cs2dp2*ps0jup(1)
        mas6x6up(6,3,ly)=cga(ly)*((cldeg+c1)*dphjup
     &                            -ps0jup(2)*c2lp1)
c
        mas6x6up(1,4,ly)=(cldeg+c1)*dphyup+c2*cs2dp2*ph1yup(1)
        mas6x6up(2,4,ly)=c2mu*(-(cldeg+c1)*(cldeg+c2)*dphyup
     &                 -c2lm1*ph0yup(1)-c4*cs2dp2*ph1yup(1)
     &                 +c2*(cldeg+c1)*ph1yup(2))
        mas6x6up(3,4,ly)=-dphyup+c2*ph1yup(2)/cldeg
        mas6x6up(4,4,ly)=c2mu*((cldeg+c2)*dphyup
     &                  -c2lm1*ph0yup(2)/cldeg
     &                  +c2*(cs2dp2*ph1yup(1)-ph1yup(2)/cldeg))
        mas6x6up(5,4,ly)=cga(ly)*cs2dp2*ps0yup(1)
        mas6x6up(6,4,ly)=-cga(ly)*(cldeg+c1)*dphyup
      else
        spb(1,3)=c1-zjup(ldeg,ly,2)
        spb(2,3)=-ci*zjup(ldeg,ly,2)
        spb(1,4)=c1
        spb(2,4)=zhup(ldeg,ly,2)
        do j=3,4
          mas6x6up(1,j,ly)=-cllp1*spb(1,j)
          mas6x6up(2,j,ly)=c2mu*cllp1*((c1-cldeg)*spb(1,j)+spb(2,j))
          mas6x6up(3,j,ly)=-(cldeg+c1)*spb(1,j)+spb(2,j)
          mas6x6up(4,j,ly)=cmu(ly)*((cxs2-c2*(cldeg**2-c1))*spb(1,j)
     &                    -c2*spb(2,j))
          mas6x6up(5,j,ly)=(0.d0,0.d0)
          mas6x6up(6,j,ly)=cga(ly)*cllp1*spb(1,j)
        enddo
      endif
c
      do j=5,6
        do i=1,6
          mas6x6up(i,j,ly)=(0.d0,0.d0)
        enddo
      enddo
      mas6x6up(5,5,ly)=c1
      mas6x6up(6,5,ly)=c2lp1
      mas6x6up(5,6,ly)=c1
c
      if(ly.eq.lylw.or.lwup.lt.0)return
c
c     for lower radius
c
      cxp2=(kp(ly)*crrlw(ly))**2
      cxs2=(ks(ly)*crrlw(ly))**2
c
      if(ksmallp(ly).or.ksmalls(ly))then
        ph0jlw(1)=spbphj(ldeg,kp(ly)*crrlw(ly))
        ph1jlw(1)=spbphj(ldeg+1,kp(ly)*crrlw(ly))
        ps0jlw(1)=spbpsj(ldeg,kp(ly)*crrlw(ly))
c
        ph0ylw(1)=spbphy(ldeg,kp(ly)*crrlw(ly))
        ph1ylw(1)=spbphy(ldeg-1,kp(ly)*crrlw(ly))
        ps0ylw(1)=spbpsy(ldeg,kp(ly)*crrlw(ly))
      endif
c
      if(ksmallp(ly))then
c
        mas6x6lw(1,1,ly)=cldeg*ph0jlw(1)-cxp2*ph1jlw(1)/c2lp3
        mas6x6lw(2,1,ly)=(-cxi*cxp2+c2mu*cllm1)*ph0jlw(1)
     &                   +c2*c2mu*cxp2*ph1jlw(1)/c2lp3
        mas6x6lw(3,1,ly)=ph0jlw(1)
        mas6x6lw(4,1,ly)=c2mu*((cldeg-c1)*ph0jlw(1)
     &                   -cxp2*ph1jlw(1)/c2lp3)
        mas6x6lw(5,1,ly)=cga(ly)*ph0jlw(1)
        mas6x6lw(6,1,ly)=cga(ly)*(cldeg+c1)*ph0jlw(1)
c
        mas6x6lw(1,2,ly)=-(cldeg+c1)*ph0ylw(1)
     &                   +cxp2*ph1ylw(1)/c2lm1
        mas6x6lw(2,2,ly)=(-cxi*cxp2+c2mu*(cldeg+c1)*(cldeg+c2))
     &                  *ph0ylw(1)-c2*c2mu*cxp2*ph1ylw(1)/c2lm1
        mas6x6lw(3,2,ly)=ph0ylw(1)
        mas6x6lw(4,2,ly)=c2mu*(-(cldeg+c2)*ph0ylw(1)
     &                         +cxp2*ph1ylw(1)/c2lm1)
        mas6x6lw(5,2,ly)=cga(ly)*ph0ylw(1)
        mas6x6lw(6,2,ly)=cga(ly)*(cldeg+c1)*ph0ylw(1)
      else
        spb(1,1)=c1-zjlw(ldeg,ly,1)
        spb(2,1)=-ci*zjlw(ldeg,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhlw(ldeg,ly,1)
        do j=1,2
          mas6x6lw(1,j,ly)=cldeg*spb(1,j)-spb(2,j)
          mas6x6lw(2,j,ly)=(-cxi*cxp2+c2mu*cllm1)*spb(1,j)
     &                  +c2*c2mu*spb(2,j)
          mas6x6lw(3,j,ly)=spb(1,j)
          mas6x6lw(4,j,ly)=c2mu*((cldeg-c1)*spb(1,j)-spb(2,j))
          mas6x6lw(5,j,ly)=cga(ly)*spb(1,j)
          mas6x6lw(6,j,ly)=cga(ly)*(cldeg+c1)*spb(1,j)
        enddo
      endif
      if(ksmalls(ly))then
        ph0jlw(2)=spbphj(ldeg,ks(ly)*crrlw(ly))
        ph1jlw(2)=spbphj(ldeg+1,ks(ly)*crrlw(ly))
        ps0jlw(2)=spbpsj(ldeg,ks(ly)*crrlw(ly))
c
        ph0ylw(2)=spbphy(ldeg,ks(ly)*crrlw(ly))
        ph1ylw(2)=spbphy(ldeg-1,ks(ly)*crrlw(ly))
        ps0ylw(2)=spbpsy(ldeg,ks(ly)*crrlw(ly))
c
        cs2dp2=(cvs(ly)/cvp(ly))**2
        dphjlw=spbdphj(ldeg,kp(ly)*crrlw(ly),ks(ly)*crrlw(ly))
     &            *(c1-(cvs(ly)/cvp(ly))**2)
        dphylw=spbdphy(ldeg,kp(ly)*crrlw(ly),ks(ly)*crrlw(ly))
     &        *(c1-(cvs(ly)/cvp(ly))**2)
c
        mas6x6lw(1,3,ly)=cldeg*dphjlw-c2*cs2dp2*ph1jlw(1)
        mas6x6lw(2,3,ly)=c2mu*(cllm1*dphjlw-c2lp3*ph0jlw(1)
     &                    +c4*cs2dp2*ph1jlw(1)+c2*cldeg*ph1jlw(2))
        mas6x6lw(3,3,ly)=dphjlw+c2*ph1jlw(2)/(cldeg+c1)
        mas6x6lw(4,3,ly)=c2mu*((cldeg-c1)*dphjlw
     &                        +c2lp3*ph0jlw(2)/(cldeg+c1)
     &                -c2*(cs2dp2*ph1jlw(1)+ph1jlw(2)/(cldeg+c1)))
        mas6x6lw(5,3,ly)=-cga(ly)*cs2dp2*ps0jlw(1)
        mas6x6lw(6,3,ly)=cga(ly)*((cldeg+c1)*dphjlw
     &                            -ps0jlw(2)*c2lp1)
c
        mas6x6lw(1,4,ly)=(cldeg+c1)*dphylw+c2*cs2dp2*ph1ylw(1)
        mas6x6lw(2,4,ly)=c2mu*(-(cldeg+c1)*(cldeg+c2)*dphylw
     &                 -c2lm1*ph0ylw(1)-c4*cs2dp2*ph1ylw(1)
     &                 +c2*(cldeg+c1)*ph1ylw(2))
        mas6x6lw(3,4,ly)=-dphylw+c2*ph1ylw(2)/cldeg
        mas6x6lw(4,4,ly)=c2mu*((cldeg+c2)*dphylw
     &                  -c2lm1*ph0ylw(2)/cldeg
     &                  +c2*(cs2dp2*ph1ylw(1)-ph1ylw(2)/cldeg))
        mas6x6lw(5,4,ly)=cga(ly)*cs2dp2*ps0ylw(1)
        mas6x6lw(6,4,ly)=-cga(ly)*(cldeg+c1)*dphylw
      else
        spb(1,3)=c1-zjlw(ldeg,ly,2)
        spb(2,3)=-ci*zjlw(ldeg,ly,2)
        spb(1,4)=c1
        spb(2,4)=zhlw(ldeg,ly,2)
        do j=3,4
          mas6x6lw(1,j,ly)=-cllp1*spb(1,j)
          mas6x6lw(2,j,ly)=c2mu*cllp1*((c1-cldeg)*spb(1,j)+spb(2,j))
          mas6x6lw(3,j,ly)=-(cldeg+c1)*spb(1,j)+spb(2,j)
          mas6x6lw(4,j,ly)=cmu(ly)*((cxs2-c2*(cldeg**2-c1))*spb(1,j)
     &                    -c2*spb(2,j))
          mas6x6lw(5,j,ly)=(0.d0,0.d0)
          mas6x6lw(6,j,ly)=cga(ly)*cllp1*spb(1,j)
        enddo
      endif
c
      do j=5,6
        do i=1,6
          mas6x6lw(i,j,ly)=(0.d0,0.d0)
        enddo
      enddo
      mas6x6lw(5,5,ly)=c1
      mas6x6lw(6,5,ly)=c2lp1
      mas6x6lw(5,6,ly)=c1
c
      if(lwup.eq.0)then
c
c       calculate inverse matrix at upper radius
c
        do j=1,6
          do i=1,6
            mas(i,j)=mas6x6up(i,j,ly)
            mas6x6inv(i,j,ly)=(0.d0,0.d0)
          enddo
          mas6x6inv(j,j,ly)=c1
        enddo
      else
c
c       calculate inverse matrix at lower radius
c
        do j=1,6
          do i=1,6
            mas(i,j)=mas6x6lw(i,j,ly)
            mas6x6inv(i,j,ly)=(0.d0,0.d0)
          enddo
          mas6x6inv(j,j,ly)=c1
        enddo
      endif
      key=0
      call cdsvd500(mas,mas6x6inv(1,1,ly),6,6,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsmat: anormal exit from cdsvd500!'
        return
      endif
c
      return
      end