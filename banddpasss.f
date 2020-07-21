      subroutine bandpass(n,f1,f2,df,nf,h)
      implicit none
      integer*4 n,nf
      real*8 f1,f2,df
      complex*16 h(nf)
c
      integer*4 l,k
      real*8 w,w1,w2,delt,arg
c
      real*8 PI
      data PI/3.14159265358979d0/
c
      if(n.le.0.or.f1.ge.f2)then
        do l=1,nf
          h(l)=(1.d0,0.d0)
        enddo
      else if(f1.le.0.d0)then
        w2=2.d0*PI*f2
        do l=1,nf
          w=2.d0*PI*dble(l-1)*df
          delt=w2-w1
          h(l)=(1.d0,0.d0)
          do k=1,n
            arg=PI*dble(2*k-1)/dble(2*n)
            h(l)=h(l)*dcmplx(0.d0,-delt)
     &          /dcmplx(w-delt*dcos(arg),-delt*dsin(arg))
          enddo
        enddo
      else
        w1=2.d0*PI*f1
        w2=2.d0*PI*f2
        do l=1,nf
          w=2.d0*PI*dble(l-1)*df
          delt=w*(w2-w1)
          h(l)=(1.d0,0.d0)
          do k=1,n
            arg=PI*dble(2*k-1)/dble(2*n)
            h(l)=h(l)*dcmplx(0.d0,-delt)
     &          /dcmplx(w*w-w1*w2-delt*dcos(arg),-delt*dsin(arg))
          enddo
        enddo
      endif
      return
      end
