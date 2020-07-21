      subroutine butterworth(n,fc,df,nf,lpf)
      implicit none
      integer*4 n,nf
      real*8 fc,df
      complex*16 lpf(nf)
c
      integer*4 l,k
      complex*16 s,cs
c
      real*8 PI
      data PI/3.14159265358979d0/
c
      if(n.le.0.or.fc.le.0.d0)then
        do l=1,nf
          lpf(l)=(1.d0,0.d0)
        enddo
      else if(n.eq.2*(n/2))then
        do l=1,nf
          s=dcmplx(0.d0,df*dble(l-1)/fc)
          lpf(l)=(1.d0,0.d0)
          do k=1,n/2
            cs=dcmplx(2.d0*dcos(PI*dble(2*k+n-1)/dble(2*n)),0.d0)
            lpf(l)=lpf(l)/(s*s-s*cs+(1.d0,0.d0))
          enddo
        enddo
      else
        do l=1,nf
          s=dcmplx(0.d0,df*dble(l-1)/fc)
          lpf(l)=(1.d0,0.d0)/(s+(1.d0,0.d0))
          do k=1,(n-1)/2
            cs=dcmplx(2.d0*dcos(PI*dble(2*k+n-1)/dble(2*n)),0.d0)
            lpf(l)=lpf(l)/(s*s-s*cs+(1.d0,0.d0))
          enddo
        enddo
      endif
      return
      end
