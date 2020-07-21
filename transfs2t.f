      subroutine transfs2t(nf,nr,df,t0,dt,dtout,fi,shift,cy0,cy,
     &                     cswap,dswap,iof,ntcutout,rname,y0,y,outfile)
      implicit none
c
      integer*4 nf,nr,iof,ntcutout
      real*8 df,t0,dt,dtout,fi
      real*8 shift(nr),dswap(4*nf),y0(nr),y(nr)
      complex*16 cy0(nf,nr),cy(nf,nr),cswap(2*nf)
      character*10 rname(nr)
      character*80 outfile
c
      real*8 PI2
      data PI2/6.28318530717959d0/
c
      integer*4 ir,lf,mf,it,itout
      real*8 t,sr,si,a,b
      complex*16 s
c
      do ir=1,nr
        do lf=1,nf
          s=dcmplx(-fi,dble(lf-1)*df)*dcmplx(PI2*shift(ir),0.d0)
          cswap(lf)=cy0(lf,ir)*cdexp(s)
        enddo
        mf=1
        do lf=2*nf,nf+2,-1
          mf=mf+1
          cswap(lf)=dconjg(cswap(mf))
        enddo
        cswap(nf+1)=(0.d0,0.d0)
c
c       convention for Fourier transform:
c       f(t)=\int F(f) exp(i2\pi f t) df
c
        call four1w(cswap,dswap,2*nf,+1)
c
        do lf=1,2*nf
          t=dble(lf-1)*dt
          cswap(lf)=cswap(lf)*dcmplx(dexp(-PI2*fi*t)*df,0.d0)
        enddo
c
        do lf=1,(ntcutout+1)/2
          t=dble(2*lf-1)*dtout
          it=1+idint(t/dt)
          if(it.ge.2*nf)then
            sr=dreal(cswap(2*nf))
          else
            b=dmod(t,dt)/dt
            a=1.d0-b
            sr=a*dreal(cswap(it))+b*dreal(cswap(it+1))
          endif
c
          t=dble(2*lf)*dtout
          it=1+idint(t/dt)
          if(it.ge.2*nf)then
            si=dreal(cswap(2*nf))
          else
            b=dmod(t,dt)/dt
            a=1.d0-b
            si=a*dreal(cswap(it))+b*dreal(cswap(it+1))
          endif
c
          cy(lf,ir)=dcmplx(sr,si)
        enddo
      enddo
c
      open(20,file=outfile,status='unknown')
      write(20,'(a,$)')' TIME       '
      do ir=1,nr-1
        write(20,'(a16,$)')rname(ir)
      enddo
      write(20,'(a16)')rname(nr)
      if(iof.eq.0)then
        do ir=1,nr
          y(ir)=0.d0
        enddo
        do it=1,ntcutout
          t=t0+dble(it-1)*dtout
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=y(ir)+dreal(cy(lf,ir))*dtout
            enddo
          else
            do ir=1,nr
              y(ir)=y(ir)+dimag(cy(lf,ir))*dtout
            enddo
          endif
        enddo
      else if(iof.eq.1)then
        do it=1,ntcutout
          t=t0+dble(it-1)*dtout
          write(20,'(f12.3,$)')t
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(cy(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(cy(lf,ir))
            enddo
          endif
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
      else if(iof.eq.2)then
        do ir=1,nr
          y0(ir)=dreal(cy(1,ir))
        enddo
        do it=1,ntcutout
          t=t0+dble(it-1)*dtout
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(cy(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(cy(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')(y(ir)-y0(ir))/dtout
          enddo
          write(20,'(E16.8)')(y(nr)-y0(nr))/dtout
          do ir=1,nr
            y0(ir)=y(ir)
          enddo
        enddo
      else
        stop ' Error in transfs2t: bad iof value!'
      endif
      close(20)
      return
      end