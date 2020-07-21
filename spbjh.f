ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Subroutine: spbjh                                                         c
c     calculates spherical Bessel functions for complex argument                c
c     zj(x)=x*j_(n+1)(x)/[j_n(x)+ix*j_(n+1)(x)], x = xa,xb                      c
c     wj(xa,xb)=ln{[j_n(xa)+i*xa*j_(n+1)(xa)]/[j_n(xb)+i*xb*j_(n+1)(xb)]}       c
c     zh(x)=x*h_(n+1)(x)/h_n(x), x=xa,xb                                        c
c     wh(xa,xb)=ln[h_n(xa)/h_n(xb)], where h_n(x) = j_n(x) + sign[im(x)]iy_n(x) c
c     where xa = ck*ra, xb = ck*rb                                              c
c     for harmonic degree nn                                                    c
c                                                                               c
c     Input:                                                                    c
c     nn = max. order                                                           c
c     xa,xb = complex arguments                                                 c
c                                                                               c
c     Output:                                                                   c
c     complex zja(n,xa), zjb(n,xb), wj(n,xa,xb)                                 c
c             zha(n,xa), zhb(n,xb), wh(n,xa,xb)                                 c
c             for n = 0, ..., nn                                                c
c                                                                               c
c     First implemention: 17 June 2007                                          c
c     Rongjiang Wang                                                            c
c     GeoForschungsZentrum Potsdam, Germany                                     c
c     Email: wang@gfz-potsdam.de                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spbjh(nn,ck,ra,rb,zja,zjb,wj,zha,zhb,wh)
      implicit none
      integer*4 nn
      real*8 ra,rb
      complex*16 ck
      complex*16 zja(0:nn),zjb(0:nn),wj(0:nn)
      complex*16 zha(0:nn),zhb(0:nn),wh(0:nn)
c
c     local memory
c
      integer*4 n,nx,nmax
      real*8 xra,xrb,xia,xib
      complex*16 xa,xb,xa2,xb2,ga,gb,gn,bda
      complex*16 zja1,zjb1,zja2,zjb2
      complex*16 cep0,cep1,cem0,cem1,cfac
      complex*16 c2xj0,c2x2j1
c
      real*8 pi,eps
      complex*16 ci,c1,c2
      data pi,eps/3.14159265358979d0,1.0d-08/
      data ci,c1,c2/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0)/
c
      xa=ck*dcmplx(ra,0.d0)
      xb=ck*dcmplx(rb,0.d0)
c
      xra=dreal(xa)
      xrb=dreal(xb)
      xia=dimag(xa)
      xib=dimag(xb)
      xa2=xa*xa
      xb2=xb*xb
c
      cep0=cdexp(dcmplx(0.d0,xra))
      cem0=c1/cep0
      cfac=dcmplx(dexp(-2.d0*dabs(xia)),0.d0)
      if(xia.ge.0.d0)then
        c2xj0=(cfac*cep0-cem0)/ci
        c2x2j1=c2xj0-(cfac*cep0+cem0)*xa
        zja(0)=c1-c1/(c1+ci+xa*(cfac*cep0+cem0)/(cfac*cep0-cem0))
      else
        c2xj0=(cep0-cfac*cem0)/ci
        c2x2j1=c2xj0-(cep0+cfac*cem0)*xa
        zja(0)=c1-c1/(c1+ci+xa*(cep0+cfac*cem0)/(cep0-cfac*cem0))
      endif
      wj(0)=c2xj0+ci*c2x2j1
c
      if(nn.le.0)goto 101
c
      zja2=(0.d0,0.d0)
      nx=10+nn/4
100   nmax=nn+nx
      gn=dcmplx(dble(2*nmax+3),0.d0)
      zja1=ci*xa2/(gn+ci*xa2)
      do n=nmax-1,nn,-1
        ga=dcmplx(dble(2*n+3),0.d0)+ci*xa2
        zja1=ci*xa2*(c1-zja1)/(ga*(c1-zja1)+ci*zja1)
      enddo
      if(cdabs(zja1-zja2).gt.eps*cdabs(zja1))then
        zja2=zja1
        nx=2*nx
        goto 100
      endif
      zja(nn)=zja1
c
      do n=nn-1,1,-1
        ga=dcmplx(dble(2*n+3),0.d0)+ci*xa2
        zja(n)=ci*xa2*(c1-zja(n+1))/(ga*(c1-zja(n+1))+ci*zja(n+1))
      enddo
c
101   continue
c
      if(rb.le.0.d0)return
c
      bda=dcmplx(dlog(rb/ra),0.d0)
c
      cep0=cdexp(dcmplx(0.d0,xrb))
      cem0=c1/cep0
      cfac=dcmplx(dexp(-2.d0*dabs(xib)),0.d0)
      if(xib.ge.0.d0)then
        c2xj0=(cfac*cep0-cem0)/ci
        c2x2j1=c2xj0-(cfac*cep0+cem0)*xb
        zjb(0)=c1-c1/(c1+ci+xb*(cfac*cep0+cem0)/(cfac*cep0-cem0))
      else
        c2xj0=(cep0-cfac*cem0)/ci
        c2x2j1=c2xj0-(cep0+cfac*cem0)*xb
        zjb(0)=c1-c1/(c1+ci+xb*(cep0+cfac*cem0)/(cep0-cfac*cem0))
      endif
      wj(0)=cdlog(wj(0)/(c2xj0+ci*c2x2j1))
     &     +bda+dcmplx(dabs(xia)-dabs(xib),0.d0)
c
      if(nn.le.0)goto 201
c
      zjb2=(0.d0,0.d0)
      nx=10+nn/4
200   nmax=nn+nx
      gn=dcmplx(dble(2*nmax+3),0.d0)
      zjb1=ci*xb2/(gn+ci*xb2)
      do n=nmax-1,nn,-1
        gb=dcmplx(dble(2*n+3),0.d0)+ci*xb2
        zjb1=ci*xb2*(c1-zjb1)/(gb*(c1-zjb1)+ci*zjb1)
      enddo
      if(cdabs(zjb1-zjb2).gt.eps*cdabs(zjb1))then
        zjb2=zjb1
        nx=2*nx
        goto 200
      endif
      zjb(nn)=zjb1
c
      do n=nn-1,1,-1
        gb=dcmplx(dble(2*n+3),0.d0)+ci*xb2
        zjb(n)=ci*xb2*(c1-zjb(n+1))/(gb*(c1-zjb(n+1))+ci*zjb(n+1))
      enddo
c
201   continue
c
      do n=1,nn
        ga=xa2+(c1-xa2+dcmplx(0.d0,dble(2*n+1)))*zja(n-1)
        gb=xb2+(c1-xb2+dcmplx(0.d0,dble(2*n+1)))*zjb(n-1)
        wj(n)=wj(n-1)+cdlog(ga/gb)+bda
      enddo
c
      if(xia+xib.ge.0.d0)then
        zha(0)=c1-ci*xa
        zhb(0)=c1-ci*xb
      else
        zha(0)=c1+ci*xa
        zhb(0)=c1+ci*xb
      endif
      wh(0)=zhb(0)-zha(0)+bda
c
      do n=1,nn
        gn=dcmplx(dble(2*n+1),0.d0)
        zha(n)=gn-xa2/zha(n-1)
        zhb(n)=gn-xb2/zhb(n-1)
        wh(n)=wh(n-1)+cdlog(zha(n-1)/zhb(n-1))+bda
      enddo
      return
      end