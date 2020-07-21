ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Function: spbdpsy                                                         c
c     calculates eq.(103) of Takeuchi & Saito related to spherical Bessel       c
c     functions for complex argument                                            c
c                                                                               c
c     Input:                                                                    c
c     n = order                                                                 c
c     x = complex argument                                                      c
c                                                                               c
c     Return:                                                                   c
c     dphy(n,xa,xb)=(phy(n,xa)-phy(n,xb))*2*(2*n-1)/(xa^2-xb^2)                 c
c     for |x/2|^2 << n                                                          c
c                                                                               c
c     First implemention: 17 June 2007                                          c
c     Rongjiang Wang                                                            c
c     GeoForschungsZentrum Potsdam, Germany                                     c
c     Email: wang@gfz-potsdam.de                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function spbdphy(n,xa,xb)
      implicit none
      integer*4 n
      complex*16 xa,xb
c
c     local memory
c
      integer*4 i
      complex*16 cxa2,cxb2,fn,gn,ca,cy
c
      real*8 eps
      data eps/1.0d-08/
c
      cxa2=(0.5d0,0.d0)*xa*xa
      cxb2=(0.5d0,0.d0)*xb*xb
c
      i=1
      ca=(1.d0,0.d0)
      gn=ca
      fn=ca
      cy=fn
100   i=i+1
      ca=(1.d0,0.d0)/dcmplx(dble((2*(n-i)+1)*i),0.d0)
      gn=gn*cxb2*ca
      fn=gn+fn*cxa2*ca
      cy=cy+fn
      if(cdabs(fn).gt.eps*cdabs(cy))goto 100
      spbdphy=cy
      return
      end