ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Function: spbphj                                                          c
c     calculates eq.(103) of Takeuchi & Saito related to spherical Bessel       c
c     functions for complex argument                                            c
c                                                                               c
c     Input:                                                                    c
c     n = order                                                                 c
c     x = complex argument                                                      c
c                                                                               c
c     Return:                                                                   c
c     phj(n,x)=j_(n)(x)*(2*n+1)!!/x^n                                           c
c     for |x/2|^2 << n                                                          c
c                                                                               c
c     First implemention: 17 June 2007                                          c
c     Rongjiang Wang                                                            c
c     GeoForschungsZentrum Potsdam, Germany                                     c
c     Email: wang@gfz-potsdam.de                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function spbphj(n,x)
      implicit none
      integer*4 n
      complex*16 x
c
c     local memory
c
      integer*4 i
      complex*16 cx2,fn,cy
c
      real*8 eps
      data eps/1.0d-08/
c
      cx2=(0.5d0,0.d0)*x*x
      i=0
      fn=(1.d0,0.d0)
      cy=(1.d0,0.d0)
100   i=i+1
      fn=-fn*cx2/dcmplx(dble((2*(n+i)+1)*i),0.d0)
      cy=cy+fn
      if(cdabs(fn).gt.eps*cdabs(cy))goto 100
c
      spbphj=cy
c
      return
      end