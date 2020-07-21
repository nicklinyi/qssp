ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Subroutine: spbh12                                                        c
c     calculates spherical Bessel functions for complex argument                c
c     zh^(1,2)(x)=x*h^(1,2)_(n+1)(x)/h^(1,2)_n(x),                              c
c     for harmonic degree nn                                                    c
c                                                                               c
c     Input:                                                                    c
c     nn = max. order                                                           c
c                                                                               c
c     Output:                                                                   c
c     complex zh^(1,2)(n,x)                                                     c
c             for n = 0, ..., nn                                                c
c                                                                               c
c     First implemention: 29 November 2016                                      c
c     Rongjiang Wang                                                            c
c     GeoForschungsZentrum Potsdam, Germany                                     c
c     Email: wang@gfz-potsdam.de                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spbh12(n,nmax,cx,zh12)
      implicit none
      integer*4 n,nmax
      complex*16 cx
      complex*16 zh12(0:nmax,2)
c
c     local memory
c
      integer*4 i
      complex*16 cx2
c
      cx2=cx*cx
c
      zh12(0,1)=(1.d0,0.d0)-(0.d0,1.d0)*cx
      zh12(0,2)=(1.d0,0.d0)+(0.d0,1.d0)*cx
c
      do i=1,n
        zh12(i,1)=dcmplx(dble(2*i+1),0.d0)-cx2/zh12(i-1,1)
        zh12(i,2)=dcmplx(dble(2*i+1),0.d0)-cx2/zh12(i-1,2)
      enddo
c
      return
      end