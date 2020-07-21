      subroutine moments(munit,strike,dip,rake,
     &                   mtt,mpp,mrr,mtp,mpr,mrt)
      implicit none
      real*8 munit,strike,dip,rake,
     &       mtt,mpp,mrr,mtp,mpr,mrt
c
      real*8 DEG2RAD
      parameter(DEG2RAD=1.745329251994328d-02)
c
      real*8 st,di,ra
      real*8 sss,sss2,ss2s,css,css2,cs2s
      real*8 ssd,ss2d,csd,cs2d,ssr,csr
c
      st=strike*DEG2RAD
      di=dip*DEG2RAD
      ra=rake*DEG2RAD
c
      sss=dsin(st)
      sss2=sss*sss
      ss2s=dsin(2.d0*st)
      css=dcos(st)
      css2=css*css
      cs2s=dcos(2.d0*st)
c
      ssd=dsin(di)
      ss2d=dsin(2.d0*di)
      csd=dcos(di)
      cs2d=dcos(2.d0*di)
c
      ssr=dsin(ra)
      csr=dcos(ra)
c
      mtt=munit*(-ssd*csr*ss2s-ss2d*ssr*sss2)
      mpp=munit*(ssd*csr*ss2s-ss2d*ssr*css2)
      mrr=-(mtt+mpp)
      mtp=munit*(-ssd*csr*cs2s-0.5d0*ss2d*ssr*ss2s)
      mpr=munit*(csd*csr*sss-cs2d*ssr*css)
      mrt=munit*(-csd*csr*css-cs2d*ssr*sss)
c
      return
      end