	subroutine disazi(rearth,lateq,loneq,latst,lonst,xnorth,yeast)
	implicit none
c
	real*8 rearth,lateq,loneq,latst,lonst,xnorth,yeast
c
	integer*4 iangle
	real*8 latb,lonb,latc,lonc,angleb,anglec
	real*8 dis,a,b,c,s,aa
c
	real*8 PI,PI2,DEGTORAD
	data PI,PI2/3.14159265358979d0,6.28318530717959d0/
	data DEGTORAD/1.745329251994328d-02/
c
c	A is north pole
c     B is station
c     C is epicentre (local cartesian coordinate origin)
c
	latb=lateq*DEGTORAD
	lonb=loneq*DEGTORAD
	latc=latst*DEGTORAD
	lonc=lonst*DEGTORAD
c
	if(lonb.lt.0.d0)lonb=lonb+PI2
	if(lonc.lt.0.d0)lonc=lonc+PI2
c
	b=0.5d0*PI-latb
	c=0.5d0*PI-latc
	if(lonc.gt.lonb)then
	  aa=lonc-lonb
	  if(aa.le.PI)then
	    iangle=1
	  else
	    aa=PI2-aa
	    iangle=-1
	  endif
	else
	  aa=lonb-lonc
	  if(aa.le.PI)then
	    iangle=-1
	  else
	    aa=PI2-aa
	    iangle=1
	  endif
	endif
	s=dcos(b)*dcos(c)+dsin(b)*dsin(c)*dcos(aa)
	a=dacos(dsign(dmin1(dabs(s),1.d0),s))
	dis=a*rearth
	if(a*b*c.eq.0.d0)then
	  angleb=0.d0
	  anglec=0.d0
	else
	  s=0.5d0*(a+b+c)
	  a=dmin1(a,s)
	  b=dmin1(b,s)
	  c=dmin1(c,s)
	  anglec=2.d0*dasin(dmin1(1.d0,
     &            dsqrt(dsin(s-a)*dsin(s-b)/(dsin(a)*dsin(b)))))
	  angleb=2.d0*dasin(dmin1(1.d0,
     &            dsqrt(dsin(s-a)*dsin(s-c)/(dsin(a)*dsin(c)))))
	  if(iangle.eq.1)then
	    angleb=PI2-angleb
	  else
	    anglec=PI2-anglec
	  endif
	endif
c
c	cartesian coordinates of the station
c
	xnorth=dis*dcos(anglec)
	yeast=dis*dsin(anglec)
c
	return
	end
