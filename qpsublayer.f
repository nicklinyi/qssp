	subroutine qpsublayer(ierr)
      use qpalloc
	implicit none
      integer*4 ierr
c
c	work space
c
	integer*4 i,j,i0,l,ly,lya,lyb,lyp,ig,di0
      real*8 f,h,dh,z,zz,slw,wvlcut,up,lw,uplw4,alfa,beta
	real*8 rrs,rrr,dvp,dvs,ro1,dro,dqp,dqs,mass,gr0
      real*8 bvn2
      logical*2 atmo,action
c
      lymax=0
      do l=1,l0
        if(vs0up(l).gt.0.d0)then
          wvlcut=vs0up(l)/fcut
        else
          wvlcut=vp0up(l)/fcut
        endif
	  h=dp0lw(l)-dp0up(l)
	  dvp=2.d0*dabs(vp0lw(l)-vp0up(l))/(vp0lw(l)+vp0up(l))
        if(vs0lw(l)+vs0up(l).gt.0.d0)then
	    dvs=2.d0*dabs(vs0lw(l)-vs0up(l))/(vs0lw(l)+vs0up(l))
        else
          dvs=0.d0
        endif
        if(dp0lw(l).le.depatmos)then
          atmo=.true.
          dro=2.d0*dabs(dlog(ro0lw(l))-dlog(ro0up(l)))
     &       /(dlog(ro0lw(l))+dlog(ro0up(l)))
        else
          atmo=.false.
          dro=2.d0*dabs(ro0lw(l)-ro0up(l))/(ro0lw(l)+ro0up(l))
        endif
        i0=1+idint(dmax1(dvp/RESOLUT,dvs/RESOLUT,dro/RESOLUT))
        i0=min0(i0,1+idint(5.d0*h/wvlcut))
        lymax=lymax+i0
      enddo
      lymax=lymax+lyadd+1
c
      allocate(ldegpsv(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: ldegpsv not allocated!'
      allocate(ldegsh(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: ldegsh not allocated!'
      allocate(nruku(0:ldegmax,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: ldegsh not allocated!'
c
      allocate(rrup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: rrup not allocated!'
      allocate(rrlw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: rrlw not allocated!'
      allocate(vpup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: vpup not allocated!'
      allocate(vplw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: vplw not allocated!'
      allocate(vsup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: vsup not allocated!'
      allocate(vslw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: vslw not allocated!'
      allocate(roup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: roup not allocated!'
      allocate(rolw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: rolw not allocated!'
      allocate(qpup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: qpup not allocated!'
      allocate(qplw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: qplw not allocated!'
      allocate(qsup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: qsup not allocated!'
      allocate(qslw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: qslw not allocated!'
c
      allocate(crrup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: crrup not allocated!'
      allocate(crrlw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: crrlw not allocated!'
c
      allocate(cvp(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvp not allocated!'
      allocate(cvpup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvpup not allocated!'
      allocate(cvplw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvplw not allocated!'
c
      allocate(cvs(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvs not allocated!'
      allocate(cvsup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvsup not allocated!'
      allocate(cvslw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cvslw not allocated!'
c
      allocate(cro(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cro not allocated!'
      allocate(croup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: croup not allocated!'
      allocate(crolw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: crolw not allocated!'
c
      allocate(cla(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cla not allocated!'
      allocate(claup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: claup not allocated!'
      allocate(clalw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: clalw not allocated!'
c
      allocate(cmu(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cmu not allocated!'
      allocate(cmuup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cmuup not allocated!'
      allocate(cmulw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cmulw not allocated!'
c
      allocate(cgr(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cgr not allocated!'
      allocate(cgrup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cgrup not allocated!'
      allocate(cgrlw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cgrlw not allocated!'
c
      allocate(cga(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cga not allocated!'
      allocate(cgaup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cgaup not allocated!'
      allocate(cgalw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cgalw not allocated!'
c
      allocate(cua(8,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cua not allocated!'
      allocate(cypnorm(6,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpsublayer: cypnorm not allocated!'
c
	ly0=0
c
      zz=0.d0
	do l=1,l0
        if(vs0up(l).gt.0.d0)then
          wvlcut=vs0up(l)/fcut
        else
          wvlcut=vp0up(l)/fcut
        endif
	  h=dp0lw(l)-dp0up(l)
	  dvp=2.d0*dabs(vp0lw(l)-vp0up(l))/(vp0lw(l)+vp0up(l))
        if(vs0lw(l)+vs0up(l).gt.0.d0)then
	    dvs=2.d0*dabs(vs0lw(l)-vs0up(l))/(vs0lw(l)+vs0up(l))
        else
          dvs=0.d0
        endif
        if(dp0lw(l).le.depatmos)then
          atmo=.true.
          dro=2.d0*dabs(dlog(ro0lw(l))-dlog(ro0up(l)))
     &       /(dlog(ro0lw(l))+dlog(ro0up(l)))
        else
          atmo=.false.
          dro=2.d0*dabs(ro0lw(l)-ro0up(l))/(ro0lw(l)+ro0up(l))
        endif
        i0=1+idint(dmax1(dvp/RESOLUT,dvs/RESOLUT,dro/RESOLUT))
        i0=min0(i0,1+idint(5.d0*h/wvlcut))
        dvp=(vp0lw(l)-vp0up(l))/h
	  dvs=(vs0lw(l)-vs0up(l))/h
        if(atmo)then
	    dro=(dlog(ro0lw(l))-dlog(ro0up(l)))/h
        else
          dro=(ro0lw(l)-ro0up(l))/h
        endif
	  dqp=(qp0lw(l)-qp0up(l))/h
	  dqs=(qs0lw(l)-qs0up(l))/h
	  dh=h/dble(i0)
	  do i=1,i0
	    ly0=ly0+1
	    if(ly0.gt.lymax-lyadd)then
	      stop ' Error in qpsublayer: lymax too small!'
	    endif
	    z=dble(i-1)*dh
	    rrup(ly0)=rratmos-(zz+z)
	    vpup(ly0)=vp0up(l)+dvp*z
	    vsup(ly0)=vs0up(l)+dvs*z
          if(atmo)then
	      roup(ly0)=ro0up(l)*dexp(dro*z)
          else
            roup(ly0)=ro0up(l)+dro*z
          endif
	    qpup(ly0)=qp0up(l)+dqp*z
	    qsup(ly0)=qs0up(l)+dqs*z
          z=z+dh
	    rrlw(ly0)=rratmos-(zz+z)
	    vplw(ly0)=vp0up(l)+dvp*z
	    vslw(ly0)=vs0up(l)+dvs*z
          if(atmo)then
	      rolw(ly0)=ro0up(l)*dexp(dro*z)
          else
            rolw(ly0)=ro0up(l)+dro*z
          endif
	    qplw(ly0)=qp0up(l)+dqp*z
	    qslw(ly0)=qs0up(l)+dqs*z
	  enddo
        zz=zz+h
      enddo
c
      deallocate(dp0,dp0up,dp0lw,vp0,vp0up,vp0lw,vs0,vs0up,vs0lw,
     &           ro0,ro0up,ro0lw,qp0,qp0up,qp0lw,qs0,qs0up,qs0lw)
c
c     add source layers
c
      do ig=1,ngrn
        rrs=rratmos-grndep(ig)
        if(rrs.lt.0.d0)then
          stop ' Wrong source depth!'
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 100
          endif
        enddo
100     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      vpup(ly+1)=vpup(ly)
	      vsup(ly+1)=vsup(ly)
	      roup(ly+1)=roup(ly)
	      qpup(ly+1)=qpup(ly)
	      qsup(ly+1)=qsup(ly)
            rrlw(ly+1)=rrlw(ly)
	      vplw(ly+1)=vplw(ly)
	      vslw(ly+1)=vslw(ly)
	      rolw(ly+1)=rolw(ly)
	      qplw(ly+1)=qplw(ly)
	      qslw(ly+1)=qslw(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
	    vplw(lys-1)=up*vpup(lys-1)+lw*vplw(lys)
	    vslw(lys-1)=up*vsup(lys-1)+lw*vslw(lys)
          if(rrs.gt.REARTH)then
            rolw(lys-1)=rolw(lys-1)
     &                 *dexp(dlog(roup(lys-1)/rolw(lys-1))*up)
          else
	      rolw(lys-1)=up*roup(lys-1)+lw*rolw(lys)
          endif
	    qplw(lys-1)=up*qpup(lys-1)+lw*qplw(lys)
	    qslw(lys-1)=up*qsup(lys-1)+lw*qslw(lys)
          rrup(lys)=rrs
	    vpup(lys)=vplw(lys-1)
	    vsup(lys)=vslw(lys-1)
	    roup(lys)=rolw(lys-1)
	    qpup(lys)=qplw(lys-1)
	    qsup(lys)=qslw(lys-1)
          ly0=ly0+1       
        endif
      enddo
c
c     add receiver layer
c
      rrr=rratmos-dpr
      if(rrr.lt.0.d0)then
        stop ' Wrong receiver depth!'
      endif
      do ly=1,ly0
        if(rrr.ge.rrlw(ly))then
          lyr=ly
          goto 200
        endif
      enddo
200   continue
      if(rrr.lt.rrup(lyr).and.rrr.gt.rrlw(lyr))then
        do ly=ly0,lyr,-1
          rrup(ly+1)=rrup(ly)
	    vpup(ly+1)=vpup(ly)
	    vsup(ly+1)=vsup(ly)
	    roup(ly+1)=roup(ly)
	    qpup(ly+1)=qpup(ly)
	    qsup(ly+1)=qsup(ly)
          rrlw(ly+1)=rrlw(ly)
	    vplw(ly+1)=vplw(ly)
	    vslw(ly+1)=vslw(ly)
	    rolw(ly+1)=rolw(ly)
	    qplw(ly+1)=qplw(ly)
	    qslw(ly+1)=qslw(ly)
        enddo
        lyr=lyr+1
        up=(rrr-rrlw(lyr))/(rrup(lyr-1)-rrlw(lyr))
        lw=1.d0-up
        rrlw(lyr-1)=rrr
	  vplw(lyr-1)=up*vpup(lyr-1)+lw*vplw(lyr)
	  vslw(lyr-1)=up*vsup(lyr-1)+lw*vslw(lyr)
        if(rrr.gt.REARTH)then
          rolw(lyr-1)=rolw(lyr-1)*dexp(dlog(roup(lyr-1)/rolw(lyr-1))*up)
        else
	    rolw(lyr-1)=up*roup(lyr-1)+lw*rolw(lyr)
        endif
	  qplw(lyr-1)=up*qpup(lyr-1)+lw*qplw(lyr)
	  qslw(lyr-1)=up*qsup(lyr-1)+lw*qslw(lyr)
        rrup(lyr)=rrr
	  vpup(lyr)=vplw(lyr-1)
	  vsup(lyr)=vslw(lyr-1)
	  roup(lyr)=rolw(lyr-1)
	  qpup(lyr)=qplw(lyr-1)
	  qsup(lyr)=qslw(lyr-1)
        ly0=ly0+1      
      endif
c
c     add min. path layer
c
      if(ipatha.eq.1)then
        rrr=rratmos-minpath
        if(rrr.lt.0.d0)then
          stop ' Wrong receiver depth!'
        endif
        do ly=1,ly0
          if(rrr.ge.rrlw(ly))then
            lya=ly
           goto 201
          endif
        enddo
201     continue
        if(rrr.lt.rrup(lya).and.rrr.gt.rrlw(lya))then
          do ly=ly0,lya,-1
            rrup(ly+1)=rrup(ly)
	      vpup(ly+1)=vpup(ly)
	      vsup(ly+1)=vsup(ly)
	      roup(ly+1)=roup(ly)
	      qpup(ly+1)=qpup(ly)
	      qsup(ly+1)=qsup(ly)
            rrlw(ly+1)=rrlw(ly)
	      vplw(ly+1)=vplw(ly)
	      vslw(ly+1)=vslw(ly)
	      rolw(ly+1)=rolw(ly)
	      qplw(ly+1)=qplw(ly)
	      qslw(ly+1)=qslw(ly)
          enddo
          lya=lya+1
          up=(rrr-rrlw(lya))/(rrup(lya-1)-rrlw(lya))
          lw=1.d0-up
          rrlw(lya-1)=rrr
	    vplw(lya-1)=up*vpup(lya-1)+lw*vplw(lya)
	    vslw(lya-1)=up*vsup(lya-1)+lw*vslw(lya)
          if(rrr.gt.REARTH)then
            rolw(lya-1)=rolw(lya-1)
     &                 *dexp(dlog(roup(lya-1)/rolw(lya-1))*up)
          else
	      rolw(lya-1)=up*roup(lya-1)+lw*rolw(lya)
          endif
	    qplw(lya-1)=up*qpup(lya-1)+lw*qplw(lya)
	    qslw(lya-1)=up*qsup(lya-1)+lw*qslw(lya)
          rrup(lya)=rrr
	    vpup(lya)=vplw(lya-1)
	    vsup(lya)=vslw(lya-1)
	    roup(lya)=rolw(lya-1)
	    qpup(lya)=qplw(lya-1)
	    qsup(lya)=qslw(lya-1)
          ly0=ly0+1
        endif
      endif
c
c     add max. path layer
c
      if(ipathb.eq.1)then
        rrr=rratmos-maxpath
        if(rrr.lt.0.d0)then
          stop ' Wrong receiver depth!'
        endif
        do ly=1,ly0
          if(rrr.ge.rrlw(ly))then
            lyb=ly
            goto 202
          endif
        enddo
202     continue
        if(rrr.lt.rrup(lyb).and.rrr.gt.rrlw(lyb))then
          do ly=ly0,lyb,-1
            rrup(ly+1)=rrup(ly)
	      vpup(ly+1)=vpup(ly)
	      vsup(ly+1)=vsup(ly)
	      roup(ly+1)=roup(ly)
	      qpup(ly+1)=qpup(ly)
	      qsup(ly+1)=qsup(ly)
            rrlw(ly+1)=rrlw(ly)
  	      vplw(ly+1)=vplw(ly)
	      vslw(ly+1)=vslw(ly)
	      rolw(ly+1)=rolw(ly)
	      qplw(ly+1)=qplw(ly)
	      qslw(ly+1)=qslw(ly)
          enddo
          lyb=lyb+1
          up=(rrr-rrlw(lyb))/(rrup(lyb-1)-rrlw(lyb))
          lw=1.d0-up
          rrlw(lyb-1)=rrr
	    vplw(lyb-1)=up*vpup(lyb-1)+lw*vplw(lyb)
	    vslw(lyb-1)=up*vsup(lyb-1)+lw*vslw(lyb)
          if(rrr.gt.REARTH)then
            rolw(lyb-1)=rolw(lyb-1)
     &                 *dexp(dlog(roup(lyb-1)/rolw(lyb-1))*up)
          else
	      rolw(lyb-1)=up*roup(lyb-1)+lw*rolw(lyb)
          endif
	    qplw(lyb-1)=up*qpup(lyb-1)+lw*qplw(lyb)
	    qslw(lyb-1)=up*qsup(lyb-1)+lw*qslw(lyb)
          rrup(lyb)=rrr
	    vpup(lyb)=vplw(lyb-1)
	    vsup(lyb)=vslw(lyb-1)
	    roup(lyb)=rolw(lyb-1)
	    qpup(lyb)=qplw(lyb-1)
	    qsup(lyb)=qslw(lyb-1)
          ly0=ly0+1      
        endif
      endif
c
c     divide thick layers
c
300   action=.false.
      do ly=1,ly0-1
        if(rrup(ly).gt.2.0d0*rrlw(ly))then
          if(ly0+1.gt.lymax)then
	      stop ' Error in qpsublayer: lymax too small!'
	    endif
          rrr=0.5d0*(rrup(ly)+rrlw(ly))
          lyp=ly
          action=.true.
          goto 301
        endif
      enddo
301   continue
      if(action)then
        do ly=ly0,lyp,-1
          rrup(ly+1)=rrup(ly)
	    vpup(ly+1)=vpup(ly)
	    vsup(ly+1)=vsup(ly)
	    roup(ly+1)=roup(ly)
	    qpup(ly+1)=qpup(ly)
	    qsup(ly+1)=qsup(ly)
          rrlw(ly+1)=rrlw(ly)
  	    vplw(ly+1)=vplw(ly)
	    vslw(ly+1)=vslw(ly)
	    rolw(ly+1)=rolw(ly)
	    qplw(ly+1)=qplw(ly)
	    qslw(ly+1)=qslw(ly)
        enddo
        lyp=lyp+1
        up=(rrr-rrlw(lyp))/(rrup(lyp-1)-rrlw(lyp))
        lw=1.d0-up
        rrlw(lyp-1)=rrr
	  vplw(lyp-1)=up*vpup(lyp-1)+lw*vplw(lyp)
	  vslw(lyp-1)=up*vsup(lyp-1)+lw*vslw(lyp)
        if(rrr.gt.REARTH)then
          rolw(lyp-1)=rolw(lyp-1)
     &               *dexp(dlog(roup(lyp-1)/rolw(lyp-1))*up)
        else
	    rolw(lyp-1)=up*roup(lyp-1)+lw*rolw(lyp)
        endif
	  qplw(lyp-1)=up*qpup(lyp-1)+lw*qplw(lyp)
	  qslw(lyp-1)=up*qsup(lyp-1)+lw*qslw(lyp)
        rrup(lyp)=rrr
	  vpup(lyp)=vplw(lyp-1)
	  vsup(lyp)=vslw(lyp-1)
	  roup(lyp)=rolw(lyp-1)
	  qpup(lyp)=qplw(lyp-1)
	  qsup(lyp)=qslw(lyp-1)
        ly0=ly0+1
        goto 300
      endif
310   continue
c
c     determine indices of main interfaces
c
      lyos=1
      lyob=1
      if(vsup(1).le.0.d0)then
        do ly=2,ly0
          if(vsup(ly).gt.0.d0.and.vslw(ly-1).le.0.d0)then
            lyob=ly
            goto 401
          endif
        enddo
401     continue
        do ly=lyob,2,-1
          if(roup(ly)-rolw(ly-1).gt.10.d0*rolw(ly-1))then
            lyos=ly
            goto 402
          endif
        enddo
402     continue
      endif
      lycm=ly0+1
      do ly=max0(2,lyob+1),ly0
        if(vsup(ly).le.0.d0.and.vslw(ly-1).gt.0.d0)then
          lycm=ly
          goto 403
        endif
      enddo
403   lycc=ly0+1
      do ly=max0(2,lycm+1),ly0
        if(vsup(ly).gt.0.d0.and.vslw(ly-1).le.0.d0)then
          lycc=ly
          goto 404
        endif
      enddo
404   continue
c
c     determine indices of receiver layer
c
      rrr=rratmos-dpr
      do ly=1,ly0
        if(rrr.ge.rrup(ly))then
          lyr=ly
          goto 500
        endif
      enddo
500   continue
c
c     determine indices of source layers
c
      do ig=1,ngrn
        rrs=rratmos-grndep(ig)
        lygrn(ig)=ly0+1
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn(ig)=ly
            goto 600
          endif
        enddo
600     continue
      enddo
c
      mass=0.d0
      do ly=ly0,1,-1
        if(ly.ge.lyos)then
          dro=(roup(ly)-rolw(ly))/(rrup(ly)-rrlw(ly))
          ro1=rolw(ly)-dro*rrlw(ly) 
          mass=mass+PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &        +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
        else
          alfa=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
          mass=mass+2.d0*PI2/alfa**3
     &        *(roup(ly)*(alfa*rrup(ly)*(alfa*rrup(ly)-2.d0)+2.d0)
     &         -rolw(ly)*(alfa*rrlw(ly)*(alfa*rrlw(ly)-2.d0)+2.d0))
        endif
        cgrup(ly)=dcmplx(BIGG*mass/rrup(ly)**2,0.d0)
      enddo
c
      do ly=1,ly0-1
        cgrlw(ly)=cgrup(ly+1)
      enddo
      cgrlw(ly0)=(0.d0,0.d0)
c
c     make Adam-Williamson condition satisfied in ocean
c
      i=0
      gr0=dreal(cgrup(1))
700   i=i+1
      do ly=lyos,lyob-1
        h=rrup(ly)-rrlw(ly)
        alfa=dreal(cgrup(ly))/vpup(ly)**2
        beta=(dreal(cgrlw(ly))/vplw(ly)**2
     &       -dreal(cgrup(ly))/vpup(ly)**2)/h
        rolw(ly)=roup(ly)*dexp(alfa*h+beta*h**2)
        if(ly.lt.lyob-1)roup(ly+1)=rolw(ly)
      enddo
      mass=0.d0
      do ly=ly0,1,-1
        if(ly.ge.lyos)then
          dro=(roup(ly)-rolw(ly))/(rrup(ly)-rrlw(ly))
          ro1=rolw(ly)-dro*rrlw(ly) 
          mass=mass+PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &        +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
        else
          alfa=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
          mass=mass+2.d0*PI2/alfa**3
     &        *(roup(ly)*(alfa*rrup(ly)*(alfa*rrup(ly)-2.d0)+2.d0)
     &         -rolw(ly)*(alfa*rrlw(ly)*(alfa*rrlw(ly)-2.d0)+2.d0))
        endif
        cgrup(ly)=dcmplx(BIGG*mass/rrup(ly)**2,0.d0)
      enddo
      do ly=1,ly0-1
        cgrlw(ly)=cgrup(ly+1)
      enddo
      cgrlw(ly0)=(0.d0,0.d0)
      if(dabs(dreal(cgrup(1))-gr0).gt.1.0d-08*dreal(cgrup(1))
     &  .and.i.le.10)then
        gr0=dreal(cgrup(1))
        goto 700
      endif
c
      freeairgrd=-dreal(cgrup(lyr))*2.d0/rrup(lyr)
      if(lyr.gt.1)then
        freeairgrd=freeairgrd+4.d0*PI*BIGG*rolw(lyr-1)
      endif
c
c     for a stable atmposheric density grandient
c
      do ly=1,lyos-1
        beta=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
        gr0=dreal(cgrup(ly))
        bvn2=-gr0*(beta+gr0/vpup(ly)**2)
        if(bvn2.lt.0.d0)vpup(ly)=dsqrt(-gr0/beta)
        gr0=dreal(cgrlw(ly))
        bvn2=-gr0*(beta+gr0/vplw(ly)**2)
        if(bvn2.lt.0.d0)vplw(ly)=dsqrt(-gr0/beta)
      enddo
      do ly=2,lyos-1
        vpup(ly)=dmax1(vpup(ly),vplw(ly-1))
        vplw(ly-1)=vpup(ly)
      enddo
c
      do ly=1,ly0
        crrup(ly)=dcmplx(rrup(ly),0.d0)
        croup(ly)=dcmplx(roup(ly),0.d0)
        claup(ly)=dcmplx(roup(ly)*(vpup(ly)**2-2.d0*vsup(ly)**2),0.d0)
        cmuup(ly)=dcmplx(roup(ly)*vsup(ly)**2,0.d0)
        cvpup(ly)=dcmplx(vpup(ly),0.d0)
        cvsup(ly)=dcmplx(vsup(ly),0.d0)
        cgaup(ly)=dcmplx(2.d0*PI2*BIGG*roup(ly),0.d0)
c
        crrlw(ly)=dcmplx(rrlw(ly),0.d0)
        crolw(ly)=dcmplx(rolw(ly),0.d0)
        clalw(ly)=dcmplx(rolw(ly)*(vplw(ly)**2-2.d0*vslw(ly)**2),0.d0)
        cvplw(ly)=dcmplx(vplw(ly),0.d0)
        cvslw(ly)=dcmplx(vslw(ly),0.d0)
        cmulw(ly)=dcmplx(rolw(ly)*vslw(ly)**2,0.d0)
        cgalw(ly)=dcmplx(2.d0*PI2*BIGG*rolw(ly),0.d0)
c
        if(rrup(ly).le.rrlw(ly))then
          cro(ly)=(0.5d0,0.d0)*(croup(ly)+crolw(ly))
        else
          if(ly.ge.lyos)then
            dro=(roup(ly)-rolw(ly))/(rrup(ly)-rrlw(ly))
            ro1=rolw(ly)-dro*rrlw(ly) 
            mass=PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*ro1
     &          *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &          +dro*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &          +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
          else
            alfa=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
            mass=2.d0*PI2/alfa**3
     &          *(roup(ly)*(alfa*rrup(ly)*(alfa*rrup(ly)-2.d0)+2.d0)
     &           -rolw(ly)*(alfa*rrlw(ly)*(alfa*rrlw(ly)-2.d0)+2.d0))
          endif
          cro(ly)=dcmplx(mass/((rrup(ly)**3-rrlw(ly)**3)
     &                        *2.d0*PI2/3.d0),0.d0)
        endif
        cla(ly)=(0.5d0,0.d0)*(claup(ly)+clalw(ly))
        cmu(ly)=(0.5d0,0.d0)*(cmuup(ly)+cmulw(ly))
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
        cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        cgr(ly)=(0.5d0,0.d0)*(cgrup(ly)+cgrlw(ly))
        cga(ly)=(0.5d0,0.d0)*(cgaup(ly)+cgalw(ly))
      enddo
c
      do ly=1,lyob-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=(1.d0,0.d0)
        cypnorm(3,ly)=cga(ly)
        cypnorm(4,ly)=cga(ly)
      enddo
      do ly=lyob,lycm-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(5,ly)=cga(ly)
        cypnorm(6,ly)=cga(ly)
      enddo
      do ly=lycm,lycc-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=(1.d0,0.d0)
        cypnorm(3,ly)=cga(ly)
        cypnorm(4,ly)=cga(ly)
      enddo
      do ly=lycc,ly0
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(5,ly)=cga(ly)
        cypnorm(6,ly)=cga(ly)
      enddo
c
	write(*,'(9a)')'  No','      R(km)','   Vp(km/s)',
     &    '   Vs(km/s)',' Ro(g/cm^3)','      Qp',
     &    '      Qs',' g(m/s^2)'
c
      rrlw(ly0)=0.d0
c
	do ly=1,ly0
	  write(*,1001)ly,rrup(ly)/1.d3,vpup(ly)/1.d3,vsup(ly)/1.d3,
     &               roup(ly)/1.d3,qpup(ly),qsup(ly),dreal(cgrup(ly))
        j=0
        do ig=1,ngrn
          if(lygrn(ig).eq.ly)then
            j=j+1
            write(*,'(a3,$)')' S '
          endif
        enddo
        if(j.eq.0)write(*,'(a3,$)')'   '
        if(lyr.eq.ly)then
          write(*,'(a3)')' R '
        else
          write(*,'(a3)')'   '
        endif
	  write(*,1002)rrlw(ly)/1.d3,vplw(ly)/1.d3,vslw(ly)/1.d3,
     &             rolw(ly)/1.d3,qplw(ly),qslw(ly),dreal(cgrlw(ly))
      enddo
c
1001	format(i4,f11.4,3f11.4,2f8.1,f9.4,$)
1002  format(f15.4,3f11.4,2f8.1,f9.4)
      return
	end
