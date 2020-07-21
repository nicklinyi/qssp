      subroutine qpgrnspec(ig)
      use qpalloc
      implicit none
      integer*4 ig
c
      integer*4 i,nn,il,istp,ly,lf,ldeg,ldeg0,ldegf,ldegup,ierr
      real*8 f,ksp,xlw,xup,dll,omi,expo,rpath,slwcut
      real*8 fac,fl,depst1,depst2,dys2,rr0a,a,b,rrs,x
      real*8 kcut1(4),kcut2(4)
      complex*16 ca,cb,cag,cs1,cs2,cs3,cs4,ct1,ct2,cll1
      complex*16 cmur,ys3d,yt1d,cg1,cg5
      complex*16 ypsv(6,4),ypsvg(6,4),ysh(2,2)
      logical*2 forruku
c
      real*8 fsimpson
c
      real*8 expos
      complex*16 c2,c3,c4
      data expos/48.d0/
      data c2,c3,c4/(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      if(ig.eq.igfirst)then
        allocate(disk(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: disk not allocated!'
c
        allocate(ksmallp(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ksmallp not allocated!'
        allocate(ksmalls(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ksmalls not allocated!'
        allocate(ksmallt(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ksmallt not allocated!'
        allocate(xp(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: xp not allocated!'
        allocate(xs(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: xs not allocated!'
        allocate(xt(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: xt not allocated!'
        allocate(kp(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: kp not allocated!'
        allocate(ks(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ks not allocated!'
        allocate(kt(lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: kt not allocated!'
        allocate(cps(6,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: cps not allocated!'
        allocate(cpt(2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: cpt not allocated!'
c
        allocate(mat2x2up(2,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2up not allocated!'
        allocate(mat2x2lw(2,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2lw not allocated!'
        allocate(mat2x2inv(2,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2inv not allocated!'
c
        allocate(mas3x3up(3,3,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas3x3up not allocated!'
        allocate(mas3x3lw(3,3,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas3x3lw not allocated!'
        allocate(mas3x3inv(3,3,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas3x3inv not allocated!'
c
        allocate(mas4x4up(4,4,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4up not allocated!'
        allocate(mas4x4lw(4,4,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4lw not allocated!'
        allocate(mas4x4inv(4,4,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4inv not allocated!'
c
        allocate(mas6x6up(6,6,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6up not allocated!'
        allocate(mas6x6lw(6,6,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6lw not allocated!'
        allocate(mas6x6inv(6,6,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6inv not allocated!'
c
        allocate(zjup(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjup not allocated!'
        allocate(zjlw(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjlw not allocated!'
        allocate(zhup(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhup not allocated!'
        allocate(zhlw(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhlw not allocated!'
        allocate(wj(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: wj not allocated!'
        allocate(wh(0:ldegmax,lymax,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: wh not allocated!'
        allocate(zh12p(0:ldegmax,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zh12p not allocated!'
        allocate(zh12sv(0:ldegmax,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zh12sv not allocated!'
        allocate(zh12sh(0:ldegmax,2,lymax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zh12sh not allocated!'
c
        allocate(zjupg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjupg not allocated!'
        allocate(zjlwg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjlwg not allocated!'
        allocate(zhupg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhupg not allocated!'
        allocate(zhlwg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhlwg not allocated!'
        allocate(wjg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: wjg not allocated!'
        allocate(whg(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: whg not allocated!'
      endif
c
c     Initiation
c
      rrs=rrup(lys)
c
      if(rr0.gt.0.d0)then
        x=dcos(rr0/rrs)
        disk(0)=1.d0
        disk(1)=(3.d0+x)/4.d0
        do ldeg=2,ldegmax
          disk(ldeg)=(dble(2*ldeg-1)*x*disk(ldeg-1)
     &               +dble(4-ldeg)*disk(ldeg-2))
     &              /dble(ldeg+3)
        enddo
        do ldeg=0,ldegmax
          disk(ldeg)=disk(ldeg)*dble(2*ldeg+1)/(4.d0*PI*rrs**2)
        enddo
      else
        do ldeg=0,ldegmax
          disk(ldeg)=dble(2*ldeg+1)/(4.d0*PI*rrs**2)
        enddo
      endif
c
      open(21,file=uspecfile(ig),form='unformatted',status='unknown')
      open(22,file=vspecfile(ig),form='unformatted',status='unknown')
      open(23,file=wspecfile(ig),form='unformatted',status='unknown')
      open(24,file=especfile(ig),form='unformatted',status='unknown')
      open(25,file=fspecfile(ig),form='unformatted',status='unknown')
      open(26,file=gspecfile(ig),form='unformatted',status='unknown')
      open(27,file=pspecfile(ig),form='unformatted',status='unknown')
      open(28,file=qspecfile(ig),form='unformatted',status='unknown')
c
      ldeg0=10+ndmax
      if(vsup(lys).gt.0.d0)then
        slwcut=1.d0/vsup(lys)
      else
        slwcut=1.d0/vpup(lys)
      endif
      if(vsup(lyr).gt.0.d0)then
        slwcut=dmin1(slwcut,1.d0/vsup(lyr))
      else
        slwcut=dmin1(slwcut,1.d0/vpup(lyr))
      endif
      if(slwmax.ge.slwcut)then
        do ldeg0=10+ndmax,ldegmin
          dll=dsqrt(dble(ldeg0)*dabs(dble(ldeg0-1)))
          expo=0.d0
          do ly=min0(lys,lyr),max0(lys,lyr)-1
            expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          enddo
          if(expo.gt.expos)goto 10
        enddo
10      continue
      endif
c
      lylwa=0
      if(ipatha.eq.1)then
        rpath=rratmos-minpath
        do ly=max0(lys,lyr),ly0
          if(rrup(ly).ge.rpath.and.rrlw(ly).lt.rpath)then
            lylwa=ly
            goto 20
          endif
        enddo
20      continue
      endif
c
      lylwb=ly0+1
      if(ipathb.eq.1)then
        rpath=rratmos-maxpath
        do ly=max0(lys,lyr,lylwa+1),ly0
          if(rrup(ly).ge.rpath.and.rrlw(ly).lt.rpath)then
            lylwb=ly
            goto 30
          endif
        enddo
30      continue
      endif
c
      omi=PI2*dble(nfcut-1)*df
      slwcut=slwmax
c
      if(lylwa.gt.lys+1)then
        fac=(1.d0+0.5d0*(rrup(lylwa)/rearth)**2)*rrup(lylwa)/rearth
        if(vsup(lylwa).gt.0.d0)then
          slwcut=dmin1(slwcut,fac/vsup(lylwa))
        else
          slwcut=dmin1(slwcut,fac/vpup(lylwa))
        endif
      endif
c
      if(idint(rearth*omi*slwcut).gt.ldegcut)then
        print *,' Warning from qpgrnspec: ldegcut may be too small!'
      endif
c
      ldegup=min0(ldegcut,max0(ldeg0,idint(rearth*omi*slwcut)))
c
      if(fgr.gt.0.d0.or.ldeggr.gt.0)then
        do ly=1,ly0
          do ldeg=0,ldegup
            nruku(ldeg,ly)=10
          enddo
        enddo
      endif
c
      write(*,*)' '
      write(*,'(a,i3,a,f7.2,a)')' ... calculate Green functions for ',
     &        ig,'. source at depth ',(grndep(ig)-depatmos)/KM2M,' km'
	write(*,'(a,i5)')'   max. harmonic degree: L_max = ',ldegup
c
      write(21)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(22)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(23)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(24)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(25)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(26)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(27)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(28)nt,ntcut,dt,nf,nfcut,df,ldegup
c
      if(.not.nogravity)then
        do ly=1,ly0
          do ldeg=0,ldegup
            nruku(ldeg,ly)=0
          enddo
        enddo
      endif
c
      do istp=1,6
        do ldeg=0,ldegmax
          ul0(ldeg,istp)=(0.d0,0.d0)
          vl0(ldeg,istp)=(0.d0,0.d0)
          wl0(ldeg,istp)=(0.d0,0.d0)
          el0(ldeg,istp)=(0.d0,0.d0)
          fl0(ldeg,istp)=(0.d0,0.d0)
          gl0(ldeg,istp)=(0.d0,0.d0)
          pl0(ldeg,istp)=(0.d0,0.d0)
          ql0(ldeg,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      do lf=1,nfcut
        f=dble(lf-1)*df
        omi=PI2*f
        comi=dcmplx(PI2*f,PI2*fi)
        comi2=comi*comi
        call qpqmodel(f)
c
        do ly=1,ly0
          kp(ly)=comi/cvp(ly)
          if(vsup(ly).gt.0.d0)then
            ks(ly)=comi/cvs(ly)
            kt(ly)=comi/cvs(ly)
          else
            ks(ly)=(0.d0,0.d0)
            kt(ly)=(0.d0,0.d0)
          endif
          if(rrlw(ly).gt.0.d0)then
            xp(ly)=cdabs(kp(ly))*rrlw(ly)
            xs(ly)=cdabs(ks(ly))*rrlw(ly)
            xt(ly)=cdabs(kt(ly))*rrlw(ly)
          else
            xp(ly)=cdabs(kp(ly))*rrup(ly)
            xs(ly)=cdabs(ks(ly))*rrup(ly)
            xt(ly)=cdabs(kt(ly))*rrup(ly)
          endif
          ksmallp(ly)=.true.
          ksmalls(ly)=.true.
          ksmallt(ly)=.true.
        enddo
c
        ldegf=min0(ldegup,max0(ldeg0,idint(rearth*omi*slwcut)))
        setzh12a=.false.
        setzh12t=.false.
c
        do ldeg=0,ldegf
          if(nogravity.or.lylwa.gt.0.or.lylwb.le.ly0.or.
     &      .not.freesurf)then
            forruku=.false.
          else
            fac=dsqrt((f/fgr)**2+(dble(ldeg)/dble(ldeggr))**2)
            forruku=fac.lt.1.d0+FLTAPER
          endif
          dll=dsqrt(dble(ldeg)*dabs(dble(ldeg-1)))
c
c         determine degree dependent starting layer number
c         of sh solution
c
          lyupt(ldeg)=ly0+1
          lylwt(ldeg)=0
          if(lys.ge.lyob.and.lys.le.min0(lycm-1,ly0).and.
     &       lyr.ge.lyob.and.lyr.le.min0(lycm-1,ly0))then
            lyupt(ldeg)=lyob
            expo=0.d0
            do ly=min0(lys,lyr)-1,lyob,-1
              ksp=omi/vsup(ly)
              if(ksp.le.0.d0)then
                expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
              else if(dll.gt.ksp*rrlw(ly))then
                xlw=ksp*rrlw(ly)
                xup=dmin1(dll,ksp*rrup(ly))
                expo=expo+fsimpson(xlw,xup,dll)
              endif
              if(expo.gt.expos)then
                lyupt(ldeg)=ly+1
                goto 101
              endif
            enddo
101         continue
c
            lylwt(ldeg)=min0(lycm,ly0)
            expo=0.d0
            do ly=max0(lys,lyr)+1,min0(lycm,ly0)-1
              ksp=omi/vsup(ly)
              if(ksp.le.0.d0)then
                expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
              else if(dll.gt.ksp*rrlw(ly))then
                xlw=ksp*rrlw(ly)
                xup=dmin1(dll,ksp*rrup(ly))
                expo=expo+fsimpson(xlw,xup,dll)
              endif
              if(expo.gt.expos)then
                lylwt(ldeg)=ly
                goto 102
              endif
            enddo
102         continue
          endif
c
          if(lys.ge.lycc.and.lyr.ge.lycc)then
            lyupt(ldeg)=lycc
            expo=0.d0
            do ly=min0(lys,lyr)-1,lycc,-1
              ksp=omi/vsup(ly)
              if(ksp.le.0.d0)then
                expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
              else if(dll.gt.ksp*rrlw(ly))then
                xlw=ksp*rrlw(ly)
                xup=dmin1(dll,ksp*rrup(ly))
                expo=expo+fsimpson(xlw,xup,dll)
              endif
              if(expo.gt.expos)then
                lyupt(ldeg)=ly+1
                goto 103
              endif
            enddo
103         continue
c
            lylwt(ldeg)=ly0
            expo=0.d0
            do ly=max0(lys,lyr)+1,ly0-1
              ksp=omi/vsup(ly)
              if(ksp.le.0.d0)then
                expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
              else if(dll.gt.ksp*rrlw(ly))then
                xlw=ksp*rrlw(ly)
                xup=dmin1(dll,ksp*rrup(ly))
                expo=expo+fsimpson(xlw,xup,dll)
              endif
              if(expo.gt.expos)then
                lylwt(ldeg)=ly
                goto 104
              endif
            enddo
104         continue
          endif
c
c         determine degree dependent starting layer number
c         of psv solution
c
          lyupp(ldeg)=1
          expo=0.d0
          do ly=min0(lys,lyr)-1,1,-1
            if(vsup(ly).gt.0.d0)then
              ksp=omi/vsup(ly)
            else
              ksp=omi/vpup(ly)
            endif
            if(ksp.le.0.d0)then
              expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
            else if(dll.gt.ksp*rrlw(ly))then
              xlw=ksp*rrlw(ly)
              xup=dmin1(dll,ksp*rrup(ly))
              expo=expo+fsimpson(xlw,xup,dll)
            endif
            if(expo.gt.expos)then
              lyupp(ldeg)=ly
              goto 201
            endif
          enddo
201       continue
c
          lylwp(ldeg)=ly0
          expo=0.d0
          do ly=max0(lys,lyr)+1,ly0-1
            if(vsup(ly).gt.0.d0)then
              ksp=omi/vsup(ly)
            else
              ksp=omi/vpup(ly)
            endif
            if(ksp.le.0.d0)then
              expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
            else if(dll.gt.ksp*rrlw(ly))then
              xlw=ksp*rrlw(ly)
              xup=dmin1(dll,ksp*rrup(ly))
              expo=expo+fsimpson(xlw,xup,dll)
            endif
            if(expo.gt.expos)then
              lylwp(ldeg)=ly+1
              goto 202
            endif
          enddo
202       continue
c
          lyups(ldeg)=lyob
          expo=0.d0
          do ly=min0(lys,lyr)-1,lyob,-1
            if(vsup(ly).gt.0.d0)then
              ksp=omi/vsup(ly)
            else
              ksp=omi/vpup(ly)
            endif
            if(ksp.le.0.d0)then
              expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
            else if(dll.gt.ksp*rrlw(ly))then
              xlw=ksp*rrlw(ly)
              xup=dmin1(dll,ksp*rrup(ly))
              expo=expo+fsimpson(xlw,xup,dll)
            endif
            if(expo.gt.expos)then
              lyups(ldeg)=ly
              goto 301
            endif
          enddo
301       continue
c
          lylws(ldeg)=ly0
          expo=0.d0
          do ly=max0(lys,lyr)+1,ly0-1
            if(vsup(ly).gt.0.d0)then
              ksp=omi/vsup(ly)
            else
              ksp=omi/vpup(ly)
            endif
            if(ksp.le.0.d0)then
              expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
            else if(dll.gt.ksp*rrlw(ly))then
              xlw=ksp*rrlw(ly)
              xup=dmin1(dll,ksp*rrup(ly))
              expo=expo+fsimpson(xlw,xup,dll)
            endif
            if(expo.gt.expos)then
              lylws(ldeg)=ly+1
              goto 302
            endif
          enddo
302       continue
        enddo
c
        if(selpsv.or.lys.lt.lyob.or.lys.ge.lycm.and.lys.lt.lycc)then
          ly=max0(lylwa,min0(lylwb,lylwp(0),lylws(0)))
          depst1=(rratmos-depatmos-rrup(ly))/KM2M
          ly=max0(lylwa,min0(lylwb,lylwp(ldegf),lylws(ldegf)))
          depst2=(rratmos-depatmos-rrup(ly))/KM2M
        else
          ly=max0(lylwa,min0(lylwb,lylwt(1)))
          depst1=(rratmos-depatmos-rrup(ly))/KM2M
          ly=max0(lylwa,min0(lylwb,lylwt(ldegf)))
          depst2=(rratmos-depatmos-rrup(ly))/KM2M
        endif
c
c       determine layer dependent max. harmonic degree
c       of sh solution
c
        if(lys.ge.lyob.and.lys.le.min0(lycm-1,ly0))then
          do ly=lyob,min0(lycm-1,ly0)
            ldegsh(ly)=1
            do ldeg=1,ldegf
              if(ly.ge.lyupt(ldeg).and.ly.le.lylwt(ldeg))then
                ldegsh(ly)=ldeg
              endif
            enddo
          enddo
        else if(lys.ge.lycc)then
          do ly=lycc,ly0
            ldegsh(ly)=1
            do ldeg=1,ldegf
              if(ly.ge.lyupt(ldeg).and.ly.le.lylwt(ldeg))then
                ldegsh(ly)=ldeg
              endif
            enddo
          enddo
        endif
c
c       determine layer dependent max. harmonic degree
c       of psv solution
c
        do ly=1,ly0
          ldegpsv(ly)=0
          do ldeg=0,ldegf
            if(ly.ge.min0(lyupp(ldeg),lyups(ldeg)).and.
     &         ly.le.max0(lylwp(ldeg),lylws(ldeg)))then
              ldegpsv(ly)=ldeg
            endif
          enddo
        enddo
c
        if(lylwa.gt.0)then
          ldegf=min0(ldegf,max0(ldegsh(lylwa),ldegpsv(lylwa)))
        endif
c
        do ldeg=0,ldegf
c
          if(.not.selpsv.or.ldeg.gt.ldegcut.or.
     &       max0(lylwp(ldeg),lylws(ldeg)).lt.max0(lylwa,lyr,lys).or.
     &       min0(lyupp(ldeg),lyups(ldeg)).gt.min0(lyr,lys))then
            do istp=1,4
              do i=1,6
                ypsv(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else if(nogravity)then
            call qppsvkern(f,ldeg,ypsv)
          else
            fac=dsqrt((f/fgr)**2+(dble(ldeg)/dble(ldeggr))**2)
            if(fac.le.1.d0)then
              call qppsvkerng(f,ldeg,ypsv)
            else if(fac.ge.1.d0+FLTAPER)then
              call qppsvkern(f,ldeg,ypsv)
            else
              call qppsvkerng(f,ldeg,ypsvg)
              call qppsvkern(f,ldeg,ypsv)
              ca=dcmplx(dsin(0.5d0*PI*(fac-1.d0)/FLTAPER)**2,0.d0)
              cb=(1.d0,0.d0)-ca
              do istp=1,4
                do i=1,6
                  ypsv(i,istp)=ca*ypsv(i,istp)+cb*ypsvg(i,istp)
                enddo
              enddo
            endif
          endif
c
          if(.not.selsh.or.ldeg.gt.ldegcut.or.ldeg.eq.0.or.
     &       lys.lt.lyob.or.lys.ge.lycm.and.lys.lt.lycc.or.
     &       lylwt(ldeg).lt.max0(lylwa,lyr,lys).or.
     &       lyupt(ldeg).gt.min0(lyr,lys))then
            do istp=1,2
              do i=1,2
                ysh(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else
            call qpshkern(f,ldeg,ysh)
          endif
c
          if(lyr.gt.1)then
            cg1=cgalw(lyr-1)
          else
            cg1=(0.d0,0.d0)
          endif
          cg5=dcmplx(dble(ldeg+1)/rrup(lyr),0.d0)
c
c         1. Vertical single force (F3=1)
c
          if(lys.lt.lyob)then
            ul0(ldeg,1)=(0.d0,0.d0)
            vl0(ldeg,1)=(0.d0,0.d0)
            wl0(ldeg,1)=(0.d0,0.d0)
            el0(ldeg,1)=(0.d0,0.d0)
            fl0(ldeg,1)=(0.d0,0.d0)
            gl0(ldeg,1)=(0.d0,0.d0)
            pl0(ldeg,1)=(0.d0,0.d0)
            ql0(ldeg,1)=(0.d0,0.d0)
          else
            cs2=dcmplx(-disk(ldeg),0.d0)
            ul0(ldeg,1)=cs2*ypsv(1,2)
            vl0(ldeg,1)=cs2*ypsv(3,2)
            wl0(ldeg,1)=(0.d0,0.d0)
            el0(ldeg,1)=cs2*ypsv(2,2)
            fl0(ldeg,1)=cs2*ypsv(4,2)
            gl0(ldeg,1)=(0.d0,0.d0)
            pl0(ldeg,1)=cs2*ypsv(5,2)
            ql0(ldeg,1)=cs2*(cg1*ypsv(1,2)-cg5*ypsv(5,2)+ypsv(6,2))
          endif
c
c         2. Explosion (M11=M22=M33=1)
c
          cs1=dcmplx( disk(ldeg)/roup(lys),0.d0)/cvpup(lys)**2
          cs2=dcmplx(-disk(ldeg)*4.d0/rrup(lys),0.d0)
     &       *(cvsup(lys)/cvpup(lys))**2
          cs4=dcmplx( disk(ldeg)*2.d0/rrup(lys),0.d0)
     &       *(cvsup(lys)/cvpup(lys))**2
          ul0(ldeg,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          vl0(ldeg,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          wl0(ldeg,2)=(0.d0,0.d0)
          el0(ldeg,2)=cs1*ypsv(2,1)+cs2*ypsv(2,2)+cs4*ypsv(2,4)
          fl0(ldeg,2)=cs1*ypsv(4,1)+cs2*ypsv(4,2)+cs4*ypsv(4,4)
          gl0(ldeg,2)=(0.d0,0.d0)
          pl0(ldeg,2)=cs1*ypsv(5,1)+cs2*ypsv(5,2)+cs4*ypsv(5,4)
          ql0(ldeg,2)=cs1*(cg1*ypsv(1,1)-cg5*ypsv(5,1)+ypsv(6,1))
     &               +cs2*(cg1*ypsv(1,2)-cg5*ypsv(5,2)+ypsv(6,2))
     &               +cs4*(cg1*ypsv(1,4)-cg5*ypsv(5,4)+ypsv(6,4))
c
c         3. CLVD (M33=1,M11=M22=-0.5)
c
          if(lys.lt.lyob)then
            ul0(ldeg,3)=(0.d0,0.d0)
            vl0(ldeg,3)=(0.d0,0.d0)
            wl0(ldeg,3)=(0.d0,0.d0)
            el0(ldeg,3)=(0.d0,0.d0)
            fl0(ldeg,3)=(0.d0,0.d0)
            gl0(ldeg,3)=(0.d0,0.d0)
            pl0(ldeg,3)=(0.d0,0.d0)
            ql0(ldeg,3)=(0.d0,0.d0)
          else
            cs1=dcmplx( disk(ldeg)/roup(lys),0.d0)/cvpup(lys)**2
            cs2=dcmplx( disk(ldeg)/rrup(lys),0.d0)
     &         *((3.d0,0.d0)-(4.d0,0.d0)*(cvsup(lys)/cvpup(lys))**2)
            cs4=-(0.5d0,0.d0)*cs2
            ul0(ldeg,3)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
            vl0(ldeg,3)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
            wl0(ldeg,3)=(0.d0,0.d0)
            el0(ldeg,3)=cs1*ypsv(2,1)+cs2*ypsv(2,2)+cs4*ypsv(2,4)
            fl0(ldeg,3)=cs1*ypsv(4,1)+cs2*ypsv(4,2)+cs4*ypsv(4,4)
            gl0(ldeg,3)=(0.d0,0.d0)
            pl0(ldeg,3)=cs1*ypsv(5,1)+cs2*ypsv(5,2)+cs4*ypsv(5,4)
            ql0(ldeg,3)=cs1*(cg1*ypsv(1,1)-cg5*ypsv(5,1)+ypsv(6,1))
     &                 +cs2*(cg1*ypsv(1,2)-cg5*ypsv(5,2)+ypsv(6,2))
     &                 +cs4*(cg1*ypsv(1,4)-cg5*ypsv(5,4)+ypsv(6,4))
          endif
c
c         4. Horizontal single force (F1=1)
c
          if(ldeg.lt.1.or.lys.lt.lyob)then
            ul0(ldeg,4)=(0.d0,0.d0)
            vl0(ldeg,4)=(0.d0,0.d0)
            wl0(ldeg,4)=(0.d0,0.d0)
            el0(ldeg,4)=(0.d0,0.d0)
            fl0(ldeg,4)=(0.d0,0.d0)
            gl0(ldeg,4)=(0.d0,0.d0)
            pl0(ldeg,4)=(0.d0,0.d0)
            ql0(ldeg,4)=(0.d0,0.d0)
          else
            cs4=dcmplx(-disk(ldeg)/(dble(ldeg)*dble(ldeg+1)),0.d0)
            ct2=cs4
            ul0(ldeg,4)=cs4*ypsv(1,4)
            vl0(ldeg,4)=cs4*ypsv(3,4)
            wl0(ldeg,4)=ct2*ysh(1,2)
            el0(ldeg,4)=cs4*ypsv(2,4)
            fl0(ldeg,4)=cs4*ypsv(4,4)
            gl0(ldeg,4)=ct2*ysh(2,2)
            pl0(ldeg,4)=cs4*ypsv(5,4)
            ql0(ldeg,4)=cs4*(cg1*ypsv(1,4)-cg5*ypsv(5,4)+ypsv(6,4))
          endif
c
c         5. Dip-slip (M13=M31=1)
c
          if(ldeg.lt.1.or.lys.lt.lyob)then
            ul0(ldeg,5)=(0.d0,0.d0)
            vl0(ldeg,5)=(0.d0,0.d0)
            wl0(ldeg,5)=(0.d0,0.d0)
            el0(ldeg,5)=(0.d0,0.d0)
            fl0(ldeg,5)=(0.d0,0.d0)
            gl0(ldeg,5)=(0.d0,0.d0)
            pl0(ldeg,5)=(0.d0,0.d0)
            ql0(ldeg,5)=(0.d0,0.d0)
          else
            ct1=dcmplx(disk(ldeg)/(dble(ldeg)*dble(ldeg+1)
     &          *roup(lys)),0.d0)/cvsup(lys)**2
            cs3=ct1
            ul0(ldeg,5)=cs3*ypsv(1,3)
            vl0(ldeg,5)=cs3*ypsv(3,3)
            wl0(ldeg,5)=ct1*ysh(1,1)
            el0(ldeg,5)=cs3*ypsv(2,3)
            fl0(ldeg,5)=cs3*ypsv(4,3)
            gl0(ldeg,5)=ct1*ysh(2,1)
            pl0(ldeg,5)=cs3*ypsv(5,3)
            ql0(ldeg,5)=cs3*(cg1*ypsv(1,3)-cg5*ypsv(5,3)+ypsv(6,3))
          endif
c
c         6. Strike-slip (M12=M21=1)
c
          if(ldeg.lt.2.or.lys.lt.lyob)then
            ul0(ldeg,6)=(0.d0,0.d0)
            vl0(ldeg,6)=(0.d0,0.d0)
            wl0(ldeg,6)=(0.d0,0.d0)
            el0(ldeg,6)=(0.d0,0.d0)
            fl0(ldeg,6)=(0.d0,0.d0)
            gl0(ldeg,6)=(0.d0,0.d0)
            pl0(ldeg,6)=(0.d0,0.d0)
            ql0(ldeg,6)=(0.d0,0.d0)
          else
            ct2=dcmplx(disk(ldeg)
     &                /(dble(ldeg)*dble(ldeg+1)*rrup(lys)),0.d0)
            cs4=-ct2
            ul0(ldeg,6)=cs4*ypsv(1,4)
            vl0(ldeg,6)=cs4*ypsv(3,4)
            wl0(ldeg,6)=ct2*ysh(1,2)
            el0(ldeg,6)=cs4*ypsv(2,4)
            fl0(ldeg,6)=cs4*ypsv(4,4)
            gl0(ldeg,6)=ct2*ysh(2,2)
            pl0(ldeg,6)=cs4*ypsv(5,4)
            ql0(ldeg,6)=cs4*(cg1*ypsv(1,4)-cg5*ypsv(5,4)+ypsv(6,4))
          endif
        enddo
c
        write(*,'(i6,a,f12.4,a,i5,a,2(f7.2,a))')lf,'.',1.0d+03*f,
     &       ' mHz: cut-off degree = ',ldegf,
     &       ', start depth = ',depst1,' - ',depst2,' km'
c
        write(21)ldegf
        write(22)ldegf
        write(23)ldegf
        write(24)ldegf
        write(25)ldegf
        write(26)ldegf
        write(27)ldegf
        write(28)ldegf
        write(21)((ul0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! Ypsv1
        write(22)((vl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! Ypsv3
        write(23)((wl0(ldeg,istp),ldeg=0,ldegf),istp=4,6)               ! Ysh1 (or Y7)
        write(24)((el0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! Ypsv2
        write(25)((fl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! Ypsv4
        write(26)((gl0(ldeg,istp),ldeg=0,ldegf),istp=4,6)               ! Ysh2 (or Y8)
        write(27)((pl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! Ypsv5
        write(28)((ql0(ldeg,istp),ldeg=0,ldegf),istp=1,6)               ! dYpsv5/dr
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
c
      if(ig.eq.iglast)then
        deallocate(disk)
        deallocate(xp,xs,xt,kp,ks,kt,cps,cpt)
        deallocate(mat2x2up,mat2x2lw,mat2x2inv,
     &             mas3x3up,mas3x3lw,mas3x3inv,
     &             mas4x4up,mas4x4lw,mas4x4inv,
     &             mas6x6up,mas6x6lw,mas6x6inv)
        deallocate(zjup,zjlw,zhup,zhlw,wj,wh,
     &             zh12p,zh12sv,zh12sh,
     &             zjupg,zjlwg,zhupg,zhlwg,wjg,whg)
      endif
      return
      end
