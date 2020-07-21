      subroutine qpgetinp(unit)
      use qpalloc
      implicit none
      integer unit
c
c     work space
c
      integer*4 i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo,sdfsel,ierr
      integer*4 ipath,isurf
      real*8 twindow,twinout,suppress,munit,sfe,sfn,sfz
      real*8 strike,dip,rake,depdif,dswap(14)
      character*80 grndir,outfile,fswap
c
c     uniform receiver depth
c     ======================
c
      lyadd=1
c
      call skipdoc(unit)
      read(unit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skipdoc(unit)
      read(unit,*)twindow,dtout
      ntcut=1+idnint(twindow/dtout)
      nt=2
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      dt=twindow/dble(nt-1)
      nf=nt/2
      df=1.d0/(dble(nt)*dt)
c
      call skipdoc(unit)
      read(unit,*)fcut
      nfcut=min0(nf,1+idnint(fcut/df))
      fcut=dble(nfcut-1)*df
      call skipdoc(unit)
      read(unit,*)slwmax
      if(slwmax.le.0.d0)then
        stop ' Error in qpgetinp: bad selection of max. slowness!'
      else
        slwmax=slwmax/KM2M
      endif
c
      call skipdoc(unit)
      read(unit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        suppress=dexp(-1.d0)
      endif
      fi=dlog(suppress)*df/PI2
      call skipdoc(unit)
      read(unit,*)ipath,minpath,maxpath
      minpath=minpath*KM2M
      maxpath=maxpath*KM2M
      call skipdoc(unit)
      read(unit,*)rearth,isurf
      rearth=rearth*KM2M
      freesurf=isurf.eq.1
      if(ipath.ne.1)then
        ipatha=0
        ipathb=0
        minpath=0.d0
        maxpath=rearth
      else
        if(minpath.lt.dpr.or.minpath.ge.dmax1(rearth,maxpath))then
          ipatha=0
          minpath=0.d0
        else
          ipatha=1
          lyadd=lyadd+1
        endif
        if(maxpath.ge.rearth.or.maxpath.le.minpath)then
          ipathb=0
          maxpath=rearth
        else
          ipathb=1
          lyadd=lyadd+1
        endif
      endif
c
c     cutoffs of spectra
c     ==================
c
      call skipdoc(unit)
      read(unit,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop ' Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      call skipdoc(unit)
      read(unit,*)i,j,ldegmin,ldegcut
      selpsv=i.eq.1
      selsh=j.eq.1
      if(.not.(selpsv.or.selsh))then
        stop ' Error in qpgetinp: none of PSV and SH is selected!'
      endif
      if(ldegcut.lt.ldegmin)ldegcut=ldegmin
      ldegmin=min0(max0(1+ndmax,ldegmin),ldegcut)
      ldegmax=ldegcut+ndmax+1
c
      allocate(plm(0:ldegmax,0:2),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: plm not allocated!'
c
      allocate(ul0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ul0 not allocated!'
      allocate(vl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: vl0 not allocated!'
      allocate(wl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: wl0 not allocated!'
      allocate(el0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: el0 not allocated!'
      allocate(fl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: fl0 not allocated!'
      allocate(gl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gl0 not allocated!'
c
      allocate(pl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: pl0 not allocated!'
      allocate(ql0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ql0 not allocated!'
c
      allocate(urlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: urlm not allocated!'
      allocate(utlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: utlm not allocated!'
      allocate(uplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: uplm not allocated!'
c
      allocate(grlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: grlm not allocated!'
      allocate(gtlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gtlm not allocated!'
      allocate(gplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gplm not allocated!'
c
      allocate(errlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: errlm not allocated!'
      allocate(ertlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ertlm not allocated!'
      allocate(erplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: erplm not allocated!'
      allocate(etrlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etrlm not allocated!'
      allocate(ett0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ett0lm not allocated!'
      allocate(ettalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ettalm not allocated!'
      allocate(ettblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ettblm not allocated!'
      allocate(etp0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etp0lm not allocated!'
      allocate(etpalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etpalm not allocated!'
      allocate(etpblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etpblm not allocated!'
      allocate(eprlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eprlm not allocated!'
      allocate(ept0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ept0lm not allocated!'
      allocate(eptalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eptalm not allocated!'
      allocate(eptblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eptblm not allocated!'
      allocate(epp0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: epp0lm not allocated!'
      allocate(eppalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eppalm not allocated!'
      allocate(eppblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eppblm not allocated!'
c
      allocate(lyupp(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyupp not allocated!'
      allocate(lyups(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyups not allocated!'
      allocate(lyupt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyupt not allocated!'
      allocate(lylwp(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylwp not allocated!'
      allocate(lylws(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylws not allocated!'
      allocate(lylwt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylwt not allocated!'
c
c     Green's function files
c     ======================
c
      call skipdoc(unit)
      read(unit,*)ngrn,rr0,grndir
      if(ngrn.le.0)then
        stop ' bad number of source depths!'
      endif
      rr0=rr0*KM2M
      lyadd=lyadd+ngrn
c
      allocate(grndep(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: grndep not allocated!'
      allocate(grnsel(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: grnsel not allocated!'
      allocate(lygrn(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lygrn not allocated!'
c
      allocate(specfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: specfile not allocated!'
      allocate(uspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: uspecfile not allocated!'
      allocate(vspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vspecfile not allocated!'
      allocate(wspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: wspecfile not allocated!'
      allocate(especfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: especfile not allocated!'
      allocate(fspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: fspecfile not allocated!'
      allocate(gspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: gspecfile not allocated!'
      allocate(pspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: pspecfile not allocated!'
      allocate(qspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qspecfile not allocated!'
c
      allocate(isg1(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: isg1 not allocated!'
      allocate(isg2(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: isg2 not allocated!'
      allocate(nsg(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: nsg not allocated!'
c
      do ig=1,ngrn
        call skipdoc(unit)
        read(unit,*)grndep(ig),specfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' bad Green function selection!'
        endif
        grndep(ig)=grndep(ig)*KM2M
      enddo
c
c     sort green function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            dswap(1)=grndep(i)
            fswap=specfile(i)
            iswap=grnsel(i)
c
            grndep(i)=grndep(j)
            specfile(i)=specfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            specfile(j)=fswap
            grnsel(j)=iswap
          endif
        enddo
      enddo
c
      do flen=80,1,-1
        if(grndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do ig=1,ngrn
        uspecfile(ig)=grndir(1:flen)//'U_'//specfile(ig)
        vspecfile(ig)=grndir(1:flen)//'V_'//specfile(ig)
        wspecfile(ig)=grndir(1:flen)//'W_'//specfile(ig)
        especfile(ig)=grndir(1:flen)//'E_'//specfile(ig)
        fspecfile(ig)=grndir(1:flen)//'F_'//specfile(ig)
        gspecfile(ig)=grndir(1:flen)//'G_'//specfile(ig)
        pspecfile(ig)=grndir(1:flen)//'P_'//specfile(ig)
        qspecfile(ig)=grndir(1:flen)//'Q_'//specfile(ig)
      enddo
c
c     multi-event source parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)ns,sdfsel
c
      allocate(sfr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfr not allocated!'
      allocate(sft(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sft not allocated!'
      allocate(sfp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfp not allocated!'
c
      allocate(mrr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrr not allocated!'
      allocate(mtt(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtt not allocated!'
      allocate(mpp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mpp not allocated!'
      allocate(mrt(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrt not allocated!'
      allocate(mpr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mpr not allocated!'
      allocate(mtp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtp not allocated!'
      allocate(lats(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lats not allocated!'
      allocate(lons(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lons not allocated!'
      allocate(deps(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: deps not allocated!'
      allocate(togs(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: togs not allocated!'
      allocate(trss(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: trss not allocated!'
c
      if(sdfsel.eq.1)then
        do is=1,ns
          call skipdoc(unit)
c
c         the six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
c
          read(unit,*)munit,mrr(is),mtt(is),mpp(is),
     &                      mrt(is),mpr(is),mtp(is),
     &                      lats(is),lons(is),deps(is),togs(is),trss(is)
          mtt(is)=mtt(is)*munit
          mpp(is)=mpp(is)*munit
          mrr(is)=mrr(is)*munit
          mtp(is)=mtp(is)*munit
          mpr(is)=mpr(is)*munit
          mrt(is)=mrt(is)*munit
          deps(is)=deps(is)*KM2M
          sfp(is)=0.d0
          sft(is)=0.d0
          sfr(is)=0.d0
        enddo
      else if(sdfsel.eq.2)then
        do is=1,ns
          call skipdoc(unit)
          read(unit,*)munit,strike,dip,rake,
     &                      lats(is),lons(is),deps(is),togs(is),trss(is)
          call moments(munit,strike,dip,rake,
     &                 mtt(is),mpp(is),mrr(is),
     &                 mtp(is),mpr(is),mrt(is))
          deps(is)=deps(is)*KM2M
          sfp(is)=0.d0
          sft(is)=0.d0
          sfr(is)=0.d0
        enddo
      else if(sdfsel.eq.3)then
        do is=1,ns
          call skipdoc(unit)
          read(unit,*)munit,sfe,sfn,sfz,
     &                      lats(is),lons(is),deps(is),togs(is),trss(is)
          sfp(is)=sfe*munit
          sft(is)=-sfn*munit
          sfr(is)=sfz*munit
          deps(is)=deps(is)*KM2M
          mtt(is)=0.d0
          mpp(is)=0.d0
          mrr(is)=0.d0
          mtp(is)=0.d0
          mpr(is)=0.d0
          mrt(is)=0.d0
        enddo
      else
        stop ' bad selection for the source data format!'
      endif
c
      togsmin=togs(1)
      do is=1,ns
        togsmin=dmin1(togsmin,togs(is))
      enddo
      do is=1,ns
        togs(is)=togs(is)-togsmin
      enddo
c
c     sort sub-events by depth
c
      do i=1,ns
        do j=i+1,ns
          if(deps(j).lt.deps(i))then
            dswap(1)=lats(i)
            dswap(2)=lons(i)
            dswap(3)=deps(i)
            dswap(4)=mtt(i)
            dswap(5)=mpp(i)
            dswap(6)=mrr(i)
            dswap(7)=mtp(i)
            dswap(8)=mpr(i)
            dswap(9)=mrt(i)
            dswap(10)=sft(i)
            dswap(11)=sfp(i)
            dswap(12)=sfr(i)
            dswap(13)=togs(i)
            dswap(14)=trss(i)
c
            lats(i)=lats(j)
            lons(i)=lons(j)
            deps(i)=deps(j)
            mtt(i)=mtt(j)
            mpp(i)=mpp(j)
            mrr(i)=mrr(j)
            mtp(i)=mtp(j)
            mpr(i)=mpr(j)
            mrt(i)=mrt(j)
            sft(i)=sft(j)
            sfp(i)=sfp(j)
            sfr(i)=sfr(j)
            togs(i)=togs(j)
            trss(i)=trss(j)
c
            lats(j)=dswap(1)
            lons(j)=dswap(2)
            deps(j)=dswap(3)
            mtt(j)=dswap(4)
            mpp(j)=dswap(5)
            mrr(j)=dswap(6)
            mtp(j)=dswap(7)
            mpr(j)=dswap(8)
            mrt(j)=dswap(9)
            sft(j)=dswap(10)
            sfp(j)=dswap(11)
            sfr(j)=dswap(12)
            togs(j)=dswap(13)
            trss(j)=dswap(14)
          endif
        enddo
      enddo
c
      isg1(1)=1
      is1=1
      do ig=1,ngrn-1
        depdif=0.5d0*(grndep(ig)+grndep(ig+1))
        isg2(ig)=isg1(ig)-1
        do is=is1,ns
          if(deps(is).lt.depdif)then
            isg2(ig)=is
          endif
        enddo
        isg1(ig+1)=isg2(ig)+1
        is1=isg1(ig+1)
      enddo
      isg2(ngrn)=ns
c
      do ig=1,ngrn
        nsg(ig)=max0(0,1+isg2(ig)-isg1(ig))
      enddo
c
c     receiver parameters
c     ===================
c
      call skipdoc(unit)
      read(unit,*)(icmp(i),i=1,11)
      call skipdoc(unit)
      read(unit,*)outfile
      call skipdoc(unit)
      read(unit,*)twinout
      ntcutout=1+idnint(twinout/dtout)
      call skipdoc(unit)
      read(unit,*)nbpf,f1corner,f2corner
      call skipdoc(unit)
      read(unit,*)slwlwcut,slwupcut
      slwlwcut=slwlwcut/KM2M
      slwupcut=slwupcut/KM2M
      if(slwupcut.gt.slwmax.or.slwlwcut.ge.slwupcut)then
        slwlwcut=0.d0
        slwupcut=slwmax
      endif
      call skipdoc(unit)
      read(unit,*)nr
c
      allocate(latr(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: latr not allocated!'
      allocate(lonr(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lonr not allocated!'
      allocate(tred(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: tred not allocated!'
      allocate(rname(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: rname not allocated!'
c
      do ir=1,nr
        call skipdoc(unit)
        read(unit,*)latr(ir),lonr(ir),rname(ir),tred(ir)
      enddo
c
      do flen=80,1,-1
        if(outfile(flen:flen).ne.' ')goto 300
      enddo
300   continue
      dispout(1)=outfile(1:flen)//'_disp_e.dat'
      dispout(2)=outfile(1:flen)//'_disp_n.dat'
      dispout(3)=outfile(1:flen)//'_disp_z.dat'
      veloout(1)=outfile(1:flen)//'_velo_e.dat'
      veloout(2)=outfile(1:flen)//'_velo_n.dat'
      veloout(3)=outfile(1:flen)//'_velo_z.dat'
      acceout(1)=outfile(1:flen)//'_acce_e.dat'
      acceout(2)=outfile(1:flen)//'_acce_n.dat'
      acceout(3)=outfile(1:flen)//'_acce_z.dat'
c
      rotaout(1)=outfile(1:flen)//'_rota_e.dat'
      rotaout(2)=outfile(1:flen)//'_rota_n.dat'
      rotaout(3)=outfile(1:flen)//'_rota_z.dat'
      rotarateout(1)=outfile(1:flen)//'_rota_rate_e.dat'
      rotarateout(2)=outfile(1:flen)//'_rota_rate_n.dat'
      rotarateout(3)=outfile(1:flen)//'_rota_rate_z.dat'
c
      strainout(1)=outfile(1:flen)//'_strain_ee.dat'
      strainout(2)=outfile(1:flen)//'_strain_en.dat'
      strainout(3)=outfile(1:flen)//'_strain_ez.dat'
      strainout(4)=outfile(1:flen)//'_strain_nn.dat'
      strainout(5)=outfile(1:flen)//'_strain_nz.dat'
      strainout(6)=outfile(1:flen)//'_strain_zz.dat'
      strainrateout(1)=outfile(1:flen)//'_strain_rate_ee.dat'
      strainrateout(2)=outfile(1:flen)//'_strain_rate_en.dat'
      strainrateout(3)=outfile(1:flen)//'_strain_rate_ez.dat'
      strainrateout(4)=outfile(1:flen)//'_strain_rate_nn.dat'
      strainrateout(5)=outfile(1:flen)//'_strain_rate_nz.dat'
      strainrateout(6)=outfile(1:flen)//'_strain_rate_zz.dat'
c
      stressout(1)=outfile(1:flen)//'_stress_ee.dat'
      stressout(2)=outfile(1:flen)//'_stress_en.dat'
      stressout(3)=outfile(1:flen)//'_stress_ez.dat'
      stressout(4)=outfile(1:flen)//'_stress_nn.dat'
      stressout(5)=outfile(1:flen)//'_stress_nz.dat'
      stressout(6)=outfile(1:flen)//'_stress_zz.dat'
      stressrateout(1)=outfile(1:flen)//'_stress_rate_ee.dat'
      stressrateout(2)=outfile(1:flen)//'_stress_rate_en.dat'
      stressrateout(3)=outfile(1:flen)//'_stress_rate_ez.dat'
      stressrateout(4)=outfile(1:flen)//'_stress_rate_nn.dat'
      stressrateout(5)=outfile(1:flen)//'_stress_rate_nz.dat'
      stressrateout(6)=outfile(1:flen)//'_stress_rate_zz.dat'
c
      gravout(1)=outfile(1:flen)//'_gravitation_e.dat'
      gravout(2)=outfile(1:flen)//'_gravitation_n.dat'
      gravout(3)=outfile(1:flen)//'_gravitation_z.dat'
c
      grmout=outfile(1:flen)//'_gravimeter.dat'
c
c     multilayered model parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)l,i
      if(i.eq.1)then
        dispersion=.true.
      else
        dispersion=.false.
      endif
c
      l0=l+1
      allocate(dp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0 not allocated!'
      allocate(vp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0 not allocated!'
      allocate(vs0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0 not allocated!'
      allocate(ro0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0 not allocated!'
      allocate(qp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0 not allocated!'
      allocate(qs0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0 not allocated!'
c
      allocate(dp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0up not allocated!'
      allocate(vp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0up not allocated!'
      allocate(vs0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0up not allocated!'
      allocate(ro0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0up not allocated!'
      allocate(qp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0up not allocated!'
      allocate(qs0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0up not allocated!'
c
      allocate(dp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0lw not allocated!'
      allocate(vp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0lw not allocated!'
      allocate(vs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0lw not allocated!'
      allocate(ro0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0lw not allocated!'
      allocate(qp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0lw not allocated!'
      allocate(qs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0lw not allocated!'
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=rearth
      do i=1,l
        call skipdoc(unit)
        read(unit,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(vp0(i).le.0.d0.or.vs0(i).lt.0.d0.or.ro0(i).le.0.d0.or.
     &     qp0(i).le.0.d0.or.qs0(i).lt.0.d0)then
          stop ' Error in qpgetinp: bad seismic parameter!'
        endif
c
        if(i.eq.1)then
          if(dp0(1).gt.0.d0)then
            stop ' Error in qpgetinp: bad start depth!'
          endif
          depatmos=-KM2M*dp0(1)
          rratmos=rearth+depatmos
        endif
c
        dp0(i)=KM2M*dp0(i)+depatmos
        vp0(i)=KM2M*vp0(i)
        vs0(i)=KM2M*vs0(i)
        ro0(i)=KM2M*ro0(i)
        if(vs0(i).gt.0.d0)qsmin=dmin1(qsmin,qs0(i))
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop ' Error in qpgetinp: bad layering of earth model!'
          endif
        endif
      enddo
c
      if(dp0(l).gt.rratmos)then
        stop ' Error in qpgetinp: earth radius larger than pre-defined!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        dp0(l)=rratmos
        vp0(l)=vp0(l-1)
        vs0(l)=vs0(l-1)
        ro0(l)=ro0(l-1)
        qp0(l)=qp0(l-1)
        qs0(l)=qs0(l-1)
      endif
c
      dpr=dpr+depatmos
      if(dpr.lt.0.d0.or.dpr.gt.dp0(l))then
        stop ' Error in qpgetinp: receiver too shallow or too deep!'
      endif
      do i=1,ngrn
        grndep(i)=grndep(i)+depatmos
        if(grndep(i).lt.0.d0.or.grndep(i).ge.dp0(l))then
          stop ' Error in qpgetinp: source too shallow or too deep!'
        endif
      enddo
      do is=1,ns
        deps(is)=deps(is)+depatmos
      enddo
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          vp0up(l0)=vp0(i-1)
          vs0up(l0)=vs0(i-1)
          ro0up(l0)=ro0(i-1)
          qp0up(l0)=qp0(i-1)
          qs0up(l0)=qs0(i-1)
c
          dp0lw(l0)=dp0(i)
          vp0lw(l0)=vp0(i)
          vs0lw(l0)=vs0(i)
          ro0lw(l0)=ro0(i)
          qp0lw(l0)=qp0(i)
          qs0lw(l0)=qs0(i)
        endif
      enddo
c
c     end of inputs
c     =============
c
      return
      end
