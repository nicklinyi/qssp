      subroutine qpwvint(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,k,id,is,ir,ig,nd,nt0,nf0,ntcut0,nfcut0,ishift
      integer*4 lf,lf1,istp,ldeg,ldegf,ldeg0
      integer*4 istat,ldegup,ldeglw,ldegneed
      real*8 depsarc,dis0,anorm,slwcut
      real*8 f,dt0,df0,rn,re,azi,bazi,bazi0
      complex*16 cll1,cp0,cp1,cp2,wavelet,muer,lamr,ksir,rrr,srt,srp
      complex*16 cfac,ca,cb,dur,dut,dup,eii,dgr,dgt,dgp
      complex*16 urdr,urdt,urdp,utdr,utdt,utdp,updr,updt,updp
      complex*16 rot(3,3),rtz(3,3),enz(3,3),swp(3,3)
      logical*2 fullwave
c
      complex*16 c1,c2,c3
      data c1,c2,c3/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      allocate(tap(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: tap not allocated!'
      allocate(ldegtap(4,ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ldegtap not allocated!'
      allocate(wvf(nf,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: wvf not allocated!'
c
      allocate(expl(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: expl not allocated!'
      allocate(clvd(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: clvd not allocated!'
      allocate(ss12(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ss12 not allocated!'
      allocate(ss11(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ss11 not allocated!'
      allocate(ds31(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ds31 not allocated!'
      allocate(ds23(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ds23 not allocated!'
c
      allocate(ue(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ue not allocated!'
      allocate(un(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: un not allocated!'
      allocate(uz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: uz not allocated!'
c
      allocate(ge(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ge not allocated!'
      allocate(gn(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: gn not allocated!'
      allocate(gz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: gz not allocated!'
      allocate(gm(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: gm not allocated!'
c
      allocate(roe(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: roe not allocated!'
      allocate(ron(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ron not allocated!'
      allocate(roz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: roz not allocated!'
c
      allocate(uee(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: uee not allocated!'
      allocate(uen(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: uen not allocated!'
      allocate(uez(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: uez not allocated!'
      allocate(unn(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: unn not allocated!'
      allocate(unz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: unz not allocated!'
      allocate(uzz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: uzz not allocated!'
c
      allocate(see(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: see not allocated!'
      allocate(sen(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: sen not allocated!'
      allocate(sez(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: sez not allocated!'
      allocate(snn(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: snn not allocated!'
      allocate(snz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: snz not allocated!'
      allocate(szz(nf,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: szz not allocated!'
c
      allocate(idr(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: idr not allocated!'
      allocate(dis(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: dis not allocated!'
      allocate(ssa(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ssa not allocated!'
      allocate(csa(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: csa not allocated!'
      allocate(ssb(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ssb not allocated!'
      allocate(ssd(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ssd not allocated!'
      allocate(csd(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: csd not allocated!'
      allocate(csb(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: csb not allocated!'
      allocate(ssf(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ssf not allocated!'
      allocate(ss2a(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: ss2a not allocated!'
      allocate(cs2a(ns,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: cs2a not allocated!'
c
      allocate(sf1(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: sf1 not allocated!'
      allocate(sf2(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: sf2 not allocated!'
      allocate(sf3(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpwvint: sf3 not allocated!'
c
      do ldeg=0,ldegmax
        tap(ldeg)=1.d0
      enddo
c
      ldeg0=10+ndmax
c
      rrr=crrup(lyr)
c
      do is=1,ns
        depsarc=PI2/dble(ldegmax)+dmax1(minpath,deps(is))/rearth
        dis0=5.d0*depsarc
c
        do ir=1,nr
          call disazi(1.d0,lats(is),lons(is),
     &                     latr(ir),lonr(ir),rn,re)
c
          dis(is,ir)=dsqrt(rn**2+re**2)
c
c         determine order of differential transform
c
          if(dis(is,ir).le.dis0)then
            idr(is,ir)=0
          else
            idr(is,ir)=min0(ndmax,idint(dlog(dis(is,ir)/dis0)))
          endif
c
          ssd(is,ir)=dcmplx(dsin(dis(is,ir)),0.d0)
          csd(is,ir)=dcmplx(dcos(dis(is,ir)),0.d0)
          ssf(is,ir)=dcmplx(2.d0*dsin(0.5d0*dis(is,ir))**2,0.d0)
c
c         azi = receiver azimuth (from south to east)
c
          if(dsqrt(re*re+rn*rn).gt.0.d0)then
            azi=datan2(re,-rn)
          else
c
c           assume southern receiver in case of 0 distance
c
            azi=datan2(0.d0,1.d0)
          endif
          ssa(is,ir)=dcmplx(dsin(azi),0.d0)
          csa(is,ir)=dcmplx(dcos(azi),0.d0)
          ss2a(is,ir)=dcmplx(dsin(2.d0*azi),0.d0)
          cs2a(is,ir)=dcmplx(dcos(2.d0*azi),0.d0)
        enddo
      enddo
c
      do ir=1,nr
        do is=1,ns
          call disazi(1.d0,latr(ir),lonr(ir),
     &                     lats(is),lons(is),rn,re)
c
c         opposite vector receiver-source
c         azmuth from south to east
c
          if(dsqrt(re*re+rn*rn).gt.0.d0)then
            bazi=datan2(-re,rn)
          else
            bazi=datan2(0.d0,-1.d0)
          endif
          ssb(is,ir)=dcmplx(dsin(bazi),0.d0)
          csb(is,ir)=dcmplx(dcos(bazi),0.d0)
        enddo
      enddo
c
      do is=1,ns
        expl(is)=dcmplx((mtt(is)+mpp(is)+mrr(is))/3.d0,0.d0)
        clvd(is)=dcmplx(mrr(is),0.d0)-expl(is)
        ss12(is)=dcmplx(mtp(is),0.d0)
        ss11(is)=dcmplx((mtt(is)-mpp(is))/2.d0,0.d0)
        ds31(is)=dcmplx(mrt(is),0.d0)
        ds23(is)=dcmplx(mpr(is),0.d0)
        sf1(is)=dcmplx(sft(is),0.d0)
        sf2(is)=dcmplx(sfp(is),0.d0)
        sf3(is)=dcmplx(sfr(is),0.d0)
      enddo
c
c     initiation
c
      do is=1,ns
        call swavelet(trss(is),df,fi,nf,wvf(1,is))
      enddo
c
      do lf=1,nf
        do ir=1,nr
          ue(lf,ir)=(0.d0,0.d0)
          un(lf,ir)=(0.d0,0.d0)
          uz(lf,ir)=(0.d0,0.d0)
c
          ge(lf,ir)=(0.d0,0.d0)
          gn(lf,ir)=(0.d0,0.d0)
          gz(lf,ir)=(0.d0,0.d0)
          gm(lf,ir)=(0.d0,0.d0)
c
          roe(lf,ir)=(0.d0,0.d0)
          ron(lf,ir)=(0.d0,0.d0)
          roz(lf,ir)=(0.d0,0.d0)
c
          uee(lf,ir)=(0.d0,0.d0)
          uen(lf,ir)=(0.d0,0.d0)
          uez(lf,ir)=(0.d0,0.d0)
          unn(lf,ir)=(0.d0,0.d0)
          unz(lf,ir)=(0.d0,0.d0)
          uzz(lf,ir)=(0.d0,0.d0)
c
          see(lf,ir)=(0.d0,0.d0)
          sen(lf,ir)=(0.d0,0.d0)
          sez(lf,ir)=(0.d0,0.d0)
          snn(lf,ir)=(0.d0,0.d0)
          snz(lf,ir)=(0.d0,0.d0)
          szz(lf,ir)=(0.d0,0.d0)
        enddo
      enddo
c
      do ig=1,ngrn
        if(nsg(ig).le.0)goto 500
        write(*,'(a)')' '
        write(*,'(a,i4,a,f5.1,a)')' processing ',1+isg2(ig)-isg1(ig),
     &    ' point source(s) at depth ',(grndep(ig)-depatmos)/KM2M,' km'
        write(*,'(a)')' open Green function data base: '
     &              //specfile(ig)(1:40)
        write(*,'(a)')' ... please wait ...'
c
        open(21,file=uspecfile(ig),form='unformatted',status='old')
        open(22,file=vspecfile(ig),form='unformatted',status='old')
        open(23,file=wspecfile(ig),form='unformatted',status='old')
        open(24,file=especfile(ig),form='unformatted',status='old')
        open(25,file=fspecfile(ig),form='unformatted',status='old')
        open(26,file=gspecfile(ig),form='unformatted',status='old')
        open(27,file=pspecfile(ig),form='unformatted',status='old')
        open(28,file=qspecfile(ig),form='unformatted',status='old')
c
        read(21)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(22)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(23)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(24)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(25)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(26)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(27)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(28)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
c
        if(ntcut0.ne.ntcut.or.dabs(dt0-dt).gt.1.0d-06*dt.or.
     &     nfcut0.lt.nfcut.or.dabs(df0-df).gt.1.0d-06*df)then
          print *,' Error in qpwvint: t/f sampling'
     &          //' inconsistent with Green functions!'
          write(*,'(a)')'               '//'  ntcut             dt'
     &                                   //'  nfcut             df'
          write(*,'(a,2(i7,f16.8))')' Current input:',
     &                               ntcut,dt,nfcut,df
          write(*,'(a,2(i7,f16.8))')'     Data base:',
     &                               ntcut0,dt0,nfcut0,df0
          stop
        endif
c
        nd=0
        do is=isg1(ig),isg2(ig)
          do ir=1,nr
            nd=max0(nd,idr(is,ir))
          enddo
        enddo
c
        if(ldegup+nd.gt.ldegmax)then
          write(*,'(a,i6,a,i6)')' Error in qpwvint: '
     &     //'max. harmonic degree required = ',ldegup+nd,
     &     ' > ldegmax defined: ',ldegmax
          stop
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
        lys=lygrn(ig)
        if(vsup(lys).gt.0.d0)then
          slwcut=1.d0/vsup(lys)
        else
          slwcut=1.d0/vpup(lys)
        endif
        fullwave=slwlwcut.le.0.d0.and.slwupcut.ge.dmin1(slwmax,slwcut)
c
        do lf=1,nfcut
          f=dble(lf-1)*df
c
          call qpqmodel(f)
          muer=cmuup(lyr)
          lamr=claup(lyr)
          ksir=lamr+c2*muer
c
          read(21)ldegf
          read(22)ldegf
          read(23)ldegf
          read(24)ldegf
          read(25)ldegf
          read(26)ldegf
          read(27)ldegf
          read(28)ldegf
c
          read(21)((ul0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !Y1
          read(22)((vl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !Y3
          read(23)((wl0(ldeg,istp),ldeg=0,ldegf),istp=4,6)              !Y7
          read(24)((el0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !Y2
          read(25)((fl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !Y4
          read(26)((gl0(ldeg,istp),ldeg=0,ldegf),istp=4,6)              !Y8
          read(27)((pl0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !Y5
          read(28)((ql0(ldeg,istp),ldeg=0,ldegf),istp=1,6)              !dY5/dr derived from Y1, Y5 and Y6
c
          ldegneed=0
          ldeglw=ldegmax
          do is=isg1(ig),isg2(ig)
            do ir=1,nr
              if(fullwave)then
                ldegtap(4,is,ir)=max0(0,ldegf-idr(is,ir)-1)
                ldegtap(3,is,ir)=ldegtap(4,is,ir)*4/5
                ldegtap(2,is,ir)=0
                ldegtap(1,is,ir)=0
              else
                ldegtap(4,is,ir)=min0(ldegf-idr(is,ir)-1,
     &                        ldeg0+idnint(rearth*PI2*f*slwupcut))
                ldegtap(3,is,ir)=ldegtap(4,is,ir)*4/5
                ldegtap(4,is,ir)=min0(ldegf-idr(is,ir)-1,
     &                                10+ldegtap(4,is,ir))
                ldegtap(2,is,ir)=min0(ldegtap(3,is,ir),
     &                           idnint(rearth*PI2*f*slwlwcut*1.25d0))
                ldegtap(1,is,ir)=ldegtap(2,is,ir)*4/5
              endif
              ldegneed=max0(ldegneed,ldegtap(4,is,ir))
              ldeglw=max0(0,min0(ldeglw,ldegtap(1,is,ir)-idr(is,ir)))
            enddo
          enddo
          ldegneed=min0(ldegneed+nd,ldegf)        
c
c         use differential filter to suppress spatial aliasing
c
c         vertical single force, explosion and clvd sources
c         ur,ut,up(=0)
c         err,ert,erp(=0),etr,ett,etp(=0),epr(=0),ept(=0),epp
c         gr,gt,gp(=0)
c
          do istp=1,3
            do ldeg=0,ldegf-1
              cll1=dcmplx(dble(ldeg*(ldeg+1)),0.d0)
              urlm(ldeg,istp,0)=ul0(ldeg,istp)                          !m=0,  1
              utlm(ldeg,istp,0)=-vl0(ldeg,istp)                         !m=1,  1
c
              grlm(ldeg,istp,0)=ql0(ldeg,istp)                          !m=0,  1
              gtlm(ldeg,istp,0)=-pl0(ldeg,istp)/rrr                     !m=1,  1
c
              errlm(ldeg,istp,0)=(el0(ldeg,istp)                        !m=0,  1
     &                     +(lamr/rrr)*(-c2*ul0(ldeg,istp)
     &                     +cll1*vl0(ldeg,istp)))/ksir
              ertlm(ldeg,istp,0)=(-ul0(ldeg,istp)                       !m=1,  1
     &                            +vl0(ldeg,istp))/rrr
c
              if(cdabs(muer).gt.0.d0)then
                srt=-fl0(ldeg,istp)
                etrlm(ldeg,istp,0)=srt/muer-ertlm(ldeg,istp,0)          !m=1,  1
              else
                etrlm(ldeg,istp,0)=ertlm(ldeg,istp,0)                   !m=1,  1
              endif
              ett0lm(ldeg,istp,0)=(ul0(ldeg,istp)
     &                        -cll1*vl0(ldeg,istp))/rrr                 !m=0,  1
              ettalm(ldeg,istp,0)=vl0(ldeg,istp)/rrr                    !m=1,  cos(t)/sin(t)
c
              epp0lm(ldeg,istp,0)= ul0(ldeg,istp)/rrr                   !m=0,  1
              eppalm(ldeg,istp,0)=-vl0(ldeg,istp)/rrr                   !m=1,  cos(t)/sin(t)
            enddo
          enddo
c
          do istp=4,6
            do ldeg=0,istp/6
              urlm(ldeg,istp,0)=(0.d0,0.d0)
              utlm(ldeg,istp,0)=(0.d0,0.d0)
              uplm(ldeg,istp,0)=(0.d0,0.d0)
              grlm(ldeg,istp,0)=(0.d0,0.d0)
              gtlm(ldeg,istp,0)=(0.d0,0.d0)
              gplm(ldeg,istp,0)=(0.d0,0.d0)
              errlm(ldeg,istp,0)=(0.d0,0.d0)
              ertlm(ldeg,istp,0)=(0.d0,0.d0)
              erplm(ldeg,istp,0)=(0.d0,0.d0)
              etrlm(ldeg,istp,0)=(0.d0,0.d0)
              ett0lm(ldeg,istp,0)=(0.d0,0.d0)
              ettalm(ldeg,istp,0)=(0.d0,0.d0)
              ettblm(ldeg,istp,0)=(0.d0,0.d0)
              etp0lm(ldeg,istp,0)=(0.d0,0.d0)
              etpalm(ldeg,istp,0)=(0.d0,0.d0)
              etpblm(ldeg,istp,0)=(0.d0,0.d0)
              eprlm(ldeg,istp,0)=(0.d0,0.d0)
              ept0lm(ldeg,istp,0)=(0.d0,0.d0)
              eptalm(ldeg,istp,0)=(0.d0,0.d0)
              eptblm(ldeg,istp,0)=(0.d0,0.d0)
              epp0lm(ldeg,istp,0)=(0.d0,0.d0)
              eppalm(ldeg,istp,0)=(0.d0,0.d0)
              eppblm(ldeg,istp,0)=(0.d0,0.d0)
            enddo
          enddo
c
c         horizontal single force and dip-slip
c
          do istp=4,5
            do ldeg=1,ldegf-1
              cll1=dcmplx(dble(ldeg*(ldeg+1)),0.d0)
              ca=dcmplx(dble(ldeg-1)**2/dble(2*ldeg-1),0.d0)
              cb=dcmplx(dble(ldeg+2)**2/dble(2*ldeg+3),0.d0)
              urlm(ldeg,istp,0)=ul0(ldeg,istp)                          !m=1,  1
              utlm(ldeg,istp,0)=ca*vl0(ldeg-1,istp)-cb*vl0(ldeg+1,istp) !m=1,  1/sin(t)
     &                       +wl0(ldeg,istp)
              uplm(ldeg,istp,0)=cb*wl0(ldeg+1,istp)-ca*wl0(ldeg-1,istp) !m=1,  1/sin(t)
     &                       -vl0(ldeg,istp)
c
              grlm(ldeg,istp,0)=ql0(ldeg,istp)                          !m=1,  1
              gtlm(ldeg,istp,0)=(ca*pl0(ldeg-1,istp)                    !m=1,  1/sin(t)
     &                          -cb*pl0(ldeg+1,istp))/rrr
              gplm(ldeg,istp,0)=-pl0(ldeg,istp)/rrr                     !m=1,  1/sin(t)
c
              errlm(ldeg,istp,0)=(el0(ldeg,istp)                        !m=1,  1
     &                       +(lamr/rrr)*(-c2*ul0(ldeg,istp)
     &                       +dble(ldeg*(ldeg+1))*vl0(ldeg,istp)))/ksir
c             ertlm(ldeg,istp,0) s. below
              erplm(ldeg,istp,0)=(-ul0(ldeg,istp)                       !m=1,  1/sin(t)
     &                          -uplm(ldeg,istp,0))/rrr
c
              if(cdabs(muer).gt.0.d0)then
                srt=(ca*fl0(ldeg-1,istp)
     &              -cb*fl0(ldeg+1,istp)+gl0(ldeg,istp))/muer
                etrlm(ldeg,istp,0)=srt+utlm(ldeg,istp,0)/rrr            !m=1,  1/sin(t)
     &             -(ca*ul0(ldeg-1,istp)-cb*ul0(ldeg+1,istp))/rrr
                ertlm(ldeg,istp,0)=srt-etrlm(ldeg,istp,0)               !m=1,  1/sin(t)
              else
                etrlm(ldeg,istp,0)=(ca*(ul0(ldeg-1,istp)                 !m=1,  1/sin(t)
     &                                 -vl0(ldeg-1,istp))
     &                             -cb*(ul0(ldeg+1,istp)
     &                                 -vl0(ldeg+1,istp)))/rrr
                ertlm(ldeg,istp,0)=etrlm(ldeg,istp,0)                   !m=1,  1/sin(t)
              endif
c
              ett0lm(ldeg,istp,0)=(ul0(ldeg,istp)                       !m=1,  1
     &                  +(c1-cll1)*vl0(ldeg,istp))/rrr
              ettalm(ldeg,istp,0)= vl0(ldeg,istp)/rrr                   !m=2,  cos(t)/sin(t)
              ettblm(ldeg,istp,0)=-wl0(ldeg,istp)/rrr                   !m=2,  1/sin(t)
c
              etp0lm(ldeg,istp,0)=-wl0(ldeg,istp)/rrr                   !m=1,  1
              etpalm(ldeg,istp,0)= vl0(ldeg,istp)/rrr                   !m=2,  1/sin(t)
              etpblm(ldeg,istp,0)=-wl0(ldeg,istp)/rrr                   !m=2,  cos(t)/sin(t)
c
              if(cdabs(muer).gt.0.d0)then
                srp=(cb*gl0(ldeg+1,istp)-ca*gl0(ldeg-1,istp)
     &             -fl0(ldeg,istp))/muer
                eprlm(ldeg,istp,0)=srp-erplm(ldeg,istp,0)               !m=1,  1/sin(t)
              else
                eprlm(ldeg,istp,0)=erplm(ldeg,istp,0)                   !m=1,  1/sin(t)
              endif
c
              ept0lm(ldeg,istp,0)=-(c1-cll1)*wl0(ldeg,istp)/rrr         !m=1,  1
              eptalm(ldeg,istp,0)=-wl0(ldeg,istp)/rrr                   !m=2,  cos(t)/sin(t)
              eptblm(ldeg,istp,0)= vl0(ldeg,istp)/rrr                   !m=2,  1/sin(t)
c
              epp0lm(ldeg,istp,0)=(ul0(ldeg,istp)-vl0(ldeg,istp))/rrr   !m=1,  1
              eppalm(ldeg,istp,0)=-vl0(ldeg,istp)/rrr                   !m=2,  cos(t)/sin(t)
              eppblm(ldeg,istp,0)= wl0(ldeg,istp)/rrr                   !m=2,  1/sin(t)
            enddo
          enddo
c
c         strike-slip
c
          do ldeg=2,ldegf-1
            cll1=dcmplx(dble(ldeg*(ldeg+1)),0.d0)
            ca=dcmplx(dble((ldeg-1)*(ldeg-2))/dble(2*ldeg-1),0.d0)
            cb=dcmplx(dble((ldeg+2)*(ldeg+3))/dble(2*ldeg+3),0.d0)
            urlm(ldeg,6,0)=ul0(ldeg,6)                                  !m=2,  1
            utlm(ldeg,6,0)= ca*vl0(ldeg-1,6)-cb*vl0(ldeg+1,6)           !m=2,  1/sin(t)
     &                     -c2*wl0(ldeg,6)
            uplm(ldeg,6,0)=-ca*wl0(ldeg-1,6)+cb*wl0(ldeg+1,6)           !m=2,  1/sin(t)
     &                     +c2*vl0(ldeg,6)
c
            grlm(ldeg,6,0)=ql0(ldeg,6)                                  !m=2,  1
            gtlm(ldeg,6,0)=(ca*pl0(ldeg-1,6)-cb*pl0(ldeg+1,6))/rrr      !m=2,  1/sin(t)
            gplm(ldeg,6,0)=c2*pl0(ldeg,6)/rrr                           !m=2,  1/sin(t)
c
            errlm(ldeg,6,0)=(el0(ldeg,6)                                !m=2,  1
     &                     +(lamr/rrr)*(-c2*ul0(ldeg,6)
     &                     +cll1*vl0(ldeg,6)))/ksir
c           ertlm(ldeg,6,0) s. below
            erplm(ldeg,6,0)=(c2*ul0(ldeg,6)-uplm(ldeg,6,0))/rrr         !m=2,  1/sin(t)
c
            if(cdabs(muer).gt.0.d0)then
              etrlm(ldeg,6,0)=ca*(fl0(ldeg-1,6)/muer                    !m=2,  1/sin(t)
     &                          +(vl0(ldeg-1,6)-ul0(ldeg-1,6))/rrr)
     &                       -cb*(fl0(ldeg+1,6)/muer
     &                          +(vl0(ldeg+1,6)-ul0(ldeg+1,6))/rrr)
     &                     -c2*(wl0(ldeg,6)/rrr+gl0(ldeg,6)/muer)
              srt=ca*fl0(ldeg-1,6)-cb*fl0(ldeg+1,6)-c2*gl0(ldeg,6)
              ertlm(ldeg,6,0)=srt/muer-etrlm(ldeg,6,0)                  !m=2,  1/sin(t)
            else
              etrlm(ldeg,6,0)=ca*(ul0(ldeg-1,6)-vl0(ldeg-1,6))/rrr      !m=2,  1/sin(t)
     &                       -cb*(ul0(ldeg+1,6)-vl0(ldeg+1,6))/rrr
              ertlm(ldeg,6,0)=etrlm(ldeg,6,0)                           !m=2  1/sin(t)
            endif
c
            ett0lm(ldeg,6,0)=(ul0(ldeg,6)                               !m=2,  1
     &                 -cll1*vl0(ldeg,6))/rrr
            ettalm(ldeg,6,0)=-utlm(ldeg,6,0)/rrr                        !m=2,  cos(t)/sin^2(t)
            ettblm(ldeg,6,0)=c2*uplm(ldeg,6,0)/rrr                      !m=2,  1/sin^2(t)
c
            etp0lm(ldeg,6,0)=c2*utlm(ldeg,6,0) /rrr                     !m=2,  1/sin^2(t)
            etpalm(ldeg,6,0)=-uplm(ldeg,6,0)/rrr                        !m=2,  cos(t)/sin^2(t)
c
            if(cdabs(muer).gt.0.d0)then
              srp=-ca*gl0(ldeg-1,6)+cb*gl0(ldeg+1,6)+c2*fl0(ldeg,6)
              eprlm(ldeg,6,0)=srp/muer-erplm(ldeg,6,0)                  !m=2,  1/sin(t)
            else
              eprlm(ldeg,6,0)=erplm(ldeg,6,0)                           !m=2,  1/sin(t)
            endif
c
            ept0lm(ldeg,6,0)=-dcmplx(dble(ldeg+2),0.d0)                 !m=2,  cos(t)/sin^2(t)
     &                      *uplm(ldeg,6,0)/rrr
            eptalm(ldeg,6,0)=dcmplx(dble(ldeg-2),0.d0)                  !m=2,  1/sin^2(t)
     &                      *uplm(ldeg-1,6,0)/rrr
c
            epp0lm(ldeg,6,0)=ul0(ldeg,6)/rrr                            !m=2,  1
            eppalm(ldeg,6,0)=utlm(ldeg,6,0)/rrr                         !m=2,  cos(t)/sin^2(t)
            eppblm(ldeg,6,0)=-c2*uplm(ldeg,6,0)/rrr                     !m=2,  1/sin^2(t)
          enddo
c
c         use differential transform to suppress spatial aliasing
c
          do id=1,nd
c
c           m = 0
c
            do istp=1,3
              urlm(0,istp,id)=urlm(0,istp,id-1)
     &                       -urlm(1,istp,id-1)/c3
              grlm(0,istp,id)=grlm(0,istp,id-1)
     &                       -grlm(1,istp,id-1)/c3
              errlm(0,istp,id)=errlm(0,istp,id-1)
     &                        -errlm(1,istp,id-1)/c3
              ett0lm(0,istp,id)=ett0lm(0,istp,id-1)
     &                         -ett0lm(1,istp,id-1)/c3
              epp0lm(0,istp,id)=epp0lm(0,istp,id-1)
     &                         -epp0lm(1,istp,id-1)/c3
            enddo
c
            do ldeg=1,ldegf-1-id
              ca=dcmplx(dble(ldeg+1)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg)/dble(2*ldeg-1),0.d0)
              do istp=1,3
                urlm(ldeg,istp,id)=urlm(ldeg,istp,id-1)
     &                         -ca*urlm(ldeg+1,istp,id-1)
     &                         -cb*urlm(ldeg-1,istp,id-1)
                grlm(ldeg,istp,id)=grlm(ldeg,istp,id-1)
     &                         -ca*grlm(ldeg+1,istp,id-1)
     &                         -cb*grlm(ldeg-1,istp,id-1)
                errlm(ldeg,istp,id)=errlm(ldeg,istp,id-1)
     &                         -ca*errlm(ldeg+1,istp,id-1)
     &                         -cb*errlm(ldeg-1,istp,id-1)
                ett0lm(ldeg,istp,id)=ett0lm(ldeg,istp,id-1)
     &                         -ca*ett0lm(ldeg+1,istp,id-1)
     &                         -cb*ett0lm(ldeg-1,istp,id-1)
                epp0lm(ldeg,istp,id)=epp0lm(ldeg,istp,id-1)
     &                         -ca*epp0lm(ldeg+1,istp,id-1)
     &                         -cb*epp0lm(ldeg-1,istp,id-1)
              enddo
            enddo
c
c           m = 1
c
            do istp=1,3
              utlm(0,istp,id)=(0.d0,0.d0)
              gtlm(0,istp,id)=(0.d0,0.d0)
              ertlm(0,istp,id)=(0.d0,0.d0)
              etrlm(0,istp,id)=(0.d0,0.d0)
              ettalm(0,istp,id)=(0.d0,0.d0)
              eppalm(0,istp,id)=(0.d0,0.d0)
            enddo
            do istp=4,5
              urlm(0,istp,id)=(0.d0,0.d0)
              utlm(0,istp,id)=(0.d0,0.d0)
              uplm(0,istp,id)=(0.d0,0.d0)
              grlm(0,istp,id)=(0.d0,0.d0)
              gtlm(0,istp,id)=(0.d0,0.d0)
              gplm(0,istp,id)=(0.d0,0.d0)
              errlm(0,istp,id)=(0.d0,0.d0)
              ertlm(0,istp,id)=(0.d0,0.d0)
              erplm(0,istp,id)=(0.d0,0.d0)
              etrlm(0,istp,id)=(0.d0,0.d0)
              ett0lm(0,istp,id)=(0.d0,0.d0)
              etp0lm(0,istp,id)=(0.d0,0.d0)
              eprlm(0,istp,id)=(0.d0,0.d0)
              ept0lm(0,istp,id)=(0.d0,0.d0)
              epp0lm(0,istp,id)=(0.d0,0.d0)
            enddo
c
            do ldeg=1,ldegf-1-id
              ca=dcmplx(dble(ldeg+2)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-1)/dble(2*ldeg-1),0.d0)
              do istp=1,3
                utlm(ldeg,istp,id)=utlm(ldeg,istp,id-1)
     &                         -ca*utlm(ldeg+1,istp,id-1)
     &                         -cb*utlm(ldeg-1,istp,id-1)
                gtlm(ldeg,istp,id)=gtlm(ldeg,istp,id-1)
     &                         -ca*gtlm(ldeg+1,istp,id-1)
     &                         -cb*gtlm(ldeg-1,istp,id-1)
                ertlm(ldeg,istp,id)=ertlm(ldeg,istp,id-1)
     &                           -ca*ertlm(ldeg+1,istp,id-1)
     &                           -cb*ertlm(ldeg-1,istp,id-1)
                etrlm(ldeg,istp,id)=etrlm(ldeg,istp,id-1)
     &                           -ca*etrlm(ldeg+1,istp,id-1)
     &                           -cb*etrlm(ldeg-1,istp,id-1)
                ettalm(ldeg,istp,id)=ettalm(ldeg,istp,id-1)
     &                           -ca*ettalm(ldeg+1,istp,id-1)
     &                           -cb*ettalm(ldeg-1,istp,id-1)
                eppalm(ldeg,istp,id)=eppalm(ldeg,istp,id-1)
     &                           -ca*eppalm(ldeg+1,istp,id-1)
     &                           -cb*eppalm(ldeg-1,istp,id-1)
              enddo
              do istp=4,5
                urlm(ldeg,istp,id)=urlm(ldeg,istp,id-1)
     &                      -ca*urlm(ldeg+1,istp,id-1)
     &                      -cb*urlm(ldeg-1,istp,id-1)
                utlm(ldeg,istp,id)=utlm(ldeg,istp,id-1)
     &                      -ca*utlm(ldeg+1,istp,id-1)
     &                      -cb*utlm(ldeg-1,istp,id-1)
                uplm(ldeg,istp,id)=uplm(ldeg,istp,id-1)
     &                      -ca*uplm(ldeg+1,istp,id-1)
     &                      -cb*uplm(ldeg-1,istp,id-1)
                grlm(ldeg,istp,id)=grlm(ldeg,istp,id-1)
     &                      -ca*grlm(ldeg+1,istp,id-1)
     &                      -cb*grlm(ldeg-1,istp,id-1)
                gtlm(ldeg,istp,id)=gtlm(ldeg,istp,id-1)
     &                      -ca*gtlm(ldeg+1,istp,id-1)
     &                      -cb*gtlm(ldeg-1,istp,id-1)
                gplm(ldeg,istp,id)=gplm(ldeg,istp,id-1)
     &                      -ca*gplm(ldeg+1,istp,id-1)
     &                      -cb*gplm(ldeg-1,istp,id-1)
                errlm(ldeg,istp,id)=errlm(ldeg,istp,id-1)
     &                       -ca*errlm(ldeg+1,istp,id-1)
     &                       -cb*errlm(ldeg-1,istp,id-1)
                ertlm(ldeg,istp,id)=ertlm(ldeg,istp,id-1)
     &                       -ca*ertlm(ldeg+1,istp,id-1)
     &                       -cb*ertlm(ldeg-1,istp,id-1)
                erplm(ldeg,istp,id)=erplm(ldeg,istp,id-1)
     &                       -ca*erplm(ldeg+1,istp,id-1)
     &                       -cb*erplm(ldeg-1,istp,id-1)
                etrlm(ldeg,istp,id)=etrlm(ldeg,istp,id-1)
     &                       -ca*etrlm(ldeg+1,istp,id-1)
     &                       -cb*etrlm(ldeg-1,istp,id-1)
                ett0lm(ldeg,istp,id)=ett0lm(ldeg,istp,id-1)
     &                        -ca*ett0lm(ldeg+1,istp,id-1)
     &                        -cb*ett0lm(ldeg-1,istp,id-1)
                etp0lm(ldeg,istp,id)=etp0lm(ldeg,istp,id-1)
     &                        -ca*etp0lm(ldeg+1,istp,id-1)
     &                        -cb*etp0lm(ldeg-1,istp,id-1)
                eprlm(ldeg,istp,id)=eprlm(ldeg,istp,id-1)
     &                       -ca*eprlm(ldeg+1,istp,id-1)
     &                       -cb*eprlm(ldeg-1,istp,id-1)
                ept0lm(ldeg,istp,id)=ept0lm(ldeg,istp,id-1)
     &                        -ca*ept0lm(ldeg+1,istp,id-1)
     &                        -cb*ept0lm(ldeg-1,istp,id-1)
                epp0lm(ldeg,istp,id)=epp0lm(ldeg,istp,id-1)
     &                        -ca*epp0lm(ldeg+1,istp,id-1)
     &                        -cb*epp0lm(ldeg-1,istp,id-1)
              enddo
            enddo
c
c           m = 2
c
            do ldeg=0,1
              do istp=4,5
                ettalm(ldeg,istp,id)=(0.d0,0.d0)
                ettblm(ldeg,istp,id)=(0.d0,0.d0)
                etpalm(ldeg,istp,id)=(0.d0,0.d0)
                etpblm(ldeg,istp,id)=(0.d0,0.d0)
                eptalm(ldeg,istp,id)=(0.d0,0.d0)
                eptblm(ldeg,istp,id)=(0.d0,0.d0)
                eppalm(ldeg,istp,id)=(0.d0,0.d0)
                eppblm(ldeg,istp,id)=(0.d0,0.d0)
              enddo
c
              urlm(ldeg,6,id)=(0.d0,0.d0)
              utlm(ldeg,6,id)=(0.d0,0.d0)
              uplm(ldeg,6,id)=(0.d0,0.d0)
              grlm(ldeg,6,id)=(0.d0,0.d0)
              gtlm(ldeg,6,id)=(0.d0,0.d0)
              gplm(ldeg,6,id)=(0.d0,0.d0)
              errlm(ldeg,6,id)=(0.d0,0.d0)
              ertlm(ldeg,6,id)=(0.d0,0.d0)
              erplm(ldeg,6,id)=(0.d0,0.d0)
              etrlm(ldeg,6,id)=(0.d0,0.d0)
              ett0lm(ldeg,6,id)=(0.d0,0.d0)
              ettalm(ldeg,6,id)=(0.d0,0.d0)
              ettblm(ldeg,6,id)=(0.d0,0.d0)
              etp0lm(ldeg,6,id)=(0.d0,0.d0)
              etpalm(ldeg,6,id)=(0.d0,0.d0)
              eprlm(ldeg,6,id)=(0.d0,0.d0)
              ept0lm(ldeg,6,id)=(0.d0,0.d0)
              eptalm(ldeg,6,id)=(0.d0,0.d0)
              epp0lm(ldeg,6,id)=(0.d0,0.d0)
              eppalm(ldeg,6,id)=(0.d0,0.d0)
              eppblm(ldeg,6,id)=(0.d0,0.d0)
            enddo
c
            do ldeg=2,ldegf-1-id
              ca=dcmplx(dble(ldeg+3)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-2)/dble(2*ldeg-1),0.d0)
              do istp=4,5
                ettalm(ldeg,istp,id)=ettalm(ldeg,istp,id-1)
     &                        -ca*ettalm(ldeg+1,istp,id-1)
     &                        -cb*ettalm(ldeg-1,istp,id-1)
                ettblm(ldeg,istp,id)=ettblm(ldeg,istp,id-1)
     &                        -ca*ettblm(ldeg+1,istp,id-1)
     &                        -cb*ettblm(ldeg-1,istp,id-1)
                etpalm(ldeg,istp,id)=etpalm(ldeg,istp,id-1)
     &                        -ca*etpalm(ldeg+1,istp,id-1)
     &                        -cb*etpalm(ldeg-1,istp,id-1)
                etpblm(ldeg,istp,id)=etpblm(ldeg,istp,id-1)
     &                        -ca*etpblm(ldeg+1,istp,id-1)
     &                        -cb*etpblm(ldeg-1,istp,id-1)
                eptalm(ldeg,istp,id)=eptalm(ldeg,istp,id-1)
     &                        -ca*eptalm(ldeg+1,istp,id-1)
     &                        -cb*eptalm(ldeg-1,istp,id-1)
                eptblm(ldeg,istp,id)=eptblm(ldeg,istp,id-1)
     &                        -ca*eptblm(ldeg+1,istp,id-1)
     &                        -cb*eptblm(ldeg-1,istp,id-1)
                eppalm(ldeg,istp,id)=eppalm(ldeg,istp,id-1)
     &                        -ca*eppalm(ldeg+1,istp,id-1)
     &                        -cb*eppalm(ldeg-1,istp,id-1)
                eppblm(ldeg,istp,id)=eppblm(ldeg,istp,id-1)
     &                        -ca*eppblm(ldeg+1,istp,id-1)
     &                        -cb*eppblm(ldeg-1,istp,id-1)
              enddo
c
              urlm(ldeg,6,id)=urlm(ldeg,6,id-1)
     &                    -ca*urlm(ldeg+1,6,id-1)
     &                    -cb*urlm(ldeg-1,6,id-1)
              utlm(ldeg,6,id)=utlm(ldeg,6,id-1)
     &                    -ca*utlm(ldeg+1,6,id-1)
     &                    -cb*utlm(ldeg-1,6,id-1)
              uplm(ldeg,6,id)=uplm(ldeg,6,id-1)
     &                    -ca*uplm(ldeg+1,6,id-1)
     &                    -cb*uplm(ldeg-1,6,id-1)
              grlm(ldeg,6,id)=grlm(ldeg,6,id-1)
     &                    -ca*grlm(ldeg+1,6,id-1)
     &                    -cb*grlm(ldeg-1,6,id-1)
              gtlm(ldeg,6,id)=gtlm(ldeg,6,id-1)
     &                    -ca*gtlm(ldeg+1,6,id-1)
     &                    -cb*gtlm(ldeg-1,6,id-1)
              gplm(ldeg,6,id)=gplm(ldeg,6,id-1)
     &                    -ca*gplm(ldeg+1,6,id-1)
     &                    -cb*gplm(ldeg-1,6,id-1)
              errlm(ldeg,6,id)=errlm(ldeg,6,id-1)
     &                     -ca*errlm(ldeg+1,6,id-1)
     &                     -cb*errlm(ldeg-1,6,id-1)
              ertlm(ldeg,6,id)=ertlm(ldeg,6,id-1)
     &                     -ca*ertlm(ldeg+1,6,id-1)
     &                     -cb*ertlm(ldeg-1,6,id-1)
              erplm(ldeg,6,id)=erplm(ldeg,6,id-1)
     &                     -ca*erplm(ldeg+1,6,id-1)
     &                     -cb*erplm(ldeg-1,6,id-1)
              etrlm(ldeg,6,id)=etrlm(ldeg,6,id-1)
     &                     -ca*etrlm(ldeg+1,6,id-1)
     &                     -cb*etrlm(ldeg-1,6,id-1)
              ett0lm(ldeg,6,id)=ett0lm(ldeg,6,id-1)
     &                      -ca*ett0lm(ldeg+1,6,id-1)
     &                      -cb*ett0lm(ldeg-1,6,id-1)
              ettalm(ldeg,6,id)=ettalm(ldeg,6,id-1)
     &                      -ca*ettalm(ldeg+1,6,id-1)
     &                      -cb*ettalm(ldeg-1,6,id-1)
              ettblm(ldeg,6,id)=ettblm(ldeg,6,id-1)
     &                      -ca*ettblm(ldeg+1,6,id-1)
     &                      -cb*ettblm(ldeg-1,6,id-1)
              etp0lm(ldeg,6,id)=etp0lm(ldeg,6,id-1)
     &                      -ca*etp0lm(ldeg+1,6,id-1)
     &                      -cb*etp0lm(ldeg-1,6,id-1)
              etpalm(ldeg,6,id)=etpalm(ldeg,6,id-1)
     &                      -ca*etpalm(ldeg+1,6,id-1)
     &                      -cb*etpalm(ldeg-1,6,id-1)
              eprlm(ldeg,6,id)=eprlm(ldeg,6,id-1)
     &                     -ca*eprlm(ldeg+1,6,id-1)
     &                     -cb*eprlm(ldeg-1,6,id-1)
              ept0lm(ldeg,6,id)=ept0lm(ldeg,6,id-1)
     &                      -ca*ept0lm(ldeg+1,6,id-1)
     &                      -cb*ept0lm(ldeg-1,6,id-1)
              eptalm(ldeg,6,id)=eptalm(ldeg,6,id-1)
     &                      -ca*eptalm(ldeg+1,6,id-1)
     &                      -cb*eptalm(ldeg-1,6,id-1)
              epp0lm(ldeg,6,id)=epp0lm(ldeg,6,id-1)
     &                      -ca*epp0lm(ldeg+1,6,id-1)
     &                      -cb*epp0lm(ldeg-1,6,id-1)
              eppalm(ldeg,6,id)=eppalm(ldeg,6,id-1)
     &                      -ca*eppalm(ldeg+1,6,id-1)
     &                      -cb*eppalm(ldeg-1,6,id-1)
              eppblm(ldeg,6,id)=eppblm(ldeg,6,id-1)
     &                      -ca*eppblm(ldeg+1,6,id-1)
     &                      -cb*eppblm(ldeg-1,6,id-1)
            enddo
          enddo
c
          do is=isg1(ig),isg2(ig)
            wavelet=wvf(lf,is)
     &             *cdexp(-dcmplx(-fi,f)*dcmplx(PI2*togs(is),0.d0))
            do ir=1,nr
              id=idr(is,ir)
c
              call taper(ldegtap(1,is,ir),ldegtap(4,is,ir),tap(0))
c
              call legendre(dis(is,ir),plm,ldegtap(4,is,ir),ldegmax)
c
              do ldeg=ldegtap(1,is,ir),ldegtap(4,is,ir)
c
                cfac=dcmplx(tap(ldeg),0.d0)*wavelet
                if(id.gt.0)cfac=cfac/ssf(is,ir)**id
c
                cp0=dcmplx(plm(ldeg,0),0.d0)*cfac
                cp1=dcmplx(plm(ldeg,1),0.d0)*cfac
                cp2=dcmplx(plm(ldeg,2),0.d0)*cfac
c
                dur=sf3(is)*urlm(ldeg,1,id)*cp0
     &             +expl(is)*urlm(ldeg,2,id)*cp0
     &             +clvd(is)*urlm(ldeg,3,id)*cp0
     &             +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &              *urlm(ldeg,4,id)*cp1*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &              *urlm(ldeg,5,id)*cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &              *urlm(ldeg,6,id)*cp2*ssd(is,ir)**2
c
                dut=sf3(is)*utlm(ldeg,1,id)*cp1*ssd(is,ir)
     &             +expl(is)*utlm(ldeg,2,id)*cp1*ssd(is,ir)
     &             +clvd(is)*utlm(ldeg,3,id)*cp1*ssd(is,ir)
     &             +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &              *utlm(ldeg,4,id)*cp1
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &              *utlm(ldeg,5,id)*cp1
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &              *utlm(ldeg,6,id)*cp2*ssd(is,ir)
c
                dup=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &              *uplm(ldeg,4,id)*cp1
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &              *uplm(ldeg,5,id)*cp1
     &             +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &              *uplm(ldeg,6,id)*cp2*ssd(is,ir)
c
                dgr=sf3(is)*grlm(ldeg,1,id)*cp0
     &             +expl(is)*grlm(ldeg,2,id)*cp0
     &             +clvd(is)*grlm(ldeg,3,id)*cp0
     &             +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &              *grlm(ldeg,4,id)*cp1*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &              *grlm(ldeg,5,id)*cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &              *grlm(ldeg,6,id)*cp2*ssd(is,ir)**2
c
                dgt=sf3(is)*gtlm(ldeg,1,id)*cp1*ssd(is,ir)
     &             +expl(is)*gtlm(ldeg,2,id)*cp1*ssd(is,ir)
     &             +clvd(is)*gtlm(ldeg,3,id)*cp1*ssd(is,ir)
     &             +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &              *gtlm(ldeg,4,id)*cp1
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &              *gtlm(ldeg,5,id)*cp1
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &              *gtlm(ldeg,6,id)*cp2*ssd(is,ir)
c
                dgp=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &              *gplm(ldeg,4,id)*cp1
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &              *gplm(ldeg,5,id)*cp1
     &             +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &              *gplm(ldeg,6,id)*cp2*ssd(is,ir)
c
                urdr=sf3(is)*errlm(ldeg,1,id)*cp0
     &              +expl(is)*errlm(ldeg,2,id)*cp0
     &              +clvd(is)*errlm(ldeg,3,id)*cp0
     &              +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &               *errlm(ldeg,4,id)*cp1*ssd(is,ir)
     &              +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &               *errlm(ldeg,5,id)*cp1*ssd(is,ir)
     &              +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &               *errlm(ldeg,6,id)*cp2*ssd(is,ir)**2
c
                urdt=sf3(is)*ertlm(ldeg,1,id)*cp1*ssd(is,ir)
     &              +expl(is)*ertlm(ldeg,2,id)*cp1*ssd(is,ir)
     &              +clvd(is)*ertlm(ldeg,3,id)*cp1*ssd(is,ir)
     &              +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &               *ertlm(ldeg,4,id)*cp1
     &              +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &               *ertlm(ldeg,5,id)*cp1
     &              +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &               *ertlm(ldeg,6,id)*cp2*ssd(is,ir)
c
                urdp=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &               *erplm(ldeg,4,id)*cp1
     &              +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &               *erplm(ldeg,5,id)*cp1
     &              +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &               *erplm(ldeg,6,id)*cp2*ssd(is,ir)
c
                utdr=sf3(is)*etrlm(ldeg,1,id)*cp1*ssd(is,ir)
     &              +expl(is)*etrlm(ldeg,2,id)*cp1*ssd(is,ir)
     &              +clvd(is)*etrlm(ldeg,3,id)*cp1*ssd(is,ir)
     &              +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &               *etrlm(ldeg,4,id)*cp1
     &              +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &               *etrlm(ldeg,5,id)*cp1
     &              +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &               *etrlm(ldeg,6,id)*cp2*ssd(is,ir)
c
                utdt=sf3(is)*(ett0lm(ldeg,1,id)*cp0
     &                       +ettalm(ldeg,1,id)*cp1*csd(is,ir))
     &              +expl(is)*(ett0lm(ldeg,2,id)*cp0
     &                        +ettalm(ldeg,2,id)*cp1*csd(is,ir))
     &              +clvd(is)*(ett0lm(ldeg,3,id)*cp0
     &                        +ettalm(ldeg,3,id)*cp1*csd(is,ir))
     &              +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &               *(ett0lm(ldeg,4,id)*cp1*ssd(is,ir)
     &                +ettalm(ldeg,4,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +ettblm(ldeg,4,id)*cp2*ssd(is,ir))
     &              +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &               *(ett0lm(ldeg,5,id)*cp1*ssd(is,ir)
     &                +ettalm(ldeg,5,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +ettblm(ldeg,5,id)*cp2*ssd(is,ir))
     &              +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &               *(ett0lm(ldeg,6,id)*cp2*ssd(is,ir)**2
     &                +ettalm(ldeg,6,id)*cp2*csd(is,ir)
     &                +ettblm(ldeg,6,id)*cp2)
c
                utdp=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &               *(etp0lm(ldeg,4,id)*cp1*ssd(is,ir)
     &                +etpalm(ldeg,4,id)*cp2*ssd(is,ir)
     &                +etpblm(ldeg,4,id)*cp2*ssd(is,ir)*csd(is,ir))
     &              +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &               *(etp0lm(ldeg,5,id)*cp1*ssd(is,ir)
     &                +etpalm(ldeg,5,id)*cp2*ssd(is,ir)
     &                +etpblm(ldeg,5,id)*cp2*ssd(is,ir)*csd(is,ir))
     &              +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &               *(etp0lm(ldeg,6,id)*cp2
     &                +etpalm(ldeg,6,id)*cp2*csd(is,ir))
c
                updr=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &               *eprlm(ldeg,4,id)*cp1
     &              +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &               *eprlm(ldeg,5,id)*cp1
     &              +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &               *eprlm(ldeg,6,id)*cp2*ssd(is,ir)
c
                updt=(sf1(is)*ssa(is,ir)-sf2(is)*csa(is,ir))
     &               *(ept0lm(ldeg,4,id)*cp1*ssd(is,ir)
     &                +eptalm(ldeg,4,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +eptblm(ldeg,4,id)*cp2*ssd(is,ir))
     &              +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &               *(ept0lm(ldeg,5,id)*cp1*ssd(is,ir)
     &                +eptalm(ldeg,5,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +eptblm(ldeg,5,id)*cp2*ssd(is,ir))
     &              +(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &               *(ept0lm(ldeg,6,id)*cp2*csd(is,ir)
     &                +eptalm(ldeg,6,id)*cp2)
c
                updp=sf3(is)*(epp0lm(ldeg,1,id)*cp0
     &                       +eppalm(ldeg,1,id)*cp1*csd(is,ir))
     &              +expl(is)*(epp0lm(ldeg,2,id)*cp0
     &                        +eppalm(ldeg,2,id)*cp1*csd(is,ir))
     &              +clvd(is)*(epp0lm(ldeg,3,id)*cp0
     &                        +eppalm(ldeg,3,id)*cp1*csd(is,ir))
     &              +(sf1(is)*csa(is,ir)+sf2(is)*ssa(is,ir))
     &               *(epp0lm(ldeg,4,id)*cp1*ssd(is,ir)
     &                +eppalm(ldeg,4,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +eppblm(ldeg,4,id)*cp2*ssd(is,ir))
     &              +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &               *(epp0lm(ldeg,5,id)*cp1*ssd(is,ir)
     &                +eppalm(ldeg,5,id)*cp2*ssd(is,ir)*csd(is,ir)
     &                +eppblm(ldeg,5,id)*cp2*ssd(is,ir))
     &              +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &               *(epp0lm(ldeg,6,id)*cp2*ssd(is,ir)**2
     &                +eppalm(ldeg,6,id)*cp2*csd(is,ir)
     &                +eppblm(ldeg,6,id)*cp2)
c
                ue(lf,ir)=ue(lf,ir)+dut*ssb(is,ir)+dup*csb(is,ir)
                un(lf,ir)=un(lf,ir)-dut*csb(is,ir)+dup*ssb(is,ir)
                uz(lf,ir)=uz(lf,ir)+dur
c
                ge(lf,ir)=ge(lf,ir)+dgt*ssb(is,ir)+dgp*csb(is,ir)
                gn(lf,ir)=gn(lf,ir)-dgt*csb(is,ir)+dgp*ssb(is,ir)
                gz(lf,ir)=gz(lf,ir)+dgr
c
                roe(lf,ir)=roe(lf,ir)+(0.5d0,0.d0)
     &              *( (urdp-updr)*ssb(is,ir)+(utdr-urdt)*csb(is,ir))
                ron(lf,ir)=ron(lf,ir)+(0.5d0,0.d0)
     &              *(-(urdp-updr)*csb(is,ir)+(utdr-urdt)*ssb(is,ir))
                roz(lf,ir)=roz(lf,ir)+(0.5d0,0.d0)*(updt-utdp)
c
                rot(1,1)=ssb(is,ir)
                rot(1,2)=csb(is,ir)
                rot(1,3)=(0.d0,0.d0)
                rot(2,1)=-csb(is,ir)
                rot(2,2)=ssb(is,ir)
                rot(2,3)=(0.d0,0.d0)
                rot(3,1)=(0.d0,0.d0)
                rot(3,2)=(0.d0,0.d0)
                rot(3,3)=(1.d0,0.d0)
c
                rtz(1,1)=utdt
                rtz(1,2)=(0.5d0,0.d0)*(utdp+updt)
                rtz(1,3)=(0.5d0,0.d0)*(utdr+urdt)
                rtz(2,1)=(0.5d0,0.d0)*(updt+utdp)
                rtz(2,2)=updp
                rtz(2,3)=(0.5d0,0.d0)*(updr+urdp)
                rtz(3,1)=(0.5d0,0.d0)*(urdt+utdr)
                rtz(3,2)=(0.5d0,0.d0)*(urdp+updr)
                rtz(3,3)=urdr
c
                do i=1,3
                  do j=1,3
                    swp(i,j)=(0.d0,0.d0)
                    do k=1,3
                      swp(i,j)=swp(i,j)+rot(i,k)*rtz(k,j)
                    enddo
                  enddo
                enddo
                do i=1,3
                  do j=1,3
                    enz(i,j)=(0.d0,0.d0)
                    do k=1,3
                      enz(i,j)=enz(i,j)+swp(i,k)*rot(j,k)
                    enddo
                  enddo
                enddo
c
                uee(lf,ir)=uee(lf,ir)+enz(1,1)
                uen(lf,ir)=uen(lf,ir)+enz(1,2)
                uez(lf,ir)=uez(lf,ir)+enz(1,3)
                unn(lf,ir)=unn(lf,ir)+enz(2,2)
                unz(lf,ir)=unz(lf,ir)+enz(2,3)
                uzz(lf,ir)=uzz(lf,ir)+enz(3,3)
c
                eii=enz(1,1)+enz(2,2)+enz(3,3)
c
                see(lf,ir)=see(lf,ir)
     &                    +lamr*eii+(2.d0,0.d0)*muer*enz(1,1)
                sen(lf,ir)=sen(lf,ir)
     &                    +(2.d0,0.d0)*muer*enz(1,2)
                sez(lf,ir)=sez(lf,ir)
     &                    +(2.d0,0.d0)*muer*enz(1,3)
                snn(lf,ir)=snn(lf,ir)
     &                    +lamr*eii+(2.d0,0.d0)*muer*enz(2,2)
                snz(lf,ir)=snz(lf,ir)
     &                    +(2.d0,0.d0)*muer*enz(2,3)
                szz(lf,ir)=szz(lf,ir)
     &                    +lamr*eii+(2.d0,0.d0)*muer*enz(3,3)
c
              enddo
            enddo
          enddo
          write(*,'(i6,a,f10.4,3(a,i5))')lf,'.',1.0d+03*f,
     &                        ' mHz: spectra read: ',ldegf,
     &                        ', used: ',ldeglw,' - ',ldegneed
        enddo
c
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        close(28)
        write(*,'(i6,a)')lf-1,' spectra read from '
     &                      //specfile(ig)(1:40)
500     continue
      enddo
c
      deallocate(plm,tap,ldegtap,wvf,
     &           expl,clvd,ss12,ss11,ds31,ds23,sf1,sf2,sf3,sfr,sft,sfp,
     &           mrr,mtt,mpp,mrt,mtp,lats,lons,deps,togs,trss,
     &           ul0,vl0,wl0,el0,fl0,gl0,pl0,ql0,
     &           urlm,utlm,uplm,errlm,ertlm,erplm,
     &           ett0lm,ettalm,ettblm,etp0lm,etpalm,etpblm,eprlm,
     &           ept0lm,eptalm,eptblm,epp0lm,eppalm,eppblm,
     &           lyupp,lyups,lyupt,lylwp,lylws,lylwt,
     &           grndep,grnsel,lygrn,
     &           specfile,uspecfile,vspecfile,wspecfile,
     &           especfile,fspecfile,gspecfile,pspecfile,qspecfile,
     &           isg1,isg2,nsg)
c
      return
      end
