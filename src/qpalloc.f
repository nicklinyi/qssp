      module qpalloc
c
c     CONSTANTS
c     =========
      integer*4 ndmax
      parameter(ndmax=2)
      real*8 PI,PI2
      parameter(PI=3.14159265358979d0,PI2=6.28318530717959d0)
      real*8 DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02,KM2M=1.0d+03)
      real*8 BIGG
      parameter(BIGG=6.6732d-11)
      real*8 RESOLUT
      parameter(RESOLUT=0.01d0)
      real*8 FLTAPER
      parameter(FLTAPER=0.2d0)
      real*8 FSBUP,FSBLW,FSBREF
      parameter(FSBUP=1.0d+02,FSBLW=2.5d-04,FSBREF=1.0d+00)
      real*8 fbvatm,fbvocean,fbvcore
      parameter(fbvatm=1.0d-06,fbvocean=0.0d+00,fbvcore=0.0d+00)
c
      logical*2 selpsv,selsh,setzh12a,setzh12t
      logical*2 nogravity,freesurf
      logical*2 dispersion
c
      integer*4 ngrn,nt,ntcut,ntcutout,nf,nfcut,nbpf
      integer*4 lyadd,ipatha,ipathb,ldeggr,ldegmin,ldegcut,ldegmax
      integer*4 nr,igfirst,iglast
      integer*4 ns,l0,lymax
      integer*4 lys,lyr,lylwa,lylwb,lyos,lyob,lycm,lycc,ly0
      integer*4 icmp(11)
c
      real*8 dt,dtout,df,fi,fcut,fgr,rratmos,depatmos
      real*8 rearth,rr0,minpath,maxpath
      real*8 slwmax,slwlwcut,slwupcut,f1corner,f2corner
      real*8 dpr,freeairgrd
      real*8 qsmin,togsmin
c
      complex*16 comi,comi2
c
      character*80 dispout(3),veloout(3),acceout(3),
     &             rotaout(3),rotarateout(3),
     &             strainout(6),strainrateout(6),
     &             stressout(6),stressrateout(6),
     &             gravout(3),grmout
c
      logical*2,allocatable:: ksmallp(:),ksmalls(:),ksmallt(:)
c
      integer*4,allocatable:: lygrn(:),grnsel(:),ldegpsv(:),ldegsh(:)
      integer*4,allocatable:: nruku(:,:)
      integer*4,allocatable:: isg1(:),isg2(:),nsg(:)
      integer*4,allocatable:: idr(:,:)
      integer*4,allocatable:: lylwp(:),lylws(:),lylwt(:),
     &                        lyupp(:),lyups(:),lyupt(:)
      integer*4,allocatable:: ldegtap(:,:,:)
c
      real*8,allocatable:: latr(:),lonr(:),tred(:),grndep(:)
      real*8,allocatable:: sfr(:),sft(:),sfp(:),
     &                     mtt(:),mpp(:),mrr(:),mtp(:),mpr(:),mrt(:),
     &                     lats(:),lons(:),deps(:),togs(:),trss(:)
      real*8,allocatable:: dis(:,:),plm(:,:)
      real*8,allocatable:: dp0(:),dp0up(:),dp0lw(:),
     &                     vp0(:),vp0up(:),vp0lw(:),
     &                     vs0(:),vs0up(:),vs0lw(:),
     &                     ro0(:),ro0up(:),ro0lw(:),
     &                     qp0(:),qp0up(:),qp0lw(:),
     &                     qs0(:),qs0up(:),qs0lw(:)
      real*8,allocatable:: rrup(:),rrlw(:),vpup(:),vplw(:),
     &                     vsup(:),vslw(:),roup(:),rolw(:),
     &                     qpup(:),qplw(:),qsup(:),qslw(:)
      real*8,allocatable:: xp(:),xs(:),xt(:)
      real*8,allocatable:: tap(:),disk(:)
c
      complex*16,allocatable:: ul0(:,:),vl0(:,:),wl0(:,:),
     &                         el0(:,:),fl0(:,:),gl0(:,:),
     &                         pl0(:,:),ql0(:,:)
      complex*16,allocatable:: urlm(:,:,:),utlm(:,:,:),uplm(:,:,:),
     &       errlm(:,:,:),ertlm(:,:,:),erplm(:,:,:),etrlm(:,:,:),
     &       ett0lm(:,:,:),ettalm(:,:,:),ettblm(:,:,:),etp0lm(:,:,:),
     &       etpalm(:,:,:),etpblm(:,:,:),eprlm(:,:,:),ept0lm(:,:,:),
     &       eptalm(:,:,:),eptblm(:,:,:),epp0lm(:,:,:),eppalm(:,:,:),
     &       eppblm(:,:,:),grlm(:,:,:),gtlm(:,:,:),gplm(:,:,:)
      complex*16,allocatable:: ue(:,:),un(:,:),uz(:,:),
     &                      roe(:,:),ron(:,:),roz(:,:),
     &                      uee(:,:),uen(:,:),uez(:,:),
     &                      unn(:,:),unz(:,:),uzz(:,:),
     &                      see(:,:),sen(:,:),sez(:,:),
     &                      snn(:,:),snz(:,:),szz(:,:),
     &                      ge(:,:),gn(:,:),gz(:,:),gm(:,:)
      complex*16,allocatable:: ssa(:,:),ss2a(:,:),csa(:,:),cs2a(:,:),
     &                  ssb(:,:),csb(:,:),ssd(:,:),csd(:,:),ssf(:,:)
      complex*16,allocatable:: crrup(:),crrlw(:),
     &           cro(:),croup(:),crolw(:),cla(:),claup(:),clalw(:),
     &           cmu(:),cmuup(:),cmulw(:),cvp(:),cvpup(:),cvplw(:),
     &           cvs(:),cvsup(:),cvslw(:),cgr(:),cgrup(:),cgrlw(:),
     &           cga(:),cgaup(:),cgalw(:)
      complex*16,allocatable:: cps(:,:),cpt(:,:),kp(:),ks(:),kt(:)
      complex*16,allocatable:: mat2x2up(:,:,:),mat2x2lw(:,:,:),
     &        mat2x2inv(:,:,:),mas3x3up(:,:,:),mas3x3lw(:,:,:),
     &        mas3x3inv(:,:,:),mas4x4up(:,:,:),mas4x4lw(:,:,:),
     &        mas4x4inv(:,:,:),mas6x6up(:,:,:),mas6x6lw(:,:,:),
     &        mas6x6inv(:,:,:)
      complex*16,allocatable:: zjup(:,:,:),zjlw(:,:,:),
     &       zhup(:,:,:),zhlw(:,:,:),wj(:,:,:),wh(:,:,:),
     &       zh12p(:,:,:),zh12sv(:,:,:),zh12sh(:,:,:)
      complex*16,allocatable:: zjupg(:),zjlwg(:),
     &           zhupg(:),zhlwg(:),wjg(:),whg(:)
      complex*16,allocatable:: cua(:,:),cypnorm(:,:)
      complex*16,allocatable:: wvf(:,:)
      complex*16,allocatable:: sf1(:),sf2(:),sf3(:),
     &                         expl(:),clvd(:),ss12(:),
     &                         ss11(:),ds31(:),ds23(:)
c
      character*80,allocatable:: specfile(:),uspecfile(:),vspecfile(:),
     &  wspecfile(:),especfile(:),fspecfile(:),gspecfile(:),
     &  pspecfile(:),qspecfile(:)
      character*10,allocatable:: rname(:)
c
      end module
