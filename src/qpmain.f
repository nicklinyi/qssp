      program qssp
      use qpalloc
      implicit none
c
c     work space
c
      integer*4 ig,ierr,runtime
      integer*4 time
      character*80 inputfile
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#    QQQQ         SSSSS        SSSSS       PPPPP     #'
      print *,'#   Q    Q       S            S            P    P    #'
      print *,'#   Q    Q        SSSS         SSSS        PPPPP     #'
      print *,'#   Q   QQ            S            S       P         #'
      print *,'#    QQQQQ       SSSSS        SSSSS        P         #'
      print *,'#                                                    #'
      print *,'#          Complete synthetic seismograms            #'
      print *,'#      (displacement/strain/stress/rotation)         #' 
      print *,'#                     based on                       #'
      print *,'#          a spherically symmetric earth model       #'
      print *,'#                                                    #'
      print *,'#                  (Version 2017)                    #'
      print *,'#   Last update (correction of errors): 2018-8-10    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#              Helmholtz Centre Potsdam              #'
      print *,'#    GFZ German Research Centre for Geosciences      #'
      print *,'#           Last modified: September 2017            #'
      print *,'#                                                    #'
      print *,'######################################################'
      print *,'                                                      '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
c     inputfile = 'Event_2016_01_31_17_38_59-prem.inp'
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpgetinp(10)
      close(10)
c
      call qpsublayer(ierr)
c
      igfirst=0
      do ig=ngrn,1,-1
        if(grnsel(ig).eq.1)igfirst=ig
      enddo
      iglast=ngrn+1
      do ig=1,ngrn
        if(grnsel(ig).eq.1)iglast=ig
      enddo
c
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          call qpgrnspec(ig)
        endif
      enddo
      call qpwvint(ierr)
      call qpfftinv(ierr)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qssp2017    #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
