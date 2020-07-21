program qssp
use qpalloc
implicit none

    !  work space

integer :: ig,ierr,runtime
integer :: time
character(len=128) :: inputfile

!     read input file file

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
! write(*,'(a,$)')' the input data file is '
! read(*,'(a)')inputfile
!First, make sure the right number of inputs have been provided
if(COMMAND_ARGUMENT_COUNT().ne.1) then
  write(*,*)'Error, one command-line argument required, stopping'
  stop
endif

call get_command_argument(1,inputfile)

runtime=time()

open(10,file=inputfile,status='old')
call qpgetinp(10)
close(10)

call qpsublayer(ierr)

igfirst=0
do ig=ngrn,1,-1
    if(grnsel(ig).eq.1) igfirst=ig
enddo
iglast=ngrn+1
do ig=1,ngrn
    if(grnsel(ig).eq.1)iglast=ig
enddo

do ig=1,ngrn
if(grnsel(ig).eq.1)then
    lys=lygrn(ig)
    call qpgrnspec(ig)
endif
enddo
call qpwvint(ierr)
call qpfftinv(ierr)

runtime=time()-runtime
write(*,'(a)')' #############################################'
write(*,'(a)')' #                                           #'
write(*,'(a)')' #      End of computations with qssp2017    #'
write(*,'(a)')' #                                           #'
write(*,'(a,i10,a)')' #       Run time: ',runtime, &
                                          ' sec            #'
write(*,'(a)')' #############################################'
stop
end
