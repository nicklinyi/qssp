      subroutine qpstart0g(ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c     for n = 0
c
      integer*4 ly,jh
      complex*16 cy(3)
c
      integer*4 i
c
      if(jh.eq.0.and.ly.ge.lys)then
        if(lylwa.gt.0.or.lylwb.le.ly0)then
          do i=1,3
            cy(i)=mas3x3up(i,1,ly)
          enddo
        else
          call qpstart0(ly,cy)
        endif
      else if(jh.eq.1.and.ly.ge.lys.or.jh.eq.2.and.ly.le.lys)then
        call qpstart0a(ly,jh,cy)
      else
        stop ' Error in qpstart0g: bad jh value!'
      endif
c
      return
      end