      subroutine qpstart6g(ldeg,ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly,jh
      complex*16 cy(6,3)
c
      integer*4 i,j
c
      if(jh.eq.0.and.ly.ge.lys)then
        if(lylwa.gt.0.or.lylwb.le.ly0)then
          do j=1,3
            do i=1,6
              cy(i,j)=mas6x6up(i,2*j-1,ly)
            enddo
          enddo
        else
          call qpstart6(ldeg,ly,cy)
        endif
      else if(jh.eq.1.and.ly.ge.lys.or.jh.eq.2.and.ly.le.lys)then
        call qpstart6a(ldeg,ly,jh,cy)
      else
        stop ' Error in qpstart6g: bad jh value!'
      endif
c
      return
      end