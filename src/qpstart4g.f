      subroutine qpstart4g(ldeg,ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly,jh
      complex*16 cy(4,2)
c
      integer*4 i,j
c
      if(jh.eq.0.and.ly.ge.lys)then
        if(lylwa.gt.0.or.lylwb.le.ly0)then
          do j=1,2
            do i=1,4
              cy(i,j)=mas4x4up(i,2*j-1,ly)
            enddo
          enddo
c
c         y1 <- y1 (normal displacement)
c         y2 <- y3 (horizontal displacement)
c         y3 <- y5 (potential)
c         y4 <- y6 (gravity)
c
          do j=1,2
            cy(2,j)=(cy(1,j)*cgrup(ly)/crrup(ly)
     &              -cy(2,j)/(croup(ly)*crrup(ly)**2)
     &              -cy(3,j))/comi2
          enddo
        else
          call qpstart4(ldeg,ly,cy)
        endif
      else if(jh.eq.1.and.ly.ge.lys.or.jh.eq.2.and.ly.le.lys)then
        call qpstart4a(ldeg,ly,jh,cy)
c
c       y1 <- y1 (normal displacement)
c       y2 <- y3 (horizontal displacement)
c       y3 <- y5 (potential)
c       y4 <- y6 (gravity)
c
        do j=1,2
          cy(2,j)=(cy(1,j)*cgrup(ly)/crrup(ly)
     &            -cy(2,j)/(croup(ly)*crrup(ly)**2)
     &            -cy(3,j))/comi2
        enddo
      else
        stop ' Error in qpstart4g: bad jh value!'
      endif
c
      return
      end