      subroutine ruku(y,i0,j0,ly,ldeg,difmat,rr1,rr2,nrr0)
      implicit none
      integer*4 i0,j0,ly,ldeg,nrr0
      real*8 rr1,rr2
      complex*16 y(i0,j0)
      external difmat
c
      integer*4 nrrmax
      parameter(nrrmax=512)
c
      integer*4 i,j,k,irr,nrr,nn
      real*8 rr,drr,epsilon
      real*8 yabs(6,3)
      complex*16 cdrr
      complex*16 y1(6,3),y2(6,3),yk(6,4),mat(6,6,3)
      save nn
c
      real*8 eps
      complex*16 c1,c2,c6
      data eps/1.0d-04/
      data c1,c2,c6/(1.d0,0.d0),(2.d0,0.d0),(6.d0,0.d0)/
c
      do j=1,j0
        do i=1,i0
          y2(i,j)=y(i,j)
        enddo
      enddo
      nrr=max0(2,nrr0/4)
100   continue
      drr=(rr2-rr1)/dble(nrr)
      cdrr=dcmplx(drr,0.d0)
      do j=1,j0
        do i=1,i0
          y1(i,j)=y2(i,j)
          y2(i,j)=y(i,j)
          yabs(i,j)=cdabs(y(i,j))
        enddo
      enddo
      do irr=1,nrr
        rr=rr1+dble(irr-1)*drr
        call difmat(ly,ldeg,rr,mat(1,1,1))
        call difmat(ly,ldeg,rr+0.5d0*drr,mat(1,1,2))
        call difmat(ly,ldeg,rr+drr,mat(1,1,3))
        do j=1,j0
          do i=1,i0
            yk(i,1)=(0.d0,0.d0)
            do k=1,i0
              yk(i,1)=yk(i,1)+cdrr*mat(i,k,1)*y2(k,j)
            enddo
          enddo
c
          do i=1,i0
            yk(i,2)=(0.d0,0.d0)
            do k=1,i0
              yk(i,2)=yk(i,2)+cdrr*mat(i,k,2)*(y2(k,j)+yk(k,1)/c2)
            enddo
          enddo
c
          do i=1,i0
            yk(i,3)=(0.d0,0.d0)
            do k=1,i0
              yk(i,3)=yk(i,3)+cdrr*mat(i,k,2)*(y2(k,j)+yk(k,2)/c2)
            enddo
          enddo
c
          do i=1,i0
            yk(i,4)=(0.d0,0.d0)
            do k=1,i0
              yk(i,4)=yk(i,4)+cdrr*mat(i,k,3)*(y2(k,j)+yk(k,3))
            enddo
          enddo
c
          do i=1,i0
            y2(i,j)=y2(i,j)+(yk(i,1)+c2*(yk(i,2)+yk(i,3))+yk(i,4))/c6
          enddo
        enddo
      enddo
c
      do j=1,j0
        do i=1,i0
          yabs(i,j)=dmax1(cdabs(y(i,j)),cdabs(y2(i,j)))
        enddo
      enddo
c
      if(nrr.lt.nrrmax)then
        do j=1,j0
          do i=1,i0
            epsilon=eps*(cdabs(y2(i,j)-y(i,j))+eps*yabs(i,j))
            if(cdabs(y2(i,j)-y1(i,j)).gt.epsilon)then
              nrr=nrr*2
              goto 100
            endif
          enddo
        enddo
      else
        print *, ' Warning in ruku: Convergence problem!'
      endif
c
      do j=1,j0
        do i=1,i0
          y(i,j)=y2(i,j)
        enddo
      enddo
c
      nrr0=nrr
c
      return
      end