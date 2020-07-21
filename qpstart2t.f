      subroutine qpstart2t(ldeg,ly,jh,cy)
      use qpalloc
      implicit none
c
c     calculate fundamental solution vectors for homogeneous sphere
c
      integer*4 ldeg,ly,jh
      complex*16 cy(2)
c
      real*8 teps
      complex*16 clm1,ca,cb
c
      real*8 eps
      data eps/1.0d-08/
c
      clm1=dcmplx(dble(ldeg-1),0.d0)
c
      if(jh.eq.1.and.ly.ge.lys)then
        cy(1)=(1.d0,0.d0)
        cy(2)=cmu(ly)*(clm1-zh12sh(ldeg,1,ly))
c
        teps=dabs(dimag(zh12sh(ldeg,1,ly)))/cdabs(zh12sh(ldeg,1,ly))
        if(teps.lt.eps)then
          ca=dcmplx((teps/eps)**2,0.d0)
          cb=(1.d0,0.d0)-ca
          cy(1)=mat2x2up(1,1,ly)
          cy(2)=ca*cy(2)*mat2x2up(1,1,ly)+cb*mat2x2up(2,1,ly)
        endif
      else if(jh.eq.2.and.ly.le.lys)then
        cy(1)=(1.d0,0.d0)
        cy(2)=cmu(ly)*(clm1-zh12sh(ldeg,2,ly))
      else
        stop ' Error in qpstart2t: bad jh value!'
      endif
      return
      end