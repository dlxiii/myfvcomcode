! $Id: pcgm.F,v 1.1.2.2 2008/04/03 17:06:14 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.2 $
#ifdef MPI_LAM
# define MPI_BCAST MPI_BCAST_
# define MPI_BARRIER MPI_BARRIER_
# define MPI_REDUCE MPI_REDUCE_
#endif

!***********************************************************************
!***********************************************************************
      subroutine pcgm(A,x,nn,iq,jq,epsilon,kmax)
use mod_prec
!-----+----------------------------------------------------------------+
! Solve Ax=b by the preconditioned conjugate gradient method, where
! A is a sparse, symmetric and positiv definite matrix, b a known 
! vector, and x the unknown vector. In the call b is x.
! The preconditioner C is a diagonal matrix stored in a vector.
! 
! Nicolai Kliem, 11 Feb. 1998
!-----+----------------------------------------------------------------+

      implicit none

      REAL(SP) :: A(*),x(*)
      integer iq(*),jq(*)

      integer nn

      REAL(SP) :: p(nn),r(nn),w(nn),b(nn),z(nn),C(nn)
      REAL(SP) :: alpha,beta,rz,rznew,rzstop,pw,epsilon

      integer i,j,l,k,n,kmax

!-----+----------------------------------------------------------------+
! Initialize
! x=0
!-----+----------------------------------------------------------------+

! The preconditioner C could be calculated once and forever in 
! subroutine STATIONARYQ4, but this recuires that the C is transfered
! through the subroutine calls, and is not implemented here to make as
! few changes as possible.
      do n=1,nn
        call calja(n,n,iq,jq,j)
        C(n)=1.0_SP/A(j)
      enddo

      do n=1,nn
        r(n)=x(n)
        b(n)=x(n)
        x(n)=0.0_SP
      enddo

      do n=1,nn
        z(n)=C(n)*r(n)
        p(n)=z(n)
      enddo

      rz=0.0_SP
      do n=1,nn
        rz=rz+r(n)*z(n)
      enddo

      k=0
      rzstop=epsilon*rz

!-----+----------------------------------------------------------------+
! Loop
!-----+----------------------------------------------------------------+

      do while (k.lt.kmax.and.rz.gt.rzstop)
!
! w=A*p
!
        do i=1,nn
           w(i)=0.0_SP
           do j=iq(i),iq(i+1)-1
              l=jq(j)
              w(i) = w(i) + A(j)*p(l)
           enddo
        enddo
!        call sprsmlt(A,iq,jq,p,w,nn)
!
! alpha=rz/(pT*w)
!
        pw=0.0_SP
        do n=1,nn
          pw=pw+p(n)*w(n)
        enddo
        alpha=rz/pw
!
! x=x+alpha*p
!
        do n=1,nn
          x(n)=x(n)+alpha*p(n)
        enddo
!
! r=r-alpha*w
!
        do n=1,nn
          r(n)=r(n)-alpha*w(n)
        enddo
!
! z=C-1*r
!
        do n=1,nn
          z(n)=C(n)*r(n)
        enddo
!
! rznew=r*z
! beta=rznew/rz
! rz=rznew
!
        rznew=0.0_SP
        do n=1,nn
          rznew=rznew+r(n)*z(n)
        enddo
        beta=rznew/rz
        rz=rznew
!
! p=z+beta*p
!
        do n=1,nn
          p(n)=z(n)+beta*p(n)
        enddo

        k=k+1
!      WRITE(*,*) 'pcgm ',k,epsilon*rz/rzstop
      enddo

      WRITE(*,*) 'pcgm ',k,epsilon*rz/rzstop

      if (k.eq.kmax) then 
!          WRITE(*,*) '# pcgm',k,epsilon*rz/rzstop
          write(*,*) 'pb in PCGM solver'
!          stop
      endif

      return
      end subroutine pcgm
!***********************************************************************

!***********************************************************************
# if defined (MULTIPROCESSOR)
!***********************************************************************
      subroutine pcgm_para(A,x,nn,iq,jq,rmask,mesg,epsilon,kmax)
use mod_prec
use lims
use mod_par
!-----+----------------------------------------------------------------+
! Solve Ax=b by the preconditioned conjugate gradient method, where
! A is a sparse, symmetric and positiv definite matrix, b a known 
! vector, and x the unknown vector. In the call b is x.
! The preconditioner C is a diagonal matrix stored in a vector.
! 
! Nicolai Kliem, 11 Feb. 1998
!-----+----------------------------------------------------------------+
use mod_par

      implicit none

      REAL(SP) ::  A(*),x(*)
      integer iq(*),jq(*)
      REAL(SP) :: rmask(*)
      logical mesg

      integer nn

      REAL(SP)    p(nn),r(nn),w(nn),b(nn),z(nn),C(nn)
      REAL(SP)    alpha,beta,rz,rznew,rzstop,pw,epsilon

      integer i,j,l,k,kmax

      REAL(SP), allocatable :: sbuf(:,:),rbuf(:,:)
      INTEGER :: IERR

!-----+----------------------------------------------------------------+
! Initialize
! x=0
!-----+----------------------------------------------------------------+

! The preconditioner C could be calculated once and forever in 
! subroutine STATIONARYQ4, but this recuires that the C is transfered
! through the subroutine calls, and is not implemented here to make as
! few changes as possible.
      do i=1,nn
        call calja(i,i,iq,jq,j)
        C(i)=1.0_SP/A(j)
      enddo

 allocate(rbuf(0:MT,4))
 rbuf(1:nn/2,1) = C(1:nn/2)
 rbuf(1:nn/2,2) = C(nn/2+1:nn)
 rbuf(1:nn/2,3) = x(1:nn/2)
 rbuf(1:nn/2,4) = x(nn/2+1:nn)
 CALL EXCHANGE(NC,MT,4,MYID,NPROCS,rbuf) !Interprocessor Exchange   !
 CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,4,MYID,NPROCS,rbuf)
 C(1:nn/2)    = rbuf(1:nn/2,1)
 C(nn/2+1:nn) = rbuf(1:nn/2,2)
 x(1:nn/2)    = rbuf(1:nn/2,3)
 x(nn/2+1:nn) = rbuf(1:nn/2,4)
 deallocate(rbuf)
! CALL EXCHANGE(NC,MT,2,MYID,NPROCS,c) !Interprocessor Exchange   !
! CALL EXCHANGE(NC,MT,2,MYID,NPROCS,x) !Interprocessor Exchange   !

      do i=1,nn
        r(i)=x(i)
        b(i)=x(i)
        x(i)=0.0_SP
      enddo

      do i=1,nn
        z(i)=C(i)*r(i)
        p(i)=z(i)
      enddo

      rz=0.0_SP
      do i=1,nn
        rz=rz+r(i)*z(i)*rmask(i)
      enddo

      call reduce(rz)

      k=0
      rzstop=epsilon*rz

!-----+----------------------------------------------------------------+
! Loop
!-----+----------------------------------------------------------------+

      do while (k.lt.kmax.and.rz.gt.rzstop)
!
! w=A*p
!
        do i=1,nn
           w(i)=0.0_SP
           do j=iq(i),iq(i+1)-1
              l=jq(j)
              w(i) = w(i) + A(j)*p(l)
!if (i-nn/2.gt.0.and.ngid(max(1,i-nn/2))==501.and.mesg) write(*,*) w(i),a(j),p(l)
           enddo
!if (i-nn/2.gt.0.and.ngid(max(1,i-nn/2))==501.and.mesg) write(*,*) i,w(i)
        enddo
!        call sprsmlt(A,iq,jq,p,w,nn)
! CALL EXCHANGE(NC,MT,2,MYID,NPROCS,w) !Interprocessor Exchange   !
 allocate(rbuf(0:MT,2))
 rbuf(1:nn/2,1) = w(1:nn/2)
 rbuf(1:nn/2,2) = w(nn/2+1:nn)
 CALL EXCHANGE(NC,MT,2,MYID,NPROCS,rbuf) !Interprocessor Exchange   !
 CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,2,MYID,NPROCS,rbuf)
 w(1:nn/2)   = rbuf(1:nn/2,1)
 w(nn/2+1:nn)= rbuf(1:nn/2,2)
 deallocate(rbuf)
!do i=1,M
! if (ngid(i)==333) write(*,*) 'cg',myid,i,k,w(i),w(i+MT)
!enddo
!
! alpha=rz/(pT*w)
!
        pw=0.0_SP
        do i=1,nn
          pw=pw+p(i)*w(i)*rmask(i)
        enddo

        call reduce(pw)
        alpha=rz/pw
!if (mesg) write(*,*) k,'pw',pw,alpha
!
! x=x+alpha*p
!
        do i=1,nn
          x(i)=x(i)+alpha*p(i)
        enddo
!
! r=r-alpha*w
!
        do i=1,nn
          r(i)=r(i)-alpha*w(i)
        enddo
!
! z=C-1*r
!
        do i=1,nn
          z(i)=C(i)*r(i)
        enddo
!
! rznew=r*z
! beta=rznew/rz
! rz=rznew
!
        rznew=0.0_SP
        do i=1,nn
          rznew=rznew+r(i)*z(i)*rmask(i)
        enddo
        call reduce(rznew)
        beta=rznew/rz
        rz=rznew
!if (mesg) write(*,*) k,'rz',rz,beta
!
! p=z+beta*p
!
        do i=1,nn
          p(i)=z(i)+beta*p(i)
        enddo

        k=k+1
!      WRITE(*,*) 'pcgm ',k,epsilon*rz/rzstop
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,IERR)

! CALL EXCHANGE(NC,MT,2,MYID,NPROCS,x) !Interprocessor Exchange   !
 allocate(rbuf(0:MT,2))
 rbuf(1:nn/2,1) = x(1:nn/2)
 rbuf(1:nn/2,2) = x(nn/2+1:nn)
 CALL EXCHANGE(NC,MT,2,MYID,NPROCS,rbuf) !Interprocessor Exchange   !
 CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,2,MYID,NPROCS,rbuf)
 x(1:nn/2)    = rbuf(1:nn/2,1)
 x(nn/2+1:nn) = rbuf(1:nn/2,2)
 deallocate(rbuf)


if (mesg) then
      WRITE(*,*) 'pcgm ',k,epsilon*rz/rzstop

      if (k.eq.kmax) then 
          WRITE(*,*) '# pcgm',k,epsilon*rz/rzstop
          write(*,*) 'pb in PCGM solver'
!          stop
      endif
endif

      return
      end subroutine pcgm_para
!***********************************************************************


!===================================================================================|

subroutine reduce(r)
use mod_par
   IMPLICIT NONE
   real(SP) :: r,rtot
   INTEGER IERR

 CALL MPI_REDUCE(r,rtot,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
 r=rtot
 call MPI_BCAST(r,1,MPI_F,0,MPI_COMM_WORLD,IERR)

 return
end subroutine reduce
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|

subroutine reduce_max(r)
use mod_par
   IMPLICIT NONE
   REAL (SP) :: r,rtot
   INTEGER IERR

 CALL MPI_REDUCE(r,rtot,1,MPI_F,MPI_MAX,0,MPI_COMM_WORLD,IERR)
 r=rtot
 call MPI_BCAST(r,1,MPI_F,0,MPI_COMM_WORLD,IERR)

 return
end subroutine reduce_max
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif
