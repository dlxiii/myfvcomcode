      subroutine icemohrcoulomb(e11,e12,e22,p,zeta,eta)
# if defined (ICE)
  use ice_variables
      implicit  none

      REAL(SP)    e11,e12,e22,p,zeta,eta,zetamax,etamax,alpha
      REAL(SP)    conv,Delta

      data    zetamax/1.0_SP/,etamax/1.0_SP/,alpha/17.0_SP/

!-----+----------------------------------------------------------------+

      conv=-e11-e22
      Delta=max(Delta0,sqrt((e11-e22)*(e11-e22) + 4.0_SP*e12*e12))

      if (conv.gt.0_SP) then
        zeta=min(zetamax,p/conv)
        eta=min(etamax,zeta*conv*tanh(alpha)/Delta)
!       eta=min(etamax,p*tanh(alpha)/Delta)
      else
        zeta=0.0_SP
        eta=0.0_SP
      endif

!-----+----------------------------------------------------------------+
      return
# endif
      end



      subroutine icenewtonian(p,zeta,eta)
# if defined (ICE)
  use ice_variables
      implicit  none

      REAL(SP)    p,zeta,eta,etamax

      etamax       = 2.5e8_SP              ! max coefficient of friction [kg/s]

!-----+----------------------------------------------------------------+

      zeta=min(0.5_SP*p/Delta0,etamax)
!      zeta=0.5_SP*p/Delta0
      eta=zeta*ellip2

!-----+----------------------------------------------------------------+
      return
# endif
      end


      subroutine iceviscplast(e11,e12,e22,p,zeta,eta)
# if defined (ICE)
  use ice_variables
      implicit  none

      REAL(SP)    e11,e12,e22,p,zeta,eta,Delta,etamax

      etamax       = 2.5e8_SP              ! max coefficient of friction [kg/s]

!-----+----------------------------------------------------------------+

      Delta=max(1e-20_SP,sqrt( (e11*e11+e22*e22)*(1.0_SP+ellip2) &
          +4.0_SP*e12*e12*ellip2 + 2.0_SP*e11*e22*(1.0_SP-ellip2) ))
      zeta=min(0.5_SP*p/Delta,etamax*p)
      eta=zeta*ellip2

!-----+----------------------------------------------------------------+
      return
# endif
      end

      subroutine icecavitating(e1,e2,p,zeta,eta)
# if defined (ICE)
  use ice_variables
  use control
      implicit  none

      REAL(SP)    e1,e2,p,eta,zeta,Delta,etamax,phi,sinphi

      etamax       = 2.5e8_SP              ! max coefficient of friction [kg/s]
      phi          = 30e0_SP              ! internal angle of friction [deg]
      sinphi     = sin (phi*deg2rad)

!-----+----------------------------------------------------------------+

          delta=max(SQRT(e1**2+e2**2),1e-12_SP)
          zeta=min(sinphi*max(0.0_SP,p)/delta,etamax)
          eta=0.0_SP

!-----+----------------------------------------------------------------+
      return
# endif
      end




      subroutine icevisco(etac,zetac)
# if defined (ICE)
  use ice_variables

      implicit  none

      integer i,ll,ii
      REAL(SP)    zeta,eta,e11,e12,e22,s11,s12,s22,e1,e2,delta
      REAL(SP)    arrcs,uel,vel,csel,pel,hel,ael,ah,tmp  !arr
      integer j,ij,k
      REAL(SP)    etac(nne),zetac(nne)
      REAL(SP)    ahi(nn)


!-----+----------------------------------------------------------------+
!  Build ice stresses

      do ll=1,nne
        csel=(cs(in(1,ll))+cs(in(2,ll))+cs(in(3,ll)))/3.0_SP
        uel=(uice(in(1,ll))+uice(in(2,ll))+uice(in(3,ll)))/3.0_SP
        vel=(vice(in(1,ll))+vice(in(2,ll))+vice(in(3,ll)))/3.0_SP
        pel=(pice(in(1,ll))+pice(in(2,ll))+pice(in(3,ll)))/3.0_SP

!        ARR = 0.5_SP*(dx(2,ll)*dy(3,ll)-dx(3,ll)*dy(2,ll))
!     ggao 11-17-2006
        arrcs=arr(ll)*csel

       if (isicec(ll)==1) then

        e11=0.0_SP
        e12=0.0_SP
        e22=0.0_SP
        do ii=1,3
          i=in(ii,ll)
          e11=e11 + uice(i)*dy(ii,ll) - vel*dx(ii,ll)*cs(i)*cs(i)/csel
          e12=e12 + vice(i)*dy(ii,ll) - uice(i)*dx(ii,ll)*cs(i) &
                                       +uel*dx(ii,ll)*cs(i)*cs(i)/csel
          e22=e22 - vice(i)*dx(ii,ll)*cs(i)
        enddo
        e11=0.50_SP*e11/arrcs
        e12=0.25_SP*e12/arrcs
        e22=0.50_SP*e22/arrcs

          e1=0.0_SP
          e2=0.0_SP
          do ii=1,3
            i=in(ii,ll)
            e1=e1 + uice(i)*dy(ii,ll)/cs(i) + vice(i)*dx(ii,ll)
            e2=e2 + vice(i)*dy(ii,ll)/cs(i) - uice(i)*dx(ii,ll)
          enddo
          e1=0.5_SP*e1/arrcs
          e2=0.5_SP*e2/arrcs

        if (rheology.eq.'viscplast') then
          call iceviscplast(e11,e12,e22,pel,zeta,eta)
        elseif (rheology.eq.'EVP') then
          call iceviscplast(e11,e12,e22,pel,zeta,eta)
        elseif (rheology.eq.'viscplast_conv') then
          call iceviscplast(e11,e12,e22,pel,zeta,eta)
        elseif (rheology.eq.'mohrcoulomb') then
          call icemohrcoulomb(e11,e12,e22,pel,zeta,eta)
        elseif (rheology.eq.'newtonian') then
          call icenewtonian(pel,zeta,eta)
        elseif (rheology.eq.'cavitating') then
          write(*,*) 'cavitating rheology not implemented'
          stop
        else
          zeta=0.0_SP
          eta=0.0_SP
        endif

        etac (ll) = eta
        zetac(ll) = zeta
       else
        etac (ll) = 0.0_SP
        zetac(ll) = 0.0_SP
       endif ! if cell is ice-covered

      enddo

!-----+----------------------------------------------------------------+
      return
# endif
      end



!-----+----------------------------------------------------------------+
! init the Gradient Solver for velocity
!-----+----------------------------------------------------------------+
subroutine calusolve(lhs,rhs,etac,zetac)
# if defined (ICE)
  use ice_variables
  use mod_par
  use all_vars

  implicit  none

  integer i,ll,ii,i1,i2,i3,ind,nb
  REAL(SP)    zeta,eta,tmp   !,arr
  integer j,ij
  REAL(SP), dimension(nne) :: etac, zetac
  complex (SP), dimension(nn) :: lhs,rhs,Fvis
  REAL(SP), dimension(2*nn) :: utmp,vtmp
  REAL(SP) :: tx,ty,px,py,fc
  REAL(SP) :: usolve(4*nn*nftrdim)   ! sparse matrix for momentum solver
  integer   iqi(2*nn+1),jqi(4*nn*nftrdim)   ! index sparse array
  REAL(SP) ::  rmask(nn,2)
  REAL(SP), allocatable, dimension(:,:) :: sbuf,rbuf

# if defined (MULTIPROCESSOR)
if (PAR) then
! rmask is used in masking out the halo nodes in the Conjugate Gradient method
! when computing the scalar product (norm, residual)

  rmask=1.0_SP

 ! the nodes along parallel bdy must be counted as a ratio
 allocate(rbuf(0:MT,2))
 rbuf(1:MT,1)=rmask(1:nn,1)
 rbuf(1:MT,2)=rmask(1:nn,2)
 CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,2,MYID,NPROCS,rbuf)
 rmask(1:nn,1)=1.0_SP/rbuf(1:MT,1)
 rmask(1:nn,2)=1.0_SP/rbuf(1:MT,2)
 deallocate(rbuf)

 ! the solid boundary velocities are zero anyway
  do i=1,nn
      if (isonb(i)==1) then 
        rmask(i,:)=0.0_SP
      endif
  enddo

 ! the halo nodes are then zeroed
  rmask(M+1:MT,:)=0.0_SP

endif
#endif

!-----+----------------------------------------------------------------+
!  Initialization

      usolve=0.0_SP
      FVis=0.0_SP*eye
!-----+----------------------------------------------------------------+
! find the index array for the uv sparse matrix

iqi(1)=1
! u-lines
do i=1,nn
nb=iaa(i+1)-iaa(i)
iqi(i+1)=iqi(i)+2*nb
ind=iqi(i)
 do j=iaa(i),iaa(i+1)-1
    jqi(ind+j-iaa(i)   )  = ja(j)
    jqi(ind+j-iaa(i)+nb)  = ja(j)+nn
 enddo
enddo
! v-lines
do i=1,nn
nb=iaa(i+1)-iaa(i)
iqi(nn+i+1)=iqi(nn+i)+2*nb
ind=iqi(nn+i)
 do j=iaa(i),iaa(i+1)-1
    jqi(ind+j-iaa(i)   )  = ja(j)
    jqi(ind+j-iaa(i)+nb)  = ja(j)+nn
 enddo
enddo
!do i=1,2*nn
! write(10,'(100(i6))') (jqi(j),j=iqi(i),iqi(i+1)-1)
!enddo
!stop

!-----+----------------------------------------------------------------+
!  Build ice stresses

!     ggao 11-17-2006

!      ARR = ART
!     ggao 11-17-2006

      do ll=1,nne
       if (isicec(ll)==1) then

         i1=in(1,ll)          
         i2=in(2,ll)          
         i3=in(3,ll)          
!         ARR = 0.5_SP*(dx(2,ll)*dy(3,ll)-dx(3,ll)*dy(2,ll))
!     ggao 11-17-2006

         do i=1,3
            ii=in(i,ll)
            do j=1,3
               ij=in(j,ll)

               px =   dy(i,ll)
               py = - dx(i,ll)*cs(ii)

               tx = (zetac(ll)+etac(ll))*dy(j,ll)
               ty = -          etac(ll) *dx(j,ll)*cs(ii)

               tmp=( tx*px+ty*py )/4.0_SP/arr(ll)/cs(ii)

               call calja(ii,ij,iqi,jqi,ind)
               usolve(ind)=usolve(ind) + tmp

! if (ngid(ii)==333.and.ll<=N.and.ii<=M.and.ij<=M) write(*,*) myid,egid(ll),ngid(ii),ngid(ij),tmp
! if (ngid(ii)==333.and.ll>N.and.ii<=M.and.ij<=M) write(*,*) myid,he_lst(ll-N),ngid(ii),ngid(ij),tmp
! if (ngid(ii)==333.and.ll<=N.and.ii>M.and.ij<=M) write(*,*) myid,egid(ll),hn_lst(ii-M),ngid(ij),tmp
! if (ngid(ii)==333.and.ll>N.and.ii>M.and.ij<=M) write(*,*) myid,he_lst(ll-N),hn_lst(ii-M),ngid(ij),tmp

! if (ngid(ii)==333.and.ll<=N.and.ii<=M.and.ij>M) write(*,*) myid,egid(ll),ngid(ii),hn_lst(ij-M),tmp
! if (ngid(ii)==333.and.ll>N.and.ii<=M.and.ij>M) write(*,*) myid,he_lst(ll-N),ngid(ii),hn_lst(ij-M),tmp
! if (ngid(ii)==333.and.ll<=N.and.ii>M.and.ij>M) write(*,*) myid,egid(ll),hn_lst(ii-M),hn_lst(ij-M),tmp
! if (ngid(ii)==333.and.ll>N.and.ii>M.and.ij>M) write(*,*) myid,he_lst(ll-N),hn_lst(ii-M),hn_lst(ij-M),tmp

               ty =             etac(ll) *dy(j,ll)
               tx = -(zetac(ll)-etac(ll))*dx(j,ll)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(ll)/cs(ii)
               call calja(ii,nn+ij,iqi,jqi,ind)
               usolve(ind)=usolve(ind) + tmp

               ty =  (zetac(ll)-etac(ll))*dy(j,ll)
               tx = -           etac(ll) *dx(j,ll)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(ll)/cs(ii)
               call calja(nn+ii,ij,iqi,jqi,ind)
               usolve(ind)=usolve(ind) + tmp

               tx =             etac(ll) *dy(j,ll)
               ty = -(zetac(ll)+etac(ll))*dx(j,ll)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(ll)/cs(ii)

               call calja(nn+ii,nn+ij,iqi,jqi,ind)
               usolve(ind)=usolve(ind) + tmp
            enddo
         enddo
       endif

      enddo

!-----+----------------------------------------------------------------+
! 1-zeroes for velocities outside the ice pack and along the solid bdy
! where no slip will be applied
! 2-add the left hand side (time tendency + friction)
!-----+----------------------------------------------------------------+
      do i=1,nn

        do j=iaa(i),iaa(i+1)-1
          ij=ja(j)
          if (isice(ij)==0 .or. isonb(ij)==1) then
             call calja(i,ij,iqi,jqi,ind)
             usolve(ind)  =0.0_SP
             call calja(nn+i,nn+ij,iqi,jqi,ind)
             usolve(ind)  =0.0_SP
             call calja(i,nn+ij,iqi,jqi,ind)
             usolve(ind)  =0.0_SP
             call calja(nn+i,ij,iqi,jqi,ind)
             usolve(ind)  =0.0_SP
             call calja(nn+i,nn+ij,iqi,jqi,ind)
             usolve(ind)  =0.0_SP
          endif
        enddo

!        call calja(i,i,iaa,ja,ind)
        if (isice(i)==0 .or. isonb(i)==1) then
          do j=iqi(i),iqi(i+1)-1
            usolve(j)=0.0_SP
          enddo
          call calja(i,i,iqi,jqi,ind)
          usolve(ind)=sv(i)

          do j=iqi(nn+i),iqi(nn+i+1)-1
            usolve(j)=0.0_SP
          enddo
          call calja(nn+i,nn+i,iqi,jqi,ind)
          usolve(ind)=sv(i)

        else
          call calja(i,i,iqi,jqi,ind)
          usolve(ind)=usolve(ind)+real(lhs(i))*sv(i)
          call calja(nn+i,nn+i,iqi,jqi,ind)
          usolve(ind)=usolve(ind)+real(lhs(i))*sv(i)
        endif
      enddo

!do i=1,M
! if (ngid(i)==333) write(*,*) myid,i,(usolve(j),j=iqi(i),iqi(i+1)-1), &
!(usolve(j),j=iqi(i+nn),iqi(i+nn+1)-1)
!enddo
!call pstop
!-----+----------------------------------------------------------------+
! iterative viscous solver
!-----+----------------------------------------------------------------+

    do i=1,nn
      fc=aimag(lhs(i))*sv(i)
      if (isice(i)==0 .or. isonb(i)==1) then
       utmp(i)   =0.0_SP
       utmp(nn+i)=0.0_SP
      else
       utmp(i)   =real (rhs(i))*sv(i) + fc * vimid(i)
       utmp(nn+i)=aimag(rhs(i))*sv(i) - fc * uimid(i)
      endif
    enddo

# if defined (MULTIPROCESSOR)
 if (PAR) then
      call pcgm_para(usolve,utmp,2*nn,iqi,jqi,rmask,msr,ice_eps,ice_kmax)
 else
      call pcgm(usolve,utmp,2*nn,iqi,jqi,ice_eps,ice_kmax)
 endif
#else
      call pcgm(usolve,utmp,2*nn,iqi,jqi,ice_eps,ice_kmax)
#endif

      uinew(1:nn)=utmp(1:nn)
      vinew(1:nn)=utmp(1+nn:2*nn)

      Fvis=-(rhs+aimag(lhs)*(vimid-eye*uimid) &
              -real (lhs)*(uinew+eye*vinew))

      rhs=rhs+Fvis
      

      return
# endif
      end



!-----+----------------------------------------------------------------+
!-----+----------------------------------------------------------------+
      subroutine icestress(lhs,rhs,etac,zetac)
# if defined (ICE)
  use ice_variables

      implicit  none

      integer i,l,ii,i1,i2,i3
      REAL(SP) :: zeta,eta,e11,e12,e22,s11,s12,s22,e1,e2,delta
      REAL(SP) :: arrcs,uel,vel,csel,pel,hel,ael,tmp !,arr
      REAL(SP) :: px,py,tx,ty,tmp1,tmp2,nx,ny,fc,lc
      integer j,ij,face
      REAL(SP), dimension(nne) :: etac, zetac
      complex (SP), dimension (nn) :: lhs, rhs, Fvis

!-----+----------------------------------------------------------------+
!  Initialization

      Fvis=eye*0.0_SP

      do l=1,nne
       if (isicec(l)==1) then

         i1=in(1,l)          
         i2=in(2,l)          
         i3=in(3,l)          
!         ARR = 0.5_SP*(dx(2,l)*dy(3,l)-dx(3,l)*dy(2,l))
!         ggao 11/13/2006
!          ARR = ART(l)
!         arr=art(l)

         do i=1,3
            ii=in(i,l)
            do j=1,3
               ij=in(j,l)

               px =   dy(i,l)
               py = - dx(i,l)*cs(ii)

               tx = (zetac(l)+etac(l))*dy(j,l)
               ty = -         etac(l) *dx(j,l)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(l)/cs(ii)
               tmp1=     tmp*uice(ij)

               ty =            etac(l) *dy(j,l)
               tx = -(zetac(l)-etac(l))*dx(j,l)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(l)/cs(ii)
               tmp1=tmp1+tmp*vice(ij)

               ty =  (zetac(l)-etac(l))*dy(j,l)
               tx = -          etac(l) *dx(j,l)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(l)/cs(ii)
               tmp2=     tmp*uice(ij)

               tx =            etac(l) *dy(j,l)
               ty = -(zetac(l)+etac(l))*dx(j,l)*cs(ii)
               tmp=( tx*px+ty*py )/4.0_SP/arr(l)/cs(ii)
               tmp2=tmp2+tmp*vice(ij)

               Fvis(ii)=Fvis(ii) - (tmp1 + eye * tmp2)

            enddo
         enddo
       endif

      enddo

rhs=rhs+Fvis/sv

!-----+----------------------------------------------------------------+
      return
# endif
      end



!-----+----------------------------------------------------------------+
!-----+----------------------------------------------------------------+
      subroutine icestress_evp(etac,zetac,rhs)
# if defined (ICE)
  use ice_variables
  use all_vars

      implicit  none
!arguments
      complex (SP) :: rhs(nn)        ! Right hand side of momentum equation
      REAL(SP), dimension(nne) :: etac,zetac
!locals
      integer i,ll,ii
      REAL(SP)    zeta,eta,e11,e12,e22,s11,s12,s22,e1,e2,delta
      REAL(SP)    arrcs,uel,vel,csel,pel,hel,ael,ah,tmp  !,arr
      integer j,ij,k,face,i1,i2
      REAL(SP)    ahi(nn)
      REAL(SP)    atmp,btmp,b1,b2,dtmp, &
                          tmp1,tmp2,px,py,tx,ty,c1dedt,nx,ny
      complex (SP) :: Fvis(nn)        ! Fviscous forces

      c1dedt=1.0_SP/ee/dtice

!-----+----------------------------------------------------------------+
!  advances ice stresses

      

      do ll=1,nne
        csel=(cs(in(1,ll))+cs(in(2,ll))+cs(in(3,ll)))/3.0_SP
        uel=(uice(in(1,ll))+uice(in(2,ll))+uice(in(3,ll)))/3.0_SP
        vel=(vice(in(1,ll))+vice(in(2,ll))+vice(in(3,ll)))/3.0_SP
        pel=(pice(in(1,ll))+pice(in(2,ll))+pice(in(3,ll)))/3.0_SP
!     ggao 11-17-2006
!       ARR = 0.5_SP*(dx(2,ll)*dy(3,ll)-dx(3,ll)*dy(2,ll))
        arrcs=arr(ll)*csel

       if (isicec(ll)==1) then

        e11=0.0_SP
        e12=0.0_SP
        e22=0.0_SP
        do ii=1,3
          i=in(ii,ll)
          e11=e11 + uice(i)*dy(ii,ll) - vel*dx(ii,ll)*cs(i)*cs(i)/csel
          e12=e12 + vice(i)*dy(ii,ll) - uice(i)*dx(ii,ll)*cs(i) &
                                      + uel*dx(ii,ll)*cs(i)*cs(i)/csel
          e22=e22 - vice(i)*dx(ii,ll)*cs(i)
        enddo
        e11=0.50_SP*e11/arrcs
        e12=0.25_SP*e12/arrcs
        e22=0.50_SP*e22/arrcs

        eta  =  etac(ll)
        zeta = zetac(ll)

        atmp = c1dedt + 0.5_SP/eta
        btmp = e12 + c1dedt * stressice(2,ll)
        stressice(2,ll) = btmp/atmp
        stressice(3,ll) = stressice(2,ll)

        btmp = (eta-zeta)/(4.0_SP*zeta*eta)
        atmp = c1dedt + 0.5_SP/eta + btmp
        dtmp = 1.0 /(atmp*atmp-btmp*btmp)

        b1 = ( e11 - pel / zeta * 0.25_SP) + c1dedt * stressice(1,ll)
        b2 = ( e22 - pel / zeta * 0.25_SP) + c1dedt * stressice(4,ll)

        stressice(1,ll) = ( b1 * atmp - b2 * btmp ) * dtmp
        stressice(4,ll) = ( b2 * atmp - b1 * btmp ) * dtmp

       else

        stressice(:,ll) = 0.0_SP

       endif

      enddo

!-----+----------------------------------------------------------------+
! advances velocities

      Fvis = eye * 0.0_SP

      do ll=1,nne
         do i=1,3
            ii=in(i,ll)
            px =   dy(i,ll)
            py = - dx(i,ll)*cs(ii)

            tx = stressice(1,ll)
            ty = stressice(2,ll)
            tmp1=( tx*px+ty*py )

            tx = stressice(3,ll)
            ty = stressice(4,ll)
            tmp2=( tx*px+ty*py )

            Fvis(ii) = Fvis(ii) - (tmp1+eye*tmp2) * 0.5_SP
         enddo
      enddo


!domain integral stuff: add a condition equivalent to d sigma_ij/dn=0

      do ii=1,nfbd
         face = sbd(ii)
         i1 = pseg(1,face)
         i2 = pseg(2,face)
         ll = tseg(1,face)
         nx =   dltyc(face)
         ny = - dltxc(face)

         tx = stressice(1,ll)
         ty = stressice(2,ll)
         tmp1=( tx*nx+ty*ny )

         tx = stressice(3,ll)
         ty = stressice(4,ll)
         tmp2=( tx*nx+ty*ny )

         Fvis(i1) = Fvis(i1) + (tmp1+eye*tmp2) * 0.5_SP ! 0.5 is due to the FE integration
         Fvis(i2) = Fvis(i2) + (tmp1+eye*tmp2) * 0.5_SP
      enddo

      Fvis = Fvis / sv

      rhs = rhs + Fvis

      return
# endif
      end



