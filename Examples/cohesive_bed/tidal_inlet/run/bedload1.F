!==========================================================================
! Given Bed Load Transports, Calculate the flux using divergence 
!   - calculate bedload transport in elements
!   - modify transport using bed slope
!   - exchange at interprocessor boundaries
!   - calculate divergence to give flux out of / into bed
!==========================================================================
  Subroutine  Calc_Bedload_Flux() 

  !use all_vars, only : niec,ntrg,ncv,nt,dltxe,dltye,nprocs,n,msr,par,wubot,wvbot,myid,art1,m
  use all_vars
  use bcs, only : iobcn,i_obc_n
  use mod_obcs, only: next_obc
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  integer  :: i1,ia,ib,ns,i,n1,n2,n3
  real(sp) :: qbl,ootaub_e,taub_e,flux
  real(sp) :: dbdx,dbdy,beta,slopex,slopey,sed_slope,max_slope,bld_dr
  real(sp), parameter :: sed_angle = 33.0  !sed "friction angle" in degrees
  real(sp), allocatable, target :: btx(:) 
  real(sp), allocatable, target :: bty(:) 
  real(sp), parameter :: taub_min = tiny(1.0_sp)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Bedload_Flux "

  !initialize working bedload transport arrays
  allocate(btx(0:nt))
  allocate(bty(0:nt))

  !main loop over sed types (cohesive?)
  do ns = 1,nsed
    sed(ns)%bflx = 0.0
    btx = 0.0
    bty = 0.0

    !loop over elements, calculate bedload transport in x (btx) and  y (bty) at each
    ! => btx (qblx, eq 28) [after loop is kg/m]
    ! => bty (qbly, eq 29) [after loop is kg/m]
    do i=1,n
  
      !magnitude (units stress/density) of bottom shear stress at elements
      taub_e = sqrt(wubot(i)**2 + wvbot(i)**2)
    
      !calculate transport using Meyer-Peter Mueller formulation
      !qbl has units of kg-m^2/s
      qbl = bedload_mpm(sed(ns)%Srho,sed(ns)%Sd50,taub_e)
      if(taub_e > taub_min)then
        ootaub_e = 1./taub_e
      else
        ootaub_e = 0.0
      endif


      !calculate the bedload transport (kg/m)
      !note the negative sign: bottom stress in FVCOM (wubot/wvbot) is stress
      !on the fluid.  So a + current will produce a - bed stress
      !we want this transport to be in the direction of the current, thus the 
      !change of sign
      btx(i) = -DTsed*qbl*wubot(i)*ootaub_e  
      bty(i) = -DTsed*qbl*wvbot(i)*ootaub_e 
  
    end do
  
    !add bed slope contribution to our horizontal fluxes
    !Eq 35, Warner etal, Comp & Geo, 2008 => qbl_slope = slopex,slopey 
    !minimize the maximum downslope, but do not truncate upslope
    !upslope is limited by Lesser formula
    !range of slopefactor (slopex)  = [.55,5.64] for ..
    !..  upslope(dhdx<0)=-infty -> downslope(dhdx>0) = infty
    !dbdx > 0 is a downslope (see bld_dr to enforce this)
    sed_slope = tan(sed_angle*3.14159/180.)
    max_slope = 0.8*sed_slope
    do i=1,n
      n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)

      bld_dr = sign(1.0_sp,btx(i))
      dbdx   = min( (awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdx)
      slopex = sed_slope/( (sed_slope-dbdx)*cos(beta)) 
      btx(i) = btx(i)*slopex

      bld_dr = sign(1.0_sp,bty(i))
      dbdy   = min( (awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdy)
      slopey = sed_slope/( (sed_slope-dbdy)*cos(beta)) 
      bty(i) = bty(i)*slopey
      
    end do
  
    ! multiply by morphological time scale (to speed up time)
!    if(morpho_model)then
       btx = btx*morpho_factor
       bty = bty*morpho_factor
!    endif
  
    ! multiply by bedload transport coefficient (fudge factor) 
     btx = btx*bedload_rate 
     bty = bty*bedload_rate 
  
    ! limit transport to available bed mass
    ! ROMS does this, but doesn't make sense, mass is linked to divergence of  
    ! transport, not horizontal fluxes , check with J. Warner
     
  
    !exchange btx/bty across interprocessor boundaries 
# if defined (MULTIPROCESSOR)
    if(par)then
      call aexchange(ec,myid,nprocs,btx,bty)
    endif
# endif
  
    !calculate the divergence to give flux (kg) on the nodes
    do i=1,ncv
       i1   = ntrg(i)
       ia   = niec(i,1)
       ib   = niec(i,2)
       flux = (-btx(i1)*dltye(i) + bty(i1)*dltxe(i))
       sed(ns)%bflx(ia) = sed(ns)%bflx(ia)-flux
       sed(ns)%bflx(ib) = sed(ns)%bflx(ib)+flux
     end do

     !update the final bedload flux (kg/m^2) on the nodes
     !this flux is + out of the bed
     !if flow dir is positive and bed_stress is increasing in x
     !this should produce a net positive (out of bed) flux
     do i=1,m
       sed(ns)%bflx(i) = sed(ns)%bflx(i)/art1(i)
     end do

     !set the flux on the open boundary
     do i=1,iobcn
       sed(ns)%bflx( i_obc_n(i) ) = sed(ns)%bflx(next_obc(i))
     end do

     !shutdown bedload if node is locally an inactive node
     do i=1,m
       sed(ns)%bflx(i) = sed(ns)%bflx(i)*bottom(i,bflag)
     end do

     !if bathy is rough - choose smooth_bedload, averaging filter
     if(bedload_smooth)then
       call n2e2d(sed(ns)%bflx,btx)
       call e2n2d(btx,sed(ns)%bflx)
     endif
       
  end do !loop over sediment types

  deallocate(btx)
  deallocate(bty)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Bedload_Flux "  

  End Subroutine Calc_Bedload_Flux
