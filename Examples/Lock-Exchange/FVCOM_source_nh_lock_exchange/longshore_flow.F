! $Id: longshore_flow.F,v 1.1.2.5 2008/04/03 16:36:02 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.5 $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETUP THE METRICS FOR LONGSHORE WIND DRIVEN FLOW
! CALCULATE BOUNDARY ANGLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine longshore_flow
  USE MOD_UTILS
  USE ALL_VARS
  USE MOD_PAR
  IMPLICIT NONE
  INTEGER :: I,J,K,nxt, prv, idx, ndx
  real(sp) :: dx, dy
  

  ! MAKE SURE THERE ARE ATLEAST TWO NODES IN THE LIST
  IF (NOBCLSF == 0 .and. PAR) THEN
     ! DO NOTHING AND RETURN
     RETURN
  ELSE IF (nobclsf < 2 .and. SERIAL) THEN
     CALL FATAL_ERROR&
          &('There are less than two long shore flow nodes',&
          &'You must change the _lsf.dat file to contain atleast two nodes!.')
  END IF

  !MAKE SURE THE FIRST NODE IN THIS LIST IS NEXT TO A SOLID BOUNDARY
  ! OR A HALO NODE (IF PAR)
  idx = ibclsf(1)
  I = 0
  DO 
     I = I + 1 
     if(I > ntsn(idx)) CALL FATAL_ERROR&
             &("THE FIRST NODE IN THE LSF LIST MUST BE NEXT TO THE SOLID BOUNDARY",&
             &"OR (IN THE PARALLEL CASE) NEXT TO THE HALO?")
     
     ndx = nbsn(idx,I)

     if(ISONB(NDX) == 1 .or. (ndx > M .and. ndx <= MT)) THEN
        exit
     end if
  END DO
     

  ! MAKE SURE ALL THE NODES IN THE LIST ARE ON THE OPEN BOUNDARY
  DO I = 1, nobclsf
     J = ibclsf(I)
     if (isonb(j) /= 2) THEN
        
        write(ipt,*) "LONG SHORE FLOW ERROR ON GLOBAL NODE NUMBER: ", NGID(J)
        CALL FATAL_ERROR&
          &("A LONG SHORE FLOW BOUNDARY NODE IS NOT ON THE OPEN BOUNDARY:", &
          &"CHECK YOUR BOUNDARY FILES!")
     END if

  END DO

  ! MAKE A LIST OF 'NEXT NODES' - used to determine angles and
  ! gradients allong the lsf boundary

  ALLOCATE(NBCLSF(nobclsf))
  DO Idx = 1, nobclsf-1
     prv= idx
     nxt= idx+1

     if(ibclsf(nxt) == ibclsf(prv)) CALL FATAL_ERROR&
          &("Two long shore flow nodes in the file list are the same!")


     I = 0
     DO 
        I = I+ 1
        if(I > ntsn(ibclsf(prv))) CALL FATAL_ERROR&
             &("Two long shore flow nodes are not next to eachother!")

        if( ibclsf(nxt) == nbsn(ibclsf(prv),I) ) exit
     END DO

     nbclsf(idx)=ibclsf(nxt)

  END DO


  ! NOW FIND THE NEXT NODE AFTER THE LAST IN THE LIST - THERE ARE ALWAYS TWO!

  idx = ibclsf(nobclsf)
  IF (nobclsf >1) THEN
     prv = ibclsf(nobclsf - 1)
  ELSE
     ! THIS CAN ONLY HAPPEN IN A PARALLEL CASE WHRE THERE IS ONLY ONE
     ! LOCAL NODE IN THE LSF BOUNDARY DUE TO A NASTY DECOMPOSITION

     ! SET PRV TO THE NEIGHBORING HALO NODE!
     I = 0
     DO 
        I = I + 1 
        if(I > ntsn(idx)) CALL FATAL_ERROR&
             &("CAN'T FIND ANOTHER OPEN BOUNDARY NODE FOR LONG SHORE FLOW?",&
             &"IF THERE IS ONLY ONE NODE IN THE DOMAIN IT MUST BE NEXT TO THE HALO!")

        ndx = nbsn(idx,I)
     
        if(isonb(ndx) == 2 .and. ndx > M ) THEN
           prv = ndx
           exit
        end if
     END DO
  END IF

  I = 0
  DO 
     I = I + 1 
     if(I > ntsn(idx)) CALL FATAL_ERROR&
             &("CAN'T FIND ANOTHER OPEN BOUNDARY NODE FOR LONG SHORE FLOW?",&
             &"The list can't contain the entire open boundary!")

     ndx = nbsn(idx,I)
     
     if(isonb(ndx) == 2 .and. ndx /= prv .and. ndx /= idx ) THEN
        nbclsf(nobclsf)= ndx
        exit
     end if
  END DO

  Allocate(WDF_ANG(nobclsf)); WDF_ANG = 0.0_sp
  Allocate(WDF_DIST(nobclsf)); WDF_DIST = 0.0_sp

  DO I = 1,nobclsf
     dx = VX(nbclsf(I)) - VX(ibclsf(I))
     dy = VY(nbclsf(I)) - VY(ibclsf(I))

     WDF_DIST(I)= sqrt(dx**2 + dy**2)
     WDF_ANG(I) =  ATAN2(dy,dx)
  END DO


  write(ipt,*) "! LONG SHORE FLOW BOUNDARY CONDITION METRICS"
  WRITE(IPT,*) "!     NodeList,   NextNode,   Angle,        Distance"
  DO I = 1, nobclsf
     
     write(ipt,*)"! ",ngid(ibclsf(I)), ngid(nbclsf(i)), wdf_ANG(I),WDF_DIST(I)
     if (nbclsf(i) > M .and. dbg_set(dbg_log)) write(ipt,*)&
          & "! The next node value shown above is in the halo and does not appear co&
          &rrectly in ngid: Don't worry about it ;-)" 

  END DO

end subroutine longshore_flow
