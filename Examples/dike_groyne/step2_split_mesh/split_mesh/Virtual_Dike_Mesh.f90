!  Virtual_Dike_Mesh.f90 
!
!  FUNCTIONS:
!  Virtual_Dike_Mesh      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Virtual_Dike_Mesh
!
!  PURPOSE:  Virtual dike mesh generation based on the original mesh without dike.
!
!****************************************************************************

Program Virtual_Dike_Mesh
USE ALL_VARS
USE MOD_THINDAM
USE MOD_SUBS
Implicit None

! Variables


! Body of Virtual_Dike_Mesh
CASENAME='tst'
FNAME1='./tst_no_dike.2dm'
FNAME2='tst_dam.2dm'
NGL= 9250
MGL= 4688

ISSpherical = .FALSE.

CALL READ_MESH
CALL TRIANGLE_GRID_EDGE
CALL READ_NODESTRING
CALL READ_THINDAM
CALL INITIALIZE_THINDAM
CALL DIKE_MESH_GENERATE


End Program Virtual_Dike_Mesh

