!!$/=============================================================================/
!!$Copyright (c) 2007, The University of Massachusetts Dartmouth 
!!$Produced at the School of Marine Science & Technology 
!!$Marine Ecosystem Dynamics Modeling group
!!$All rights reserved.
!!$
!!$The FVCOM Offline Lagrangian Model has been developed by the joint UMASSD-WHOI
!!$research team.   For details of authorship and attribution of credit please see
!!$the FVCOM technical manual or contact the MEDM group.
!!$
!!$ 
!!$This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu/ The
!!$full copyright notice is contained in the file COPYRIGHT located in the root
!!$directory of the FVCOM code. This original header must be maintained in all
!!$distributed versions.
!!$
!!$THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!!$ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
!!$IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
!!$ARE  DISCLAIMED.  
!!$
!!$/-----------------------------------------------------------------------------/
!!$CVS VERSION INFORMATION
!!$$Id: $
!!$$Name: $
!!$$Revision: $
!!$/=============================================================================/
!
! VERSION 2.0
! created by Martin Huret
! modified by J. Churchill

!==============================================================================|
!   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
!==============================================================================|

MODULE MOD_LAG

  USE MOD_PREC
  IMPLICIT NONE
  SAVE

  TYPE LAG_OBJ          
     INTEGER                         :: NP_OUT     !!NUMBER OF PARTICLES OUTSIDE DOMAIN
     INTEGER,  POINTER, DIMENSION(:) :: ITAG       !!LABEL FOR THE PARTICLE 
     INTEGER,  POINTER, DIMENSION(:) :: HOST       !!ELEMENT CONTAINING PARTICLE
     INTEGER,  POINTER, DIMENSION(:) :: FOUND      !!HOST ELEMENT IS FOUND       
     INTEGER,  POINTER, DIMENSION(:) :: INDOMAIN   !!PARTICLE IS IN THE DOMAIN 

     INTEGER,  POINTER, DIMENSION(:) :: INWATER    !!PARTICLE IS WATER (ADDED JHC 07/07) 
	   
     INTEGER,  POINTER, DIMENSION(:) :: SBOUND     !!HOST ELEMENT HAS A SOLID BOUNDARY NODE   
     REAL(SP), POINTER, DIMENSION(:) :: XP         !!X POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: YP         !!Y POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: ZP         !!SIGMA POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: XPLST      !!LAST X POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: YPLST      !!LAST Y POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: ZPLST      !!LAST SIGMA POSITION OF PARTICLE 
!  ZPTLST ADDED BY JHC 3/07
     REAL(SP), POINTER, DIMENSION(:) :: ZPTLST     !!LAST Z POSITION OF PARTICLE 
!  ZPIN ADDED BY JHC 4/07
     REAL(SP), POINTER, DIMENSION(:) :: ZPIN       !!INPUT VERTICAL POSITION OF PARTICLE 

     REAL(SP), POINTER, DIMENSION(:) :: XPT        !!X ABSOLUTE POSITION OF PARTICLE
     REAL(SP), POINTER, DIMENSION(:) :: YPT        !!Y ABSOLUTE POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: ZPT        !!Z POSITION OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: HP         !!BATHYMETRY AT PARTICLE POSITION
     REAL(SP), POINTER, DIMENSION(:) :: EP         !!FREE SURFACE HEIGHT AT PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: UP         !!U VELOCITY OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: VP         !!V VELOCITY OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: WP         !!W VELOCITY OF PARTICLE 
     REAL(SP), POINTER, DIMENSION(:) :: TEMP       !!TEMPERATURE AT PARTICLE POSITION
     REAL(SP), POINTER, DIMENSION(:) :: SAL        !!SALINITY AT PARTICLE POSITION
  END TYPE LAG_OBJ

  TYPE(LAG_OBJ)  ::  LAG
  INTEGER        ::  IELAG,ISLAG,INSTP
  INTEGER        ::  TDRIFT
  INTEGER        ::  DTOUT
  INTEGER        ::  ITOUT        !!ITERATION FOR OUTPUT
  CHARACTER(LEN=80)  LAGINI       !!PARTICLES STARTING LOCATIONS

  LOGICAL            P_SIGMA      !!INITIAL PARTICLE POSITION COORDINATE SYSTEM 
  LOGICAL            OUT_SIGMA    !!OUT PARTICLE POSITION COORDINATE SYSTEM 
  LOGICAL            F_DEPTH      !!FIXED DEPTH TRACKING OPTION 

  !-- RANDOM WALK
  INTEGER    ::  IRW
  REAL(SP)   ::  DHOR
  REAL(SP)   ::  DTRW

  !--File Specifiers ------------------------------------------------------------!
  CHARACTER(LEN=80) :: CASENAME    !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
  CHARACTER(LEN=80) :: GEOAREA     !!LETTER ACRONYM SPECIFYING ZONE IDENTITY (MAX 80 CHARS)
  CHARACTER(LEN=80) :: OUTDIR      !!PARENT OUTPUT DIRECTORY
  CHARACTER(LEN=80) :: INPDIR      !!MAIN   INPUT DIRECTORY   
  CHARACTER(LEN=80) :: INFOFILE    !!

  !--File Unit Specifiers -------------------------------------------------------!
  INTEGER  :: IOPAR,IPT,INLAG

  !--Parameters Controlling Time/Time Stepping-----------------------------------!
  INTEGER  :: YEARLAG,MONTHLAG,DAYLAG,HOURLAG !! Starting date of the tracking
  INTEGER  :: IINT           !!TIME STEP ITERATION NUMBER
  REAL(SP) :: DTI

  ! -- Parameters for specifying start time and input file name (JHC 3/07)

  CHARACTER(LEN=4) :: YEAR0
  REAL(SP) timein
  INTEGER JD19950
END MODULE MOD_LAG

!==============================================================================|
!   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
!==============================================================================|

MODULE LIMS
  USE MOD_PREC
  IMPLICIT NONE
  SAVE

  INTEGER N                  !!NUMBER OF ELEMENTS 
  INTEGER M                  !!NUMBER OF NODES
  INTEGER NDRFT              !!NUMBER OF LAGRANGIAN TRACKING PARTICLES
  INTEGER KB                 !!NUMBER OF SIGMA LEVELS
  INTEGER KBM1               !!NUMBER OF SIGMA LEVELS-1
  INTEGER KBM2               !!NUMBER OF SIGMA LEVELS-2
  INTEGER NE                 !!NUMBER OF UNIQUE EDGES
  INTEGER MX_NBR_ELEM        !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE

END MODULE LIMS

!==============================================================================|
!  VARIABLES                                                                   |
!==============================================================================|

MODULE ALL_VARS

  USE MOD_PREC
  USE LIMS
  IMPLICIT NONE
  SAVE

  !--Constants-------------------------------------------------------------------!

  REAL(SP), PARAMETER :: GRAV      = 9.81_SP
  REAL(SP), PARAMETER :: PI        = 3.141592653_SP
  REAL(SP), PARAMETER :: PI2       = 6.283185307_SP
  REAL(SP), PARAMETER :: ZERO      = 0.0_SP 
  REAL(SP), PARAMETER :: ONE_THIRD = 1.0_SP/3.0_SP 
  REAL(SP), PARAMETER :: REARTH    = 6371.0E03_SP   !!Earth Radius in Meters
  REAL(SP), PARAMETER :: DEG2RAD   = PI2/360.0_SP   !!Radians/Degree
  REAL(SP), PARAMETER :: TPI       = DEG2RAD*REARTH !TPI=pi*rearth/180.=3.14159265/180.0*6371.*1000.

  !--------------------------Grid Metrics---------------------------------------------!

  REAL(SP)              :: VXMIN,VYMIN,VXMAX,VYMAX
  REAL(SP), ALLOCATABLE :: XC(:)               !!X-COORD AT FACE CENTER 
  REAL(SP), ALLOCATABLE :: YC(:)               !!Y-COORD AT FACE CENTER
  REAL(SP), ALLOCATABLE :: VX(:)               !!X-COORD AT GRID POINT
  REAL(SP), ALLOCATABLE :: VY(:)               !!Y-COORD AT GRID POINT

  !----------------Node, Boundary Condition, and Control Volume-----------------------!

  INTEGER, ALLOCATABLE :: NV(:,:)             !!NODE NUMBERING FOR ELEMENTS
  INTEGER, ALLOCATABLE :: NBE(:,:)            !!INDICES OF ELMNT NEIGHBORS
  INTEGER, ALLOCATABLE :: NTVE(:)         
  INTEGER, ALLOCATABLE :: ISONB(:)            !!NODE MARKER = 0,1,2   
  INTEGER, ALLOCATABLE :: ISBCE(:)     
  INTEGER, ALLOCATABLE :: NBVE(:,:)
  INTEGER, ALLOCATABLE :: NBVT(:,:)

  !----------------2-d arrays for the vertical(sigma) coordinate -------------------------------!

  REAL(SP), ALLOCATABLE :: Z(:,:)                    !!SIGMA COORDINATE VALUE 
  REAL(SP), ALLOCATABLE :: ZZ(:,:)                   !!INTRA LEVEL SIGMA VALUE
  REAL(SP), ALLOCATABLE :: DZ(:,:)                   !!DELTA-SIGMA VALUE
  REAL(SP), ALLOCATABLE :: DZZ(:,:)                  !!DELTA OF INTRA LEVEL SIGMA 
  REAL(SP), ALLOCATABLE :: Z1(:,:)                   !!SIGMA COORDINATE VALUE 
  REAL(SP), ALLOCATABLE :: ZZ1(:,:)                  !!INTRA LEVEL SIGMA VALUE
  REAL(SP), ALLOCATABLE :: DZ1(:,:)                  !!DELTA-SIGMA VALUE
  REAL(SP), ALLOCATABLE :: DZZ1(:,:)                 !!DELTA OF INTRA LEVEL SIGMA 

  !---------------2-d flow variable arrays at nodes----------------------------------!

  REAL(SP), ALLOCATABLE :: H(:)            !!BATHYMETRIC DEPTH   
  REAL(SP), ALLOCATABLE :: D(:)            !!CURRENT DEPTH   
  REAL(SP), ALLOCATABLE :: EL(:)           !!CURRENT SURFACE ELEVATION
  REAL(SP), ALLOCATABLE :: ET(:)           !!SURFACE ELEVATION AT PREVIOUS TIME STEP

  !---------------- internal mode   arrays-(element based)----------------------------!

  REAL(SP), ALLOCATABLE :: U(:,:)         !X-VELOCITY
  REAL(SP), ALLOCATABLE :: V(:,:)         !Y-VELOCITY
  REAL(SP), ALLOCATABLE :: W(:,:)         !VERTICAL VELOCITY IN SIGMA SYSTEM
  REAL(SP), ALLOCATABLE :: WW(:,:)        !Z-VELOCITY
  REAL(SP), ALLOCATABLE :: UT(:,:)        !X-VELOCITY FROM PREVIOUS TIMESTEP
  REAL(SP), ALLOCATABLE :: VT(:,:)        !Y-VELOCITY FROM PREVIOUS TIMESTEP
  REAL(SP), ALLOCATABLE :: WT(:,:)        !VELOCITY SIGMA FROM PREVIOUS TIMESTEP
  REAL(SP), ALLOCATABLE :: WWT(:,:)       !Z-VELOCITY FROM PREVIOUS TIMESTEP
  REAL(SP), ALLOCATABLE :: RHO(:,:)       !DENSITY AT ELEMENTS

  !-----------------------3d variable arrays-(node based)-----------------------------!

  REAL(SP), ALLOCATABLE :: T1(:,:)         !!TEMPERATURE AT NODES               
  REAL(SP), ALLOCATABLE :: S1(:,:)         !!SALINITY AT NODES               
  REAL(SP), ALLOCATABLE :: RHO1(:,:)       !!DENSITY AT NODES               
  REAL(SP), ALLOCATABLE :: TT1(:,:)        !!TEMPERATURE FROM PREVIOUS TIME
  REAL(SP), ALLOCATABLE :: ST1(:,:)        !!SALINITY FROM PREVIOUS TIME 
  REAL(SP), ALLOCATABLE :: WTS(:,:)        !!VERTICAL VELOCITY IN SIGMA SYSTEM       
  REAL(SP), ALLOCATABLE :: KH(:,:)         !!TURBULENT DIFFUSIVITY

  !------------shape coefficient arrays and control volume metrics--------------------!

  REAL(SP), ALLOCATABLE :: A1U(:,:)      
  REAL(SP), ALLOCATABLE :: A2U(:,:)     
  REAL(SP), ALLOCATABLE :: AWX(:,:)   
  REAL(SP), ALLOCATABLE :: AWY(:,:)  
  REAL(SP), ALLOCATABLE :: AW0(:,:) 

END MODULE ALL_VARS

