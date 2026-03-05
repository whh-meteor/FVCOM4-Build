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
!==============================================================================|
!    Allocate and Initialize Most Arrays                                       !
!==============================================================================|

SUBROUTINE ALLOC_VARS

  !==============================================================================!
  USE ALL_VARS
  IMPLICIT NONE
  INTEGER :: NCT 
  !==============================================================================!

  NCT = N*3

  !==============================================================================!
  !  ALLOCATE:                                                                   !
  !==============================================================================!

  !--------------------------Grid Metrics---------------------------------------------!

  ALLOCATE(XC(0:N))            ;XC   = ZERO   !!X-COORD AT FACE CENTER 
  ALLOCATE(YC(0:N))            ;YC   = ZERO   !!Y-COORD AT FACE CENTER
  ALLOCATE(VX(0:M))            ;VX   = ZERO   !!X-COORD AT GRID POINT
  ALLOCATE(VY(0:M))            ;VY   = ZERO   !!Y-COORD AT GRID POINT

  !----------------Node, Boundary Condition, and Control Volume-----------------------!

  ALLOCATE(NV(0:N,4))           ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
  ALLOCATE(NBE(0:N,3))          ;NBE      = 0  !!INDICES OF ELEMENT NEIGHBORS
  ALLOCATE(NTVE(0:M))           ;NTVE     = 0 
  ALLOCATE(ISONB(0:M))          ;ISONB    = 0  !!NODE MARKER = 0,1,2
  ALLOCATE(ISBCE(0:N))          ;ISBCE    = 0 

  !----------------2-d arrays for the vertical coordinate -------------------------------!

  ALLOCATE(Z(0:M,KB))           ; Z      = ZERO    !!SIGMA COORDINATE VALUE 
  ALLOCATE(ZZ(0:M,KB))          ; ZZ     = ZERO    !!INTRA LEVEL SIGMA VALUE
  ALLOCATE(DZ(0:M,KB))          ; DZ     = ZERO    !!DELTA-SIGMA VALUE
  ALLOCATE(DZZ(0:M,KB))         ; DZZ    = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
  ALLOCATE(Z1(0:N,KB))          ; Z1     = ZERO    !!SIGMA COORDINATE VALUE 
  ALLOCATE(ZZ1(0:N,KB))         ; ZZ1    = ZERO    !!INTRA LEVEL SIGMA VALUE
  ALLOCATE(DZ1(0:N,KB))         ; DZ1    = ZERO    !!DELTA-SIGMA VALUE
  ALLOCATE(DZZ1(0:N,KB))        ; DZZ1   = ZERO    !!DELTA OF INTRA LEVEL SIGMA 

  !---------------2-d flow variable arrays at nodes----------------------------------!

  ALLOCATE(H(0:M))             ;H    = ZERO       !!BATHYMETRIC DEPTH   
  ALLOCATE(D(0:M))             ;D    = ZERO       !!DEPTH   
  ALLOCATE(EL(0:M))            ;EL   = ZERO       !!SURFACE ELEVATION
  ALLOCATE(ET(0:M))           ;ET  = ZERO       !!SURFACE ELEVATION PREVIOUS TIMESTEP

  !---------------- internal mode   arrays-(element based)----------------------------!

  ALLOCATE(U(0:N,KB))          ;U     = ZERO   !!X-VELOCITY
  ALLOCATE(V(0:N,KB))          ;V     = ZERO   !!Y-VELOCITY
  ALLOCATE(W(0:N,KB))          ;W     = ZERO   !!VERTICAL VELOCITY IN SIGMA SYSTEM
  ALLOCATE(WW(0:N,KB))         ;WW    = ZERO   !!Z-VELOCITY
  ALLOCATE(UT(0:N,KB))         ;UT    = ZERO   !!X-VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(VT(0:N,KB))         ;VT    = ZERO   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(WT(0:N,KB))         ;WT    = ZERO   !!VERTICAL VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(WWT(0:N,KB))        ;WWT   = ZERO   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
!  ALLOCATE(KH(0:N,KB))         ;KH    = ZERO   !!TURBULENT QUANTITY

  !-----------------------3d variable arrays-(node based)-----------------------------!

  ALLOCATE(T1(0:M,KB))         ;T1     = ZERO  !!TEMPERATURE AT NODES               
  ALLOCATE(S1(0:M,KB))         ;S1     = ZERO  !!SALINITY AT NODES               
  ALLOCATE(RHO1(0:M,KB))       ;RHO1   = ZERO  !!DENSITY AT NODES               
  ALLOCATE(TT1(0:M,KB))        ;TT1    = ZERO  !!TEMPERATURE FROM PREVIOUS TIME
  ALLOCATE(ST1(0:M,KB))        ;ST1    = ZERO  !!SALINITY FROM PREVIOUS TIME 
  ALLOCATE(WTS(0:M,KB))        ;WTS    = ZERO  !!VERTICAL VELOCITY IN SIGMA SYSTEM      
  ALLOCATE(KH(0:M,KB))         ;KH    = ZERO   !!TURBULENT QUANTITY

  !------------shape coefficient arrays and control volume metrics--------------------!

  ALLOCATE(A1U(0:N,4))         ;A1U   = ZERO
  ALLOCATE(A2U(0:N,4))         ;A2U   = ZERO 
  ALLOCATE(AWX(0:N,3))         ;AWX   = ZERO 
  ALLOCATE(AWY(0:N,3))         ;AWY   = ZERO 
  ALLOCATE(AW0(0:N,3))         ;AW0   = ZERO 

  RETURN
END SUBROUTINE ALLOC_VARS
!==============================================================================|
