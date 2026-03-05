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
!   Input Parameters Which Control the Model Run                               |
!==============================================================================|

SUBROUTINE DATA_RUN            

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_LAG
  USE MOD_INP
  IMPLICIT NONE
  INTEGER  ISCAN
  CHARACTER(LEN=120) :: FNAME, ISTR
  INTEGER I

  !==============================================================================|
  !   READ IN VARIABLES AND SET VALUES                                           |
  !==============================================================================|

  FNAME = "./"//TRIM(CASENAME)//"_run.dat"

  !------------------------------------------------------------------------------|
  !     "INFO FILE"   !!
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(FNAME,"INFOFILE",CVAL = INFOFILE)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INFOFILE: ',ISCAN
     STOP
  END IF
  !-----------------OPEN RUNTIME INFO FILE---------------------------------------!
  IPT = 71
  IF(TRIM(INFOFILE) /= "screen")THEN
     OPEN(IPT, FILE=TRIM(INFOFILE))
  ELSE
     IPT = 6
  END IF
  !------------------------------------------------------------------------------|
  !   TIME STEP (DTI in second) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"DTI",FSCAL = DTI)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DTI: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT TIME STEP OF FLOW FIELDS (INSTP in second) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"INSTP",ISCAL = INSTP)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INSTP: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   OUTPUT TIME STEP (DTOUT in hour) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"DTOUT",ISCAL = DTOUT)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DTOUT: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT TIME STEP OF FLOW FIELDS (TDRIFT in hour) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"TDRIFT",ISCAL = TDRIFT)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING TDRIFT: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT YEAR OF RUN (YEARLAG) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"YEARLAG",ISCAL = YEARLAG)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING YEARLAG: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT MONTH OF RUN (MONTHLAG) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"MONTHLAG",ISCAL = MONTHLAG)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING YEARLAG: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT DAY OF RUN (DAYLAG) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"DAYLAG",ISCAL = DAYLAG)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DAYLAG: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   INPUT HOUR OF RUN (HOURLAG) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"HOURLAG",ISCAL = HOURLAG)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING HOURLAG: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !   "P_SIGMA" TURNS ON INPUT VERTICAL LOCATION OF PARTICLES IN SIGMA 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"P_SIGMA",LVAL = P_SIGMA)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING P_SIGMA: ',ISCAN
     IF(ISCAN == -2)THEN
        WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     STOP
  END IF
  !------------------------------------------------------------------------------|
  !   "OUT_SIGMA" TURNS ON OUTPUT VERTICAL LOCATION OF PARTICLES IN SIGMA 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"OUT_SIGMA",LVAL = OUT_SIGMA)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OUT_SIGMA: ',ISCAN
     IF(ISCAN == -2)THEN
        WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     STOP
  END IF
  !------------------------------------------------------------------------------|
  !   "F_DEPTH" KEEP SAME Z DEPTH ALONG THE TRACKING
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"F_DEPTH",LVAL = F_DEPTH)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING F_DEPTH: ',ISCAN
     IF(ISCAN == -2)THEN
        WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     STOP
  END IF
  !------------------------------------------------------------------------------|
  !  RANDOM WALK CHOICE (IRW)
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"IRW",ISCAL = IRW)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING IRW: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  ! HORIZONTAL DIFFUSION COEFFICIENT (DHOR) 
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"DHOR",FSCAL = DHOR)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DHOR: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !  RANDOM WALK TIME STEP (DTRW)
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(TRIM(FNAME),"DTRW",FSCAL = DTRW)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING DTRW: ',ISCAN
     STOP 
  END IF
  !------------------------------------------------------------------------------|
  !     "GEOAREA"   !!DIRECTORY FOR INPUT FILES            
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(FNAME,"GEOAREA",CVAL = GEOAREA)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING GEOAREA: ',ISCAN
     STOP 
  END IF
  I = LEN_TRIM(GEOAREA)
  IF(GEOAREA(I:I) == "/") GEOAREA(I:I) = " "
  !------------------------------------------------------------------------------|
  !     "INPDIR"   !!DIRECTORY FOR INPUT FILES            
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(FNAME,"INPDIR",CVAL = INPDIR)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING INPDIR: ',ISCAN
     STOP 
  END IF
  I = LEN_TRIM(INPDIR)
  IF(INPDIR(I:I) == "/") INPDIR(I:I) = " "
  !------------------------------------------------------------------------------|
  !     "LAGINI"   !!DIRECTORY FOR INPUT FILES            
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(FNAME,"LAGINI",CVAL = LAGINI)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING LAGINI: ',ISCAN
     STOP 
  END IF
  I = LEN_TRIM(LAGINI)
  IF(LAGINI(I:I) == "/") LAGINI(I:I) = " "
  !------------------------------------------------------------------------------|
  !     "OUTDIR"   !!
  !------------------------------------------------------------------------------|
  ISCAN = SCAN_FILE(FNAME,"OUTDIR",CVAL = OUTDIR)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OUTDIR: ',ISCAN
     STOP 
  END IF
  I = LEN_TRIM(OUTDIR)
  IF(OUTDIR(I:I) == "/") OUTDIR(I:I) = " "
  !==============================================================================|
  !            SET THE UNIT VALUES FOR INPUT OUTPUT FILES                        !
  !==============================================================================|
  IOPAR=11
  INLAG=13

  RETURN
END SUBROUTINE DATA_RUN
!------------------------------------------------------------------------------|
