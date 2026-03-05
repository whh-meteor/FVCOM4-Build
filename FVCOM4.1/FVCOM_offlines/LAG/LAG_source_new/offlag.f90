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
PROGRAM PARTICLE_TRAJ

  !==============================================================================!
  !                LAGRANGIAN PARTICLE TRACKING OFF-LINE PROGRAM                 !
  !                                                                              !
  !   CASENAME FOR RUNFILE SHOULD BE ENTERED ON COMMAND LINE                     ! 
  !      IF RUNFILE IS LAG_RUN.DAT THEN USE : ptraj lag for EXAMPLE              !
  !   OUTPUT TYPE (NETCDF OR ASCII) IS CONTROLED IN THE MAKEFILE                 !
  !   EVERYTHING ELSE IS GIVEN IN THE RUNFILE                                    !
  !                                                                              !
  !------------------------------------------------------------------------------!
  ! SET_LAG             :   READ IN INITIAL PARTICLE LOCATIONS                   !
  !                     :   DETERMINE INITIAL HOST ELEMENTS                      !
  ! LAG_UPDATE          :   UPDATE PARTICLE LOCATION AND WRITE DATA              !
  ! TRAJECT             :   USE 4 STAGE RK SCHEME TO TRACK PARTICLE              !
  ! RAND_WALK           :   APPLY THE RANDOM WALK PROCESS                        !
  ! INTERP_V            :   INTERPOLATE VELOCITY AT PARTICLE POSITION            !
  ! INTERP_ELH          :   INTERPOLATE ETA/H AT PARTICLE POSITION               !
  ! INTERP_KH           :   INTERPOLATE KH AND ITS DERIVATIVE AT PART. POSITION  !
  ! FHE_ROBUST          :   FIND ELEMENT CONTAINING PARTICLE (ROBUST)            !
  ! FHE_QUICK           :   FIND ELEMENT CONTAINING PARTICLE (QUICK)             !
  ! ISINTRIANGLE        :   TRUE IF PARTICLE IS IN GIVEN TRIANGLE                !
  ! -----------------------------------------------------------------            ! 
  !==============================================================================!
!
! VERSION 3.0

!  04/2010 IN THIS VERSION THE LAGRANGIAN PARTICLE TRACKING ON SPHERICAL COORDINATES
!          

!  10/29/2009 the input file has been modified to match the standard FVCOM monthly output from Netcdf------------
!             it allows you run for continuour months or year, it detect leap year too. 
!             Please always check the ncdio.F to check the vertical omega velocity output, 
!             in different netcdf : it maybe called omega or wts, it maybe output from 
!             cell(N) or node(M) 

!  04/07 ************** IN THIS VERSION THE CONSTANT DEPTH PARTICLES ARE *************
!        ************** AT A FIXED DISTANCED FROM THE MOVING SEA SURFACE *************
!        **************        AND NOT THE MEAN SEA SURFACE 		 *************

  USE MOD_LAG
  USE LIMS
  IMPLICIT NONE

  !------------------------------------------------------------------------------!
  CHARACTER(LEN=100) :: TEMPSTR, FNAME, NCFILE
  CHARACTER(LEN= 4)  :: YEAR,MONTH_ch
  LOGICAL FEXIST
  
  
  CALL GETARG(1,TEMPSTR)
  IF(LEN_TRIM(TEMPSTR) == 0)THEN
     PRINT*, 'PLEASE PROVIDE CASENAME ON COMMAND LINE'
     PRINT*, 'STOPPING...'
     STOP
  END IF
  CASENAME = ADJUSTL(TEMPSTR)
  FNAME = "./"//TRIM(CASENAME)//"_run.dat"
  INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
  IF(.NOT.FEXIST)THEN
     PRINT*, 'FILE ',FNAME,' DOES NOT EXIST'
     PRINT*, 'STOPPING...'
     STOP
  END IF

  !--READ PARAMETERS CONTROLLING MODEL RUN
  CALL DATA_RUN

  !--SELECT NETCDF FILE TO READ FOR GRID CONSTRUCTION
  WRITE(YEAR,'(I4.4)') YEARLAG
  WRITE(YEAR0,'(I4.4)') YEARLAG  ! initial setup year
  WRITE(MONTH_ch,'(I4.4)') MONTHLAG 

  NCFILE ="./"//TRIM(INPDIR)//"/"//YEAR0//"/"//TRIM(GEOAREA)//"_"//MONTH_ch//".nc"
 
PRINT * , nCFILE

  !--DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
  CALL NCD_READ_GRID(NCFILE)
  WRITE(IPT,*)  '!  # OF NODES            :',M
  WRITE(IPT,*)  '!  # OF ELEMENTS         :',N
  WRITE(IPT,*)  '!  # OF SIGMA LAYERS     :',KBM1
  WRITE(IPT,*) 
  WRITE(IPT,*)  '!      MESH READING      :    FINISHED'

  !--ALLOCATE VARIABLES
  CALL ALLOC_VARS

  !--SET UP GRID METRICS 
  CALL NCD_READ_SHAPE(NCFILE)
  WRITE(IPT,*)  '!  GRID METRICS READING  :    FINISHED'
  CALL TRIANGLE_GRID_EDGE
  WRITE(IPT,*)  '!  SET DOMAIN            :    COMPLETED'

  !--PARTICLE TRACKING
  CALL SET_LAG
  
  CALL LAG_UPDATE

END PROGRAM PARTICLE_TRAJ

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE SET_LAG
  !------------------------------------------------------------------------------|
  !  READ IN LAGRANGIAN CONTROL PARAMETERS AND INITIAL LAGRANGIAN POSITIONS      |
  !------------------------------------------------------------------------------|

  USE MOD_LAG
  USE ALL_VARS
  IMPLICIT NONE

  !------------------------------------------------------------------------------!
  LOGICAL FEXIST
  CHARACTER(LEN=100) :: FNAME,NCFILE
  CHARACTER(LEN=4)   :: TMPSTR
  INTEGER I
  LOGICAL ALL_FOUND
  REAL(SP) :: LAG_TIME
  REAL(SP), ALLOCATABLE, DIMENSION(:) :: ZTEMP

  CHARACTER(LEN=100) :: FILENUMBER,PREFF
  CHARACTER(LEN=4)   :: YEAR,MONTH_ch
  INTEGER  :: HOUR,DAY
  REAL(SP), DIMENSION(0:N,KB) :: UNC       !!U VELOCITY FIELD 
  REAL(SP), DIMENSION(0:N,KB) :: VNC       !!V VELOCITY FIELD
  REAL(SP), DIMENSION(0:N,KB) :: WNC       !!W VELOCITY FIELD 
  REAL(SP), DIMENSION(0:M,KB) :: KHNC      !!KH  FIELD
  REAL(SP), DIMENSION(0:M)    :: ELNC      !!FREE SURFACE HEIGHT FIELD
  INTEGER JD   !FUNCTION TO COMPUTE JULIAN DAY (JHC 4/07)
  
  !------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------|
  !  READ INITIAL POSITION FILE                                                  |
  !------------------------------------------------------------------------------|


  !----Check for Existence of File-----------------------------------------------|
  FNAME ="./"//TRIM(INPDIR)//"/"//TRIM(LAGINI)//".dat"
  INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
  IF(.NOT.FEXIST)THEN
     WRITE(*,*)'LAGRANGIAN PARTICLE INITIAL POSITION FILE: '
     WRITE(*,*) FNAME,' DOES NOT EXIST'
     WRITE(*,*)'HALTING.....'
     STOP
  END IF


  !----Read Particle Positions---------------------------------------------------|

  OPEN(UNIT=INLAG,FILE=FNAME,FORM='FORMATTED')
  READ(INLAG,*)NDRFT


 ! MOFIFIED BY JHC 04/07 SO THAT THE PARTICLE STARTING VERTICAL POSITION IS INPUT AS 
 ! LAG%ZPIN (RATHER THAN ZTEMP).  LAG%ZPIN IS SAVED FOR LATER USE (e.g. in the fixed
 ! depth option procedure, LAG%ZPIN is used in keeping the particle at its initial 
 ! distance from the sea surface).
 ! 
  ALLOCATE(LAG%ITAG(NDRFT),LAG%XPT(NDRFT),LAG%YPT(NDRFT),ZTEMP(NDRFT),LAG%ZPIN(NDRFT))

  DO I=1,NDRFT
     READ(INLAG,*) LAG%ITAG(I),LAG%XPT(I),LAG%YPT(I),LAG%ZPIN(I)
  END DO
# if defined (SPHERICAL)  
  WHERE(LAG%XPT < 0.0_SP) LAG%XPT=360.0_SP+LAG%XPT
# endif  
  CLOSE(INLAG)
  ZTEMP = LAG%ZPIN
  ! END MODIFICATION

  !----Allocate Lagrangian Tracking Arrays---------------------------------------|
  ALLOCATE( LAG%XP(NDRFT),LAG%YP(NDRFT),LAG%ZP(NDRFT))
  ALLOCATE( LAG%ZPT(NDRFT)   ) ; LAG%ZPT   = 0
  ALLOCATE( LAG%HP(NDRFT)   ) ; LAG%HP   = 0
  ALLOCATE( LAG%EP(NDRFT)   ) ; LAG%EP   = 0
  ALLOCATE( LAG%FOUND(NDRFT)   ) ; LAG%FOUND   = 0
  ALLOCATE( LAG%HOST(NDRFT) ) ; LAG%HOST = 0
  ALLOCATE( LAG%SBOUND(NDRFT) ) ; LAG%SBOUND = 0
  ALLOCATE( LAG%INDOMAIN(NDRFT) ) ; LAG%INDOMAIN = 1
  ALLOCATE( LAG%TEMP(NDRFT),LAG%SAL(NDRFT))
  ALLOCATE( LAG%UP(NDRFT),LAG%VP(NDRFT),LAG%WP(NDRFT))

! modified by JHC 07/07 to include INWATER in the Tracking Arrays
  ALLOCATE(LAG%INWATER(NDRFT)) ; LAG%INWATER = 1

  !----Shift Coordinates to Model Coordinate System------------------------------|
  LAG%XP = LAG%XPT - VXMIN
  LAG%YP = LAG%YPT - VYMIN

  !----Determine Element Containing Each Particle--------------------------------|
  CALL FHE_ROBUST(NDRFT,LAG%XP,LAG%YP,LAG%HOST,LAG%FOUND,LAG%SBOUND, &
       LAG%INDOMAIN,LAG%INWATER)
  WHERE(LAG%FOUND == 0)
     LAG%INDOMAIN = 0
  END WHERE

  !--Convert Coordinates Of Initial Particle Position
  HOUR=HOURLAG
  WRITE(YEAR0,'(i4.4)') YEARLAG

  WRITE(FILENUMBER,'(I4.4)') monthlag !pengfei
  PREFF ="./"//TRIM(INPDIR)//"/"//YEAR0//"/"//trim(GEOAREA)//"_"
  
  NCFILE = TRIM(PREFF)//TRIM(FILENUMBER)//'.nc'
  WRITE(*,*) '1 ', NCFILE 

! ***** MODIFIED BY JHC on 04/07 TO INCLUDE TIME OF DATA READ (timein)******
!       WHERE timein IS IN SEC PAST SOME START TIME

  CALL NCD_READ(NCFILE,UNC,VNC,WNC,KHNC,ELNC,timein,HOUR)
!       Modified by JHC 07/07 to acquire INWATER 
  CALL INTERP_ELH(NDRFT,LAG%HOST,LAG%INDOMAIN,LAG%SBOUND, &
       LAG%XP,LAG%YP,H,ELNC,LAG%HP,LAG%EP,0,LAG%INWATER)

  IF(P_SIGMA)THEN
     LAG%ZP   = ZTEMP   ! sigma should be negative
     LAG%ZPT  = LAG%ZP*(LAG%HP + LAG%EP) + LAG%EP
  ELSE 
!  ******  MODIFIED BY JCH ON 04/07 SO THAT THE INITIAL CARTESIAN Z IS RELATIVE TO THE MOVING SEA SURFACE  
!     LAG%ZPT  = -ZTEMP  ! in case z are input in >0 values (otherwise +ZTEMP)
     LAG%ZPT  = -ZTEMP + LAG%EP	! ADD SEA SURFACE ELEVATION TO INPUT Z COORDINATE OF PARTICLE
     LAG%ZP  = (LAG%ZPT-LAG%EP)/(LAG%HP+LAG%EP)
!   
  ENDIF


  !----Store Current Position Array in Previous Time Level Array-----------------|
  !ALLOCATE(LAG%XPLST(NDRFT),LAG%YPLST(NDRFT),LAG%ZPLST(NDRFT))
   ALLOCATE(LAG%XPLST(NDRFT),LAG%YPLST(NDRFT),LAG%ZPLST(NDRFT),LAG%ZPTLST(NDRFT))
  ! 3/07  JHC ADDED LAG%ZPTLST TO STORE LAST CARTESIAN Z POSITION OF PARTICLES

  LAG%XPLST = LAG%XP
  LAG%YPLST = LAG%YP
  LAG%ZPLST = LAG%ZP;

  LAG%ZPTLST = LAG%ZPT
  !  3/07 JC ADDED TO SET ZPTLST TO THE LAST GEOGR. BASED VALUE OF Z

  !----Count Number of Particles Outside Domain----------------------------------|
  LAG%NP_OUT = NDRFT - SUM(LAG%INDOMAIN)

  !--Calculate Particle Tracking Begin/End Iterations----------------------------|

  ISLAG = (JD(YEARLAG,MONTHLAG,DAYLAG) - JD(YEARLAG,1,1))*24+HOURLAG+1

!  IELAG=ISLAG+TDRIFT
  IELAG=ISLAG+TDRIFT-1 ! 

  !------------------------------------------------------------------------------|
  !  PRINT STATISTICS ON LAGRANGIAN TRACKING TO OUTPUT                           |
  !------------------------------------------------------------------------------|

  WRITE(IPT,*) '!'
  WRITE(IPT,*) '!                LAGRANGIAN TRACKING INFO'
  WRITE(IPT,*) '!'
  WRITE(IPT,*) '!  # LAG POINTS          :',NDRFT              
  WRITE(IPT,*) '!  # LAG START ITERATION :',ISLAG              
  WRITE(IPT,*) '!  # LAG END ITERATION :',IELAG              
  IF(LAG%NP_OUT /= 0)THEN
     WRITE(IPT,*) '!  # PTS OUTSIDE DOMAIN  :',LAG%NP_OUT         
  ELSE
     WRITE(IPT,*) '!  # PTS OUTSIDE DOMAIN  :  NONE'    
  END IF

  !------------------------------------------------------------------------------|
  !  OPEN UP OUTPUT FILES                                                        |
  !------------------------------------------------------------------------------|
  
  ! **ddd** open up diagnostics file as unit 20
  FNAME='./'//TRIM(OUTDIR)//'/'//TRIM(GEOAREA)//'_'//TRIM(CASENAME)//'diag.dat'

  OPEN(UNIT=20,FILE=FNAME,STATUS = 'UNKNOWN')



  !--WRITE INITIAL PARTICLE POSITIONS AND VELOCITIES TO FILE
  LAG_TIME = FLOAT(ISLAG)
  IF(OUT_SIGMA)THEN
     ZTEMP = LAG%ZP
  ELSE
     ZTEMP = LAG%ZPT
  ENDIF

!          MODIFIED BY JHC 04/07 TO INCLUDE SURFACE ELEVATION AND BOTTOM DEPTH IN OUTPUT NCD FILE
!          Modified by JHC 07/07 to include the inwater variable in the output NCD file
#  if defined (OUT_NETCDF)
  ITOUT=1
  FNAME='./'//TRIM(OUTDIR)//'/'//TRIM(GEOAREA)//'_'//TRIM(CASENAME)//'.nc'
!  CALL NCD_WRITE(FNAME,NDRFT,LAG_TIME,LAG%ITAG,LAG%INDOMAIN, &
!       LAG%XPT,LAG%YPT,ZTEMP,LAG%UP,LAG%VP,LAG%WP,ITOUT)
  CALL NCD_WRITE(FNAME,NDRFT,LAG_TIME,LAG%ITAG,LAG%INDOMAIN, &
       LAG%XPT,LAG%YPT,ZTEMP,LAG%UP,LAG%VP,LAG%WP,LAG%EP,LAG%HP,LAG%INWATER,ITOUT)

#  else

  WRITE(TMPSTR,'(i4.4)') ISLAG
  FNAME='./'//TRIM(OUTDIR)//'/'//TRIM(GEOAREA)//'_'//TRIM(CASENAME)//''//TRIM(TMPSTR)//'_out.dat'
  OPEN(IOPAR,FILE=FNAME,STATUS='UNKNOWN') 
  DO I=1, NDRFT
     WRITE(IOPAR,11) LAG_TIME,LAG%ITAG(I),LAG%XPT(I),LAG%YPT(I),ZTEMP(I)
  ENDDO
#  endif

!   END MODIFICATION 

  DEALLOCATE(ZTEMP)


11 FORMAT(1F7.2,I4,2F12.3,F10.3)

END SUBROUTINE SET_LAG

!==============================================================================!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE LAG_UPDATE
  !==============================================================================|
  !  UPDATE PARTICLE POSITIONS, CALCULATE SCALAR FIELDS AND PARTICLE VELOCITIES  |
  !==============================================================================|

  USE MOD_LAG
  USE ALL_VARS
  IMPLICIT NONE
  
  !------------------------------------------------------------------------------!
  REAL(SP), DIMENSION(0:N,KB) :: UNC1       !!U VELOCITY FIELD START OF HOUR
  REAL(SP), DIMENSION(0:N,KB) :: VNC1       !!V VELOCITY FIELD START OF HOUR
  REAL(SP), DIMENSION(0:N,KB) :: WNC1       !!W VELOCITY FIELD START OF HOUR
  REAL(SP), DIMENSION(0:M,KB) :: KHNC1      !!KH  FIELD START OF HOUR
  REAL(SP), DIMENSION(0:M)    :: ELNC1      !!FREE SURFACE HEIGHT FIELD START OF HOUR
  REAL(SP), DIMENSION(0:N,KB) :: UNC2       !!U VELOCITY FIELD END OF HOUR
  REAL(SP), DIMENSION(0:N,KB) :: VNC2       !!V VELOCITY FIELD END OF HOUR
  REAL(SP), DIMENSION(0:N,KB) :: WNC2       !!W VELOCITY FIELD END OF HOUR
  REAL(SP), DIMENSION(0:M,KB) :: KHNC2      !!KH  FIELD END OF HOUR
  REAL(SP), DIMENSION(0:M)    :: ELNC2      !!FREE SURFACE HEIGHT FIELD END OF HOUR
  REAL(SP), DIMENSION(NDRFT) :: ZTEMP
  INTEGER  :: I,NH,I1,I2,IT,INT2,HOUR,DAY,month,yeartime
  REAL(SP) :: TMP1,TMP2
  REAL(SP) :: LAG_TIME
  CHARACTER(LEN=100) :: FNAME,NCFILE,PREFF,FILENUMBER
  CHARACTER(LEN=4)   :: YEAR

  REAL(SP), DIMENSION(2) :: TDUM               !! VARIABLES TO KEEP TRACK OF CPU TIME
  REAL(SP)   :: DELTAT,TOTALT
  REAL(SP)   :: ETIME,DTIME
  INTEGER JD     ! FUNCTION TO COMPUTE JULIAN DAY (JHC 4/07)
  INTEGER hourfeb ! check leap year
  integer,dimension(12) :: month_hour !pengfei
  month_hour(1)=744
  month_hour(2)=672 ! initial
  month_hour(3)=744
  month_hour(4)=720
  month_hour(5)=744
  month_hour(6)=720
  month_hour(7)=744
  month_hour(8)=744
  month_hour(9)=720
  month_hour(10)=744
  month_hour(11)=720
  month_hour(12)=744  


  !--FIRST HOUR READING OF VELOCITY FIELDS FROM NETCDF FILE

  HOUR=HOURLAG+(daylag-1)*24
  MONTH=MONTHLAG !pengfei
!  DAY=DAYLAG
  yeartime=yearlag
  WRITE(YEAR,'(i4.4)') YEARLAG
  month_hour(2)=hourfeb(yeartime) ! check leap year
  WRITE(FILENUMBER,'(I4.4)') month !pengfei
  PREFF ="./"//TRIM(INPDIR)//"/"//YEAR//"/"//trim(GEOAREA)//"_"


  NCFILE = TRIM(PREFF)//TRIM(FILENUMBER)//'.nc'

  WRITE(*,*) NCFILE


! ***** JHC 04/07 INCLUDES TIME OF DATA READ (timein)******
!       WHERE timein IS IN SEC PAST SOME START TIME

  CALL NCD_READ(NCFILE,UNC1,VNC1,WNC1,KHNC1,ELNC1,timein,HOUR)
  HOUR=HOUR+1
!  IF (HOUR > 23) THEN
!     HOUR=0   ! CHANGED 4/07 BY JCH (FROM HOUR=1) SO THAT NCD_READ READS THE FIRST HOUR OF EVERY
!     DAY=DAY+1     ! FILE AND A FULL 24 HOURS ARE READ EACH DAY.
!  ENDIF
  IF (HOUR > month_hour(month)-1) THEN !pengfei
     HOUR=0   
     month=month+1     ! FILE AND A FULL 24 HOURS ARE READ EACH DAY.
     if (month>12) then
        yeartime= yeartime+1 
        month =1
        hour=0
        month_hour(2)=hourfeb(year) !check if leap year
     end if
  ENDIF


  UT=UNC1
  VT=VNC1
  WWT=WNC1
  ET=ELNC1

  !------------------------------------------------------------------------------|
  !   LOOP OVER THE TRACKING PERIOD                                              |
  !------------------------------------------------------------------------------|
  IINT=0

  DO NH=ISLAG,IELAG
    DELTAT = DTIME(TDUM)
    TOTALT = ETIME(TDUM)

     WRITE(IPT,*)NH-ISLAG+1,'/',IELAG-ISLAG+1, 'finished (hours)'
     WRITE(20,*) NH-ISLAG+1,'/',IELAG-ISLAG+1, 'finished (hours)'
		 WRITE(20,22) DELTAT,TOTALT
22       FORMAT(' DELTA-T(CPU)=',F8.2,'; TOTAL-T(CPU)=',F8.2)

     !--HOURLY READING OF VELOCITY FIELDS IN NETCDF FILE




! ***** JHC 04/07 INCLUDES TIME OF DATA READ (timein)******
!       WHERE timein IS IN SEC PAST SOME START TIME
  WRITE(YEAR,'(i4.4)') YEARTIME
  month_hour(2)=hourfeb(yeartime) ! check leap year
  WRITE(FILENUMBER,'(I4.4)') month !pengfei
  PREFF ="./"//TRIM(INPDIR)//"/"//YEAR//"/"//trim(GEOAREA)//"_"




     NCFILE = TRIM(PREFF)//TRIM(FILENUMBER)//'.nc'
     CALL NCD_READ(NCFILE,UNC2,VNC2,WNC2,KHNC2,ELNC2,timein,HOUR)
     write(*,*)NCFILE

     HOUR=HOUR+1
!     IF (HOUR > 23) THEN
!        HOUR=0   ! CHANGED 4/07 BY JCH (FROM HOUR=1) SO THAT NCD_READ READS THE FIRST HOUR OF EVERY
!        DAY=DAY+1    !FILE AND A FULL 24 HOURS ARE READ EACH DAY.
!     ENDIF
  IF (HOUR > month_hour(month)-1) THEN !pengfei
     HOUR=0   
     month=month+1     ! FILE AND A FULL 24 HOURS ARE READ EACH DAY.
     if (month>12) then
        yeartime= yeartime+1
        month =1
        hour=0
        month_hour(2)=hourfeb(year) !check if leap year
     end if
 
 ENDIF

     !------------------------------------------------------------------------------|
     !   LOOP WITHIN ONE HOUR
     !------------------------------------------------------------------------------|
     I1=1
     I2=INT(INSTP/DTI)  ! caution of time step here

     DO IT=I1,I2

        IINT=IINT+1

        !--Time interpolation of physical fields
        TMP2=FLOAT(IT-I1+1)/FLOAT(I2-I1+1)
        TMP1=FLOAT(IT-I2)/FLOAT(I1-1-I2)

        U=TMP1*UNC1+TMP2*UNC2
        V=TMP1*VNC1+TMP2*VNC2
        WW=TMP1*WNC1+TMP2*WNC2
        KH=TMP1*KHNC1+TMP2*KHNC2
        EL=TMP1*ELNC1+TMP2*ELNC2

        !--Transformation from w to omega
		! -- COMMENTED OUT FOR NOW
!        CALL WOMEGA(DTI,WWT,WT,UT,VT,ET,EL)
        WT=WWT
        W=WW

        !--Particle is moving
!
!       MODIFIED BY JHC 04/07 TO ACQUIRE SURFACE ELEVATION AND DEPTH AT THE PARTICLE LOCATIONS
!       Modified by JHC 07/07 to acquire INWATER 
 
!        CALL TRAJECT(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,LAG%ZPT,LAG%HOST,&
!             LAG%INDOMAIN,LAG%SBOUND,UT,U,VT,V,WT,W,H,ET,EL)
        CALL TRAJECT(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,LAG%ZPT,LAG%HOST,&
             LAG%INDOMAIN,LAG%SBOUND,UT,U,VT,V,WT,W,H,ET,EL,LAG%EP,LAG%HP,LAG%INWATER)
        !--Particle is randomly moving
!        IF (IRW >= 1) CALL RAND_WALK(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,&
!             LAG%ZPT,LAG%HOST,LAG%INDOMAIN,LAG%SBOUND,KH,H,EL) 
        IF (IRW >= 1) CALL RAND_WALK(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,&
             LAG%ZPT,LAG%HOST,LAG%INDOMAIN,LAG%SBOUND,KH,H,EL,LAG%EP,LAG%HP,LAG%INWATER) 

!       END MODIFICATION

        !--WRITE DATA TO FILE 
        INT2=DTOUT*3600./DTI    
        IF(MOD(IINT,INT(INT2)) == 0)THEN

           !--CALCULATE PARTICLE VELOCITIES
!           LAG%UP =   100.0_SP*(LAG%XP - LAG%XPLST)/(DTI*DTOUT*3600.)
!           LAG%VP =   100.0_SP*(LAG%YP - LAG%YPLST)/(DTI*DTOUT*3600.)
!           LAG%WP =  1000.0_SP*(LAG%ZP - LAG%ZPLST)/(DTI*DTOUT*3600.)
!  3/07 Modified by JC (what was DTI doing in the denominator, and why use sigma z to cal vet vel?)

#           if defined (SPHERICAL)
            LAG%UP =   100.0_SP*(LAG%XP - LAG%XPLST)*TPI*COS(DEG2RAD*(LAG%YP+LAG%YPLST)*0.5)/(DTOUT*3600.)
            LAG%VP =   100.0_SP*(LAG%YP - LAG%YPLST)*TPI/(DTOUT*3600.)
#           else
            LAG%UP =   100.0_SP*(LAG%XP - LAG%XPLST)/(DTOUT*3600.)
            LAG%VP =   100.0_SP*(LAG%YP - LAG%YPLST)/(DTOUT*3600.)
#           endif
            LAG%WP =  1000.0_SP*(LAG%ZPT - LAG%ZPTLST)/(DTOUT*3600.)


           !--WRITE PARTICLE POSITION AND VELOCITIES TO FILE        
           LAG_TIME=(FLOAT(IINT)*DTI)/(3600.)+FLOAT(ISLAG)
           
!  ****ADDED BY JC TO CHECK THAT TIME WRITTEN TO THE OUTPUT FILE (LAG_TIME) IS 1 HOUR BEHIND
!                    THE TIME OF THE MOST RECENT INPUT FILE READ (timein ; FROM NCD_READ) *******           
		   
!		   IF(ABS(LAG_TIME+1-(timein/3600.))>0.00001)THEN
!		        WRITE(*,*) ' **DIFFERENCE BETWEEN OUTPUT AND INPUT TIMES**'
!				WRITE(*,*) LAG_TIME+1,timein/3600.
!		   ENDIF


           LAG%XPT = LAG%XP + VXMIN
           LAG%YPT = LAG%YP + VYMIN
!           print *, LAG%YPT,LAG%YP, VYMIN !should be mark

           IF(OUT_SIGMA)THEN
              ZTEMP = LAG%ZP
           ELSE
              ZTEMP = LAG%ZPT
           ENDIF
!          MODIFIED BY JHC 04/07 TO INCLUDE SURFACE ELEVATION AND BOTTOM DEPTH IN OUTPUT NCD FILE
!          Modified by JHC 07/07 to include the inwater variable in the output NCD file
#  if defined (OUT_NETCDF)
           ITOUT=ITOUT+1
           FNAME='./'//TRIM(OUTDIR)//'/'//TRIM(GEOAREA)//'_'//TRIM(CASENAME)//'.nc'
  
!           CALL NCD_WRITE(FNAME,NDRFT,LAG_TIME,LAG%ITAG,LAG%INDOMAIN, &
!                LAG%XPT,LAG%YPT,ZTEMP,LAG%UP,LAG%VP,LAG%WP,ITOUT)
           CALL NCD_WRITE(FNAME,NDRFT,LAG_TIME,LAG%ITAG,LAG%INDOMAIN, &
                LAG%XPT,LAG%YPT,ZTEMP,LAG%UP,LAG%VP,LAG%WP,LAG%EP,LAG%HP,LAG%INWATER,ITOUT)

#    else
           DO I=1, NDRFT
              WRITE(IOPAR,11) LAG_TIME,LAG%ITAG(I),LAG%XPT(I),LAG%YPT(I), &
                               ZTEMP(I),LAG%ZP(I),LAG%ZPT(I)
           ENDDO
#  endif

!          END MODIFICATION


           LAG%XPLST = LAG%XP
           LAG%YPLST = LAG%YP
		   LAG%ZPLST = LAG%ZP

           LAG%ZPTLST = LAG%ZPT
!    3/07  JC ADDED THE ABOVE TO SET ZPTLST TO LAST GEOGR BASED Z        

        END IF

        !--TIME STEP UPDATE OF VELOCITY FIELDS
        UT=U
        VT=V
        WWT=W
        ET=EL

     ENDDO

     !--HOURLY UPDATE OF PHYSICAL FIELDS
     UNC1=UNC2
     VNC1=VNC2
     WNC1=WNC2
     ELNC1=ELNC2
     KHNC1=KHNC2

  ENDDO

11 FORMAT(1F7.2,I4,2F12.3,F10.3,2F12.3)
10 FORMAT(10F12.4)

  RETURN
END SUBROUTINE LAG_UPDATE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

!SUBROUTINE TRAJECT(NPTS,DELTAT,PDXN,PDYN,PDZN,PDZNT,HOST,INDOMAIN, & 
!      SBOUND,U1,U2,V1,V2,W1,W2,HL,EL1,EL2)

!  MODIFIED BY JHC 04/07 TO PASS THE SURFACE ELEVATION AND BOTTOM DEPTH AT THE 
!  PARTICLE POSITION
!  Modified by JHC 07/07 to pass INWATER to calling program
SUBROUTINE TRAJECT(NPTS,DELTAT,PDXN,PDYN,PDZN,PDZNT,HOST,INDOMAIN, & 
      SBOUND,U1,U2,V1,V2,W1,W2,HL,EL1,EL2,EP,HP,INWATER)
  !==============================================================================|
  !  INTEGRATE PARTICLE POSITION FROM X0 TO XN USING VELOCITY FIELDS AT TIME     |
  !  T0 AND TIME TN                                                              |
  !  NPTS:     NUMBER OF PARTICLES					         |
  !  DELTAT:   TIME STEP (USUALLY DTI)				                 !
  !  PDXN:     X POSITION AT TIME T0, RETURNS X POSITION AT TIME TN              |
  !  PDYN:     Y POSITION AT TIME T0, RETURNS Y POSITION AT TIME TN              |
  !  PDZN:     SIGMA POSITION AT TIME T0, RETURNS SIGMA POSITION AT TIME TN      |
  !  PDZNT:    Z POSITION AT TIME TN                                             |
  !  HOST:     CURRENT SET OF ELEMENTS CONTAINING PARTICLES                      |
  !  INDOMAIN: CURRENT STATUS OF PARTICLES, UPDATED FOR PARTICLE POSITION AT TN  |
  !  SBOUND:   HOST ON SOLIDE BOUNDARY                                           |
  !  U1/V1/W1: VELOCITY FIELD (U,V,OMEGA) AT TIME T0                             |
  !  U2/V2/W2: VELOCITY FIELD (U,V,OMEGA) AT TIME TN                             |
  !  HL:       BATHYMETRY                                                        | 
  !  EL1,EL2:  FREE SURFACE ELEVATION AT TIME T0 and TN                          !
  !  EP:       FREE SUR. ELEV. AT PARTICLE AT TIME T0 (INPUT) AND TN (RETURN)    !
  !  HP:       BOTTOM DEPTH AT PARTICLE AT TIME T0 (INPUT) AND TN (RETURN)       ! 
  !==============================================================================|

  USE LIMS
  USE MOD_LAG
  USE ALL_VARS

  IMPLICIT NONE

  !------------------------------------------------------------------------------!
  INTEGER,                     INTENT(IN)     :: NPTS
  REAL(SP),                    INTENT(IN)     :: DELTAT
  REAL(SP), DIMENSION(NPTS),   INTENT(INOUT)  :: PDXN,PDYN,PDZN,PDZNT
  INTEGER , DIMENSION(NPTS),   INTENT(INOUT)  :: HOST 
  INTEGER , DIMENSION(NPTS),   INTENT(INOUT)  :: INDOMAIN,SBOUND
  REAL(SP), DIMENSION(0:N,KB), INTENT(IN)     :: U1,U2,V1,V2,W1,W2
  REAL(SP), DIMENSION(0:M),    INTENT(IN)     :: HL,EL1,EL2
  !------------------------------------------------------------------------------|
  !  NS   : NUMBER OF STAGES IN EXPLICIT RUNGA-KUTTA                             |
  !  CHIX : STAGE FUNCTION EVALUATION FOR X-VELOCITY                             | 
  !  PDX  : STAGE PARTICLE X-POSITION                                            |
  !  UL   : STAGE U VELOCITY                                                     |
  !  EPS  : PARAMETER DEFINING DEPTH OF DRY ELEMENT                              |
  !  DMAX : MAXIMIM SIGMA DEPTH                                                  |
  !  A_RK : ERK COEFFICIENTS (A)                                                 |
  !  B_RK : ERK COEFFICIENTS (B)                                                 |
  !  C_RK : ERK_COEFFICIENTS (C)                                                 |
  !------------------------------------------------------------------------------|
  REAL(SP),DIMENSION(NPTS)               :: PDXT,PDYT,PDZT
  INTEGER, DIMENSION(NPTS)               :: INWATER,INWATERI
  INTEGER, DIMENSION(NPTS)               :: FOUND
  INTEGER                                :: NS
  INTEGER,  PARAMETER                    :: MSTAGE = 4
  REAL(SP), DIMENSION(0:N,KB)            :: UL,VL,WL
  REAL(SP), DIMENSION(0:M)               :: ELL
  REAL(SP), DIMENSION(NPTS,0:MSTAGE)     :: CHIX,CHIY,CHIZ 
  REAL(SP), DIMENSION(NPTS)              :: PDX,PDY,PDZ
  REAL(SP), DIMENSION(NPTS)              :: UP,VP,WP,HP,EP
  REAL(SP), PARAMETER                    :: EPS  = 1.0E-5
  REAL(SP), PARAMETER, DIMENSION(MSTAGE) :: A_RK = (/0.0_DP,0.5_DP,0.5_DP,1.0_DP/) 
  REAL(SP), PARAMETER, DIMENSION(MSTAGE) :: B_RK = (/1.0_DP/6.0_DP,1.0_DP/3.0_DP, &
       1.0_DP/3.0_DP,1.0_DP/6.0_DP/) 
  REAL(SP), PARAMETER, DIMENSION(MSTAGE) :: C_RK = (/0.0_DP,0.5_DP,0.5_DP,1.0_DP/) 
  LOGICAL ALL_FOUND
  !------------------------------------------------------------------------------|

  !--Initialize Stage Functional Evaluations
  CHIX = 0.0_SP  
  CHIY = 0.0_SP  
  CHIZ = 0.0_SP  

  !--Loop over RK Stages 
  DO NS=1,MSTAGE

     !!Particle Position at Stage N (x,y,sigma)
#    if defined (SPHERICAL)
!!     PDX(:)  = PDXN(:)  + A_RK(NS)*DELTAT*CHIX(:,NS-1)/(TPI*COS(DEG2RAD*PDYN(:))+1.0E-8)
     PDY(:)  = PDYN(:)  + A_RK(NS)*DELTAT*CHIY(:,NS-1)/TPI
     PDX(:)  = PDXN(:)  + A_RK(NS)*DELTAT*CHIX(:,NS-1)/(TPI*COS(DEG2RAD*(PDYN(:)+PDY(:))*0.5_SP)+1.0E-10)
#    else
     PDX(:)  = PDXN(:)  + A_RK(NS)*DELTAT*CHIX(:,NS-1)
     PDY(:)  = PDYN(:)  + A_RK(NS)*DELTAT*CHIY(:,NS-1)
#    endif
     PDZ(:)  = PDZN(:)  + A_RK(NS)*DELTAT*CHIZ(:,NS-1)
     
!     if(ns == 4)write(100,'(2f20.12,2f10.6)') PDX(1),PDXN(1),CHIX(1,NS-1),CHIY(1,NS-1)

#    if defined (SPHERICAL)
     WHERE(PDX < 0.0_SP)
       PDX = PDX + 360.0_SP
     ELSEWHERE(PDX > 360.0_SP)
       PDX = PDX - 360.0_SP
     END WHERE
             
     WHERE(PDY > 90.0_SP)
       PDY = 180.0_SP - PDY
     ELSEWHERE(PDY < -90.0_SP)
       PDY = -180.0_SP - PDY
     END WHERE
#    endif
             
     !!Adjust Sigma Position to Reflect Off Bottom (Mirroring)
     PDZ = MAX(PDZ,-(2.0+PDZ)) 

     !!Adjust Sigma Position to Remain Below Free Surface
     PDZ = MIN(PDZ,0.0_SP)

     !!Calculate Velocity Field for Stage N Using C_RK Coefficients
     UL  = (1.0_SP-C_RK(NS))*U1 + C_RK(NS)*U2 
     VL  = (1.0_SP-C_RK(NS))*V1 + C_RK(NS)*V2 
     WL  = (1.0_SP-C_RK(NS))*W1 + C_RK(NS)*W2 
     ELL  = (1.0_SP-C_RK(NS))*EL1 + C_RK(NS)*EL2 

     !!Evaluate Velocity (u,v,w) at Stage NS Particle Position
	 !       Modified by JHC 07/07 to acquire INWATER 
     CALL INTERP_V(NPTS,HOST,INDOMAIN,SBOUND,PDX,PDY,PDZ,UL,VL,WL,UP,VP,WP,INWATER)

     !!Evaluate EL/H at Stage NS Particle Position
     CALL INTERP_ELH(NPTS,HOST,INDOMAIN,SBOUND,PDX,PDY,HL,ELL,HP,EP,0,INWATER)

     CHIX(:,NS) = UP(:)
     CHIY(:,NS) = VP(:)
     CHIZ(:,NS) = WP(:)/(HP(:)+EP(:))    !delta_sigma/deltaT = ww/D

     !!Limit vertical motion in very shallow water
     WHERE( (HP + EP) < EPS)
        CHIZ(:,NS) = 0.0_SP
     END WHERE

  END DO

  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  PDXT(:)  = PDXN(:)
  PDYT(:)  = PDYN(:)
  PDZT(:)  = PDZN(:)
  DO NS=1,MSTAGE
#    if defined (SPHERICAL)
!!     PDXT(:) = PDXT(:) + DELTAT*CHIX(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))/(TPI*COS(DEG2RAD*PDYT(:))+1.0E-8)
     PDYT(:) = PDYT(:) + DELTAT*CHIY(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))/TPI
     PDXT(:) = PDXT(:) + DELTAT*CHIX(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))/(TPI*COS(DEG2RAD*(PDYT(:)+PDYN(:))*0.5)+1.0E-10)
#    else
     PDXT(:) = PDXT(:) + DELTAT*CHIX(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))
     PDYT(:) = PDYT(:) + DELTAT*CHIY(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))
#    endif
     PDZT(:) = PDZT(:) + DELTAT*CHIZ(:,NS)*B_RK(NS)*FLOAT(INDOMAIN(:))

#    if defined (SPHERICAL)
     WHERE(PDXT < 0.0_SP)
       PDXT = PDXT + 360.0_SP
     ELSEWHERE(PDXT > 360.0_SP)
       PDXT = PDXT - 360.0_SP
     END WHERE
             
     WHERE(PDYT > 90.0_SP)
       PDYT = 180.0_SP - PDYT
     ELSEWHERE(PDYT < -90.0_SP)
       PDYT = -180.0_SP - PDYT
     END WHERE
#    endif
             
  END DO

  !--Evaluate Temporary Location
  FOUND=0
  INWATER=1
  CALL FHE_QUICK(NPTS,PDXT,PDYT,HOST,FOUND,SBOUND,INDOMAIN,ALL_FOUND)
  IF(.NOT.ALL_FOUND)THEN
      CALL FHE_ROBUST(NPTS,PDXT,PDYT,HOST,FOUND,SBOUND,INDOMAIN,INWATER)
  ENDIF

  !--Update Only Particle Still in Water
  PDXN(:)  = PDXN(:)*(1.0_SP-FLOAT(INWATER(:)))+PDXT(:)*FLOAT(INWATER(:))
  PDYN(:)  = PDYN(:)*(1.0_SP-FLOAT(INWATER(:)))+PDYT(:)*FLOAT(INWATER(:)) 
  PDZN(:)  = PDZT(:)

  !--Adjust Depth of Updated Particle Positions----------------------------------!
  PDZN = MAX(PDZN,-(2.0+PDZN))                 !Reflect off Bottom
  PDZN = MIN(PDZN,0.0_SP)                      !Don t Pierce Free Surface

  !--Evaluate Bathymetry and Free Surface Height at Updated Particle Position----!
  CALL INTERP_ELH(NPTS,HOST,INDOMAIN,SBOUND,PDXN,PDYN,HL,ELL,HP,EP,1,INWATERI)


  !--Sigma adjustment if fixed depth tracking------------------------------------!
  IF(F_DEPTH)THEN
!
!    *********************** MODIFIED BY JHC 04/07 SO THAT **********************
!    ***********************  FIXED DEPTH TRACKING KEEPS   **********************
!    ***********************  PARTICLES AT A FIXED DEPTH   **********************
!    ***********************  FROM THE MOVING SEA SURFACE  **********************
!     PDZN = (PDZNT-EP)/(HP+EP)
      PDZN = (-LAG%ZPIN)/(HP+EP)  !  THIS IS REALLY PDZN = ((-LAG%ZPIN+EP) - EP)/(HP+EP)
	                          !  WHERE ZPIN IS THE SPECIFIED FIXED DEPTH RELATIVE TO THE SS

     PDZN = MAX(PDZN,-1.0_SP) ! Depth can change though if particle goes into shallower areas
  ENDIF

  !--Calculate Particle Location in Cartesian Vertical Coordinate----------------!
  PDZNT = PDZN*(HP+EP)  + EP
  
  RETURN
END SUBROUTINE TRAJECT

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  MODIFIED BY JHC 04/07 TO PASS THE SURFACE ELEVATION AND BOTTOM DEPTH AT THE 
!  PARTICLE POSITION
!  Modified by JHC 07/07 to pass INWATER to calling program
!  Bug found by Mark D. Rowe NOAA 11-13-2014  ---   DO I=1,N => DO I=1,M

!SUBROUTINE RAND_WALK(NPTS,DELTAT,PDXN,PDYN,PDZN,PDZNT,HOST,INDOMAIN,SBOUND,KHL,HL,ELL)
SUBROUTINE RAND_WALK(NPTS,DELTAT,PDXN,PDYN,PDZN,PDZNT,HOST,INDOMAIN,SBOUND,KHL, &
                       HL,ELL,EP,HP,INWATER)

  !==============================================================================|
  !  VERTICAl AND HORIZONTAL RANDOM WALK FOR THE PARTCICLE TRACKING              |
  !  USING VISSER'S METHOD FOR THE VERTICAL                                      |
  !  NPTS:     NUMBER OF PARTICLES					         |
  !  DELTAT:   TIME STEP (USUALLY DTI)				                 |
  !  PDXN:     X POSITION AT TIME T0, RETURNS X POSITION AT TIME TN              |
  !  PDYN:     Y POSITION AT TIME T0, RETURNS Y POSITION AT TIME TN              |
  !  PDZN:     SIGMA POSITION AT TIME T0, RETURNS SIGMA POSITION AT TIME TN      |
  !  PDZNT:    Z POSITION AT TIME TN                                             |
  !  HOST:     CURRENT SET OF ELEMENTS CONTAINING PARTICLES                      |
  !  INDOMAIN: CURRENT STATUS OF PARTICLES, UPDATED FOR PARTICLE POSITION AT TN  |
  !  SBOUND:   HOST ON SOLIDE BOUNDARY                                           |
  !  HL:       BATHYMETRY                                                        | 
  !  ELL:      FREE SURFACE ELEVATION AT TIME TN                                 | 
  !  EP:       FREE SUR. ELEV. AT PARTICLE AT TIME T0 (INPUT) AND TN (RETURN)    !
  !  HP:       BOTTOM DEPTH AT PARTICLE AT TIME T0 (INPUT) AND TN (RETURN)       !                          | 
  !==============================================================================|

  USE ALL_VARS
  USE MOD_LAG
  IMPLICIT NONE 

  !------------------------------------------------------------------------------!
  INTEGER, DIMENSION(NPTS)                 :: FOUND
  INTEGER,   INTENT(IN)                    :: NPTS
  REAL(SP),  INTENT(IN)                    :: DELTAT
  REAL(SP), DIMENSION(NPTS), INTENT(INOUT) :: PDXN,PDYN,PDZN,PDZNT
  INTEGER , DIMENSION(NPTS), INTENT(INOUT) :: HOST,SBOUND
  INTEGER , DIMENSION(NPTS), INTENT(INOUT) :: INDOMAIN
  REAL(SP), DIMENSION(0:M,KB), INTENT(IN)  :: KHL
  REAL(SP), DIMENSION(0:M),    INTENT(IN)  :: HL,ELL
  !------------------------------------------------------------------------|
  INTEGER , DIMENSION(NPTS)    :: INWATER
  INTEGER , DIMENSION(NPTS)    :: PDXT,PDYT,HDIFFX,HDIFFY
  REAL(SP), DIMENSION(NPTS)    :: HP,EP
  REAL(SP), DIMENSION(NPTS)    :: KHR,DKHR
  REAL(SP), DIMENSION(KB)      :: KHSP1,KHSP2,ZSP
  REAL(SP), DIMENSION(0:M,KB)  :: DKHSP,KHSMTH1,KHSMTH2 
  REAL(SP)                     :: ZDIFF,WDIFF   
  INTEGER             :: RD,ITEND,I,IT,IK,IP
  REAL(SP)            :: GASDEV,RAN1,RAND,LB,LT,LBS,LTS
  REAL(SP), PARAMETER :: VAR=ONE_THIRD
  REAL(SP), PARAMETER :: AC=1.0_SP/6.0_SP
  REAL(SP), PARAMETER :: BIG=1.0E30
  LOGICAL ALL_FOUND
  !------------------------------------------------------------------------|

  !==============================================================================|
  !    Spline is recommended by Ross and Sharples (2004), to create 
  !    a continuous and differentiable diffusivity profile, in order
  !    to meet the time step criterion : DT<MIN(1/K'')
  !
  !    Before the splines are computed, the vertical diffusivity
  !    is smoothed with a three point filter [1/6 2/3 1/6] to
  !    remove spurious negative diffusivities in the resulting
  !    spline interpolations. 
  !==============================================================================|
  
  IF(IRW >= 2) THEN    

     IF (KB <= 3) THEN
       PRINT*,'float interp code cant handle'
       PRINT*,'kb.le.3'
       CALL EXIT(1)
     ENDIF
  !--as k increases, sigma decreases, which is not what the cubic spline routine wants.
!     DO IK=1,KB
!        ZSP(KB-IK+1)=Z(IK)
!     ENDDO
     KHSMTH2(:,1)=KHL(:,1) 
     KHSMTH2(:,KB)=KHL(:,KB-1)
     DO IK=2,KB-1
        KHSMTH2(:,IK)=(KHL(:,IK-1)+KHL(:,IK))/2.
     ENDDO
     KHSMTH1(:,1)=KHSMTH2(:,1)
     KHSMTH1(:,KB)=KHSMTH2(:,KB)
     DO IK=2,KB-1
        KHSMTH1(:,IK)=AC*KHSMTH2(:,IK-1)+(1-2*AC)*KHSMTH2(:,IK)+AC*KHSMTH2(:,IK+1)
     ENDDO
!MDR 11-13-2014     DO I=1,N
     DO I=1,M            !KH is on nodes -- Mark D. Rowe NOAA
        DO IK=1,KB
           KHSP1(KB-IK+1)=KHSMTH1(I,IK)
           ZSP(KB-IK+1)=Z1(I,IK)
        ENDDO  
        CALL spline(ZSP,KHSP1,KB,BIG,BIG,KHSP2)
        DO IK=1,KB
           DKHSP(I,IK)=KHSP2(KB-IK+1) 
        ENDDO
     ENDDO
     
  ENDIF

  CALL INTERP_ELH(NPTS,HOST,INDOMAIN,SBOUND,PDXN,PDYN,HL,ELL,HP,EP,0,INWATER)

  RD=-IINT
  ITEND=INT(DELTAT/DTRW)
  DO IT=1,ITEND

     !--Horizontal Random Walk (Not validated)----------------------------------------|
     IF(IRW == 1.OR.IRW == 3)THEN
        DO IP=1,NPTS
           HDIFFX(IP) = GASDEV(RD)*SQRT(2.0*DHOR*DTRW)
           HDIFFY(IP) = GASDEV(RD)*SQRT(2.0*DHOR*DTRW)
        ENDDO
        PDXT = PDXN+HDIFFX*FLOAT(INDOMAIN)
        PDYT = PDYN+HDIFFY*FLOAT(INDOMAIN)
        
        !--Evaluate Temporary Location
        FOUND=0
        INWATER=0
        CALL FHE_QUICK(NPTS,PDXT,PDYT,HOST,FOUND,SBOUND,INDOMAIN,ALL_FOUND)
        IF(.NOT.ALL_FOUND)THEN
           CALL FHE_ROBUST(NPTS,PDXT,PDYT,HOST,FOUND,SBOUND,INDOMAIN,INWATER)
        ENDIF
        
        !--Update Only Particle Still in Water
        PDXN = PDXN*(1.0_SP-FLOAT(INWATER))+PDXT*(FLOAT(INWATER))
        PDYN = PDYN*(1.0_SP-FLOAT(INWATER))+PDYT*(FLOAT(INWATER))

        CALL INTERP_ELH(NPTS,HOST,INDOMAIN,SBOUND,PDXN,PDYN,HL,ELL,HP,EP,1,INWATER)
  
     ENDIF
     
     !--Vertical Random Walk----------------------------------------------------------|
     IF(IRW >= 2) THEN 
        CALL INTERP_KH(NPTS,HOST,INDOMAIN,PDXN,PDYN,PDZN,HP,EP,KHSMTH1,DKHSP,KHR,DKHR)

        DO IP=1,NPTS

           IF(KHR(IP) < 0.0) THEN
              KHR(IP) = 0.0
              DKHR(IP) = 0.0
           ENDIF

           !--Uniform distribution between -1 and 1
           RAND=(RAN1(RD)*2.0)-1.0
           WDIFF = DKHR(IP)*DTRW+RAND*SQRT(2/VAR*KHR(IP)*DTRW)

           !--Need displacement in sigma
           ZDIFF = WDIFF/(HP(IP)+EP(IP)) 
           PDZN(IP) = PDZN(IP)+ZDIFF

           ! Hereafter either either random mixed layer or simple reflection-----------------|

           !--Creating a random mixed layer for the boundary, as in Ross and Sharples,
           !  to avoid accumulation at the boundaries. LBS and LTS in sigma.
           !  LB and LT have to be larger than the distance to the boundary where 
           !  the probability density is not equal to 1
           !  LB (m) random mixed layer for the reflecting boundary issue
           !  LT (m) random mixed layer for the reflecting boundary issue
           LB=1.0_SP
           LT=1.0_SP
           LTS=LT/(HP(IP)+EP(IP))
           LBS=LB/(HP(IP)+EP(IP))
           IF (PDZN(IP) > (-LTS)) PDZN(IP)=-RAN1(RD)*LTS
           IF (PDZN(IP) < (-1+LBS)) PDZN(IP)=-1+RAN1(RD)*LBS

           !--As the spline is natural (K'=0) at the boundaries, no accumulation 
           !  should occur with simple reflection 
           !   PDZN(I) = MIN(PDZN(I),-PDZN(I))
           !   PDZN(I) = MAX(PDZN(I),-(2.0+PDZN(I)))
           IF (PDZN(IP) > 0.0.or.PDZN(IP) < -1.0) THEN
              PRINT*, ' Sigma is out [0;-1]', PDZN(IP)
              STOP
           ENDIF
       ENDDO   
    ENDIF

  ENDDO

  !--Calculate Particle Location in Cartesian Vertical Coordinate----------------!
  PDZNT = PDZN*(HP+EP)  + EP 

  RETURN
END SUBROUTINE RAND_WALK

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  Modified by JHC 07/07 to pass INWATER to calling program


SUBROUTINE INTERP_V(NPTS,HOST,INDOMAIN,SBOUND,XP,YP,ZP,UIN,VIN,WIN,UP,VP,WP,INWATER) 
  !==============================================================================|
  !  Given a set of pts (xp,yp,zp) of size npts, obtain a linear interpolation   |
  !  of the provided velocity field (UIN,VIN,WIN) at these points                |
  !                                                                              |
  !  RETURNS:                                                                    |
  !     UP,VP,WP (Velocity Field at x,y,z)                                       |
  !                                                                              |
  !==============================================================================|

  USE ALL_VARS
  IMPLICIT NONE

  !------------------------------------------------------------------------------!
  INTEGER,  INTENT(IN)                       :: NPTS
  INTEGER,  INTENT(INOUT), DIMENSION(NPTS)   :: HOST,INDOMAIN,SBOUND
  REAL(SP), INTENT(IN), DIMENSION(NPTS)      :: XP,YP,ZP
  REAL(SP), INTENT(IN), DIMENSION(0:N,1:KB)  :: UIN,VIN,WIN
  REAL(SP), INTENT(OUT), DIMENSION(NPTS)     :: UP,VP,WP
  !------------------------------------------------------------------------------!
  INTEGER, DIMENSION(NPTS) :: FOUND,INWATER
  INTEGER  :: NP,I,E1,E2,E3,K1,K2,K
  REAL(SP) :: W01,WX1,WY1,W02,WX2,WY2
  REAL(SP) :: X0C,Y0C 
  REAL(SP) :: DUDX,DUDY,DVDX,DVDY,DWDX,DWDY
  REAL(SP) :: UE01,UE02,VE01,VE02,WE01,WE02
  REAL(SP) :: ZF1,ZF2,DDZ1
  LOGICAL     ALL_FOUND
  REAL(DP) XTMP,XTMP1	    
  !==============================================================================|

  !==============================================================================!
  !  DETERMINE ELEMENT CONTAINING POINT (XP,YP,ZP)                               !
  !==============================================================================!

  INWATER = 1
  FOUND = 0
  CALL FHE_QUICK(NPTS,XP,YP,HOST,FOUND,SBOUND,INDOMAIN,ALL_FOUND)
  IF(.NOT. ALL_FOUND)THEN
      CALL FHE_ROBUST(NPTS,XP,YP,HOST,FOUND,SBOUND,INDOMAIN,INWATER)
  ENDIF

  !===============================================================================!
  !  DETERMINE VELOCITY, BATHYMETRY, AND FREE SURFACE HEIGHT AT (XP,YP,ZP)        !
  !===============================================================================!

  UP = 0.0_SP
  VP = 0.0_SP
  WP = 0.0_SP

  DO NP=1,NPTS
     IF(INDOMAIN(NP) == 0 .OR. INWATER(NP) == 0)CYCLE !!Particle not in Domain
     I  = HOST(NP)
     E1  = NBE(I,1)
     E2  = NBE(I,2)
     E3  = NBE(I,3)

#    if defined (SPHERICAL)
     XTMP  = XP(NP)*TPI-XC(I)*TPI
     XTMP1 = XP(NP)-XC(I)
     IF(XTMP1 > 180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP = 360.0_SP*TPI+XTMP
     END IF		
     X0C =XTMP*COS(DEG2RAD*(YP(NP)+YC(I))*0.5)
     Y0C =(YP(NP)-YC(I))*TPI
#    else
     X0C = XP(NP) - XC(I)
     Y0C = YP(NP) - YC(I)
#    endif
     !----Determine Sigma Layers Above and Below Particle (For U/V/W Velocities)-------!
     IF(ZP(NP) > ZZ1(I,1))THEN     !!Particle Near Surface
        K1  = 1
        K2  = 1
        ZF1 = 1.0_SP 
        ZF2 = 0.0_SP
     ELSE IF(ZP(NP) < ZZ1(I,KBM1)) THEN !!Particle Near Bottom
        K1 = KBM1
        K2 = KBM1
        ZF1 = 1.0_SP
        ZF2 = 0.0_SP
     ELSE
        DO K=2,KBM1
          K1 = K-1
	  K2 = K
	  IF(ZZ1(I,K) < ZP(NP))EXIT
	END DO
	    
        DDZ1 = 0.5*(DZ1(I,K1)+DZ1(I,K2))
        ZF1 = (ZP(NP)-ZZ1(I,K2))/DDZ1
        ZF2 = (ZZ1(I,K1)-ZP(NP))/DDZ1

!        K1  = INT( (ZZ1(I,1)-ZP(NP))/DZ1(I,1) ) + 1
!        K2  = K1 + 1
!        ZF1 = (ZP(NP)-ZZ1(I,K2))/DZ1(I,1)
!       ZF2 = (ZZ1(I,K1)-ZP(NP))/DZ1(I,1)
     END IF

     !----Linear Interpolation of Particle Velocity in Sigma Level Above Particle-----!
     K = K1 

     DUDX = A1U(I,1)*UIN(I,K)+A1U(I,2)*UIN(E1,K)+A1U(I,3)*UIN(E2,K)+A1U(I,4)*UIN(E3,K)
     DUDY = A2U(I,1)*UIN(I,K)+A2U(I,2)*UIN(E1,K)+A2U(I,3)*UIN(E2,K)+A2U(I,4)*UIN(E3,K)
     DVDX = A1U(I,1)*VIN(I,K)+A1U(I,2)*VIN(E1,K)+A1U(I,3)*VIN(E2,K)+A1U(I,4)*VIN(E3,K)
     DVDY = A2U(I,1)*VIN(I,K)+A2U(I,2)*VIN(E1,K)+A2U(I,3)*VIN(E2,K)+A2U(I,4)*VIN(E3,K)
     DWDX = A1U(I,1)*WIN(I,K)+A1U(I,2)*WIN(E1,K)+A1U(I,3)*WIN(E2,K)+A1U(I,4)*WIN(E3,K)
     DWDY = A2U(I,1)*WIN(I,K)+A2U(I,2)*WIN(E1,K)+A2U(I,3)*WIN(E2,K)+A2U(I,4)*WIN(E3,K)
     UE01 = UIN(I,K) + DUDX*X0C + DUDY*Y0C
     VE01 = VIN(I,K) + DVDX*X0C + DVDY*Y0C
     WE01 = WIN(I,K) + DWDX*X0C + DWDY*Y0C

     !----Linear Interpolation of Particle Position in Sigma Level Below Particle-----!
     K = K2 

     DUDX = A1U(I,1)*UIN(I,K)+A1U(I,2)*UIN(E1,K)+A1U(I,3)*UIN(E2,K)+A1U(I,4)*UIN(E3,K)
     DUDY = A2U(I,1)*UIN(I,K)+A2U(I,2)*UIN(E1,K)+A2U(I,3)*UIN(E2,K)+A2U(I,4)*UIN(E3,K)
     DVDX = A1U(I,1)*VIN(I,K)+A1U(I,2)*VIN(E1,K)+A1U(I,3)*VIN(E2,K)+A1U(I,4)*VIN(E3,K)
     DVDY = A2U(I,1)*VIN(I,K)+A2U(I,2)*VIN(E1,K)+A2U(I,3)*VIN(E2,K)+A2U(I,4)*VIN(E3,K)
     DWDX = A1U(I,1)*WIN(I,K)+A1U(I,2)*WIN(E1,K)+A1U(I,3)*WIN(E2,K)+A1U(I,4)*WIN(E3,K)
     DWDY = A2U(I,1)*WIN(I,K)+A2U(I,2)*WIN(E1,K)+A2U(I,3)*WIN(E2,K)+A2U(I,4)*WIN(E3,K)
     UE02 = UIN(I,K) + DUDX*X0C + DUDY*Y0C
     VE02 = VIN(I,K) + DVDX*X0C + DVDY*Y0C
     WE02 = WIN(I,K) + DWDX*X0C + DWDY*Y0C

     !----Interpolate Particle Velocity Between Two Sigma Layers----------------------!

     UP(NP) = UE01*ZF1 + UE02*ZF2
     VP(NP) = VE01*ZF1 + VE02*ZF2 
     WP(NP) = WE01*ZF1 + WE02*ZF2 

  END DO !!LOOP OVER PARTICLES


  RETURN
END SUBROUTINE INTERP_V

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  Modified by JHC 07/07 to pass INWATER to calling program

SUBROUTINE INTERP_ELH(NPTS,HOST,INDOMAIN,SBOUND,XP,YP,HIN,EIN,HP,EP,IFC,INWATER) 
  !==============================================================================|
  !  Given a set of pts (xp,yp) of size npts, obtain a linear interpolation      |
  !  of the provided free surface/bathymetry (HIN,EIN) at these points           |
  !                                                                              |
  !  RETURNS:                                                                    |
  !     HP(bathymetry at xp,yp) and EP (Free Surface Elevation at xp,yp)         |
  !                                                                              |
  !  FLAGS:                                                                      |
  !     IFC  = 0 IF HOST IS KNOWN TO HAVE CORRECT HOST ELEMENTS FOR CURRENT   |
  !              PARTICLE POSITIONS                                              |
  !     IFC  = 1 IF HOST SHOULD BE UPDATED (ADDS COMPUTATIONAL WORK)          |
  !              PARTICLE POSITIONS                                              |
  !                                                                              |
  !==============================================================================|

  USE ALL_VARS
  IMPLICIT NONE
  !------------------------------------------------------------------------------!
  INTEGER,  INTENT(IN)                     :: NPTS,IFC
  INTEGER,  INTENT(INOUT), DIMENSION(NPTS) :: HOST,INDOMAIN,SBOUND
  REAL(SP), INTENT(IN), DIMENSION(NPTS)    :: XP,YP
  REAL(SP), INTENT(IN), DIMENSION(0:M)     :: HIN,EIN
  REAL(SP), INTENT(OUT), DIMENSION(NPTS)   :: HP,EP
  !------------------------------------------------------------------------------!
  INTEGER, DIMENSION(NPTS) :: FOUND,INWATER
  INTEGER  :: IP,I,N1,N2,N3,K1,K2,K
  REAL(SP) :: H0,HX,HY,E0,EX,EY
  REAL(SP) :: X0C,Y0C 
  LOGICAL ALL_FOUND
  REAL(DP) XTMP,XTMP1	    
  !------------------------------------------------------------------------------!

  !===============================================================================!
  !  DETERMINE ELEMENT CONTAINING POINT (XP,YP)                                   !
  !===============================================================================!

  IF(IFC == 1)THEN
     INWATER=1
     FOUND  = 0
     CALL FHE_QUICK(NPTS,XP,YP,HOST,FOUND,SBOUND,INDOMAIN,ALL_FOUND)
     IF(.NOT. ALL_FOUND)THEN
        CALL FHE_ROBUST(NPTS,XP,YP,HOST,FOUND,SBOUND,INDOMAIN,INWATER)
     ENDIF
  END IF

  !===============================================================================!
  !  LINEARLY INTERPOLATE FREE SURFACE HEIGHT AND BATHYMETRY                      !
  !===============================================================================!
  HP = 1.0_SP
  EP = 1.0_SP

  DO IP=1,NPTS
     IF(INDOMAIN(IP) == 0)CYCLE
     I  = HOST(IP)
     N1  = NV(I,1)
     N2  = NV(I,2)
     N3  = NV(I,3)

#    if defined (SPHERICAL)
     XTMP  = XP(IP)*TPI-XC(I)*TPI
     XTMP1 = XP(IP)-XC(I)
     IF(XTMP1 > 180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP = 360.0_SP*TPI+XTMP
     END IF		
     X0C =XTMP*COS(DEG2RAD*(YP(IP)+YC(I))*0.5)
     Y0C =(YP(IP)-YC(I))*TPI
#    else
     X0C = XP(IP) - XC(I)
     Y0C = YP(IP) - YC(I)
#    endif

     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0(I,1)*HIN(N1)+AW0(I,2)*HIN(N2)+AW0(I,3)*HIN(N3)
     HX = AWX(I,1)*HIN(N1)+AWX(I,2)*HIN(N2)+AWX(I,3)*HIN(N3)
     HY = AWY(I,1)*HIN(N1)+AWY(I,2)*HIN(N2)+AWY(I,3)*HIN(N3)
     HP(IP) = H0 + HX*X0C + HY*Y0C

     !----Linear Interpolation of Free Surface Height---------------------------------!
     E0 = AW0(I,1)*EIN(N1)+AW0(I,2)*EIN(N2)+AW0(I,3)*EIN(N3)
     EX = AWX(I,1)*EIN(N1)+AWX(I,2)*EIN(N2)+AWX(I,3)*EIN(N3)
     EY = AWY(I,1)*EIN(N1)+AWY(I,2)*EIN(N2)+AWY(I,3)*EIN(N3)
     EP(IP) = E0 + EX*X0C + EY*Y0C

  END DO !!LOOP OVER PARTICLES


  RETURN
END SUBROUTINE INTERP_ELH

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE INTERP_KH(NPTS,HOST,INDOMAIN,XP,YP,ZP,HP,EP,ZKH,DZKH,KHOUT,DKHOUT)
  !==============================================================================|
  !  Given a set of pts (xp,yp,zp) of size npts, obtain a spline interpolation   |
  !  on the vertical with the provided eddy diffusivity (ZKH) and its derivative |
  !  (DZKH) at grid point points, then linear interpolation on the horizontal    |
  !                                                                              |
  !  RETURNS:                                                                    |
  !             Both dkh/dz (DKHOUT) and kh (KHOUT)                              |
  !                                                                              |
  !==============================================================================|

  USE ALL_VARS
  USE MOD_LAG
  IMPLICIT NONE
  !------------------------------------------------------------------------------!
  INTEGER,  INTENT(IN)                      :: NPTS
  REAL(SP), INTENT(IN), DIMENSION(NPTS)     :: XP,YP,ZP,HP,EP
  INTEGER,  INTENT(IN), DIMENSION(NPTS)     :: HOST,INDOMAIN
  REAL(SP), INTENT(OUT),  DIMENSION(NPTS)   :: DKHOUT, KHOUT
  REAL(SP), INTENT(IN) ,DIMENSION(0:M,KB)   :: DZKH,ZKH
  !---------------------------------------------------------------------------!
  REAL(SP) :: X0C,Y0C,COF1,COF2,COF3
  INTEGER  :: J1,J2,J3,I,IP
  REAL(SP) :: DKHR1,DKHR2,DKHR3,DKHR4
  REAL(SP) :: KHR1,KHR2,KHR3,KHR4
  REAL(SP) :: DDKHR1,DDKHR2,DDKHR3,DDKHR4
  REAL(SP) :: DKHTMP,DZP
  INTEGER, DIMENSION(0:KB+1) :: NZRINDX
  REAL(SP) :: HK,HK2,AK,AK2,AK3,BK,BK2,BK3
  INTEGER  :: KLO,KHI,NZR
  
  REAL(SP) :: Z_TMP(KB)
  REAL(DP) XTMP,XTMP1	    
  !------------------------------------------------------------------------------!

  !===============================================================================!
  !  INTERPOLATE EDDY DIFFUSIVITY AND ITS DERIVATIVE                              !
  !===============================================================================!

  KHOUT  = 0.0_SP
  DKHOUT = 0.0_SP

  NZRINDX(0)=1
  DO I=1,KB
     NZRINDX(I)=I
  ENDDO
  NZRINDX(KB+1)=KB

  DO IP=1,NPTS
     IF(INDOMAIN(IP) == 0)CYCLE  

     I  = HOST(IP)
     J1  = NV(I,1)
     J2  = NV(I,2)
     J3  = NV(I,3)

#    if defined (SPHERICAL)
     XTMP  = XP(IP)*TPI-XC(I)*TPI
     XTMP1 = XP(IP)-XC(I)
     IF(XTMP1 > 180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP = 360.0_SP*TPI+XTMP
     END IF		
     X0C =XTMP*COS(DEG2RAD*(YP(IP)+YC(I))*0.5)
     Y0C =(YP(IP)-YC(I))*TPI
#    else
     X0C = XP(IP) - XC(I)
     Y0C = YP(IP) - YC(I)
#    endif

!     I=HOST(IP)
!     X0C=XP(IP)-XC(I)
!     Y0C=YP(IP)-YC(I)
!     J1=NBE(I,1)
!     J2=NBE(I,2)
!     J3=NBE(I,3)

     !--DERIVATIVE OF THE DIFFUSION ---------------------------------------

     !--find vertical location
     NZR = -ZP(IP)*(KB-1)+1 !guess value for hunt 
     Z_TMP(1:KB)=Z1(I,1:KB)          
     CALL HUNT(Z_TMP,KB,ZP(IP),NZR)     
     
     KHI=NZRINDX(NZR)
     KLO=NZRINDX(NZR+1)

     !--as k in z(k) increases, sigma decreases.
     HK=Z1(I,KHI)-Z1(I,KLO)
     AK=(Z1(I,KHI)-ZP(IP))/HK
     AK2=AK**2
     BK=(ZP(IP)-Z1(I,KLO))/HK
     BK2=BK**2

 !    DKHR1=(-ZKH(I,KLO)+ZKH(I,KHI))/HK+((-3*AK2+1)*DZKH(I,KLO)+(3*BK2-1)*DZKH(I,KHI))*HK/6
 !    DKHR2=(-ZKH(J1,KLO)+ZKH(J1,KHI))/HK+((-3*AK2+1)*DZKH(J1,KLO)+(3*BK2-1)*DZKH(J1,KHI))*HK/6
 !    DKHR3=(-ZKH(J2,KLO)+ZKH(J2,KHI))/HK+((-3*AK2+1)*DZKH(J2,KLO)+(3*BK2-1)*DZKH(J2,KHI))*HK/6
 !    DKHR4=(-ZKH(J3,KLO)+ZKH(J3,KHI))/HK+((-3*AK2+1)*DZKH(J3,KLO)+(3*BK2-1)*DZKH(J3,KHI))*HK/6

 !    COF1=A1U(I,1)*DKHR1 +A1U(I,2)*DKHR2+A1U(I,3)*DKHR3 +A1U(I,4)*DKHR4
 !    COF2=A2U(I,1)*DKHR1 +A2U(I,2)*DKHR2+A2U(I,3)*DKHR3 +A2U(I,4)*DKHR4

 !    DKHTMP=DKHR1+COF1*X0C+COF2*Y0C

     DKHR1=(-ZKH(J1,KLO)+ZKH(J1,KHI))/HK+((-3*AK2+1)*DZKH(J1,KLO)+(3*BK2-1)*DZKH(J1,KHI))*HK/6
     DKHR2=(-ZKH(J2,KLO)+ZKH(J2,KHI))/HK+((-3*AK2+1)*DZKH(J2,KLO)+(3*BK2-1)*DZKH(J2,KHI))*HK/6
     DKHR3=(-ZKH(J3,KLO)+ZKH(J3,KHI))/HK+((-3*AK2+1)*DZKH(J3,KLO)+(3*BK2-1)*DZKH(J3,KHI))*HK/6

     COF1=AW0(I,1)*DKHR1 +AW0(I,2)*DKHR2+AW0(I,3)*DKHR3
     COF2=AWX(I,1)*DKHR1 +AWX(I,2)*DKHR2+AWX(I,3)*DKHR3
     COF3=AWY(I,1)*DKHR1 +AWY(I,2)*DKHR2+AWY(I,3)*DKHR3

     DKHTMP=COF1+COF2*X0C+COF3*Y0C
     !--want the answer to be dkh/dz not dkh/dsigma
     DKHOUT(IP)=DKHTMP/(HP(IP)+EP(IP))

     !--SECOND DERIVATIVE used for time step calculation regarding the time step criterion of Visser
     !     IF (IINT == 1) ZKHZZ=1000.0
     !     DDKHR1=AK*DZKH(I,KLO)+BK*DZKH(I,KHI)
     !     DDKHR2=AK*DZKH(J1,KLO)+BK*DZKH(J1,KHI)
     !     DDKHR3=AK*DZKH(J2,KLO)+BK*DZKH(J2,KHI)
     !     DDKHR4=AK*DZKH(J3,KLO)+BK*DZKH(J3,KHI)

     !     COF1=A1U(I,1)*DDKHR1 +A1U(I,2)*DDKHR2+A1U(I,3)*DDKHR3 +A1U(I,4)*DDKHR4
     !     COF2=A2U(I,1)*DDKHR1 +A2U(I,2)*DDKHR2+A2U(I,3)*DDKHR3 +A2U(I,4)*DDKHR4
     !     
     !     DDKHRTMP=DDKHR1+COF1*X0C+COF2*Y0C
     !     DDKHRTMP=1.0_SP/ABS(DDKHRTMP/(HP(IP)+EP(IP))**2)
     !     IF (DDKHRTMP < ZKHZZ) THEN
     !     ZKHZZ=DDKHRTMP
     !     PRINT*, ZKHZZ, IP,(HP(IP)+EP(IP)), NZR, J1,J2,J3
     !     ENDIF

     !--DIFFUSION ITSELF -----------------------------------------

     !--find z in grid again, as per visser, but in sigma
     DZP = ZP(IP)+0.5*DKHOUT(IP)*DTRW/(HP(IP)+EP(IP))
     !--After adding 0.5*dkhtmp*dtrw, new z can be out of [0;-1]
     DZP= MIN(DZP,ZERO)
     DZP = MAX(DZP,-1.0_SP) 

     !--find vertical location
     CALL HUNT(Z_TMP,KB,DZP,NZR)
     KHI=NZRINDX(NZR)
     KLO=NZRINDX(NZR+1)       

     HK=Z1(I,KHI)-Z1(I,KLO)
     HK2=HK**2
     AK=(Z1(I,KHI)-DZP)/HK
     AK3=AK**3
     BK=(DZP-Z1(I,KLO))/HK
     BK3=BK**3

!     KHR1=AK*(ZKH(I,KLO))+BK*(ZKH(I,KHI))+((AK3-AK)*(DZKH(I,KLO))+(BK3-BK)*(DZKH(I,KHI)))*HK2/6
!     KHR2=AK*(ZKH(J1,KLO))+BK*(ZKH(J1,KHI))+((AK3-AK)*(DZKH(J1,KLO))+(BK3-BK)*(DZKH(J1,KHI)))*HK2/6
!     KHR3=AK*(ZKH(J2,KLO))+BK*(ZKH(J2,KHI))+((AK3-AK)*(DZKH(J2,KLO))+(BK3-BK)*(DZKH(J2,KHI)))*HK2/6
!     KHR4=AK*(ZKH(J3,KLO))+BK*(ZKH(J3,KHI))+((AK3-AK)*(DZKH(J3,KLO))+(BK3-BK)*(DZKH(J3,KHI)))*HK2/6

!     COF1=A1U(I,1)*KHR1 +A1U(I,2)*KHR2+A1U(I,3)*KHR3 +A1U(I,4)*KHR4
!     COF2=A2U(I,1)*KHR1 +A2U(I,2)*KHR2+A2U(I,3)*KHR3 +A2U(I,4)*KHR4

!     KHOUT(IP)=KHR1+COF1*X0C+COF2*Y0C

     KHR1=AK*(ZKH(J1,KLO))+BK*(ZKH(J1,KHI))+((AK3-AK)*(DZKH(J1,KLO))+(BK3-BK)*(DZKH(J1,KHI)))*HK2/6
     KHR2=AK*(ZKH(J2,KLO))+BK*(ZKH(J2,KHI))+((AK3-AK)*(DZKH(J2,KLO))+(BK3-BK)*(DZKH(J2,KHI)))*HK2/6
     KHR3=AK*(ZKH(J3,KLO))+BK*(ZKH(J3,KHI))+((AK3-AK)*(DZKH(J3,KLO))+(BK3-BK)*(DZKH(J3,KHI)))*HK2/6

     COF1=AW0(I,1)*KHR1 +AW0(I,2)*KHR2+AW0(I,3)*KHR3
     COF2=AWX(I,1)*KHR1 +AWX(I,2)*KHR2+AWX(I,3)*KHR3
     COF3=AWY(I,1)*KHR1 +AWY(I,2)*KHR2+AWY(I,3)*KHR3

     KHOUT(IP)=COF1+COF2*X0C+COF3*Y0C

  END DO


  RETURN
END SUBROUTINE INTERP_KH

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|

SUBROUTINE FHE_QUICK(NPTS,XP,YP,HOST,IFOUND,SBOUND,INDOMAIN,ALL_FOUND)
  !==============================================================================|
  !  Determine Which Element A List of Particles Reside in By Searching          |
  !  Neighboring Elements.  Updates HOST component of Lagrangian Particle        |  
  !  Type and updates logical array "ELEM_FOUND" flagging whether the host       |
  !  Has been found							         |
  !==============================================================================|

  USE ALL_VARS
  IMPLICIT NONE
  !------------------------------------------------------------------------------!
  INTEGER,                   INTENT(IN)    :: NPTS
  REAL(SP), DIMENSION(NPTS), INTENT(IN)    :: XP,YP
  INTEGER,  DIMENSION(NPTS), INTENT(IN)    :: INDOMAIN
  INTEGER,  DIMENSION(NPTS), INTENT(INOUT) :: HOST,IFOUND,SBOUND
  LOGICAL, INTENT(OUT) ::  ALL_FOUND
  !------------------------------------------------------------------------------!
  INTEGER I,J,K,ILAST,INEY,NCHECK
  REAL(SP), DIMENSION(3) :: XLAST,YLAST,XNEY,YNEY
  REAL(SP) :: XLAG,YLAG
  LOGICAL ::  ISINTRIANGLE
  !==============================================================================|

  ALL_FOUND = .FALSE.

  DO I=1,NPTS
     IF(INDOMAIN(I) == 0)CYCLE
     XLAG  = XP(I) 
     YLAG  = YP(I) 
     ILAST = HOST(I) 
     XLAST = VX(NV(ILAST,1:3))
     YLAST = VY(NV(ILAST,1:3))
     IF(ISINTRIANGLE(XLAST,YLAST,XLAG,YLAG))THEN      !!PARTICLE REMAINS IN ELEMENT
        IFOUND(I) = 1 
     ELSE                                             !!CHECK NEIGHBORS           
        outer: DO J=1,3
           NCHECK = NV(ILAST,J)
           DO K=1,NTVE(NV(ILAST,J))
              INEY = NBVE(NCHECK,K) 
              XNEY = VX(NV(INEY,1:3))
              YNEY = VY(NV(INEY,1:3))
              IF(ISINTRIANGLE(XNEY,YNEY,XLAG,YLAG))THEN
                 IFOUND(I)  = 1 
                 HOST(I) = INEY
                 IF(ISONB(NV(INEY,1)) == 1.OR.ISONB(NV(INEY,2)) == 1.OR.ISONB(NV(INEY,3)) == 1)THEN
                    SBOUND(I)=1
                 ELSE 
                    SBOUND(I)=0
                 ENDIF
                 EXIT outer
              END IF
           END DO
        END DO outer
     END IF
  END DO

  IF(SUM(IFOUND) == SUM(INDOMAIN)) ALL_FOUND = .TRUE.

  RETURN
END SUBROUTINE FHE_QUICK

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE FHE_ROBUST(NPTS,XP,YP,HOST,IFOUND,SBOUND,INDOMAIN,INWATER)

!*********** Modified by Pengfei 04/07 ******************************************
  !==============================================================================|
  !  Find Home Element For Points (X,Y)                                          |
  !  Search Nearest Element to Progressively Further Elements. Updates Lagrangian| 
  !  component "host" and marks Lagrangian component "ifound" with 1 if          |
  !  found.  Returns logical variable "all_found" if all lagrangian variables    |
  !  have a known host element.  The host element may have been found prior to   |
  !  entry in this routine.                                                      |
  !==============================================================================|

  USE ALL_VARS
  IMPLICIT NONE
  !------------------------------------------------------------------------------!
  INTEGER,                   INTENT(IN)    :: NPTS
  REAL(SP), DIMENSION(NPTS), INTENT(IN)    :: XP,YP 
  INTEGER,  DIMENSION(NPTS), INTENT(INOUT) :: HOST,IFOUND,INDOMAIN,SBOUND,INWATER
  !------------------------------------------------------------------------------|
  INTEGER I,MIN_LOC
  REAL(SP), DIMENSION(1:N,1) :: RADLIST
  REAL(SP), DIMENSION(3) :: XTRI,YTRI
  REAL(SP) :: XLAG,YLAG,RADLAST
  INTEGER  :: LOCIJ(2)
  LOGICAL  :: ISINTRIANGLE 
  REAL(SP) :: X0C,Y0C 
  REAL(DP) XTMP,XTMP1	  
  INTEGER :: J  
!**********************pengfei******************************
  INTEGER:: PENG


!  DO I=1,NPTS
!     IF(IFOUND(I) == 1 .OR. INDOMAIN(I) == 0)CYCLE
!     XLAG  = XP(I) 
!     YLAG  = YP(I) 
!     XTRI    = VX(NV(MIN_LOC,1:3)) 
!     YTRI    = VY(NV(MIN_LOC,1:3)) 
!     IF(ISINTRIANGLE(XTRI,YTRI,XLAG,YLAG))THEN
!       IFOUND(I)  = 1 
!       HOST(I) = MIN_LOC
!       IF(ISONB(NV(HOST(I),1)) == 1.OR.ISONB(NV(HOST(I),2)) == 1.OR.ISONB(NV(HOST(I),3)) == 1)THEN
!          SBOUND(I)=1
!       ELSE
!          SBOUND(I)=0
!       ENDIF
!     END IF
!  END DO
!*********************************************************
  !==============================================================================|


  DO I=1,NPTS
    IF(IFOUND(I) == 1 .OR. INDOMAIN(I) == 0)CYCLE
    XLAG  = XP(I) 
    YLAG  = YP(I) 

#   if defined (SPHERICAL)
    DO J=1,N
     XTMP  = XC(J)*TPI-XLAG*TPI
     XTMP1 = XC(J)-XLAG
!     write(100,*) j,xtmp1,xc(j),xlag
     IF(XTMP1 > 180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP = 360.0_SP*TPI+XTMP
     END IF		
     X0C =XTMP*COS(DEG2RAD*(YC(J)+YLAG)*0.5)
     Y0C =(YC(J)-YLAG)*TPI

     RADLIST(J,1) = SQRT(X0C**2 + Y0C**2)
    END DO 
#   else
    RADLIST(1:N,1) = SQRT((XC(1:N)-XLAG)**2 + (YC(1:N)-YLAG)**2)
#   endif

    RADLAST = 0.0_SP
!!**********************pengfei********************
     peng=0
!!*********************************************
     in:  DO WHILE(.TRUE.)
!!*****************pengfei************changed by Pengfei 04/07 ************************************
!!              to only check eight closest elements before determining that particle is at boundary
        peng=peng+1
!***************************************
        LOCIJ   = MINLOC(RADLIST,RADLIST>RADLAST)
        MIN_LOC = LOCIJ(1)
        IF(MIN_LOC == 0 .OR. peng>=8) THEN
!           print *, min_loc ,peng
!           WRITE(*,*) 'particle',i,'hit solid or open boundary'
!           WRITE(*,*) 'location is: ' ,xlag+vxmin,ylag+vymin,xlag,ylag,vxmin,vymin

           EXIT in
!****************************end pengfei***************************************
        ELSEIF (MIN_LOC<0) THEN
           WRITE(20,*) 'MIN_LOC IS NEGATIVE'
        END IF
        XTRI    = VX(NV(MIN_LOC,1:3)) 
        YTRI    = VY(NV(MIN_LOC,1:3)) 
        RADLAST = RADLIST(MIN_LOC,1)
        IF(ISINTRIANGLE(XTRI,YTRI,XLAG,YLAG))THEN
           IFOUND(I)  = 1 
           HOST(I) = MIN_LOC
           IF(ISONB(NV(HOST(I),1)) == 1.OR.ISONB(NV(HOST(I),2)) == 1.OR.ISONB(NV(HOST(I),3)) == 1)THEN
              SBOUND(I)=1
           ELSE
              SBOUND(I)=0
           ENDIF
           EXIT in 
        END IF
        RADLAST = RADLIST(MIN_LOC,1)
     END DO in
  END DO

  WHERE(IFOUND == 0.AND.SBOUND == 0)
     INDOMAIN = 0
  END WHERE
  WHERE(IFOUND == 0.AND.SBOUND == 1)
     INWATER = 0
  END WHERE

  RETURN
END SUBROUTINE FHE_ROBUST

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!==============================================================================|

LOGICAL FUNCTION ISINTRIANGLE(XT,YT,X0,Y0) 
  !==============================================================================|
  !  Determine if Point (X0,Y0) Is In triangle defined by nodes (XT(3),YT(3))    |
  !  Using Algorithm Used for Scene Rendering in Computer Graphics               |
  !==============================================================================|
  !------------------------------------------------------------------------------|

  USE MOD_PREC
  IMPLICIT NONE
  !------------------------------------------------------------------------------!
  REAL(SP), INTENT(IN) :: X0,Y0
  REAL(SP), INTENT(IN) :: XT(3),YT(3)
  REAL(SP) :: F1,F2,F3
  !------------------------------------------------------------------------------|

  ISINTRIANGLE = .FALSE.  

  F1 = (Y0-YT(1))*(XT(2)-XT(1)) - (X0-XT(1))*(YT(2)-YT(1))
  F2 = (Y0-YT(3))*(XT(1)-XT(3)) - (X0-XT(3))*(YT(1)-YT(3))
  F3 = (Y0-YT(2))*(XT(3)-XT(2)) - (X0-XT(2))*(YT(3)-YT(2))
  IF(F1*F3 >= 0.0_SP .AND. F3*F2 >= 0.0_SP) ISINTRIANGLE = .TRUE.

  RETURN
END FUNCTION ISINTRIANGLE

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
SUBROUTINE WOMEGA(DELTAT,WWIN,WOUT,UIN,VIN,ELIN,ETIN)     
  !============================================================================|
  !  Compute Omega Velocity from cartesian                                     |
  !============================================================================|
  !============================================================================|
  USE MOD_PREC
  USE ALL_VARS
  IMPLICIT NONE
  !----------------------------------------------------------------------------|
  REAL(SP)                                   :: DELTAT
  REAL(SP), INTENT(IN),DIMENSION(0:N,KB)     :: WWIN,UIN,VIN
  REAL(SP), INTENT(INOUT),DIMENSION(0:N,KB)  :: WOUT
  REAL(SP), INTENT(IN), DIMENSION(0:M)       :: ELIN,ETIN
  !----------------------------------------------------------------------------|
  REAL(SP) :: DDDX,DDDY,DEDX,DEDY,ETF1AA,WW1,WW2,ET1
  INTEGER  :: I,K,J1,J2,J3
  !============================================================================|

  D(:) = H(:) + ELIN(:)

  DO I=1,N

        J1=NV(I,1)
        J2=NV(I,2)
        J3=NV(I,3)

        DDDX=AWX(I,1)*D(J1) + AWX(I,2)*D(J2) + AWX(I,3)*D(J3)
        DDDY=AWY(I,1)*D(J1) + AWY(I,2)*D(J2) + AWY(I,3)*D(J3)
        DEDX=AWX(I,1)*ELIN(J1) + AWX(I,2)*ELIN(J2) + AWX(I,3)*ELIN(J3)
        DEDY=AWY(I,1)*ELIN(J1) + AWY(I,2)*ELIN(J2) + AWY(I,3)*ELIN(J3)
        ETF1AA=ONE_THIRD*(ELIN(NV(I,1))+ELIN(NV(I,2))+ELIN(NV(I,3)))
        ET1=ONE_THIRD*(ETIN(NV(I,1))+ETIN(NV(I,2))+ETIN(NV(I,3)))

        WOUT(I,1) = ZERO
        DO K=1,KBM1
           WW2=(ZZ1(I,K)+1.)*(ETF1AA-ET1)/DELTAT
           WW1=WWIN(I,K)-WW2
           WOUT(I,K+1)=2.0_SP*(WW1-U(I,K)*(Z1(I,K)*DDDX+DEDX)-V(I,K)*(Z1(I,K)*DDDY+DEDY))-WOUT(I,K) 
        END DO
        
        !--for conservation issue
        DO K=1,KBM1
           WOUT(I,K) = WOUT(I,K) + WOUT(I,KB)/KBM1
        END DO
        WOUT(I,KB) = ZERO
  
  END DO

  RETURN
END SUBROUTINE WOMEGA
!==============================================================================|

!********************ADDED BY JHC 3/07 *********************

      INTEGER FUNCTION JD(YEAR,MONTH,DAY)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
!
      INTEGER YEAR,MONTH,DAY,I,J,K
!
      I= YEAR
      J= MONTH
      K= DAY
!
      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12) &
         /12-3*((I+4900+(J-14)/12)/100)/4
!
      RETURN
      END FUNCTION JD

!

!---------------------------check day of Feb by pengfei-------------------------

      INTEGER FUNCTION hourfeb(YEAR)
!
!
      INTEGER YEAR,I,J,K
      hourfeb=28*24
      if (mod(year,100)==0 ) then
          if ( mod(year,400)==0) then
          hourfeb=29*24
          print *, 'FYI, leap year!'
          end if
      else
         if (mod(year,4)==0) then
            hourfeb=29*24
            print *, "FYI, leap year!"
         end if
      end if 
      RETURN
      END FUNCTION hourfeb



