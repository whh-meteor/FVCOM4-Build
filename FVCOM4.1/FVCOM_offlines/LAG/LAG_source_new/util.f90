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
!==============================================================================!
!      From Numerical Recipes
!==============================================================================!

SUBROUTINE hunt(xx,n,x,jlo)
  ! from numerical recipies vol 2 
  USE MOD_PREC

  INTEGER, INTENT(INOUT) :: jlo
  INTEGER , INTENT(IN)   :: n
  REAL(SP) , INTENT(IN)  :: xx(n)
  REAL(SP), INTENT(IN)   :: x

  INTEGER   :: inc,jhi,jm
  LOGICAL   :: ascnd

  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1) then
     if (x.eq.xx(n)) jlo=n-1
     if(x.eq.xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3

END SUBROUTINE hunt

!==============================================================================|

FUNCTION gasdev(idum)

  USE MOD_PREC

  INTEGER , INTENT(INOUT) :: idum

  REAL(SP) :: gasdev
  INTEGER  :: iset
  REAL(SP) :: fac,gset,rsq,v1,v2,ran1
  SAVE iset,gset
  DATA iset/0/

  if (iset.eq.0) then
1    v1=2.*ran1(idum)-1.
     v2=2.*ran1(idum)-1.
     rsq=v1**2+v2**2
     if(rsq.ge.1..or.rsq.eq.0.)goto 1
     fac=sqrt(-2.*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
  else
     gasdev=gset
     iset=0
  endif
  return
END FUNCTION gasdev

!==============================================================================|

FUNCTION ran1(idum)

  USE MOD_PREC
  ! from numerical recipies vol 2 
  INTEGER , INTENT(INOUT) :: idum

  INTEGER   ::IA,IM,IQ,IR,NTAB,NDIV
  REAL(SP)  :: ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER   :: j,k,iv(NTAB),iy
  SAVE iv,iy
  DATA iv /NTAB*0/, iy /0/

  if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
     end do
     iy=iv(1)
  endif
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1=min(AM*iy,RNMX)
  return
END FUNCTION ran1

!==============================================================================|

SUBROUTINE spline(x,y,n2,yp1,ypn,y2)

  USE MOD_PREC
  !     From numerical recipies vol 2, but modfied so that nmax=50
  INTEGER  :: n2
  REAL(SP), INTENT(OUT), DIMENSION(n2) :: y2
  REAL(SP), INTENT(IN)                 :: yp1,ypn
  REAL(SP), INTENT(IN), DIMENSION(n2)  :: x,y

  INTEGER  :: NMAX
  PARAMETER (NMAX=50)
  INTEGER  :: i,k
  REAL(SP) :: p,qn,sig,un,u(NMAX)
 
 if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n2-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))--sig*u(i-1))/p
  END DO
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n2)-x(n2-1)))*(ypn-(y(n2)-y(n2-1))/(x(n2)-x(n2-1)))
  endif
  y2(n2)=(un-qn*u(n2-1))/(qn*y2(n2-1)+1.)
  do k=n2-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  END DO
  return
END SUBROUTINE spline

!==============================================================================|
