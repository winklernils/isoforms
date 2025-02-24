!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!        continuation 2 species BVP
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!--------- ---- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
  
  DOUBLE PRECISION s,x,cA,cB,Z,V,lam,cinf,cr,Kr,chir,s0

  s = U(1)
  x = U(2)
  cA = U(3)
  cB = U(4)

  Z = PAR(1)
  V = PAR(2)
  lam = PAR(3) * 100.0
  cinf = PAR(4)
  cr = PAR(5)
  Kr = PAR(6)
  chir = PAR(7)
  s0 = PAR(8)

  
  F(1) = x
  F(2) = s/Z + s0/Z - lam/Z * (cA + cB / chir)
  F(3) = cA * (cinf-(cA+cB))* (cB+(cA-cinf)*Kr) *(V-x) / (cinf*cinf*Kr)
  F(4) = cB * (cinf-(cA+cB))* (cB-cinf+cA*Kr)   *(V-x) / (cinf*cinf*Kr)
 
  
  
  !-- Provide user-supplied derivatives of functions with respect to parameters  
  IF(IJAC.EQ.0)RETURN 
  
  
  IF(IJAC.EQ.1)RETURN 
  ! F(1) is independent of parameters, DFDP(1,*) = 0 which is already default value
  

END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T) 
!--------- ----- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  PAR(1) = 1.0d0
  PAR(2) = 0.0d0
  PAR(3) = 0.d0/100
  PAR(4) = 1.5d0
  PAR(5) = 4.5d0
  PAR(6) = 0.2d0
  PAR(7) = 0.05d0/0.18d0
  PAR(8) = PAR(3)*100 * (1 + 1/ (PAR(5)*PAR(7)) )

  U(1) = 0.0d0
  U(2) = 0.0d0
  U(3) = 1.0d0
  U(4) = 1.0d0 / PAR(5)

END SUBROUTINE STPNT

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!--------- ---- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)


  ! REDUCED SINCE ONE B.C. SATISFIED AUTOMATICALLY
  FB(1)=U0(1) 
  FB(2)=U1(1) 
  FB(3)=U0(2)-PAR(2)
  FB(4)=U1(2)-PAR(2)
  
  !-- Provide user-supplied derivatives of boundary conditions with respect to B.C.s  
  IF(IJAC.EQ.0)RETURN 
  
  DBC(1,1) = 1.0 ! DFB(1)/DU0(1)
  DBC(2,5) = 1.0 ! DFB(2)/DU1(1)
  DBC(3,2) = 1.0 ! DFB(3)/DU0(2)
  DBC(4,6) = 1.0 ! DFB(4)/DU1(2)

  IF(IJAC.EQ.1)RETURN 
  ! parameter derivatives
  DBC(3,9) = -PAR(3) ! DFB(3)/DPAR(1)
  DBC(3,11) = -PAR(1)
  
  DBC(4,9) = -PAR(3)
  DBC(4,11) = -PAR(1)
  
END SUBROUTINE BCND

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!--------- ----

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

  FI(1) = U(3) - 1.0 
  FI(2) = U(4) - 1.0/PAR(5)

  IF(IJAC.EQ.0)RETURN 
  
  DINT(1,3) = 0.d0 ! DFI(1)/DU(3)
  DINT(2,4) = 0.d0
  
  IF(IJAC.EQ.1)RETURN
  ! parameter derivatives
  DINT(1,5) = 0.d0 ! DFI(1)/DPAR(5)
  DINT(2,5) = -1.d0/(PAR(5)*PAR(5))
  DINT(2,10) = 0.d0 

END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

SUBROUTINE PVLS
END SUBROUTINE PVLS
