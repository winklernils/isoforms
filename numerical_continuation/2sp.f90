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
  
  DOUBLE PRECISION s,x,cA,cB,L,Z,V,kap,cinf,cr,Kr,PA,PB

  s = U(1)
  x = U(2)
  cA = U(3)
  cB = U(4)

  L = PAR(1)
  Z = PAR(2)
  V = PAR(3)
  kap = PAR(4) * 1000.0
  cinf = PAR(5)
  cr = PAR(6)
  Kr = PAR(7)
  PA = PAR(8)
  PB = PAR(9)

  
  F(1) = x
  F(2) = L*L/Z* s - L*L/Z * kap * (L - 1.0) - L*L/Z * kap  * (PA*cA + PB*cB)
  F(3) = cA * (cinf-(cA+cB))* (cB+(cA-cinf)*Kr) *(L*V-x) / (cinf*cinf*Kr)
  F(4) = cB * (cinf-(cA+cB))* (cB-cinf+cA*Kr)   *(L*V-x) / (cinf*cinf*Kr)
 
  
  
  !-- Provide user-supplied derivatives of functions with respect to parameters  
  IF(IJAC.EQ.0)RETURN 
  
  DFDU(1,1) = 0.0
  DFDU(1,2) = 1.0
  DFDU(1,3) = 0.0
  DFDU(1,4) = 0.0
  
  DFDU(2,1) = L*L/Z
  DFDU(2,2) = 0.0
  DFDU(2,3) = - L*L/Z * kap * PA
  DFDU(2,4) = - L*L/Z * kap * PB
  
  DFDU(3,1) = 0.0
  DFDU(3,2) = - cA * (cinf-(cA+cB))* (cB+(cA-cinf)*Kr) / (cinf*cinf*Kr)
  DFDU(3,3) = -(3.d0*cA*cA *Kr + 2.d0*cA*(cB + cB*Kr - 2.d0*cinf*Kr)+(cB-cinf) * (cB-cinf*Kr))*(L*V-x) / (cinf*cinf*Kr)
  DFDU(3,4) = -cA*(2.d0*cB+cA*(1.d0+Kr)-cinf*(1.d0 + Kr))*(L*V - x) / (cinf*cinf*Kr)
  
  DFDU(4,1) = 0.0
  DFDU(4,2) = -cB * (cinf-(cA+cB))* (cB-cinf+cA*Kr)/ (cinf*cinf*Kr)
  DFDU(4,3) = -cB*(2.0*cA*Kr + cB*(1.0+Kr)-cinf * (1.0+Kr)) * (L*V - x) / (cinf*cinf*Kr)
  DFDU(4,4) = -(3.0*cB*cB - 4.0*cB*cinf + cinf*cinf + cA*cA * Kr + cA*(2.0*cB-cinf) * (1.0+Kr))*(L*V - x) / (cinf*cinf*Kr)
  
  
  IF(IJAC.EQ.1)RETURN 
  ! F(1) is independent of parameters, DFDP(1,*) = 0 which is already default value
  DFDP(2,1) = 2.0*L/Z*S - kap/Z*(3.0*L*L - 2.0*L) - 2.0*L/Z*kap * (PA*cA + PB*cB)
  DFDP(2,2) = -1.0/(Z*Z) * (L*L*s - L*L* kap * (L - 1.0) - L*L* kap * (PA*cA + PB*cB))
  DFDP(2,3) = 0.0
  DFDP(2,4) = - L*L/Z * 1000.0 * (L - 1.0) - L*L/Z * 1000.0 * (PA*cA + PB*cB)
  DFDP(2,5) = 0.0
  DFDP(2,6) = 0.0
  DFDP(2,7) = 0.0
  DFDP(2,8) = - L*L/Z * kap * cA
  DFDP(2,9) = - L*L/Z * kap * cB
  
  DFDP(3,1) = cA * (cinf-(cA+cB))* (cB+(cA-cinf)*Kr)*V / (cinf*cinf*Kr)
  DFDP(3,2) = 0.0
  DFDP(3,3) = cA * (cinf-(cA+cB))* (cB+(cA-cinf)*Kr) * L / (cinf*cinf*Kr)
  DFDP(3,4) = 0.0
  DFDP(3,5) = (cA*(2.0*cB*(cA+cB)-cB*cinf+ 2.0*cA*(cA+cB)*Kr-(2.0*cA+cB)*cinf*Kr)*(L*V-x))/(cinf*cinf*cinf * Kr)
  DFDP(3,6) = 0.0
  DFDP(3,7) = (cA*cB*(cA+cB-cinf)*(L*V-x))/(cinf*cinf * Kr*Kr)
  DFDP(3,8) = 0.0
  DFDP(3,9) = 0.0
  
  DFDP(4,1) = cB * (cinf-(cA+cB))* (cB-cinf+cA*Kr)*V / (cinf*cinf*Kr)
  DFDP(4,2) = 0.0
  DFDP(4,3) = cB * (cinf-(cA+cB))* (cB-cinf+cA*Kr)*L / (cinf*cinf*Kr)
  DFDP(4,4) = 0.0
  DFDP(4,5) = (cB * (2*cB*(cB-cinf)+2.0*cA*cA *Kr+cA*(2.0*cB-cinf)*(1.0+Kr))*(L*V-x))/(cinf*cinf*cinf *Kr)
  DFDP(4,6) = 0.0
  DFDP(4,7) = (cB*(cB-cinf)*(cA+cB-cinf)*(L*V-x))/(cinf*cinf * Kr*Kr)
  DFDP(4,8) = 0.0
  DFDP(4,9) = 0.0
  

END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T) 
!--------- ----- 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  
  PAR(2) = 1.25d0
  PAR(3) = 0.d0
  PAR(4) = 0.1/1000.0
  PAR(5) = 4.0d0
  PAR(6) = 1.5d0
  PAR(7) = 0.2d0
  PAR(8) = 0.05d0
  PAR(9) = 0.18d0
  PAR(1) = ( 1.d0 + SQRT(1.d0 - 4.d0*PAR(8) - 4.d0*PAR(9)/PAR(6)) )/2.d0

  U(1) = 0.0d0
  U(2) = 0.0d0
  U(3) = 1.0d0 / PAR(1)
  U(4) = 1.0d0 / (PAR(1) * PAR(6))

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
  FB(3)=U0(2)-PAR(1)*PAR(3)
  FB(4)=U1(2)-PAR(1)*PAR(3)
  
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

  FI(1) = PAR(1)*U(3) - 1.0 
  FI(2) = PAR(1)*U(4) - 1.0/PAR(6)

  IF(IJAC.EQ.0)RETURN 
  
  DINT(1,3) = PAR(1) ! DFI(1)/DU(3)
  DINT(2,4) = PAR(1)
  
  IF(IJAC.EQ.1)RETURN
  ! parameter derivatives
  DINT(1,5) = U(3) ! DFI(1)/DPAR(1)
  DINT(2,5) = U(4)
  DINT(2,10) = 1.0/(PAR(6)*PAR(6)) 

END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

SUBROUTINE PVLS
END SUBROUTINE PVLS
