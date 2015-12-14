SUBROUTINE tranem (mcsid, ngrid, r, icomp, u, rc)
!*****
!     COMPUTES A STRESS TRANSFORMATION MATRIX U FOR TRIANGLES AND QUADS.
!     INPUTS
!        MCSID  ID OF COORDINATE SYSTEM REFERENCED ON MAT1,2 DATA CARD.
!        NGRID  3 FOR TRIANGLES, 4 FOR QUADS.
!        R      ARRAY OF BASIC LOCATIONS OF ELEMENT GRID PTS (3,NGRID).
!     OUTPUTS
!        ICOMP  1 (IF MAT X-AXIS IS USED) OR 2 (IF Y-AXIS IS USED).
!        U      ARRAY (3X3) FOR TRANSFORMATION, STORED BY ROW.
!        RC     BASIC LOCATION COORDINATES OF ELEMENT CENTER.
!     REQUIREMENTS
!        SUBROUTINE PRETRS MUST SET UP FOR TRANSS. SEE P.M. PAGE 3.4-66
!*****
 
 INTEGER, INTENT(IN)                      :: mcsid
 INTEGER, INTENT(IN)                      :: ngrid
 REAL, INTENT(IN)                         :: r(9)
 INTEGER, INTENT(OUT)                     :: icomp
 REAL, INTENT(IN OUT)                     :: u(9)
 REAL, INTENT(OUT)                        :: rc(3)
 INTEGER :: ecpt(4),subnam(2)
 
 
 
 
 REAL :: rcent(4)
 
 EQUIVALENCE  (rcent(1), ecpt(1))
 
 DATA subnam /4HTRAN,2HEM/
 
!-----------------------------------------------------------------------
 
 IF(ngrid /= 3 .AND. ngrid /= 4 )  CALL mesage(-61,0,subnam)
!*****
!     FIND THE UNIT NORMAL OF THE ELEMENT
!*****
 i = 3*(ngrid-3)
 vn1 = (r(8)-r(2))*(r(i+9)-r(6))-(r(9)-r(3))*(r(i+8)-r(5))
 vn2 = (r(9)-r(3))*(r(i+7)-r(4))-(r(7)-r(1))*(r(i+9)-r(6))
 vn3 = (r(7)-r(1))*(r(i+8)-r(5))-(r(8)-r(2))*(r(i+7)-r(4))
 temp = SQRT(vn1**2+vn2**2+vn3**2)
 IF(temp <= 0.0) CALL mesage(-61,0,subnam)
 vn1 = vn1 / temp
 vn2 = vn2 / temp
 vn3 = vn3 / temp
!*****
!     GET THE UNIT VECTORS OF MCSID AT ELEM CENTER. PUT IN U TEMPORARILY
!*****
 grds = ngrid
 DO  ic=1,3
   sum = 0.0
   DO  ig=1,ngrid
     k = 3*ig + ic-3
     sum = sum +r(k)
   END DO
   rcent(ic+1) = sum / grds
   rc(ic) = rcent(ic+1)
 END DO
 ecpt(1) = mcsid
 CALL transs(ecpt,u)
!*****
!     SELECT FIRST OR SECOND VECTOR TO PROJECT FOR ELEM-MAT X-AXIS
!*****
 vndotm=vn1*u(1)+vn2*u(4)+vn3*u(7)
 IF( vndotm**2 > 0.4) GO TO 30
 icomp = 1
 vm1 = u(1)
 vm2 = u(4)
 vm3 = u(7)
 GO TO 40
 30 CONTINUE
 icomp = 2
 vm1 = u(2)
 vm2 = u(5)
 vm3 = u(8)
 vndotm = vn1*vm1+vn2*vm2+vn3*vm3
 40 CONTINUE
!*****
!     FIND COSINE AND SINE OF ANGLE
!*****
 ve1 = r(4) - r(1)
 ve2 = r(5) - r(2)
 ve3 = r(6) - r(3)
 c = ve1*(vm1-vndotm*vn1) + ve2*(vm2-vndotm*vn2)  &
     + ve3*(vm3-vndotm*vn3)
 s = ve1*(vm2*vn3-vm3*vn2) + ve2*(vm3*vn1-vm1*vn3)  &
     + ve3*(vm1*vn2-vm2*vn1)
 temp = SQRT(c*c+s*s)
 IF(temp <= 0.0) CALL mesage(-61,0,subnam)
 c = c/temp
 s = s/temp
!*****
!     FILL IN THE U MATRIX, ROW STORED.
!*****
 u(1) = c*c
 u(4) = s*s
 u(7) = -c*s
 u(2) = u(4)
 u(5) = u(1)
 u(8) = -u(7)
 u(3) = 2.0*u(8)
 u(6) = -u(3)
 u(9) = u(1)-u(4)
 
 RETURN
 
END SUBROUTINE tranem
