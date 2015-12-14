SUBROUTINE selas1(iarg)
!*****
! THIS ROUTINE IS PHASE I OF STRESS DATA RECOVERY FOR THE ELAS ELEMENTS.
 
 
 
!*****
 
 
 
!              E C P T - S  F O R  E L A S  E L E M E N T S
 
 
 
!                  TYPE             TYPE           TYPE           TYPE
!         CELAS1           CELAS2         CELAS3         CELAS4
! ECPT(1) IELID     I      IELID     I    IELID      I   IELID      I
! ECPT(2) IGP1      I      K         R    IS1        I   K          R
! ECPT(3) IGP2      I      IGP1      I    IS2        I   IS1        I
! ECPT(4) IC1       I      IGP2      I    K          R   IS2        I
! ECPT(5) IC2       I      IC1       I    GSUBE      R
! ECPT(6) K         R      IC2       I    S          R
! ECPT(7) GSUBE     R      GSUBE     R
! ECPT(8) S         R      S         R
 
 
 
 
 INTEGER, INTENT(IN OUT)                  :: iarg
 DIMENSION iecpt(6)
 
! SDR2 PHASE I INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x5/ ecpt(100),  &
     jelid              ,isilno(2)  &
     ,                  stiff              ,scoeff ,                  dummy(120)
 
 
 
 EQUIVALENCE (iecpt(1),ecpt(1)) ,(scoeff,icoeff)
 
! BUILD UP OUTPUT BLOCK DEPENDING UPON WHICH ELEMENT TYPE, ELAS1, ELAS2,
! ELAS3 OR ELAS4, IS BEING WORKED ON.
 
 SELECT CASE ( iarg )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
 END SELECT
 
! ELAS1
 
 10 isilno(1) = iecpt(2) + iecpt(4)
 isilno(2) = iecpt(3) + iecpt(5)
 IF (iecpt(4) > 0) isilno(1) = isilno(1) - 1
 IF (iecpt(5) > 0) isilno(2) = isilno(2) - 1
 stiff  = ecpt(6)
 scoeff = ecpt(8)
 GO TO 50
 
! ELAS2
 
 20 isilno(1) = iecpt(3) + iecpt(5)
 isilno(2) = iecpt(4) + iecpt(6)
 IF (iecpt(5) > 0) isilno(1) = isilno(1) - 1
 IF (iecpt(6) > 0) isilno(2) = isilno(2) - 1
 stiff  = ecpt(2)
 scoeff = ecpt(8)
 GO TO 50
 
! ELAS3
 
 30 isilno(1) = iecpt(2)
 isilno(2) = iecpt(3)
 stiff  = ecpt(4)
 scoeff = ecpt(6)
 GO TO 50
 
! ELAS4
 
 40 isilno(1) = iecpt(3)
 isilno(2) = iecpt(4)
 stiff  = ecpt(2)
 icoeff = -1
 
! STORE ELEMENT ID.
 
 50 jelid = iecpt(1)
 RETURN
END SUBROUTINE selas1
