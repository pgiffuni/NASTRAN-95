SUBROUTINE kelas (ijklmn)
!*****
! THIS ROUTINE COMPUTES THE ELEMENT STIFFNESS AND STIFFNESS DAMPING
! 1 X 1 MATRICES FOR ELEMENTS ELAS1, ELAS2, ELAS3, ELAS4.
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
 
 
 
 
 INTEGER, INTENT(IN)                      :: ijklmn
 DOUBLE PRECISION :: ke
 
 
 
 DIMENSION iecpt(5)
 
 
 
 COMMON   /system/ isys
 
! SMA1 I/O PARAMETERS
 
 COMMON   /sma1io/ ifcstm             ,ifmpt  &
     ,                  ifdit              ,idum1  &
     ,                  ifecpt             ,igecpt  &
     ,                  ifgpct             ,iggpct  &
     ,                  ifgei              ,iggei  &
     ,                  ifkgg              ,igkgg  &
     ,                  if4gg              ,ig4gg  &
     ,                  ifgpst             ,iggpst  &
     ,                  inrw               ,outrw  &
     ,                  clsnrw             ,clsrw  &
     ,                  neor               ,eor  &
     ,                  mcbkgg(7)          ,mcb4gg(7)
 
! SMA1 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON   /sma1bk/ icstm              ,ncstm  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6k              ,n6x6k  &
     ,                  i6x64              ,n6x64
 
! SMA1 PROGRAM CONTROL PARAMETERS
 
 COMMON   /sma1cl/ iopt4              ,k4ggsw  &
     ,                  npvt               ,left  &
     ,                  frowic             ,lrowic  &
     ,                  nrowsc             ,tnrows  &
     ,                  jmax               ,nlinks  &
     ,                  link(10)           ,idetck  &
     ,                  dodet              ,nogo
 
! ECPT COMMON BLOCK
 
 COMMON   /sma1et/ ecpt(100)
 
 
 
 EQUIVALENCE (iecpt(1),ecpt(1))
 
 
 
 DATA iscalr /0/
 
 
 
 iarg = ijklmn
 
! MAKE THE ECPT-S FOR ALL ELAS ELEMENTS LOOK EXACTLY LIKE THE ECPT FOR
! ELAS1
 
 SELECT CASE ( iarg )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 10
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
 END SELECT
 
! ELAS2
 
 10 SAVE = ecpt(2)
 DO  i = 3,6
   iecpt(i-1) = iecpt(i)
 END DO
 ecpt(6) = SAVE
 GO TO 50
 
! ELAS3
 
 30 ecpt(7) = ecpt(5)
 ecpt(6)  = ecpt(4)
 iecpt(4) = 1
 iecpt(5) = 1
 GO TO 50
 
! ELAS4
 
 40 ecpt(6)  = ecpt(2)
 iecpt(2) = iecpt(3)
 iecpt(3) = iecpt(4)
 iecpt(4) = 1
 iecpt(5) = 1
 
! DETERMINE WHICH POINT IS THE PIVOT POINT AND SET APPROPRIATE POINTERS
 
 50 ind = 2
 IF (iecpt(2) == npvt) GO TO 60
 IF (iecpt(3) /= npvt) RETURN
 ipvt  = 3
 ipdof = 5
 inpvt = 2
 inpdof = 4
 IF (iecpt(2) == 0) ind = 1
 GO TO 80
 
! CHECK TO SEE IF BOTH POINTS MATCH THE PIVOT POINT.
 
 60 IF (iecpt(3) /= npvt) GO TO 70
 IF (iscalr == 0) GO TO 65
 iscalr = 0
 RETURN
 65 iscalr = 1
 ind = 4
 70 ipvt   = 2
 ipdof  = 4
 inpvt  = 3
 inpdof = 5
 IF (iecpt(3) == 0) ind = 1
 80 IF (iecpt(ipdof)  <= 0) iecpt(ipdof)  = 1
 IF (iecpt(inpdof) <= 0) iecpt(inpdof) = 1
 
! II AND JJ ARE THE ROW AND COLUMN INDICES OF THE MATRIX INTO WHICH THE
! SPRING AND SPRING DAMPING CONSTANTS WILL BE ADDED.
 
 ii = iecpt(ipvt)  + iecpt(ipdof)  - 1
 jj = iecpt(inpvt) + iecpt(inpdof) - 1
 ke = ecpt(6)
 INDEX = 6
 ifile = ifkgg
 85 ASSIGN 100 TO iretrn
 i = ii
 j = ii
 90 CALL sma1b (ke,j,i,ifile,0.0D0)
 IF (ind == 1) GO TO 130
 GO TO iretrn, (100,110,120,130)
 100 ASSIGN 110 TO iretrn
 ke = - ke
 j  =   jj
 GO TO 90
 110 IF (ind /= 4) GO TO 130
 ASSIGN 120 TO iretrn
 ke = ecpt(6)
 i = jj
 GO TO 90
 120 ASSIGN 130 TO iretrn
 ke = -ke
 j  = ii
 GO TO 90
 130 IF (INDEX == 7)  RETURN
 IF (iopt4 == 0  .OR.  iarg == 4) RETURN
 
! IF G SUB E IS NON-ZERO, SET PARAMETERS FOR K4GG INSERTION.
 
 IF (ecpt(7) == 0.0) RETURN
 k4ggsw = 1
 ifile = if4gg
 ke = ecpt(7) * ecpt(6)
 INDEX = 7
 GO TO 85
END SUBROUTINE kelas
