SUBROUTINE fndset (gpid,x,ibuf,n)
     
!     GPID = GRID POINT TABLE FOR THIS SET
 
!     N  = 0 INPUT
!     FNDSET READS THE COORDINATES OF THE GRID POINTS IN THIS SET.
!     IF THE GRID POINT TABLE VALUE IS ZERO THE CORRESPONDING GRID
!     POINT IS NOT USED IN THIS SET AND ITS VALUES SKIPPED, OTHERWISE
!     THE XYZ COORDINDATE VALUES ARE READ FROM BGPDT AND PACKED INTO
!     X SPACE. TOTALLY THERE ARE NGPSET GRID DATA SAVED IN X.
!     CORE NEEDED FOR X = 3*NGPSET (PROVIDED BY CALLING ROUTINE)
 
!     N  = 1 INPUT/OUTPUT
!     FNDSET POSITIONS THE STRESS FILE TO THE SUBCASE/VALUE LAST
!     PROCESSED
 
 
 INTEGER, INTENT(IN OUT)                  :: gpid(1)
 REAL, INTENT(IN OUT)                     :: x(3,1)
 INTEGER, INTENT(IN OUT)                  :: ibuf
 INTEGER, INTENT(OUT)                     :: n
 INTEGER :: bgpdt,oes1,rew,subc
 REAL :: u(3)
 COMMON /BLANK / ngp,skp11(4),ngpset,skp12(4),skp21(4),bgpdt, skp22(8),oes1
 COMMON /names / nirew,inprew,skpn1(2),rew,norew
 COMMON /xxparm/ skpp(211),subc,flag,DATA
 EQUIVALENCE     (u(1),insub)
 DATA    twopi / 0.0 /
 
 IF (n /= 0) GO TO 30
 CALL gopen (bgpdt,gpid(ibuf),inprew)
 j = 1
 DO  i = 1,ngp
   IF (gpid(i) /= 0) GO TO 10
   CALL fread (bgpdt,0,-4,0)
   CYCLE
   10 CALL fread (bgpdt,0,-1,0)
   CALL fread (bgpdt,x(1,j),3,0)
   j = j + 1
 END DO
 CALL CLOSE (bgpdt,rew)
 GO TO 110
 
!     POSITION OES1
 
 30 IF (twopi < 6.2) twopi = 8.0*ATAN(1.0)
 CALL gopen (oes1,gpid(ibuf),inprew)
 40 CALL READ  (*90,*90,oes1,j,1,0,i)
 CALL fread (oes1,0,-2,0)
 CALL fread (oes1,u,3,0)
 IF (subc /= insub) GO TO 70
 IF (flag-1.0 < 0.0) THEN
   GO TO   100
 ELSE IF (flag-1.0 == 0.0) THEN
   GO TO    60
 END IF
 50 j = j/10
 
!     REAL EIGENVALUE ANALYSIS - CONVERT TO FREQUENCY
 
 IF (j == 2) u(3) = SQRT(ABS(u(3)))/twopi
 IF (DATA-u(3) > 1.0E-6) GO TO 70
 GO TO 100
 60 IF (DATA-u(2) == 0.0) THEN
   GO TO   100
 END IF
 
!     WRONG CASE
 
 70 CALL fwdrec (*90,oes1)
 CALL fwdrec (*90,oes1)
 GO TO 40
 90 n = n + 1
 100 CALL bckrec (oes1)
 CALL CLOSE  (oes1,norew)
 
 110 RETURN
END SUBROUTINE fndset
