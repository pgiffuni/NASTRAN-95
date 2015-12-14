SUBROUTINE sd2rhd (istyp,isetup)
     
!     THIS ROUTINE WRITES HEADING FOR PRECISION CHECK IN SDR2E.
!     WORDS 1,2,6 AND 7 PRESET BY CALLING ROUTINE.
!     ISETUP.NE.0 FIRST CALL.
 
 
 INTEGER, INTENT(OUT)                     :: istyp(7)
 INTEGER, INTENT(IN OUT)                  :: isetup
 INTEGER :: branch, ldmd(8)
 COMMON /sdr2x4/ dummy(50),branch
 COMMON /system/ isysb    ,nout
 EQUIVALENCE     (istyp6,rstyp6), (istyp7,rstyp7)
 
 DATA    ldmd  / 4HLOAD, 4HMODE, 4H, fr, 4HEQ.=, 4H, ei, 4HGEN=,  &
     4H, ti, 4HME =  /
 
 IF (isetup == 0) GO TO 1510
 SELECT CASE ( branch )
   CASE (    1)
     GO TO 1501
   CASE (    2)
     GO TO 1503
   CASE (    3)
     GO TO 1501
   CASE (    4)
     GO TO 1501
   CASE (    5)
     GO TO 1503
   CASE (    6)
     GO TO 1505
   CASE (    7)
     GO TO 1501
   CASE (    8)
     GO TO 1507
   CASE (    9)
     GO TO 1507
   CASE (   10)
     GO TO 1501
 END SELECT
 
!     STATICS
 
 1501 n1 = 3
 istyp(3) = ldmd(1)
 GO TO 1510
 
!     EIGR,FREQ
 
 1503 n1 = 6
 istyp(3) = ldmd(2)
 istyp(4) = ldmd(3)
 istyp(5) = ldmd(4)
 GO TO 1510
 
!     TRANSIENT
 
 1505 n1 = 6
 istyp(3) = ldmd(1)
 istyp(4) = ldmd(7)
 istyp(5) = ldmd(8)
 GO TO 1510
 
!     BUCKLING, COMPLEX EIGENVALUE
 
 1507 n1 = 6
 istyp(3) = ldmd(2)
 istyp(4) = ldmd(5)
 istyp(5) = ldmd(6)
 IF (branch == 9) n1 = 7
 
 1510 CALL page2 (3)
 istyp6 = istyp(6)
 istyp7 = istyp(7)
 IF (n1 == 3) WRITE(nout,1512) (istyp(i),i=1,n1)
 IF (n1 == 6) WRITE(nout,1512) (istyp(i),i=1,5),rstyp6
 IF (n1 == 7) WRITE(nout,1512) (istyp(i),i=1,5),rstyp6,rstyp7
 1512 FORMAT (1H0,5X,45HE l e m e n t   p r e c i s i o n   c h e c k,  &
     /4X,32HSIGNIFICANT digits for subcase =,i7,1H,,i7,3H = ,3A4, 1P,2E15.6)
 RETURN
END SUBROUTINE sd2rhd
