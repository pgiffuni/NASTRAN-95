SUBROUTINE stpax3 ( again)
     
 LOGICAL, INTENT(IN OUT)                  :: again
 INTEGER :: iforce(25), istres(100), elemid, iblock(62,14)  &
     ,       iclock(62,14),cangle
 
 REAL :: savef(75), saves(75)
 
 
 
 COMMON /sdr2x7/ dum(100),stress(100),force(25)  &
     ,               skip(729),BLOCK(62,14),clock(62,14)
 
! SCRATCH BLOCK
 COMMON /sdr2x8/ disp(59),nangle,elemid,unu(94),kangle,klemid
 COMMON /isave / isavef(75),isaves(75)
 
 COMMON /sdr2de/ dum5(33), ipart
 
 COMMON /sdr2x4/ dum4(51),ktype
 
 
 EQUIVALENCE ( istres(1), stress(1)), ( iforce(1), force(1))  &
     ,            (iblock(1,1), BLOCK(1,1))  &
     ,            (iclock(1,1),clock(1,1)),(isavef(1),savef(1))  &
     ,            (isaves(1),saves(1)),(nangle,cangle)
 
 IF ( again ) GO TO 10
 again = .true.
 kangle = 0
 10 nangle = kangle
 elemid = klemid
 nangle = nangle + 1
 kangle = nangle
 
 
!  BRANCH TO INSERT STRESSES AND FORCES INTO FORCE AND STRESS OR
!                                            SAVEF AND SAVES
 
!    KTYPE=1 - REAL OUTPUT FROM BLOCK IS TRANSFERED TO CLOCK, THEN
!              STORED IN FORCE AND STRESS, NOTHING IN SAVEF AND SAVES
!    KTYPE=2 - COMPLEX OUTPUT
!    IPART=1 - IMAGINARY PART OF COMPLEX OUTPUT FROM BLOCK, STORED
!              IN SAVEF AND SAVES
!    IPART=2 - REAL PART OF COMPLEX OUTPUT FROM CLOCK STORED IN
!              FORCE AND STRESS
 
 IF(ktype == 2) GO TO 19
 DO  i=1,62
   DO  j=1,14
     clock(i,j) = BLOCK(i,j)
   END DO
 END DO
 19 CONTINUE
 
!  OUTPUT FORCES FOR THIS ANGLE
 iforce(1)=elemid
 force(2) = clock(1,cangle)
 DO  i=1,16
   force(2+i) = clock(46+i,cangle)
 END DO
 
! OUTPUT STRESSES
 istres  (1) = elemid
 stress(2) = clock(1,cangle)
 DO  i=1,45
   stress(2+i) = clock(i+1,cangle)
 END DO
 
 IF(ktype == 2) GO TO 40
 IF(cangle == 14) GO TO 100
 IF(iclock(1,cangle+1) == 1) GO TO 100
 GO TO 70
 
 40 CONTINUE
 
!  OUTPUT FORCES FOR THIS ANGLE
 isavef(1)=elemid
 savef(2) = BLOCK (1,nangle)
 DO  i=1,16
   savef (2+i) = BLOCK (46+i, nangle)
 END DO
 
! OUTPUT STRESSES
 isaves(1) = elemid
 saves(2) = BLOCK(1,nangle)
 DO  i=1,45
   saves(2+i) = BLOCK(i+1,nangle)
 END DO
 
 IF (nangle == 14) GO TO 100
 IF (iblock(1,nangle+1) == 1)  GO TO 100
 70 CONTINUE
 
 RETURN
 100 again = .false.
 RETURN
END SUBROUTINE stpax3
