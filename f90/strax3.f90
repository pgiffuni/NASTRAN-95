SUBROUTINE strax3 ( again)
     
 
 LOGICAL, INTENT(IN OUT)                  :: again
 INTEGER :: iforce(25), istres(100), elemid, iblock(22,14)  &
     ,       iclock(22,14),cangle
 REAL :: savef(75), saves(75)
 
 
 
 COMMON /sdr2x7/ dum(100),stress(100),force(25)  &
     ,               skip(25),BLOCK(22,14),clock(22,14)
 
! SCRATCH BLOCK
 COMMON  /sdr2x8/ disp(9), eforc(9), estres(9)  &
     ,                harm, n, sinphi, conphi, nphi, nangle  &
     ,                elemid, unu(123), nelhar
 
 COMMON /isave / isavef(75),isaves(75)
 
 COMMON /sdr2de/ dum5(33), ipart
 
 COMMON /sdr2x4/ dum4(51),ktype
 
 EQUIVALENCE ( istres(1), stress(1)), ( iforce(1), force(1))  &
     , (iblock(1,1), BLOCK(1,1)), (iclock(1,1),clock(1,1))  &
     , (isavef(1),savef(1)),(isaves(1),saves(1)) , (nangle,cangle)
 
 IF ( again ) GO TO 10
 again = .true.
 nangle = 0
 10 nangle = nangle + 1
 
 
!  BRANCH TO INSERT STRESSES AND FORCES INTO FORCE AND STRESS OR
!                                            SAVEF AND SAVES
 
!    KTYPE=1 - REAL OUTPUT FROM BLOCK IS TRANSFERED TO CLOCK, THEN
!            STORED IN FORCE AND STRESS, NOTHING IN SAVEF AND SAVES
!    KTYPE=2 - COMPLEX OUTPUT
!    IPART=1 - IMAGINARY PART OF COMPLEX OUTPUT FROM BLOCK, STORED
!              IN SAVEF AND SAVES
!    IPART=2 - REAL PART OF COMPLEX OUTPUT FROM CLOCK STORED IN
!              FORCE AND STRESS
 
 IF (ktype == 2) GO TO 19
 DO  i=1,22
   DO  j=1,14
     clock(i,j) = BLOCK(i,j)
   END DO
 END DO
 19 CONTINUE
 
!  OUTPUT FORCES FOR THIS ANGLE
 
 iforce(1)=elemid
 force(2) = clock(1,cangle)
 j = 2
 DO  i=1,9
   j = j + 1
   force(j) = clock(i+7,cangle)
   
!  OUTPUT CHARGES
   IF((i /= 3).AND.(i /= 6).AND.(i /= 9)) CYCLE
   j = j + 1
   k=19+i/3
   force(j) = clock(k,cangle)
 END DO
 
! OUTPUT STRESSES
 istres  (1) = elemid
 stress(2) = clock(1,cangle)
 DO  i = 1, 6
   stress(2+i) = clock(i+1,cangle)
 END DO
 
!  OUTPUT FLUXES
 DO  i=1,3
   stress(i+8) = clock(i+16,cangle)
 END DO
 
 IF(ktype == 2) GO TO 40
 IF(cangle == 14) GO TO 100
 IF(iclock(1,cangle+1) == 1) GO TO 100
 GO TO 90
 
 40 CONTINUE
 
!  OUTPUT FORCES FOR THIS ANGLE
 
 isavef(1)=elemid
 savef(2) = BLOCK (1,nangle)
 j = 2
 DO  i=1,9
   j = j + 1
   savef(j) = BLOCK(i+7,nangle)
   
!  OUTPUT CHARGES
   IF((i /= 3).AND.(i /= 6).AND.(i /= 9)) CYCLE
   j = j + 1
   k=19+i/3
   savef(j) = BLOCK(k,nangle)
 END DO
 
! OUTPUT STRESSES
 isaves(1) = elemid
 saves(2) = BLOCK (1, nangle)
 DO  i=1,6
   saves(2+i) = BLOCK ( i+1,nangle)
 END DO
 
!  OUTPUT FLUXES
 DO  i=1,3
   saves(i+8) = BLOCK(i+16,nangle)
 END DO
 
 IF (nangle == 14) GO TO 100
 IF (iblock (1,nangle+1) == 1)  GO TO 100
 90 CONTINUE
 
 RETURN
 100 again = .false.
 RETURN
END SUBROUTINE strax3
