SUBROUTINE dsmg2
!*****
! THIS MODULE PERFORMS THE FOLLOWING MATRIX OPERATIONS...
 
!     KBAA  = KAA + BETA * KDAA
!     KBFS  = KFS + BETA * KDFS
!     KBSS  = KSS + BETA * KDSS
!     PBL   = BETA * PL
!     PBS   = BETA * PS
!     YBS   = BETA * YS
!     UBOOV = BETA * UOOV
 
! THE VALUE OF BETA LIES IN THE DIFFERENTIAL STIFFNESS COEFFICIENT SET
! NO. SPECIFIED BY THE INPUT PARAMETER DSCSET.  THE PARTICULAR VALUE OF
! BETA TO BE USED ON ANY PASS THROUGH THIS MODULE IS THE NDSKIP(TH)
! VALUE IN THE DSCSET SET.  IREPTD IS SET EQUAL TO -1 AFTER THE LAST
! BETA IN THE SET HAS BEEN ENCOUNTERED, THEREBY TERMINATING THE
! DIFFERENTIAL STIFFNESS RIGID FORMAT DMAP LOOP BEGINNING AT THE
! LABEL DSLOOP AND ENDING AT THE REPT DSLOOP,100$ STATEMENT.
 
 
!  DMAP CALL...
 
 
!     DSMG2    MPT,KAA,KDAA,KFS,KDFS,KSS,KDSS,PL,PS,YS,UOOV/KBAA,KBFS,
!              KBSS,PBL,PBS,YBS,UBOOV/V,N,NDSKIP/V,N,REPEATD/
!              V,N,DSCOSET/ $
!*****
 INTEGER :: mpt                ,setno  &
     ,                  kaa                ,kdaa  &
     ,                  kfs                ,kdfs  &
     ,                  kss                ,kdss  &
     ,                  pl                 ,ps  &
     ,                  ys                 ,uoov  &
     ,                  kbaa               ,kbfs  &
     ,                  kbss               ,pbl  &
     ,                  pbs                ,ybs  &
     ,                  uboov              ,sysbuf  &
     ,                  buffr1             ,buffr2  &
     ,                  file1              ,file2  &
     ,                  file3              ,clsrw  &
     ,                  dsnos              ,dscset
 
 
 
 DIMENSION mcb(7)             ,dsnos(2)  &
     ,                  NAME(2)            ,BLOCK(11)  &
     ,                  iblock(11)
 
 
 
 EQUIVALENCE (BLOCK(1),iblock(1))  &
     ,                  (betasp,ibeta)
 
! MODULE PARAMETERS
 
 COMMON /BLANK/ ndskip             ,ireptd  &
     ,                  dscset
 
! SYSTEM COMMON
 
 COMMON   /system/ sysbuf
 
 
 
 COMMON   /zzzzzz/ z(1)
 
 
 
 DATA NAME(1)/4HDSMG/    ,NAME(2)/4H2   /
 DATA mpt  &
     ,                  kaa                ,kdaa  &
     ,                  kfs                ,kdfs  &
     ,                  kss                ,kdss  &
     ,                  pl                 ,ps  &
     ,                  ys                 ,uoov  &
     ,                  kbaa               ,kbfs  &
     ,                  kbss               ,pbl  &
     ,                  pbs                ,ybs ,                  uboov /  &
     101,102,103,104,105,106,107,108,109,110, 111,201,202,203,204,205,206,207 /
 DATA clsrw/1/
 DATA neor,dsnos(1),dsnos(2) / 0,53,10/
 
 
 
 izmax  = korsz (z)
 buffr1 = izmax  - sysbuf
 buffr2 = buffr1 - sysbuf
 left   = buffr2 - 1
 
! TURN DIFFERENTIAL STIFFNESS LOOPING FLAG ON AND INCREMENT THE INDEX
! OF BETA.  NOTE THAT NDSKIP MUST BE SET TO ZERO IN THE MODULE
! PROPERTIES TABLE.
 
 ireptd = 1
 ndskip = ndskip + 1
 
! CALL LOCATE TO FIND THE RECORD OF THE MPT WHERE THE DSFACT CARDS ARE.
! THIS IS DONE ONLY IF A D.S. COEFFICIENT SET NO. IS GIVEN.
 
 IF (dscset /= (-1)) GO TO 5
 
! THERE IS NO LOOPING.  TURN LOOPING INDICATOR OFF.
! SEE COMMENTS ABOVE FORTRAN STATEMENT NO. 70 RE THE 4 SCALAR MULTIPLI-
! CATIONS WHEN DSCSET = -1.
 
 ireptd = -1
 GO TO 165
 5 CALL preloc(*1030,z(buffr1),mpt)
 CALL locate(*1035,z(buffr1),dsnos,idummy)
 
 
 
 10 CALL READ(*1040,*1050,mpt,setno,1,neor,idummy)
 IF (setno == dscset) GO TO 30
 
! READ ONE WORD AT A TIME UNTIL A -1 (END OF SET INDICATOR) IS READ.
 
 DO  i = 1,32000
   CALL READ(*1060,*1070,mpt,j,1,neor,idummy)
   IF (j == (-1)) GO TO 10
 END DO
 CALL mesage (-30,84,1)
 
! TEST TO SEE IF WORDS MUST BE SKIPPED.
 
 30 IF (ndskip == 1) GO TO 40
 
! SKIP NDSKIP - 1 WORDS
 
 CALL READ(*1080,*1090,mpt,0,-(ndskip-1),neor,idummy)
 
! READ THE VALUE OF BETA
 
 40 CALL READ(*1100,*1110,mpt,betasp,1,neor,idummy)
 IF (ibeta == (-1)) CALL mesage (-30,84,2)
 
! IF THE NEXT WORD IS A -1, WE HAVE READ THE LAST BETA.  HENCE SET
! IREPTD = -1
 
 CALL READ(*1120,*1130,mpt,j,1,neor,iflag)
 IF (j == (-1)) ireptd = -1
 CALL CLOSE (mpt,clsrw)
 
! PERFORM THE 4 SCALAR MULTIPLICATIONS.  N.B.---IF DSCSET = -1, THAT IS,
! ONLY ONE BETA WILL BE USED AND THAT HAS AN ASSUMED VALUE OF 1.0, IT IS
! ASSUMED THAT EQUIVALENCES HAVE BEEN MADE BETWEEN PL AND PBL, PS AND
! PBS, YS AND YBS, AND UOOV AND UBOOV.
 
 ind = 0
 file1 = pl
 file2 = pbl
 70 mcb(1) = file1
 CALL rdtrl (mcb)
 
! A FATAL ERROR OCCURS IF PL IS PURGED.
 
 IF (mcb(1) < 0  .AND.  ind == 0) CALL mesage (-1,file1,NAME)
 
! IF THE INPUT FILE IS NOT PURGED AND IS NOT PL, SKIP THE OPERATION.
 
 IF (mcb(1) < 0) GO TO 130
 
! THE INPUT FILE IS NOT PURGED. IF THE OUTPUT FILE IS PURGED, A FATAL
! ERROR OCCURS.
 
 mcb(1) = file2
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) CALL mesage (-1,file2,NAME)
 iblock(1) = 1
 BLOCK (2) = betasp
 BLOCK (3) = 0.0
 BLOCK (4) = 0.0
 BLOCK (5) = 0.0
 BLOCK (6) = 0.0
 iblock(7) = 1
 BLOCK (8) = 0.0
 BLOCK (9) = 0.0
 BLOCK(10) = 0.0
 BLOCK(11) = 0.0
 CALL ssg2c (file1,0,file2,0,BLOCK(1))
 130 ind = ind + 1
 SELECT CASE ( ind )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 160
   CASE (    4)
     GO TO 170
 END SELECT
 140 file1= ps
 file2 = pbs
 GO TO 70
 150 file1 = ys
 file2 = ybs
 GO TO 70
 160 file1 = uoov
 file2 = uboov
 GO TO 70
 
! PERFORM MATRIX ADDITIONS
 
 165 betasp= 1.0
 170 file1 = kaa
 file2 = kdaa
 file3 = kbaa
 ind   = 0
 180 iundef= 0
 mcb(1)= file1
 CALL rdtrl (mcb(1))
 IF (mcb(1) < 0)  IF (ind) 190,190,200
 GO TO 210
 190 CALL mesage (-1,file1,NAME)
 200 iundef = 1
 210 mcb(1) = file2
 CALL rdtrl (mcb(1))
 IF (mcb(1) < 0)  IF (iundef) 260,260,220
 IF (iundef == 1) CALL mesage (-30,84,15)
 iblock(1) = 1
 BLOCK (2) = 1.0
 BLOCK (3) = 0.0
 BLOCK (4) = 0.0
 BLOCK (5) = 0.0
 BLOCK (6) = 0.0
 iblock(7) = 1
 BLOCK (8) = betasp
 BLOCK (9) = 0.0
 BLOCK(10) = 0.0
 BLOCK(11) = 0.0
 CALL ssg2c (file1,file2,file3,0,BLOCK(1))
 220 ind = ind + 1
 SELECT CASE ( ind )
   CASE (    1)
     GO TO 230
   CASE (    2)
     GO TO 240
   CASE (    3)
     GO TO 250
 END SELECT
 230 file1 = kfs
 file2 = kdfs
 file3 = kbfs
 GO TO 180
 240 file1 = kss
 file2 = kdss
 file3 = kbss
 GO TO 180
 250 RETURN
 260 CALL mesage (-1,file2,NAME )
 1030 i = 3
 GO TO 1099
 1035 i = 4
 GO TO 1099
 1040 i = 5
 GO TO 1099
 1050 i = 6
 GO TO 1099
 1060 i = 7
 GO TO 1099
 1070 i = 8
 GO TO 1099
 1080 i = 9
 GO TO 1099
 1090 i = 10
 GO TO 1099
 1100 i = 11
 GO TO 1099
 1110 i = 12
 GO TO 1099
 1120 i = 13
 GO TO 1099
 1130 i = 14
 1099 CALL mesage (-30,84,i)
 RETURN
END SUBROUTINE dsmg2
