SUBROUTINE mred1d
     
!     THIS SUBROUTINE GENERATES THE EEDX DATA BLOCK USING THE EED DATA
!     BLOCK FORMAT FROM THE EIGR OR EIGC AND EIGP BULK DATA FOR THE
!     MRED1 MODULE.
 
!     INPUT  DATA
!     GINO - DYNAMICS - EIGC DATA
!                       EIGP DATA
!                       EIGR DATA
 
!     OUTPUT DATA
!     GINO - EEDX     - EIGC DATA
!                       EIGP DATA
!                       EIGR DATA
 
!     PARAMETERS
!     INPUT  - DNAMIC - DYNAMICS DATA BLOCK INPUT FILE NUMBER
!              GBUF1  - GINO BUFFER
!              EEDX   - EEDX DATA BLOCK OUTPUT FILE NUMBER
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              IEIG   - EIGENVALUE EXTRACTION SET IDENTIFICATION NUMBER
!     OUTPUT - DRY    - MODULE OPERATION FLAG
!     OTHERS - EIGTYP - EIG CARD TYPE PROCESSING FLAG
!                     = 1, PROCESS EIGC DATA
!                     = 2, PROCESS EIGP DATA
!                     = 3, PROCESS EIGR DATA
!              EIGCP  - EIGC AND EIGP DATA ERROR FLAG
!                     = 0, NO EIGC, EIGP DATA - NO ERROR
!                     = 1, EIGC DATA ONLY - NO ERROR
!                     = 2, EIGP DATA ONLY - ERROR
!                     = 3, EIGC AND EIGP DATA - NO ERROR
!              EIGTRL - EEDX TRAILER
!              EIGCPR - DUMMY EIG(C,P,R) ARRAY
!              EIG    - ARRAY OF EIG(C,P,R) CARD TYPES AND HEADER
!                       INFORMATION
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              NWDS2R - NUMBER OF EIG(C,P,R) WORDS TO READ ON DYNAMIC
!                       DATA FILE
 
 EXTERNAL        orf
 LOGICAL :: usrmod
 INTEGER :: orf,oldnam,dry,TYPE,gbuf1,gbuf2,z,dnamic,  &
     eig(3,3),eigcpr(3),eedx,eigtrl(7),eigtyp,eigcp
 DIMENSION       modnam(2),letr(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / oldnam(2),dry,idum1(3),TYPE(2),idum5,gbuf1,gbuf2,  &
     idum2(3),korlen,idum7(4),ieig,idum3(6),korbgn, idum6(12),usrmod
 COMMON /zzzzzz/ z(1)
 COMMON /two   / itwo(32)
 COMMON /system/ idum4,iprntr
 DATA    dnamic, eig,eedx/103,207,2,0,257,4,0,307,3,0,202/
 DATA    modnam, letr /4HMRED,4H1D  ,1HC,1HP,1HR/
 DATA    komplx, kreal/4HCOMP,4HREAL/
 
!     OPEN DYNAMICS, EEDX DATA BLOCKS
 
 IF (dry == -2) RETURN
 IF (usrmod) GO TO 175
 CALL preloc (*180,z(gbuf1),dnamic)
 CALL gopen (eedx,z(gbuf2),1)
 
!     SET PROCESSING FLAGS
 
 eigtyp = 0
 eigcp  = 0
 eigtrl(1) = eedx
 DO  i = 2,7
   eigtrl(i) = 0
 END DO
 
!     INCREMENT EIG PROCESSING FLAG
!     EIGTYP .EQ. 1, PROCESS EIGC DATA
!     EIGTYP .EQ. 2, PROCESS EIGP DATA
!     EIGTYP .EQ. 3, PROCESS EIGR DATA
 
 20 eigtyp = eigtyp + 1
 IF (eigtyp == 4) GO TO 170
 
!     SELECT EIG MODE
 
 IF (TYPE(1) == kreal  .AND. eigtyp < 3) GO TO 20
 IF (TYPE(1) == komplx .AND. eigtyp == 3) GO TO 20
 DO  i = 1,3
   eigcpr(i) = eig(i,eigtyp)
 END DO
 
!     LOCATE EIG(C,P,R) DATA CARD
 
 CALL locate (*20,z(gbuf1),eigcpr,itest)
 
!     SET UP EEDX DATA RECORD
 
 DO  i = 1,3
   z(korbgn+i-1) = eigcpr(i)
 END DO
 
!     FIND CORRECT EIG(C,P,R) DATA CARD
 
 SELECT CASE ( eigtyp )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 70
 END SELECT
 50 nwds2r = 10
 GO TO 80
 60 nwds2r = 4
 GO TO 80
 70 nwds2r = 18
 80 CALL READ (*190,*200,dnamic,z(korbgn+3),nwds2r,0,nowdsr)
 IF (z(korbgn+3) == ieig) GO TO 100
 SELECT CASE ( eigtyp )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 80
 END SELECT
 
!     READ REST OF EIGC DATA
 
 90 CALL READ (*190,*200,dnamic,z(korbgn+3),7,0,nowdsr)
 IF (z(korbgn+3) == -1) GO TO 80
 GO TO 90
 
!     SELECT EIG PROCESSING MODE
 
 100 SELECT CASE ( eigtyp )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 140
   CASE (    3)
     GO TO 150
 END SELECT
 
!     WRITE EIGC DATA ONTO EEDX DATA BLOCK
 
 110 CALL WRITE (eedx,z(korbgn),13,0)
 eigtrl(2) = orf(eigtrl(2),16384)
 eigcp = eigcp + 1
 120 CALL READ (*190,*200,dnamic,z(korbgn),7,0,nowdsr)
 IF (z(korbgn) == -1) GO TO 130
 CALL WRITE (eedx,z(korbgn),7,0)
 GO TO 120
 130 CALL WRITE (eedx,z(korbgn),7,1)
 GO TO 20
 
!     WRITE EIGP DATA ONTO EEDX DATA BLOCK
 
 140 CALL WRITE (eedx,z(korbgn),7,1)
 eigcp = eigcp + 2
 eigtrl(2) = orf(eigtrl(2),4096)
 GO TO 20
 
!     WRITE EIGR DATA ONTO EEDX DATA BLOCK
 
 150 CALL WRITE (eedx,z(korbgn),21,1)
 eigtrl(2) = orf(eigtrl(2),8192)
 GO TO 20
 
!     CLOSE DYNAMICS, EEDX DATA BLOCKS
 
 170 CALL CLOSE (dnamic,1)
 CALL CLOSE (eedx,1)
 
!     TEST FOR EIG CARD ERRORS
 
 IF (eigtrl(2) == 0) GO TO 230
 IF (eigcp == 2) GO TO 240
 
!     WRITE EEDX DATA BLOCK TRAILER
 
 CALL wrttrl (eigtrl)
 175 CONTINUE
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 180 imsg = -1
 GO TO 220
 190 imsg = -2
 IF (eigtyp == 2) GO TO 20
 GO TO 210
 200 imsg = -3
 IF (eigtyp == 2) GO TO 20
 210 WRITE (iprntr,900) ufm,letr(eigtyp),ieig,oldnam
 220 CALL sofcls
 CALL mesage (imsg,dnamic,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 230 WRITE (iprntr,901) ufm,ieig,oldnam
 GO TO 250
 240 WRITE (iprntr,902) ufm,ieig,oldnam
 250 dry = -2
 RETURN
 
 900 FORMAT (a23,' 6627, NO EIG',a1,' DATA CARD ',  &
     'SPECIFIED FOR SET ID',i9,', SUBSTRUCTURE ',2A4,1H.)
 901 FORMAT (a23,' 6628, NO EIGC OR EIGR CARD SPECIFIED FOR SET ID',i9,  &
     ', SUBSTRUCTURE ',2A4,1H.)
 902 FORMAT (a23,' 6629, NO EIGC DATA CARD SPECIFHIED WITH EIGP DATA ',  &
     'CARD SET ID',i9,', SUBSTRUCTURE ',2A4,1H.)
 
END SUBROUTINE mred1d
