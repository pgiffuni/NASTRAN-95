SUBROUTINE mred2a
     
!     THIS SUBROUTINE PARTITIONS THE STIFFNESS MATRIX INTO BOUNDARY AND
!     INTERIOR POINTS AND THEN SAVES THE PARTITIONING VECTOR ON THE SOF
!     AS THE UPRT ITEM FOR THE MRED2 MODULE.
 
!     INPUT DATA
!     GINO - USETMR   - USET TABLE FOR REDUCED SUBSTRUCTURE
!            KAA      - SUBSTRUCTURE STIFFNESS MATRIX
 
!     OUTPUT DATA
!     GINO - KBB      - KBB PARTITION MATRIX
!            KIB      - KIB PARTITION MATRIX
!            KII      - KII PARTITION MATRIX
!     SOF  - UPRT     - PARTITION VECTOR FOR ORIGINAL SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT  - GBUF   - GINO BUFFER
!              INFILE - INPUT FILE NUMBERS
!              ISCR   - SCRATCH FILE NUMBERS
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!     OTHERS - USETMR - USETMR INPUT FILE NUMBER
!              KAA    - KAA INPUT FILE NUMBER
!              KBB    - KBB OUTPUT FILE NUMBER
!              KIB    - KIB OUTPUT FILE NUMBER
!              KII    - KII OUTPUT FILE NUMBER
!              UPRT   - KAA PARTITION VECTOR FILE NUMBER
 
 LOGICAL :: bounds
 INTEGER :: dry,gbuf1,oldnam,usrmod,z,un,ub,ui,fuset,  &
     usetmr,uprt,eqst,modnam(2),itrlr(7)
 COMMON /BLANK / idum1,dry,idum6,gbuf1,idum2(5),infile(12),  &
     idum3(6),iscr(10),korlen,korbgn,oldnam(2), idum7(6),usrmod,idum9,bounds
 COMMON /zzzzzz/ z(1)
 COMMON /bitpos/ idum4(9),un,idum5(10),ub,ui
 COMMON /patx  / lcore,nsub(3),fuset
 COMMON /system/ idum8,iprntr
 EQUIVALENCE     (eqst,infile(4)),(usetmr,infile(5)),  &
     (kaa,infile(6)),(kbb,iscr(1)),(kib,iscr(2)), (kii,iscr(3)),(uprt,iscr(5))
 DATA    modnam/ 4HMRED,4H2A  /
 DATA    item  / 4HUPRT/
 
!     LOCATE PARTITIONING VECTOR
 
 IF (dry == -2) GO TO 100
 IF (bounds) GO TO 10
 lcore = korlen
 fuset = usetmr
 CALL calcv (uprt,un,ui,ub,z(korbgn))
 GO TO 20
 10 CALL mtrxi (uprt,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 30
 itrlr(1) = eqst
 CALL rdtrl (itrlr)
 nsub(1) = itrlr(6)
 nsub(2) = itrlr(7)
 
!     PARTITION STIFFNESS MATRIX
 
!                  **         **
!                  *     .     *
!        **   **   * KBB . KBI *
!        *     *   *     .     *
!        * KAA * = *...........*
!        *     *   *     .     *
!        **   **   * KIB . KII *
!                  *     .     *
!                  **         **
 
 20 CONTINUE
 CALL gmprtn (kaa,kii,0,kib,kbb,uprt,uprt,nsub(1),nsub(2), z(korbgn),korlen)
 
!     SAVE PARTITIONING VECTOR
 
 IF (bounds) GO TO 25
 CALL mtrxo (uprt,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 30
 25 CONTINUE
 GO TO 100
 
!     PROCESS MODULE FATAL ERRORS
 
 30 SELECT CASE ( itest )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 45
   CASE (    3)
     GO TO 50
   CASE (    4)
     GO TO 55
   CASE (    5)
     GO TO 60
   CASE (    6)
     GO TO 80
 END SELECT
 40 imsg = -9
 GO TO 90
 45 imsg = -11
 GO TO 90
 50 imsg = -1
 GO TO 70
 55 imsg = -2
 GO TO 70
 60 imsg = -3
 70 CALL smsg (imsg,item,oldnam)
 GO TO 100
 80 imsg = -10
 90 dry = -2
 CALL smsg1 (imsg,item,oldnam,modnam)
 100 RETURN
 
END SUBROUTINE mred2a
