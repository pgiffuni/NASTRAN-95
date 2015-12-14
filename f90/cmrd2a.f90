SUBROUTINE cmrd2a
     
!     THIS SUBROUTINE PARTITIONS THE STIFFNESS MATRIX INTO BOUNDARY AND
!     INTERIOR POINTS AND THEN SAVES THE PARTITIONING VECTOR ON THE SOF
!     AS THE UPRT ITEM FOR THE CMRED2 MODULE.
 
!     INPUT DATA
!     GINO - USETMR - USET TABLE FOR REDUCED SUBSTRUCTURE
!            KAA    - SUBSTRUCTURE STIFFNESS MATRIX
 
!     OUTPUT DATA
!     GINO - KBB  - KBB PARTITION MATRIX
!            KIB  - KIB PARTITION MATRIX
!            KII  - KII PARTITION MATRIX
!     SOF  - UPRT - PARTITION VECTOR FOR ORIGINAL SUBSTRUCTURE
 
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
!              KBI    - KBI OUTPUT FILE NUMBER
!              KII    - KII OUTPUT FILE NUMBER
!              UPRT   - KAA PARTITION VECTOR FILE NUMBER
 
 INTEGER :: dry,gbuf1,otfile,oldnam,z,un,ub,ui,fuset,usetmr, uprt,modnam(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / idum1,dry,idum6,gbuf1,idum2(5),infile(11),  &
     otfile(6),iscr(11),korlen,korbgn,oldnam(2)
 COMMON /zzzzzz/ z(1)
 COMMON /bitpos/ idum4(9),un,idum5(10),ub,ui
 COMMON /patx  / lcore,nsub(3),fuset
 COMMON /system/ idum3,iprntr
 EQUIVALENCE     (usetmr,infile(6)),(kaa,infile(7)),  &
     (kbb,iscr(1)),(kib,iscr(2)),(kii,iscr(4)), (kbi,iscr(3)),(uprt,iscr(5))
 DATA    modnam/ 4HCMRD,4H2A  /
 DATA    item  / 4HUPRT       /
 
!     SET UP PARTITIONING VECTOR
 
 IF (dry == -2) RETURN
 lcore = korlen
 fuset = usetmr
 CALL calcv(uprt,un,ui,ub,z(korbgn))
 
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
 
 CALL gmprtn(kaa,kii,kbi,kib,kbb,uprt,uprt,nsub(1),nsub(2), z(korbgn),korlen)
 
!     SAVE PARTITIONING VECTOR
 
 CALL mtrxo(uprt,oldnam,  item,0,itest)
 IF (itest /= 3) GO TO 30
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 30 SELECT CASE ( itest )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 40
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 50
   CASE (    5)
     GO TO 60
   CASE (    6)
     GO TO 80
 END SELECT
 40 WRITE (iprntr,900) ufm,modnam,item,oldnam
 dry = -2
 RETURN
 50 imsg = -2
 GO TO 70
 60 imsg = -3
 70 CALL smsg(imsg,item,oldnam)
 RETURN
 
 80 WRITE (iprntr,901) ufm,modnam,item,oldnam
 dry = -2
 RETURN
 
 900 FORMAT (a23,' 3211, MODULE ',2A4,8H - item ,a4,  &
     ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
 901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR ',  &
     'I/O OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURGED.')
 
END SUBROUTINE cmrd2a
