SUBROUTINE mred2p (nus,nuf,n2)
     
!     THIS SUBROUTINE OUTPUTS THE HAB MATRIX TO THE SOF AS THE HORG ITEM
!     FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN)                      :: nus
 INTEGER, INTENT(IN)                      :: nuf
 INTEGER, INTENT(IN)                      :: n2
 INTEGER :: dry,gbuf1,otfile,z,typin,typout,hab
 DIMENSION       rz(1),modnam(2),itrlr1(7)
 COMMON /BLANK / idum1,dry,idum2,gbuf1,idum3(17),otfile(6),  &
     iscr(10),korlen,korbgn,oldnam(2),idum4(12),nmodes
 COMMON /zzzzzz/ z(1)
 COMMON /packx / typin,typout,irow,nrow,incr
 COMMON /system/ idum5,iprntr
 EQUIVALENCE     (hab,iscr(2)),(rz(1),z(1))
 DATA    modnam/ 4HMRED,4H2P  /
 DATA    item  / 4HHORG/
 
!     FORM HAB MATRIX
 
!        **   **   **     **
!        *     *   *   .   *
!        * HAB * = * I . 0 *
!        *     *   *   .   *
!        **   **   **     **
 
 IF (dry == -2) GO TO 160
 kolmns = nus + nuf + n2
 IF (n2 == 0) kolmns = kolmns + (nmodes - nuf)
 typin  = 1
 typout = 1
 irow   = 1
 nrow   = nus + nuf
 incr   = 1
 iform  = 2
 CALL makmcb (itrlr1,hab,nrow,iform,typin)
 CALL gopen (hab,z(gbuf1),1)
 DO  i = 1,kolmns
   DO  j = 1,nrow
     rz(korbgn+j-1) = 0.0
     IF (i > nus+nuf) CYCLE
     IF (j == i) rz(korbgn+j-1) = 1.0
   END DO
   CALL pack (z(korbgn),hab,itrlr1)
 END DO
 CALL CLOSE (hab,1)
 CALL wrttrl (itrlr1)
 
!     STORE HAB MATRIX AS HORG ON SOF
 
 CALL mtrxo (hab,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 70
 GO TO 160
 
!     PROCESS MODULE FATAL ERRORS
 
 70 SELECT CASE ( itest )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 90
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 110
   CASE (    5)
     GO TO 120
   CASE (    6)
     GO TO 140
 END SELECT
 80 imsg = -9
 GO TO 150
 90 imsg = -11
 GO TO 150
 100 imsg = -1
 GO TO 130
 110 imsg = -2
 GO TO 130
 120 imsg = -3
 130 CALL smsg (imsg,item,oldnam)
 GO TO 160
 140 imsg = -10
 150 dry = -2
 CALL smsg1 (imsg,item,oldnam,modnam)
 160 RETURN
END SUBROUTINE mred2p
