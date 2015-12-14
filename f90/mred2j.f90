SUBROUTINE mred2j (nuf,n2)
     
!     THIS SUBROUTINE PARTITIONS THE PHISS MATRIX FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN)                      :: nuf
 INTEGER, INTENT(OUT)                     :: n2
 INTEGER :: dry,gbuf1,otfile,typin,typout,phiss,rprtn, itrlr1(7),modnam(2)
 COMMON /BLANK / idum1,dry,idum4,gbuf1,idum2(5),infile(12),  &
     otfile(6),iscr(10),korlen,korbgn,idum3(14),nmodes
 COMMON /zzzzzz/ rz(1)
 COMMON /packx / typin,typout,irow,nrow,incr
 EQUIVALENCE     (phiss,infile(3)),(phiss1,iscr(8)),  &
     (phiss2,iscr(9)) ,(rprtn,iscr(10))
 DATA    modnam/ 4HMRED,4H2J  /
 
!     SET UP PARTITIONING VECTOR
 
 IF (dry == -2) RETURN
 typin  = 1
 typout = 1
 irow   = 1
 incr   = 1
 
!     COMMENTS FROM G.CHAN/UNISYS    4/92
!     ORIGINALLY AT THIS POINT, THE FOLLOWING DO 20 LOOP IS IN ERROR
!       1. KOLUMN AND J ARE NOT DEFINED
!       2. NROW AND ITRLR1 ARE ALSO NOT YET DEFINED
 
!     MY BEST GUESS IS THE NEXT 10 LINES THAT FOLLOW
 
 ifile = phiss
 itrlr1(1) = phiss
 CALL rdtrl (itrlr1)
 IF (itrlr1(1) < 0) GO TO 30
 kolumn = itrlr1(2)
 nrow   = itrlr1(3)
 DO  i = 1,kolumn
   rz(korbgn+i-1) = 0.0
   IF (i > nuf) rz(korbgn+i-1) = 1.0
 END DO
 
 iform = 7
 CALL makmcb (itrlr1,rprtn,nrow,iform,itrlr1(5))
 CALL gopen  (rprtn,rz(gbuf1),1)
 CALL pack   (rz(korbgn),rprtn,itrlr1)
 CALL CLOSE  (rprtn,1)
 CALL wrttrl (itrlr1)
 
!     PARTITION PHISS MATRIX
 
!        **     **   **               **
!        *       *   *        .        *
!        * PHISS * = * PHISS1 . PHISS2 *
!        *       *   *        .        *
!        **     **   **               **
 
 itrlr1(1) = phiss
 CALL rdtrl (itrlr1)
 n2 = nmodes - nuf
 CALL gmprtn (phiss,phiss1,0,phiss2,0,rprtn,0,nuf,n2,rz(korbgn), korlen)
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 30 imsg = -1
 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
END SUBROUTINE mred2j
