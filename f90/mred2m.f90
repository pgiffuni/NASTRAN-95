SUBROUTINE mred2m (nuf,n2,nus)
     
!     THIS SUBROUTINE FORMS THE HK MATRIX FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN OUT)                  :: nuf
 INTEGER, INTENT(IN OUT)                  :: n2
 INTEGER, INTENT(IN)                      :: nus
 INTEGER :: dry,gbuf1,otfile,z,typin,typout,hkpg,phi12i,hk, rprtn,cprtn
 DIMENSION       itrlr1(7),modnam(2),itrlr2(7),isub(4),rz(1)
 COMMON /BLANK / idum1,dry,idum4,gbuf1,idum2(17),otfile(6),  &
     iscr(10),korlen,korbgn,idum3(14),nmodes
 COMMON /zzzzzz/ z(1)
 COMMON /packx / typin,typout,irow,nrow,incr
 EQUIVALENCE     (hk,iscr(2)),(ident,iscr(8)),(hkpg,iscr(3)),  &
     (phi12i,iscr(8)),(cprtn,iscr(9)),(rprtn,iscr(9))
 EQUIVALENCE     (rz(1),z(1))
 DATA    modnam/ 4HMRED,4H2M  /
 
!     FORM HK MATRIX
 
!        **  **   **             **
!        *    *   *      .        *
!        * HK * = * HKPG . PHI12I *
!        *    *   *      .        *
!        **  **   **             **
 
 IF (dry == -2) RETURN
 IF (nuf ==  0) GO TO 30
 itrlr1(1) = hkpg
 CALL rdtrl (itrlr1)
 itrlr2(1) = phi12i
 CALL rdtrl (itrlr2)
 incr    = 1
 typin   = 1
 typout  = 1
 irow    = 1
 nrow    = itrlr1(3) + itrlr2(3)
 isub(1) = itrlr1(3)
 isub(2) = itrlr2(3)
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > itrlr1(3)) rz(korbgn+i-1) = 1.0
 END DO
 iform = 7
 CALL makmcb (itrlr2,rprtn,nrow,iform,typin)
 CALL gopen (rprtn,z(gbuf1),1)
 CALL pack (z(korbgn),rprtn,itrlr2)
 CALL CLOSE (rprtn,1)
 CALL wrttrl (itrlr2)
 itype = 2
 CALL gmmerg (hk,hkpg,0,phi12i,0,rprtn,0,isub,itype,z(korbgn), korlen)
 RETURN
 
!     NO UF POINTS
 
!        **  **   **     **
!        *    *   *   .   *
!        * HK * = * 0 . I *
!        *    *   *   .   *
!        **  **   **     **
 
 30 typin  = 1
 typout = 1
 irow   = 1
 nrow   = nmodes
 incr   = 1
 iform  = 8
 IF (korbgn+nmodes >= korlen) GO TO 100
 
!     GENERATE IDENTITY MATRIX
 
 CALL makmcb (itrlr2,ident,nmodes,iform,typin)
 CALL gopen  (ident,z(gbuf1),1)
 DO  i = 1,nmodes
   DO  j = 1,nmodes
     rz(korbgn+j-1) = 0.0
     IF (j == i) rz(korbgn+j-1) = 1.0
   END DO
   CALL pack (z(korbgn),ident,itrlr2)
 END DO
 CALL CLOSE (ident,1)
 CALL wrttrl (itrlr2)
 
!     GENERATE ROW PARTITIONING VECTOR
 
 nrow = nus + nmodes
 IF (korbgn+nrow >= korlen) GO TO 100
 j = nrow
 DO  i = 1,j
   rz(korbgn+i-1) = 0.0
   IF (i > nus) rz(korbgn+i-1) = 1.0
 END DO
 iform = 7
 CALL makmcb (itrlr2,rprtn,nrow,iform,typin)
 CALL gopen (rprtn,z(gbuf1),1)
 CALL pack (z(korbgn),rprtn,itrlr2)
 CALL CLOSE (rprtn,1)
 CALL wrttrl (itrlr2)
 
!     FORM HK MATRIX
 
 isub(1) = nus
 isub(2) = nmodes
 itype   = 2
 CALL gmmerg (hk,0,0,ident,0,rprtn,0,isub,itype,z(korbgn),korlen)
 RETURN
 
!     PROCESS SYSTEM ERRORS
 
 100 imsg  =-8
 ifile = 0
 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
END SUBROUTINE mred2m
