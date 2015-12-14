SUBROUTINE mred2o (nus)
     
!     THIS SUBROUTINE FORMS THE M MATRIX FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN OUT)                  :: nus
 INTEGER :: dry,gbuf1,gbuf2,sbuf1,sbuf2,sbuf3,otfile,z,  &
     typin,typout,prec,typea,typeb,hk,gs,rprtn,hm, gszero,dblkor
 DOUBLE PRECISION :: dz
 DIMENSION       itrlr1(7),itrlr2(7),modnam(2),isub(4), rz(1),BLOCK(11),dz(1)
 COMMON /BLANK / idum1,dry,idum2,gbuf1,gbuf2,idum3,sbuf1,sbuf2,  &
     sbuf3,infile(12),otfile(6),iscr(10),korlen,korbgn,  &
     idum5(14),nmodes,idum4(2),lstzwd
 COMMON /zzzzzz/ z(1)
 COMMON /mpy3tl/ itrlra(7),itrlrb(7),itrlre(7),itrlrc(7),jscr(3),  &
     lkore,icode,prec,dummy(13)
 COMMON /packx / typin,typout,irow,nrow,incr
 EQUIVALENCE     (lamamr,infile(2)),(gs,iscr(7)),(dz(1),z(1)),  &
     (hk,iscr(2)),(kmw2,iscr(5)),(hm,iscr(9)),  &
     (gszero,iscr(10)),(m,iscr(10)),(rprtn,iscr(8)),  &
     (rz(1),z(1)),(typea,BLOCK(1)),(typeb,BLOCK(7))
 DATA    modnam/ 4HMRED,4H2O  /
 
!     FORM HM MATRIX
 
!        **  **   **  **   **          **
!        *    *   *    *   *    .   .   *
!        * HM * = * HK * + * GS . 0 . 0 *
!        *    *   *    *   *    .   .   *
!        **  **   **  **   **          **
 
 IF (dry == -2) RETURN
 IF (nus ==  0) GO TO 60
 
!     GENERATE ROW PARTITION VECTOR
 
 itrlr1(1) = hk
 CALL rdtrl (itrlr1)
 itrlr2(1) = gs
 CALL rdtrl (itrlr2)
 typin = 1
 typout= 1
 irow  = 1
 nrow  = itrlr1(2)
 incr  = 1
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > itrlr2(2)) rz(korbgn+i-1) = 1.0
 END DO
 iform = 7
 CALL makmcb (itrlr2,rprtn,nrow,iform,typin)
 CALL gopen  (rprtn,z(gbuf1),1)
 CALL pack   (z(korbgn),rprtn,itrlr2)
 CALL CLOSE  (rprtn,1)
 CALL wrttrl (itrlr2)
 
!     MERGE GS, ZERO MATRICES
 
 isub(1) = itrlr2(2)
 isub(2) = itrlr1(2) - itrlr2(2)
 CALL gmmerg (gszero,gs,0,0,0,rprtn,0,isub,itrlr2(5),z(korbgn), korlen)
 
!     FORM HM MATRIX
 
 itrlr2(1) = gszero
 CALL rdtrl (itrlr2)
 DO  i = 1,11
   BLOCK(i) = 0.0
 END DO
 BLOCK(2) = 1.0
 BLOCK(8) = 1.0
 typea = itrlr1(5)
 typeb = itrlr2(5)
 iop   = 1
 CALL sofcls
 CALL ssg2c (hk,gszero,hm,iop,BLOCK)
 GO TO 70
 
!     IF NO US POINTS
 
!        **  **   **  **
!        *    *   *    *
!        * HM * = * HK *
!        *    *   *    *
!        **  **   **  **
 
 60 hm = hk
 CALL sofcls
 
!     FORM KMW2 = M  MATRIX
!                  I
 
 70 ifile = lamamr
 CALL gopen  (lamamr,z(gbuf1),0)
 CALL fwdrec (*140,lamamr)
 iform = 3
 itype = 1
 CALL makmcb (itrlr1,kmw2,nmodes,iform,itype)
 typin = 1
 typout= 1
 irow  = 1
 nrow  = nmodes
 incr  = 1
 CALL gopen (kmw2,z(gbuf2),1)
 DO  i = 1,nmodes
   CALL READ (*130,*140,lamamr,z(korbgn),7,0,nwdsrd)
   DO  j = 1,nmodes
     rz(korbgn+7+j-1) = 0.0
     IF (j == i) rz(korbgn+7+j-1) = rz(korbgn+5)
   END DO
   CALL pack  (z(korbgn+7),kmw2,itrlr1)
 END DO
 CALL CLOSE (lamamr,1)
 CALL CLOSE (kmw2,1)
 CALL wrttrl (itrlr1)
 
!     FORM M MATRIX
 
!                      T
!                **  ** **     ** **  **
!        ** **   *    * * .     * *    *
!        *   *   *    * *  .    * *    *
!        * M * = * HM * *   M   * * HM *     WHERE M = M
!        *   *   *    * *    .  * *    *                I
!        ** **   *    * *     . * *    *
!                **  ** **     ** **  **
 
 itrlr1(1) = hm
 itrlr2(1) = kmw2
 CALL rdtrl (itrlr1)
 CALL rdtrl (itrlr2)
 DO  i = 1,7
   itrlra(i) = itrlr1(i)
   itrlrb(i) = itrlr2(i)
   itrlre(i) = 0
 END DO
 iprc = 1
 ityp = 0
 IF ((itrlra(5) == 2) .OR. (itrlra(5) == 4)) iprc = 2
 IF ((itrlrb(5) == 2) .OR. (itrlrb(5) == 4)) iprc = 2
 IF (itrlra(5) >= 3) ityp = 2
 IF (itrlrb(5) >= 3) ityp = 2
 itype = iprc + ityp
 iform = 6
 CALL makmcb (itrlrc,m,itrlr1(3),iform,itype)
 jscr(1) = iscr(7)
 jscr(2) = iscr(8)
 jscr(3) = iscr(6)
 icode   = 0
 prec    = 0
 dblkor  = (korbgn/2) + 1
 lkore   = lstzwd - (2*dblkor - 1)
 CALL mpy3dr (dz(dblkor))
 CALL wrttrl (itrlrc)
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 130 imsg = -2
 GO TO 150
 140 imsg = -3
 150 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
END SUBROUTINE mred2o
