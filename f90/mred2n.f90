SUBROUTINE mred2n
     
!     THIS SUBROUTINE CALCULATES THE K MATRIX FOR THE MRED2 MODULE.
 
 INTEGER :: dry,gbuf1,gbuf2,sbuf1,sbuf2,sbuf3,otfile,prec,  &
     typin,typout,hk,dblkor
 DOUBLE PRECISION :: dz
 DIMENSION       itrlr1(7),itrlr2(7),modnam(2),rz(1),dz(1)
 COMMON /BLANK / idum1,dry,idum2,gbuf1,gbuf2,idum3,sbuf1,sbuf2,  &
     sbuf3,infile(12),otfile(6),iscr(10),korlen,korbgn,  &
     idum4(14),nmodes,idum6(2),lstzwd
 COMMON /zzzzzz/ z(1)
 COMMON /condas/ idum5(4),forpi2
 COMMON /mpy3tl/ itrlra(7),itrlrb(7),itrlre(7),itrlrc(7),jscr(3),  &
     lkore,icode,prec,dummy(13)
 COMMON /packx / typin,typout,irow,nrow,incr
 EQUIVALENCE     (lamamr,infile(2)),(rz(1),z(1)),(dz(1),z(1)),  &
     (hk,iscr(2)),(kmw2,iscr(5)),(k,iscr(3))
 DATA    modnam/ 4HMRED,4H2N  /
 
!                    2
!     FORM KMW2 = M W  MATRIX
!                  I I
 
 IF (dry == -2) RETURN
 IF (korbgn+7+nmodes >= korlen) GO TO 75
 ifile = lamamr
 CALL gopen  (lamamr,z(gbuf1),0)
 CALL fwdrec (*70,lamamr)
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
   CALL READ (*60,*70,lamamr,z(korbgn),7,0,nwdsrd)
   DO  j = 1,nmodes
     rz(korbgn+7+j-1) = 0.0
     IF (j == i) rz(korbgn+7+j-1) = forpi2*rz(korbgn+5)*(rz(korbgn+4)**2)
   END DO
   CALL pack (z(korbgn+7),kmw2,itrlr1)
 END DO
 CALL CLOSE (lamamr,1)
 CALL CLOSE (kmw2,1)
 CALL wrttrl (itrlr1)
 
!     FORM K MATRIX
 
!                      T
!                       **      **
!        ** **   **  ** * .    0 * **  **
!        *   *   *    * *  .     * *    *                  2
!        * K * = * HK * *   K    * * HK *     WHERE K = M W
!        *   *   *    * *    .   * *    *                I I
!        ** **   **  ** * 0   .  * **  **
!                       **      **
 
 itrlr2(1) = hk
 CALL rdtrl (itrlr2)
 DO  i = 1,7
   itrlra(i) = itrlr2(i)
   itrlrb(i) = itrlr1(i)
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
 CALL makmcb (itrlrc,k,itrlr2(3),iform,itype)
 jscr(1) = iscr(8)
 jscr(2) = iscr(6)
 jscr(3) = iscr(2)
 icode   = 0
 prec    = 0
 dblkor  = (korbgn/2) + 1
 lkore   = lstzwd - (2*dblkor - 1)
 CALL sofcls
 CALL mpy3dr (dz(dblkor))
 CALL wrttrl (itrlrc)
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 60 imsg = -2
 GO TO 80
 70 imsg = -3
 GO TO 80
 75 imsg = -8
 ifile = 0
 80 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
END SUBROUTINE mred2n
