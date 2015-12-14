SUBROUTINE mred2i (kode,nuf,n2)
     
!     THIS SUBROUTINE COMPUTES THE GS MATRIX FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN OUT)                  :: kode
 INTEGER, INTENT(IN OUT)                  :: nuf
 INTEGER, INTENT(IN OUT)                  :: n2
 INTEGER :: dry,gbuf1,otfile,z,typin,typout,typinu,  &
     qsmrow,qsmcol,gs,qsm,dblkor,sglkor,qsmtyp
 DOUBLE PRECISION :: dz
 DIMENSION        modnam(2),itrlr1(7),rz(1),dz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  idum1,dry,idum2,gbuf1,idum3(5),infile(12),  &
     otfile(6),iscr(10),korlen,korbgn,idum6(14), nmodes,modlen
 COMMON /zzzzzz/  z(1)
 COMMON /system/  idum4,iprntr
 COMMON /condas/  idum5(4),forpi2
 COMMON /packx /  typin,typout,irow,nrow,incr
 COMMON /unpakx/  typinu,irowu,nrowu,incru
 EQUIVALENCE      (lamamr,infile(2)),(qsm,infile(12)),  &
     (gs,iscr(7)),(rz(1),z(1)),(dz(1),z(1))
 DATA    modnam/  4HMRED,4H2I  /, idiag / 3 /
 
!     TEST OPERATION MODE
 
 IF (dry == -2) RETURN
 
!     FORM GS MATRIX
 
!                 **     **
!                 *       *        T
!        **  **   * .   0 * **   **
!        *    *   *  .    * *     *                 2
!        * GS * =-*  1/K  * * QSM *    WHERE K = M W
!        *    *   *    .  * *     *               I I
!        **  **   * 0   . * **   **
!                 *       *
!                 **     **
 
 itrlr1(1) = qsm
 CALL rdtrl (itrlr1)
 IF (itrlr1(1) < 0) GO TO 100
 qsmrow = itrlr1(2)
 qsmcol = itrlr1(3)
 
!                        2
!     FORM K = 1.0 / (M W )
!                      I I
 
 IF (korbgn+7+qsmrow >= korlen) GO TO 80
 ifile = lamamr
 CALL gopen (lamamr,z(gbuf1),0)
 CALL fwdrec (*70,lamamr)
 nmodes = 0
 10 CALL READ (*60,*15,lamamr,z(korbgn),7,0,nwdsrd)
 rz(korbgn+7+nmodes) = 1.0/(forpi2*rz(korbgn+5)*(rz(korbgn+4)**2))
 nmodes = nmodes + 1
 IF (korbgn+7+nmodes >= korlen) GO TO 80
 GO TO 10
 15 CALL CLOSE (lamamr,1)
 IF (nmodes /= qsmrow) GO TO 110
 modlen = nmodes
 
!     READ QSM INTO CORE
 
 kore   = korbgn
 korbgn = korbgn + 7 + itrlr1(2)
 IF (korbgn+qsmrow*(qsmcol+1) >= korlen) GO TO 80
 typinu = itrlr1(5)
 irowu  = 1
 nrowu  = itrlr1(3)
 incru  = 1
 qsmtyp = itrlr1(5)
 dblkor = korbgn/2 + 1
 sglkor = 2*dblkor - 1
 CALL gopen (qsm,z(gbuf1),0)
 IF (qsmtyp == 2) GO TO 26
 locqsm = sglkor
 DO  i = 1,qsmrow
   CALL unpack (*20,qsm,rz(sglkor))
   GO TO 25
   20 DO  j = 1,qsmcol
     rz(sglkor+j-1) = 0.0E0
   END DO
   25 sglkor = sglkor + itrlr1(3)
 END DO
 korbgn = sglkor
 GO TO 30
 26 locqsm = dblkor
 DO  i = 1,qsmrow
   CALL unpack (*27,qsm,dz(dblkor))
   GO TO 29
   27 DO  j = 1,qsmcol
     dz(dblkor+j-1) = 0.0D0
   END DO
   29 dblkor = dblkor + itrlr1(3)
 END DO
 korbgn = dblkor
 30 CALL CLOSE (qsm,1)
 
!     FORM GS MATRIX
 
 typin  = itrlr1(5)
 typout = itrlr1(5)
 irow   = 1
 nrow   = qsmrow
 incr   = 1
 CALL makmcb (itrlr1,gs,qsmrow,idiag,typin)
 dblkor = korbgn/2 + 1
 sglkor = 2*dblkor - 1
 CALL gopen (gs,z(gbuf1),1)
 DO  i = 1,qsmcol
   DO  j = 1,qsmrow
     k = 3*(j-1)
     IF (qsmtyp == 2) GO TO 32
     rz(sglkor+j-1) = rz(kore+7+j-1)*rz(locqsm+k)
     CYCLE
     32 CONTINUE
     dz(dblkor+j-1) = rz(kore+7+j-1)*dz(locqsm+k)
   END DO
   CALL pack (dz(dblkor),gs,itrlr1)
 END DO
 korbgn = kore
 CALL CLOSE (gs,1)
 CALL wrttrl (itrlr1)
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 60 imsg = -2
 GO TO 90
 70 imsg = -3
 GO TO 90
 80 imsg = -8
 ifile = 0
 90 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 100 WRITE (iprntr,900) ufm
 GO TO 120
 110 WRITE (iprntr,901) ufm,qsmrow,qsmcol,nmodes
 120 dry = -2
 RETURN
 
 900 FORMAT (a23,' 6638, IN MODULE MREDUCE WITH USERMODE=2, THE ',  &
     'CONSTRAINT FORCES MATRIX (QSM) CANNOT BE PURGED.')
 901 FORMAT (a23,' 6634, IN MODULE MREDUCE WITH USERMODE=2, THE ',  &
     'CONSTRAINT FORCES MATRIX (',i3,3H x ,i3,1H), /30X,  &
     'IS INCOMPATABLE WITH THE NUMBER OF MODES (',i3,2H).)
 
END SUBROUTINE mred2i
