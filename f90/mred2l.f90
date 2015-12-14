SUBROUTINE mred2l (nuf,n2,nus,ufbits)
     
!     THIS SUBROUTINE PREFORMS PRELIMINARY CALCULATIONS AND MERGES OF
!     THE HK MATRIX FOR THE MRED2 MODULE.
 
 
 INTEGER, INTENT(IN)                      :: nuf
 INTEGER, INTENT(IN)                      :: n2
 INTEGER, INTENT(IN)                      :: nus
 INTEGER, INTENT(IN OUT)                  :: ufbits
 EXTERNAL        andf
 INTEGER :: dry,gbuf1,sbuf1,sbuf2,sbuf3,otfile,t,signab,signc,  &
     prec,scr,typin,typout,andf,rows, phiss1, phiss,dblkor,sglkor
 DOUBLE PRECISION :: dz
 DIMENSION       itrlr1(7),itrlr2(7),modnam(2),isub(4),rz(1),dz(1)
 COMMON /BLANK / idum1,dry,idum2,gbuf1,idum3(2),sbuf1,sbuf2,sbuf3,  &
     infile(12),otfile(6),iscr(10),korlen,korbgn, idum4(14),modpts,idum5(2),lstzwd
 COMMON /zzzzzz/ z(1)
 COMMON /mpyadx/ itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nz,t,  &
     signab,signc,prec,scr
 COMMON /packx / typin,typout,irow,nrow,incr
 EQUIVALENCE     (gs,iscr(7)),(phiss1,iscr(8)),(phiss2,iscr(9)),  &
     (ident,iscr(5)),(phiss,iscr(6)),(phigs,iscr(2)),  &
     (phis12,iscr(2)),(phi12i,iscr(8)),(rprtn,iscr(5)),  &
     (cprtn,iscr(10)),(rz(1),z(1)),(dz(1),z(1))
 DATA    modnam/ 4HMRED,4H2L  /
 
!                  -1
!     COMPUTE PHISS1
 
 IF (dry == -2) RETURN
 CALL sofcls
 ifile = phiss1
 itrlr1(1) = phiss1
 CALL rdtrl (itrlr1)
 CALL gopen (phiss1,z(gbuf1),0)
 kolumn = itrlr1(2)
 rows  = itrlr1(3)
 itest = kolumn * rows
 IF ((korbgn+itest+(3*kolumn)) >= korlen) GO TO 190
 kore = 0
 dblkor = (korbgn/2) + 1
 sglkor = (2*dblkor) - 1
 IF (itrlr1(5) == 2) GO TO 15
 DO  i = 1,kolumn
   CALL READ (*170,*180,phiss1,z(sglkor+kore),rows,0,nwdsrd)
   kore  = kore + rows
 END DO
 icore = ((sglkor+itest)/2) + 1
 GO TO 19
 15 DO  i = 1,kolumn
   CALL READ (*170,*180,phiss1,dz(dblkor+kore),rows,0,nwdsrd)
   kore  = kore + rows
 END DO
 icore = dblkor + itest
 19 CALL CLOSE (phiss1,1)
 invert = 0
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (rows,dz(dblkor),kolumn,b,invert,determ,ising, dz(icore))
 IF (ising == 2) GO TO 192
 kore  = 0
 incr  = 1
 typin = 1
 typout= 1
 irow  = 1
 nrow  = kolumn
 CALL makmcb (itrlr2,phiss1,nrow,itrlr1(4),itrlr1(5))
 CALL gopen (phiss1,z(gbuf1),1)
 DO  i = 1,rows
   IF (itrlr1(5) == 2) GO TO 25
   CALL pack (rz(sglkor+kore),phiss1,itrlr2)
   GO TO 29
   25 CALL pack (dz(dblkor+kore),phiss1,itrlr2)
   29 kore = kore + rows
 END DO
 CALL CLOSE (phiss1,1)
 
!     COMPUTE PHIGS
 
!                               -1
!        **     **    **      **  **     ** **  **
!        *       *    *        *  *       * *    *
!        * PHIGS * = -* PHISS1 *  * PHISS * * GS *
!        *       *    *        *  *       * *    *
!        **     **    **      **  **     ** **  **
 
 itrlr1(1) = phiss1
 itrlr2(1) = phiss
 CALL rdtrl (itrlr1)
 CALL rdtrl (itrlr2)
 icol = itrlr2(3)
 DO  i = 1,7
   itrlra(i) = itrlr1(i)
   itrlrb(i) = itrlr2(i)
   itrlrc(i) = 0
 END DO
 CALL makmcb (itrlrd,phissi,itrlr2(3),itrlr2(4),itrlr2(5))
 t = 0
 signab = -1
 signc = 1
 prec = 0
 scr = iscr(10)
 nz = lstzwd - ((2*dblkor)-1)
 CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 itrlr1(1) = gs
 CALL rdtrl (itrlr1)
 DO  i = 1,7
   itrlra(i) = itrlrd(i)
   itrlrb(i) = itrlr1(i)
 END DO
 CALL makmcb (itrlrd,phigs,itrlr1(3),itrlr1(4),itrlr1(5))
 signab = 1
 CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 
!     FORM HKPG MATRIX
 
!                   **               **
!                   *       .         *
!                   *       .      -1 *
!        **    **   * PHIGS . PHISS   *
!        *      *   *       .         *
!        * HKPG * = *.................*
!        *      *   *       .         *
!        **    **   *   0   .    0    *
!                   *       .         *
!                   **               **
 
 nrow = nuf + n2
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > n2) rz(korbgn+i-1) = 1.0
 END DO
 typin  = 1
 typout = 1
 irow   = 1
 incr   = 1
 iform  = 7
 CALL makmcb (itrlr1,cprtn,nrow,iform,typin)
 CALL gopen (cprtn,z(gbuf1),1)
 CALL pack (z(korbgn),cprtn,itrlr1)
 CALL CLOSE (cprtn,1)
 CALL wrttrl (itrlr1)
 nrow  = nus + nuf
 ifile = usetmr
 CALL gopen (usetmr,z(gbuf1),0)
 DO  i = 1,nrow
   CALL READ (*170,*180,usetmr,z(korbgn),1,0,nwdsrd)
   rz(korbgn+i) = 0.0
   IF (andf(z(korbgn),ufbits) /= 0) rz(korbgn+i) = 1.0
 END DO
 CALL CLOSE (usetmr,1)
 nrow = rows
 rows = nuf + n2
 CALL makmcb (itrlr2,rprtn,nrow,iform,typin)
 CALL gopen (rprtn,z(gbuf1),1)
 CALL pack (z(korbgn+1),rprtn,itrlr2)
 CALL CLOSE (rprtn,1)
 CALL wrttrl (itrlr2)
 isub(1) = nuf
 isub(2) = n2
 isub(3) = nus
 isub(4) = n2
 itype   = 2
 CALL gmmerg (hkpg,phigsh,0,phiss1,0,rprtn,cprtn,isub,itype, z(korbgn),korlen)
 
!     COMPUTE PHIS12
 
!                                -1
!        **      **    **      **  **      **
!        *        *    *        *  *        *
!        * PHIS12 * = -* PHISS1 *  * PHISS2 *
!        *        *    *        *  *        *
!        **      **    **      **  **      **
 
 itrlr1(1) = phiss1
 itrlr2(1) = phiss2
 CALL rdtrl (itrlr1)
 CALL rdtrl (itrlr2)
 modpts = itrlr1(3) + itrlr2(3)
 DO  i = 1,7
   itrlra(i) = itrlr1(i)
   itrlrb(i) = itrlr2(i)
 END DO
 CALL makmcb (itrlrd,phis12,itrlr2(3),itrlr2(4),itrlr2(5))
 signab = -1
 CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 
!     GENERATE IDENTITY MATRIX
 
 nrow = itrlrd(3)
 CALL makmcb (itrlr1,ident,nrow,itrlrd(4),itrlrd(5))
 CALL gopen (ident,z(gbuf1),1)
 DO  i = 1,nrow
   DO  j = 1,nrow
     rz(korbgn+j-1) = 0.0
     IF (j == i) rz(korbgn+j-1) = 1.0
   END DO
   CALL pack (z(korbgn),ident,itrlr1)
 END DO
 CALL CLOSE (ident,1)
 CALL wrttrl (itrlr1)
 
!     GENERATE PHI12I MATRIX
 
!                     **      **
!                     *        *
!        **      **   * PHIS12 *
!        *        *   *        *
!        * PHI12I * = *........*
!        *        *   *        *
!        **      **   *   I    *
!                     *        *
!                     **      **
 
 itrlr1(1) = phis12
 CALL rdtrl (itrlr1)
 isub(3) = itrlr1(3)
 isub(4) = nrow
 nrow = itrlr1(2) + nrow
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > itrlr1(2)) rz(korbgn+i-1) = 1.0
 END DO
 incr = 1
 CALL makmcb (itrlr2,cprtn,nrow,iform,typin)
 CALL gopen (cprtn,z(gbuf1),1)
 CALL pack (z(korbgn),rprtn,itrlr2)
 CALL CLOSE (cprtn,1)
 CALL wrttrl (itrlr2)
 CALL gmmerg (phi12i,phis12,ident,0,0,0,cprtn,isub,itype, z(korbgn),korlen)
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 170 imsg = -2
 GO TO 200
 180 imsg = -3
 GO TO 200
 190 imsg = -8
 GO TO 194
 192 imsg = -37
 194 ifile = 0
 200 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
END SUBROUTINE mred2l
