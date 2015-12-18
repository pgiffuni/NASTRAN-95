SUBROUTINE xrecps (inew,iold)
     
!     ******************************************************************
!     * ATTENTION CDC 6600 SET-UPS ** THESE ENTRY POINTS MAY BE        *
!     * SEPARATED EACH ENTRY MAY BE MADE A SUBROUTINE (EXCEPT /CRDFLG/ *
!     * AND /INTEXT/ WHICH USE COMMON CODE)  DUPE THE SPECIFICATION    *
!     * STMTS FOR EACH SUB                                             *
!     ******************************************************************
 
 
 INTEGER, INTENT(IN)                      :: inew
 INTEGER, INTENT(OUT)                     :: iold
 INTEGER :: kbmsk1(8),sft(4),nrecps(2),con(38),mk(4),c10c(7),extab(37),&
            outtap
 
!     ENTRY XFADJ (BF,SD,KK)
!     * XFADJ ADJUSTS 4 CHARACTER FIELDS, LEFT OR RIGHT, 2 OR 4 FIELDS
!       AT A TIME - IF FIELDS CONTAIN ONLY INTEGERS 0 THRU 9, SHIFT IS
!       RIGHT, OTHERWISE SHIFT IS LEFT  / BF= ADDR OF LEFT MOST FIELD /
!       SD= 0 SINGLE (2 FIELDS), 1 DOUBLE (4 FIELDS).  THIS ROUTINE
!       DETERMINES ONLY TYPE OF SHIFT NEEDED, SHIFTING IS DONE BY XFADJ1
!       KK IS RETURNED EQUAL TO 0 FOR INTEGER, 1 FOR NON-INTEGER
 
 INTEGER :: bf(1)
 
!     ENTRY XBCDBI (BA)
!     * XBCDBI CONVERTS 2, 4 CHARACTER BCD INTEGER FIELDS (RIGHT
!       ADJUSTED IN THE LEFT MOST 4 CHAR) INTO A SINGLE FORTRAN BINARY
!       INTEGER (RIGHT ADJUSTED IN THE WORD IN THE RIGHT FIELD)
!       BA= ADDR OF LEFT FIELD
 
 INTEGER :: ba(2)
 
!     ENTRY XPRETY (BFF)
!     * ROUTINE PRETTIES UP SORT OUTPUT BY LEFT ADJUSTING ALL FIELDS
 
 INTEGER :: bff(2)
 
!     ENTRY CRDFLG (CARD)
!     * ROUTINE SETS CARD TYPE FLAGS IN RESTART TABLES
!       CONVERTS TO EXTERNAL CODE FIRST
!       IF CARD TYPE IS PARAM, SET FLAG FOR PARAM NAME (FIELD 2)
 
 INTEGER :: card(4)
 
!     ENTRY EXTINT (EXTWRD)
!     * ROUTINE CONVERTS FROM EXTERNAL MACHINE DEPENDENT CHARACTER CODES
!       TO AN INTERNAL MACHINE INDEPENDENT INTEGER
 
 INTEGER :: extwrd(1)
 
!     ENTRY INTEXT (INTWRD)
!     * ROUTINE CONVERTS FROM INTERNAL MACHINE INDEPENDENT INTEGERS TO
!       AN EXTERNAL MACHINE DEPENDENT CHARACTER CODE
 
 INTEGER :: intwrd(2)
 
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: dec
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /system/ b,outtap,d1(6),nlpp,d2(2),lcnt,d3(26), nbpc,nbpw,ncpw
 COMMON /xsrtcm/ bimsk1(6),bimsk2(5),bimsk3(4),bimsk4(4),bimsk5(2),  &
     bimsk6,bkmsk1(8),bkmsk2,shifts(4),  &
     icon1,icon2,star,plus,dollar,starl,slash,sftm, mask,BLANK,mka,is,mbit4
 COMMON /two   / itwo(32)
 COMMON /ifpx0 / lbd,lcc,ibits(1)
 COMMON /ifpx1 / num,icards(2)
 EQUIVALENCE     (sft(1),shifts(1)),(mk(1),bimsk3(1)),  &
     (sft1,shifts(2)),(extab(1),con(1))
 DATA itape4/304/,nrecps/4HXREC,4HPS  /
 DATA con/4H    ,4H   0,4H   1,4H   2,4H   3,4H   4,4H   5,4H   6,  &
     4H   7,4H   8,4H   9,4H   a,4H   b,4H   c,4H   d,4H   e,4H   f,  &
     4H   g,4H   h,4H   i,4H   j,4H   k,4H   l,4H   m,4H   n,4H   o,  &
     4H   p,4H   q,4H   r,4H   s,4H   t,4H   u,4H   v,4H   w,4H   x,  &
     4H   y,4H   z,4H              /
 DATA c10c/10,100,1000,10000,100000,1000000,10000000/
 DATA par1,par2/4HPARA,4HM       /
 DATA kpret1,kpret2/4H.   ,4H0.0 /
 
 DATA kbmsk1 / 4H0000, 4H000$, 4H00$$, 4H0$$$,  &
     4H$$$ , 4H$$  , 4H$   , 4H            /
 DATA istr   , istrl , ipls  , idollr, islsh , izero /  &
     4H   * , 4H*   , 4H+   , 4H$   , 4H/   , 4H0   /
 
 
!     THE ARRAYS IN /XSRTCM/ WILL BE SET BY INITO AS FOLLOWS
 
!                           VAX
!                    CDC    IBM   UNIVAC
!        SHIFTS(1) =  0      0      0
!        SHIFTS(2) =  6      8      9
!        SHIFTS(3) = 12     16     18
!        SHIFTS(4) = 18     24     27
!        SFTM      = 36      0      0
 
!                      ----------- BYTE --------------
!                      1ST   2ND   3RD   4TH   5TH,...
!        BIMSK1(1) = / 777 / 777 / 777 / 000 / 00..        CDC USES /77/
!        BIMSK1(2) = / 777 / 777 / 000 / 000 / 00..        INSTEAD OF
!        BIMSK1(3) = / 777 / 000 / 000 / 000 / 00..        /777/ IN A
!        BIMSK1(4) = / 000 / 000 / 000 / 777 / 00..        BYTE
!        BIMSK1(5) = / 000 / 000 / 777 / 777 / 00..
!        BIMSK1(6) = / 000 / 777 / 777 / 777 / 00..
 
!        BIMSK2(1) = / 777 / 777 / 777 / 777 / 77.. (FOR CDC ONLY)
!                  = / 377 / 777 / 777 / 777 / 00.. (FOR IBM,VAX,UNIVAC)
!        BIMSK2(2) = / 777 / 777 / 777 / 000 / 77..
!        BIMSK2(3) = / 777 / 777 / 000 / 000 / 77..
!        BIMSK2(4) = / 777 / 000 / 000 / 000 / 77..
!        BIMSK2(5) = / 000 / 000 / 000 / 000 / 77..
 
!        BIMSK3(1) = / 777 / 000 / 000 / 000 / 00..
!        BIMSK3(2) = / 000 / 777 / 000 / 000 / 00..
!        BIMSK3(3) = / 000 / 000 / 777 / 000 / 00..
!        BIMSK3(4) = / 000 / 000 / 000 / 777 / 00..
 
!        BIMSK4(1) = / 000 / 777 / 777 / 777 / 77..
!        BIMSK4(2) = / 777 / 000 / 777 / 777 / 77..
!        BIMSK4(3) = / 777 / 777 / 000 / 777 / 77..
!        BIMSK4(4) = / 777 / 777 / 777 / 000 / 77..
 
!        BIMSK5(1) = / 377 / 777 / 777 / 777 / 00..
!        BIMSK5(2) = / 377 / 777 / 777 / 000 / 00..
!        BIMSK6    = / 000 / 000 / 000 / 000 / 77..
 
!        IS        = / 400 / 000 / 000 / 000 / 77..
!        MKA       = / 000 / 000 / 000 / 777 / 77..
!        MASK      = 4TH OR 10TH BYTE IS /777/, ZERO FILLED
!        BLANK     = 4TH OR 10TH BYTE IS BLANK, ZERO FILLED
 
!     ARRAY BKMSK1 IS SAME AS KBMSK1 EXCEPT THAT THE DOLLARS ARE
!     REPLACED BY BINARY ZEROS
!     SIMILARY, THE BLANKS IN ISTR,ISTRL,IPLS,IDOLLR,ISLSH, AND ARRAY
!     CON ARE ALSO REPLACED BY BINARY ZEROS.
!     ICON1 AND ICON2 ARE LEFT ADJUSTED CON(1) AND CON(2), ZERO FILLED.
 
!     THIS ROUTINE POSITIONS ITAPE4 TO THE PROPER CONTINUATION RECORD
 
 dec  = mach == 5 .OR. mach == 6 .OR. mach == 21
 IF (inew /= 1) GO TO 10
 CALL REWIND (itape4)
 iold = 2
 RETURN
 
 10 idif = inew - iold
 IF (idif < 0) THEN
   GO TO    50
 ELSE IF (idif == 0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 20 iold = inew + 1
 RETURN
 
 30 DO  i = 1,idif
   CALL fwdrec (*65,itape4)
 END DO
 GO TO 20
 50 idif = IABS(idif)
 DO  i = 1,idif
   CALL bckrec (itape4)
 END DO
 GO TO 20
 65 WRITE  (outtap,66) sfm
 66 FORMAT (a25,' 217, ILLEGAL EOF ON ITAPE4.')
 CALL mesage (-37,0,nrecps)
 RETURN
 
!     INITIALIZES BCD CONSTANTS FOR USE WITHIN SORT
 
 ENTRY initco
!     ============
 
!     INITIALIZE (CREATE) BINARY CHARACTER MASKS
 
 dec       = mach == 5 .OR. mach == 6 .OR. mach == 21
 shifts(1) = 0
 shifts(2) = nbpc
 shifts(3) = nbpc*2
 shifts(4) = nbpc*3
 mbits     = complf(0)
 sftm      = (ncpw-4)*nbpc
 mbit4     = lshift(mbits,sftm)
 bimsk1(1) = lshift(mbit4,nbpc)
 bimsk1(2) = lshift(bimsk1(1),nbpc)
 bimsk1(3) = lshift(bimsk1(2),nbpc)
 bimsk1(4) = rshift(bimsk1(3),nbpc*3)
 bimsk1(5) = rshift(bimsk1(2),nbpc*2)
 bimsk1(6) = rshift(bimsk1(1),nbpc)
 bimsk2(1) = mbits
 bimsk2(2) = complf(bimsk1(4))
 bimsk2(3) = complf(bimsk1(5))
 bimsk2(4) = complf(bimsk1(6))
 bimsk2(5) = rshift(mbits,nbpc*4)
 bimsk3(4) = bimsk1(4)
 bimsk3(3) = lshift(bimsk3(4),nbpc)
 bimsk3(2) = lshift(bimsk3(3),nbpc)
 bimsk3(1) = bimsk1(3)
 bimsk4(1) = complf(bimsk3(1))
 bimsk4(2) = complf(bimsk3(2))
 bimsk4(3) = complf(bimsk3(3))
 bimsk4(4) = complf(bimsk3(4))
 bimsk5(1) = rshift(bimsk2(1),1)
 bimsk5(2) = rshift(lshift(bimsk2(2),1),1)
 bimsk6    = bimsk2(5)
 IF (mach == 2 .OR. dec) bimsk2(1) = bimsk5(1)
 
!     NEXT CARD FOR UNIVAC ASCII VERSION ONLY (NOT FORTRAN 5)
 
 IF (mach == 3) bimsk2(1) = bimsk5(1)
 mask  = rshift(bimsk3(4),sftm)
 BLANK = rshift(kbmsk1(8),(3*nbpc+sftm))
 is    = complf(bimsk5(1))
 mka   = orf(bimsk3(4),bimsk6)
 
!     INITIALIZE THE BCD BLANK DATA
 
 IF (dec) GO TO 92
 
!     IBM, CDC, UNIVAC
 
 bkmsk1(1) = kbmsk1(1)
 bkmsk1(2) = andf(kbmsk1(2),bimsk2(2))
 bkmsk1(3) = andf(kbmsk1(3),bimsk2(3))
 bkmsk1(4) = andf(kbmsk1(4),bimsk2(4))
 bkmsk1(5) = andf(kbmsk1(5),orf(bimsk1(4),bimsk6))
 bkmsk1(6) = andf(kbmsk1(6),orf(bimsk1(5),bimsk6))
 bkmsk1(7) = andf(kbmsk1(7),orf(bimsk1(6),bimsk6))
 bkmsk1(8) = kbmsk1(8)
 bkmsk2    = andf(bkmsk1(1),bimsk6)
 star      = andf(istr  ,orf(bimsk1(4),bimsk6))
 plus      = andf(ipls  ,bimsk2(4))
 dollar    = andf(idollr,bimsk2(4))
 starl     = andf(istrl ,bimsk2(4))
 slash     = andf(islsh ,bimsk2(4))
 DO  i = 1,38
   con(i) = andf(con(i),bimsk3(4))
 END DO
 icon1  = lshift(con(1),sft(4)-1)
 icon2  = lshift(con(2),sft(4)-1)
 RETURN
 
!     VAX
 
 92 bkmsk2    = 0
 bkmsk1(1) = kbmsk1(1)
 bkmsk1(2) = khrfn3(bkmsk2,kbmsk1(2),-1,1)
 bkmsk1(3) = khrfn3(bkmsk2,kbmsk1(3),-2,1)
 bkmsk1(4) = khrfn3(bkmsk2,kbmsk1(4),-3,1)
 bkmsk1(5) = khrfn3(bkmsk2,kbmsk1(5),-3,0)
 bkmsk1(6) = khrfn3(bkmsk2,kbmsk1(6),-2,0)
 bkmsk1(7) = khrfn3(bkmsk2,kbmsk1(7),-1,0)
 bkmsk1(8) = kbmsk1(8)
 star      = khrfn1(bkmsk2,4,istr  ,4)
 plus      = khrfn1(bkmsk2,1,ipls  ,1)
 dollar    = khrfn1(bkmsk2,1,idollr,1)
 starl     = khrfn1(bkmsk2,1,istrl ,1)
 slash     = khrfn1(bkmsk2,1,islsh ,1)
 DO  i = 1,38
   con(i) = khrfn1(bkmsk2,4,con(i),4)
 END DO
 icon1  = rshift(khrfn1(bkmsk2,1,con(1),4),1)
 icon2  = rshift(khrfn1(bkmsk2,1,con(2),4),1)
 RETURN
 
 
 ENTRY xfadj (bf,sd,kk)
!     ======================
 
!     DATA SFT /0,6,12,18/
!     DATA MK  /O770000000000,O007700000000,O000077000000,O000000770000/
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 ii  = 2
 IF (sd == 1) ii = 4
 DO  i = 1,ii
   bfi = bf(i)
   DO  j = 1,4
     ji  = 5 - j
     IF (.NOT.dec) test = rshift(andf(bfi,mk(j)),sft(ji))
     IF (     dec) test = khrfn1(bkmsk2,4,bfi,j)
     DO  k = 1,11
       IF (test == con(k)) GO TO 200
     END DO
     
!     CHARACTER NON-INTEGER
     
     CALL xfadj1 (bf,lshift,sd)
     kk = 1
     RETURN
     
     200 IF (k == 1) CYCLE
     
!     CHARACTER INTEGER
     
     CALL xfadj1 (bf,rshift,sd)
     kk = 0
     RETURN
     
   END DO
 END DO
 
!     ALL FIELDS BLANK
 
 kk = 0
 RETURN
 
 
 ENTRY xbcdbi (ba)
!     =================
 
!     DATA SFT1/6/,SFTM/12/,MASK/O77/,BLANK/O60/
 
!     IF MACHINE IS VAX-11/780, ORDER OF CHARACTERS IN A WORD IS REVERSE
!     OF THAT ON OTHER MACHINES.  THE CHARACTER ORDER MUST THEREFORE BE
!     REVERSED BEFORE DECODING TO AN INTEGER VALUE.
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 IF (.NOT.dec) GO TO 430
 DO  iba = 1,2
   itemp = 0
   DO  ivax = 1,4
     jtemp = rshift(ba(iba),8*(ivax-1))
     jtemp = andf(mask,jtemp)
     jtemp = lshift(jtemp,8*(4-ivax))
     itemp = orf(itemp,jtemp)
   END DO
   ba(iba) = itemp
 END DO
 
 430 CONTINUE
 ba(1) = rshift(ba(1),sftm)
 ba(2) = rshift(ba(2),sftm)
 ivar  = andf(ba(2),mask)
 IF (ivar /= BLANK) GO TO 490
 ba(2) = 0
 RETURN
 
 490 IF (mach == 4) ivar = ivar - 27
 ivar  = andf(ivar,15)
 DO  i = 1,3
   ba(2) = rshift(ba(2),sft1)
   ICHAR = andf(ba(2),mask)
   IF (mach == 4) ICHAR = ICHAR - 27
   ivar  = ivar + c10c(i)*andf(15,ICHAR)
 END DO
 ICHAR = andf(ba(1),mask)
 IF (mach == 4) ICHAR = ICHAR - 27
 ivar  = ivar + c10c(4)*andf(15,ICHAR)
 DO  i = 5,7
   ba(1) = rshift(ba(1),sft1)
   ICHAR = andf(ba(1),mask)
   IF (mach == 4) ICHAR = ICHAR - 27
   ivar  = ivar + c10c(i)*andf(15,ICHAR)
 END DO
 ba(2) = ivar
 RETURN
 
 
 ENTRY xprety (bff)
!     ==================
 
!     DATA  MKA/O000000777777/, STAR/4H000*/
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 IF (.NOT.dec) itst = andf(mka,bff(2))
 IF (     dec) itst = khrfn1(bkmsk2,4,bff(2),4)
 IF (itst == star) GO TO 610
 DO  i = 3,17,2
   IF (bff(i) == bkmsk1(8) .AND. bff(i+1) == bkmsk1(8)) CYCLE
   CALL xfadj1 (bff(i),lshift,0)
   IF (bff(i) == kpret1) bff(i) = kpret2
   IF (bff(i) /= bkmsk1(8)) CYCLE
   IF (.NOT.dec) bff(i) = orf(rshift(bff(i),sft(2)),bkmsk1(4))
   IF (     dec) bff(i) = khrfn3(izero,bff(i),1,0)
 END DO
 RETURN
 
 610 DO  i = 3,15,4
   IF (bff(i) == bkmsk1(8) .AND. bff(i+1) == bkmsk1(8) .AND.  &
       bff(i+2) == bkmsk1(8) .AND. bff(i+3) == bkmsk1(8)) CYCLE
   CALL xfadj1 (bff(i),lshift,1)
   IF (bff(i) == kpret1) bff(i) = kpret2
   IF (bff(i) /= bkmsk1(8)) CYCLE
   IF (.NOT.dec) bff(i) = orf(rshift(bff(i),sft(2)),bkmsk1(4))
   IF (     dec) bff(i) = khrfn3(izero,bff(i),1,0)
 END DO
 RETURN
 
 
 ENTRY crdflg (card)
!     ===================
 
 dec    = mach == 5 .OR. mach == 6 .OR. mach == 21
 inwrdi = card(1)
 kard2  = card(2)
 kbrn   = -1
 ASSIGN 640 TO iret
 GO TO 770
 640 IF (.NOT.dec) kard2 = orf(andf(bimsk1(1),kard2),bkmsk1(5))
 IF (     dec) kard2 = khrfn1(kard2,4,bkmsk1(8),4)
 IF (kard1 /= par1 .OR. kard2 /= par2) GO TO 645
 kard1 = card(3)
 kard2 = card(4)
 645 lmt   = num* 2
 DO  i = 1,lmt,2
   IF (kard1 == icards(i) .AND. kard2 == icards(i+1)) GO TO 660
 END DO
 RETURN
 
 660 j = i/2
 icycl = (j/31) + 1
 ipos  = MOD(j,31) + 2
 ibits(icycl) = orf(ibits(icycl),itwo(ipos))
 RETURN
 
 
 ENTRY extint (extwrd)
!     =====================
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 DO  i = 1,2
   exwrdi = extwrd(i)
   DO  j = 1,4
     ji = 5 - j
     sftji = sft(ji)
     IF (.NOT.dec) test = rshift(andf(exwrdi,mk(j)),sftji)
     IF (     dec) test = khrfn1(bkmsk2,4,exwrdi,j)
     DO  k = 1,37
       IF (test == extab(k)) GO TO 720
     END DO
     k = 1
     EXIT
     720 IF (.NOT.dec)  &
         exwrdi = orf(andf(exwrdi,bimsk4(j)),lshift(k,sftji+sftm))
     IF (dec) exwrdi = khrfn1(exwrdi,j,k,-1)
     IF (k == 1) EXIT
   END DO
   740 extwrd(i) = exwrdi
   IF (k == 1) RETURN
 END DO
 RETURN
 
 
 ENTRY intext (intwrd)
!     =====================
 
 dec    = mach == 5 .OR. mach == 6 .OR. mach == 21
 ASSIGN 800 TO iret
 inwrdi = intwrd(1)
 kbrn   = 0
 770 DO  j = 1,4
   ji     = 5 - j
   sftji  = sft(ji)
   IF (.NOT.dec) test = rshift(andf(inwrdi,mk(j)),sftji+sftm)
   IF (     dec) test = khrfn1(bkmsk2,-1,inwrdi,j)
   IF (test > 37) EXIT
   IF (.NOT.dec)  &
       inwrdi = orf(andf(inwrdi,bimsk4(j)),lshift(extab(test),sftji))
   IF (dec) inwrdi = khrfn1(inwrdi,j,extab(test),4)
   IF (test == 1) EXIT
 END DO
 781 IF (kbrn < 0) THEN
   GO TO   782
 ELSE IF (kbrn == 0) THEN
   GO TO   784
 ELSE
   GO TO   786
 END IF
 782 kard1  = inwrdi
 inwrdi = card(2)
 kbrn   = +2
 GO TO 810
 784 intwrd(1) = inwrdi
 inwrdi = intwrd(2)
 kbrn = +1
 GO TO 810
 786 IF (kbrn == 1) GO TO 788
 kard2 = inwrdi
 GO TO 790
 788 intwrd(2) = inwrdi
 790 GO TO iret, (800,640)
 800 RETURN
 
 810 IF (test == 1 .OR. test > 37) GO TO 790
 GO TO 770
 
END SUBROUTINE xrecps
