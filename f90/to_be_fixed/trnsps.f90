SUBROUTINE trnsps (z,iz)
     
!     MATRIX TRANSPOSE ROUTINE REPLACING NASTRAN ORIGINAL TRNSP, WHICH
!     IS AT LEAST 2 TO 4 TIMES SLOWER (COMPARISON DONE ON VAX), AND
!     USING UP TO 8 SCRATCH FILES
 
!     WITH BOTH IN-CORE AND OUT-OF-CORE LOGICS
!     (USE TRANSP FOR IN-CORE MATRIX TRANSPOSE)
 
!     IF DGFLAG = -123457890 (SET BY DTRANP), AND INPUT IS A UPPER OR
!     LOWER TRIANGULAR MATRIX, THE DIAGONAL ELEMENTS ARE REPLACED BY
!     UNITY (1.0)
 
!     CALLER MUST SUPPLY A SCRATCH FILE ISCR, IF MATRIX TO BE TRANSPOSED
!     IS SQUARE, RECTANGULAR, LOWER, AND UPPER TRIAGULAR (FORM 1,2,4,5).
 
!     THIS ROUTINE SETS UP THE OUTPUT MATRIX TRAILER WORDS IN NAMEAT
!     (FILEAT) BUT IT DOES NOT CALL WRTTRL TO WRITE THEM OUT
 
!     WRITTEN BY G.CHAN/UNISYS  12/91
 
 
 REAL, INTENT(OUT)                        :: z(6)
 INTEGER, INTENT(IN OUT)                  :: iz(2)
 LOGICAL :: debug
 INTEGER :: sysbuf,base,FILE,dgflag,filea(7),fileat(7)
 DIMENSION  a(2),nam(2)
 DOUBLE PRECISION :: da
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  dgflag
 COMMON /trnspx/  namea, ncola, nrowa, iforma,itypa, ia(2),  &
     nameat,ncolat,nrowat,iforat,itypat,iat(2), lcore,nscr,iscr
 COMMON /system/  sysbuf,nout
 COMMON /packx /  iotyp,iotypa,ip,jp,incr
 COMMON /unpakx/  iotyp1,iu,ju,incr1
 COMMON /TYPE  /  rc(2),iwords(4)
 COMMON /names /  rd,rdrew,wrt,wrtrew,clsrew
 EQUIVALENCE      (filea(1),namea),(fileat(1),nameat),(a(1),da)
 DATA    nam   /  4HTRNS,4HPS    /, debug / .false. /
 
 CALL sswtch (19,i)
 IF (i == 1) debug = .true.
 last   = 1
 ntype  = iotypa
 IF (ntype == 3) ntype = 2
 ibuf1  = lcore - sysbuf
 ibuf   = ibuf1 - sysbuf
 nz     = ibuf  - 1
 imhere = 10
 IF (nz <= 0) GO TO 820
 nrec   = 0
 FILE   = namea
 IF (iforma > 2 .OR. ncola == 1) CALL OPEN (*800,namea,z(ibuf1),rdrew)
 DO  i = 2,7
   fileat(i) = filea(i)
 END DO
 IF (debug) WRITE (nout,20) fileat
 20 FORMAT (' TRNSPS/@5 BEFORE TRANSPOSE, TRAIL-AT =',7I8)
 SELECT CASE ( iforma )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 600
   CASE (    5)
     GO TO 600
   CASE (    6)
     GO TO 500
   CASE (    7)
     GO TO 730
   CASE (    8)
     GO TO 550
 END SELECT
 
!     SQUARE AND RECTANGULAR MATRICES
!     ===============================
 
 30 IF (ncola == 1) GO TO 580
 nrowat = ncola
 ncolat = 0
 iat(1) = 0
 iat(2) = 0
 ip   = 1
 jp   = nrowat
 incr = 1
 nwd  = iwords(itypa)
 nwd1 = nwd - 1
 nwds = ncola*nwd
 IF (nrec /= 0) GO TO 40
 irat = MIN0(MAX0((lcore/100000+4)*ncola/nrowa,3),10)
 iend = (ibuf1-1-nwds)/irat
 iend = MAX0(iend,5000)
 iend1= iend + 1
 CALL unpscr (filea,iscr,z,ibuf1,ibuf,iend,0,1)
 nrec = filea(4)/10
 40 FILE = iscr
 CALL OPEN (*800,iscr,z(ibuf1),rdrew)
 j    = filea(6) - iend*irat
 IF (j > 0) GO TO 200
 
!     ENTIRE FILEA (FROM ISCR FILE) FITS INTO CORE
 
 IF (debug) WRITE (nout,50) uim
 50 FORMAT (a29,', MATRIX TRANSPOSE WAS PORCESSED BY THE NEW TRNSP ',  &
     'IN-CORE METHOD')
 CALL fwdrec (*810,iscr)
 ll = nwds + 1
 DO  i = 1,nrec
   CALL READ (*810,*60,iscr,z(ll),iend1,1,k)
   imhere = 60
   GO TO 820
   ll = ll + k
 END DO
 CALL CLOSE (iscr,clsrew)
 
 FILE = nameat
 CALL OPEN (*800,nameat,z(ibuf1),wrtrew)
 CALL fname (nameat,a(1))
 CALL WRITE (nameat,a(1),2,1)
 
 DO  k = 1,nrowa
   DO   j = 1,nwds
     z(j) = 0.0
   END DO
   base = nwds + 2
   IF (nwd-2 < 0) THEN
     GO TO   130
   ELSE IF (nwd-2 == 0) THEN
     GO TO   110
   END IF
   80 DO  i = 1,ncola
     ii = iz(base-1)
     jj = iz(base  )
     IF (k < ii .OR. k > jj) GO TO 100
     kx = (k-ii)*nwd + base
     lx = (i- 1)*nwd
     DO  j = 1,nwd
       z(j+lx) = z(j+kx)
     END DO
     100 base = base + (jj-ii+1)*nwd + 2
   END DO
   GO TO 150
   110 DO  i = 1,ncola
     ii = iz(base-1)
     jj = iz(base  )
     IF (k < ii .OR. k > jj) GO TO 120
     kx = (k-ii)*2 + base
     lx = (i- 1)*2
     z(lx+1) = z(kx+1)
     z(lx+2) = z(kx+2)
     120 base = base + (jj-ii+2)*2
   END DO
   GO TO 150
   130 DO  i = 1,ncola
     ii = iz(base-1)
     jj = iz(base  )
     IF (k < ii .OR. k > jj) GO TO 140
     kx = k - ii + base
     z(i) = z(kx+1)
     140 base = base + jj - ii + 3
   END DO
   150 CALL pack (z(1),nameat,nameat)
 END DO
 GO TO 450
 
!     ENTIRE FILEA CAN NOT FIT INTO CORE
 
!     OPEN CORE ALLOCATION -             N1    N2              NZ
!                                        /     /  <-- IEND --> /
!     +----------------------------------+-----+---------------+---+---+
!      /          OPEN CORE               /     /                GINO
!     I1                                 I2    I3               BUFFERS
 
!      Z(I1)... Z(N1) FOR TRANSPOSED OUTPUT MATRIX NAMEAT
!     IZ(I2)...IZ(N2) IS A (3 x NREC) TABLE, (MIN, MAX, COLUMN COUNTER)
!               CONTROLLING DATA TRANSFER FROM SCRATCH FILE ISCR.
!      Z(I3)... Z(NZ) FOR INPUT MATRIX NAMEA COMING FROM ISCR
 
!     NOTE - THE RATIO OF (N1-I1)/(NZ-I3), WHICH IS IRAT, IS A FUNCTION
!            OF OPEN CORE SIZE, AND THE MATRIX COLUMN AND ROW SIZES.
!            IRAT IS LIMITED TO 10:1
!     NCPP = NO. OF COULMNS PER PASS, OF THE TRANSPOSE MATRIX NAMEAT
 
!     THE TERMS 'ROW' AND 'COLUMN' ARE LOOSELY DEFINED IN COMMENT LINES
 
 200 n2   = nz - iend
 i3   = n2 + 1
 n1   = n2 - 3*nrec
 i2   = n1 + 1
 ncpp = n1/nwds
 ncp7 = ncpp*7
 npas = (ncola+ncpp-1)/ncpp
 IF (.NOT.debug .AND. j > 3*nz) GO TO 230
 WRITE  (nout,210) uim,npas,j
 210 FORMAT (a29,', MATRIX TRANSPOSE WAS PROCESSED BY THE NEW TRNSP ',  &
     'OUT-OF-CORE METHOD WITH',i5,' NO. OF PASSES', /5X,  &
     '(FOR MAXIMUM EFFECIENCY, THE IN-CORE METHOD COULD BE ',  &
     'ACTIVATED WITH',i9,' ADDITIONAL OPEN CORE WORDS)')
 WRITE  (nout,220) n1,iend,irat,ncpp,npas,nrec
 220 FORMAT (/5X,'OPEN CORE -',i9,' WORDS USED FOR TRANSPOSE OUTPUT ',  &
     'MATRIX, AND',i8,' WORDS FOR INPUT MATRIX (',i2,'/1 RATIO)'  &
     ,      /5X,'NO. OF COLUMNS PER PASS =',i5,',  NO. OF PASSES =',i6,  &
     ',  INPUT MATRIX REWRITTEN IN',i4,' RECORDS')
 230 FILE = nameat
 CALL OPEN  (*800,nameat,z(ibuf),wrtrew)
 CALL fname (nameat,a(1))
 CALL WRITE (nameat,a(1),2,1)
 DO  mm = i2,n2,3
   iz(mm  ) = nrowa
   iz(mm+1) = 0
 END DO
 CALL tmtogo (t1)
 
!     OUTER KB-KE LOOP
 
!     MAP DATA INTO TRANSPOSE OUTPUT MATRIX SPACE, Z(I1)...Z(N1), BY
!     PASSES. EACH PASS RANGES FROM KB THRU KE COLUMNS
 
 FILE = iscr
 ke = 0
 250 kb = ke + 1
 ke = ke + ncpp
 IF (ke > nrowa) ke = nrowa
 IF (ke /=  ncp7) GO TO 270
 IF (debug) WRITE (nout,260) (iz(j),j=i2,n2)
 260 FORMAT ('  IZ(I2...N2) =',18I6, /,(15X,18I6))
 CALL tmtogo (t2)
 t1 = (t1-t2)*0.143
 t1 = t1*FLOAT(npas)
 IF (t1 > t2) GO TO 880
 270 CALL REWIND (iscr)
 CALL fwdrec (*810,iscr)
 kbe = (ke-kb+1)*nwds
 DO  j = 1,kbe
   z(j) = 0.0
 END DO
 mm = n1 - 3
 ll = 0
 base = 2
 
!     MIDDLE I-LOOP
 
!     LOAD DATA FROM ISCR/NAMEA INTO Z(I3)...Z(NZ) WHEN NEEDED.
!     AND RUN THRU EACH ROW OF MATRIX NAMEA IN THIS LOOP
 
 i  = 0
 300 i  = i + 1
 IF (i > ncola) GO TO 430
 IF (base < ll) GO TO 340
 mm = mm + 3
 IF (kb == 1) GO TO 320
 
!     IF NOT FIRST PASS, CHECK KB AND KE AGAINST MIN/MAX TABLE IN IZ(I2)
!     THRU IZ(N2). IF THEY ARE OUTSIDE RANGE, SKIP NEXT DATA RECORD FROM
!     ISCR FILE AND UPDATE COLUMN COUNTER I
 
 IF (.NOT.(kb > iz(mm+2) .OR. ke < iz(mm+1))) GO TO 320
 CALL fwdrec (*810,iscr)
 i  = iz(mm+3)
 GO TO 300
 320 CALL READ (*810,*330,iscr,z(i3),iend1,1,ll)
 imhere = 160
 GO TO 820
 330 ll = n2 + ll
 base = n2 + 2
 340 ii = iz(base-1)
 jj = iz(base  )
 IF (kb > 1) GO TO 350
 
!     DURING FIRST PASS, SAVE MIN-II, MAX-JJ, AND COLUMN I IN IZ(MM)
!     TABLE. MM RUNS FROM I2 THRU N2.
 
 IF (ii < iz(mm+1)) iz(mm+1) = ii
 IF (jj > iz(mm+2)) iz(mm+2) = jj
 iz(mm+3) = i
 
 350 iikb = MAX0(ii,kb)
 jjke = MIN0(jj,ke)
 IF (jjke < iikb) GO TO 420
 
!     INNER K-LOOP
 
!     RUN THRU THE IIKB-JJKE ELEMENTS FOR EACH ROW OF MATRIX NAMEA,
 
!     KK = (IIKB-KB)*NWDS
!     LX = (I-1)*NWD + KK + 1
!     KK = BASE -  II*NWD + 1
 
 lx = (i-1)*nwd + (iikb-kb)*nwds + 1
 kx = (iikb-ii)*nwd + base + 1
 IF (nwd-2 < 0) THEN
   GO TO   360
 ELSE IF (nwd-2 == 0) THEN
   GO TO   380
 ELSE
   GO TO   400
 END IF
 360 DO  k = iikb,jjke
   z(lx) = z(kx)
   kx = kx + 1
   lx = lx + nwds
 END DO
 GO TO 420
 380 DO  k = iikb,jjke
   z(lx  ) = z(kx  )
   z(lx+1) = z(kx+1)
   kx = kx + 2
   lx = lx + nwds
 END DO
 GO TO 420
 400 DO  k = iikb,jjke
   z(lx  ) = z(kx  )
   z(lx+1) = z(kx+1)
   z(lx+2) = z(kx+2)
   z(lx+3) = z(kx+3)
   kx = kx + 4
   lx = lx + nwds
 END DO
 
!     END OF INNER K-LOOP
 
!     ADJUST BASE FOR ANOTHER ROW OF MATRIX NAMEA
 
 420 base = base + (jj-ii+1)*nwd + 2
 GO TO 300
 
!     END OF MIDDLE I-LOOP
 
!     PACK THE KB THRU KE COLUMNS OF THE TRANSPOSE MATRIX NAMEAT OUT
 
 430 DO  j = 1,kbe,nwds
   CALL pack (z(j),nameat,nameat)
 END DO
 
 IF (ke < nrowa) GO TO 250
 CALL CLOSE (iscr,1)
 
!     END OF OUTTER KB-KE LOOP, AND
!     END OF SQUARE AND RECTANGULAR MATRIX TRNASPOSE
 
!     OPEN AND CLOSE SCRATCH FILE AGAIN TO PHYSICALLY DELETE THE FILE.
!     MATRIX TRAILER WILL BE WRITTEN OUT BY DTRANP
 
 450 CALL CLOSE (nameat,clsrew)
 CALL gopen (iscr,z(ibuf1),wrtrew)
 CALL CLOSE (iscr,clsrew)
 GO TO 900
 
!     SYMMETRIC MATRIX
!     ================
 
 500 IF (ncola == nrowa) GO TO 520
 CALL fname (namea,a)
 WRITE  (nout,510) uwm,a,ncola,nrowa
 510 FORMAT (a25,' FROM TRNSP, ',2A4,' MATRIX (',i7,4H by ,i7,  &
     ') IS NOT SYMMETRIC NOR SQUARE ', /5X, 'IT WILL BE TREATED AS RECTANGULAR')
 CALL CLOSE (namea,clsrew)
 GO TO 30
 520 FILE   = nameat
 CALL OPEN (*800,nameat,z(ibuf),wrtrew)
 CALL cpyfil (namea,nameat,z(1),nz,k)
 CALL CLOSE (nameat,clsrew)
 CALL CLOSE (namea, clsrew)
 IF (debug) WRITE (nout,525) fileat
 525 FORMAT (' TRNSPS/@525 AFTER TRANSPOSE, TRAIL-AT =',7I8)
 GO TO 900
 
!     DIAGONAL MATRIX
!     ===============
!     DIAGONAL MATRIX (IFORMA=3) IS A ONE-COLUMN MATRIX. (1xN)
 
!     THE MATRIX AT RIGHT IS SQUARE (IFORMA=1),      1.  0.  0.
!     OR RECTANGULAR (IFORMA=2), AND IS NOT          0.  2.  0.
!     DIAGONAL (IFORMA=3) IN NASTRAN TERMINOLOGY     0.  0.  1.
 
 530 GO TO 520
 
!     IDENTITY MATRIX
!     ===============
!     SIMILAR TO DIAGONAL MATRIX, INDENTITY MATRIX (IFORMA = 8) IS ALSO
!     IN ONE-COLUMN MATRIX FORM
 
!     ALSO, THE IDENTITY MATRIX MAY EXIST ONLY IN THE MATRIX TRAILER.
!     IT DOES NOT PHYSICALLY EXIST.
 
 
 550 CALL READ (*900,*900,namea,z(1),1,1,j)
 CALL bckrec (namea)
 GO TO 520
 
!     ONE-COLUMN (1xN) RECTANGUALR MATRIX
!     ===================================
!     TRANSPOSE IS A ROW VECTOR, FORM=7. THE TRAILER REMAINS 1xN.
 
 580 IF (ncola /= 1) GO TO 860
 iforat = 8
 GO TO 520
 
!     UPPER OR LOWER TRIANGULAR MATRICES
!     ==================================
 
!     TRANSPOSE OF UPPER TRIANGULAR MATRIX IS THE LOWER TRIANG. MATRIX
!     AND VISE VERSA
 
!     (IS THIS HOW THE UPPER OR LOWER TRIANGULAR MATRIX WRITTEN? <==?
 
!     NO! IT IS NOT. WE STOP TRNSP SENDING THESE MATRICES OVER HERE.
!     BESIDE, THE LOGIC OF WRITING THE MATRIX BACKWARD HERE IS NOT
!     CORRECT. WE HAVE NOT ACCOMPLISHED THE TRANSPOSE OF THE ORIGINAL
!     MATRIX YET. ALSO, WE SHOULD WRITE THE TRANSPOSE MATRIX OUT BY
!     STRINGS, OR PACK THE MATRIX OUT)
 
 600 imhere = 600
 n1   = -37
 IF (n1 == -37) GO TO 830
 CALL gopen (iscr,z(ibuf),wrtrew)
 CALL skprec (namea,ncola)
 nwd  = iwords(itypa)
 irat = 3
 iend = (ibuf-1-nwd*ncola)/irat
 iend1= iend + 1
 isum = 0
 DO  i = 1,ncola
   iu   = 0
   CALL unpack (*830,namea,z(3))
   iz(1) = iu
   iz(2) = ju
   ll   = (ju-iu+1)*nwd + 2
   isum = isum + ll
   IF (isum <= iend) GO TO 610
   nrec = nrec + 1
   CALL WRITE (iscr,0,0,1)
   isum = ll
   610 IF (dgflag /= -123457890) GO TO 710
   IF (iforma == 5) GO TO 660
   SELECT CASE ( itypa )
     CASE (    1)
       GO TO 620
     CASE (    2)
       GO TO 630
     CASE (    3)
       GO TO 640
     CASE (    4)
       GO TO 650
   END SELECT
   620 z(3) = 1.0
   GO TO 710
   630 da = 1.0D+0
   z(3) = a(1)
   z(4) = a(2)
   GO TO 710
   640 z(4) = 0.0
   GO TO 620
   650 z(5) = 0.0
   z(6) = 0.0
   GO TO 630
   660 SELECT CASE ( itypa )
     CASE (    1)
       GO TO 670
     CASE (    2)
       GO TO 680
     CASE (    3)
       GO TO 690
     CASE (    4)
       GO TO 700
   END SELECT
   670 z(ju+2) = 1.0
   GO TO 710
   680 da = 1.0D+0
   z(ju*2+1) = a(1)
   z(ju*2+2) = a(2)
   GO TO 710
   690 z(ju*2+1) = 1.0
   z(ju*2+2) = 0.0
   GO TO 710
   700 j  = ju*4 - 3
   da = 1.0D+0
   z(j+1) = a(1)
   z(j+2) = a(2)
   z(j+3) = 0.0
   z(j+4) = 0.0
   710 CALL WRITE (iscr,z(1),ll,0)
   CALL bckrec (namea)
   CALL bckrec (namea)
 END DO
 nrec = nrec + 1
 CALL WRITE (iscr,0,0,1)
 CALL CLOSE (namea,clsrew)
 CALL CLOSE (iscr ,clsrew)
 itypat = itypa
 IF (iforma == 4) iforat = 5
 IF (iforma == 5) iforat = 4
 iat(1) = ia(1)
 iat(2) = ia(2)
 dgflag = 0
 filea(4) = nrec*10
 filea(6) = isum
 GO TO 30
 
!     ROW VECTOR (IFORMA=7, 1xN)
!     ==========================
 
!     A ROW VECTOR IS A ROW OF MATRIX ELEMENTS STORED IN COLUMN FORMAT
!     WITH TRAILER 1xN (NOT Nx1). THEREFORE THE TRANSPOSE OF ROW VECTOR
!     (IFORMA=7) IS A COLUMN VECTOR, WHICH IS RECTANG. (IFORAT=2).
!     THE TRAILER REMAINS UNCHANGED
 
 730 IF (ncola /= 1) GO TO 860
 iforat = 2
 GO TO 520
 
!     ERROR MESSAGES
 
 800 IF (iforma == 8) GO TO 900
 n1 = -1
 GO TO 850
 810 n1 = -2
 GO TO 850
 820 n1 = -8
 830 WRITE  (nout,840) imhere
 840 FORMAT (/5X,'IMHERE =',i5)
 850 CALL mesage (n1,FILE,nam)
 860 CALL fname (namea,a)
 WRITE  (nout,870) ufm,a,iforma,ncola,nrowa
 870 FORMAT (a23,' FROM TRNSPS, INPUT MATRIX ',2A4,' IS NOT SUITABLE ',  &
     'FOR MATRIX TRANSPOSE.', /5X,'FORM, COLUMN, ROW =',3I6)
 CALL mesage (-37,namea,nam)
 880 WRITE  (nout,890) ufm,t1
 890 FORMAT (a23,', INSUFFICIENT TIME REMAINING FOR MATRIX TRANSPOSE',  &
     /5X,'ESTIMATED TIME NEEDED (FOR TRANSPOSE ALONE) =',i9, ' CPU SECONDS')
 CALL mesage (-37,0,nam)
 
 900 RETURN
END SUBROUTINE trnsps
