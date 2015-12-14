SUBROUTINE mpydri (a,da,b,db,c,dc)
     
!     SPECIAL MPYAD PERFORMS THE MATRIX OPERATION
!        (+/-)A   *B (+/-)C = D   OR
!        (+/-)A(T)*B (+/-)C = D
 
!     WHERE A, OR B IS , OR BOTH ARE, DIAGONAL, ROW VECTOR, OR IDENTITY
!     MATRIX.  MATRIX C CAN BE PURGED
 
!     THIS ROUITNE DOES NOT HANDEL A-TRANSPOSE, WHILE B IS DIAGNOL, ROW
!     VECTOR, OR IDENTIY MASTRIX. ONLY EXCEPTION IS A IS TRULY (Nx1).
 
!     NOTE -
!     1. IN NASTRAN GINO, THE TRAILER 2ND AND 3RD WORDS FOR A ROW-VECTOR
!        IS (1xM), AND THE DIAGONAL MATRIX IS ALSO (1xM)
!     2. THE ROW-VECTOR AND DIAGONAL MATRIX ARE PACKED IN ONE RECORD.
!        AND THUS, THEY REQUIRE SPECIAL ATTENTION DEALING WITH THE FILEB
!        WHILE FILEA IS ALREADY A ROW-VECTOR, OR A DIAGONAL MATRIX
 
!     WRITTEN BY G.CHAN/UNISYS,  1/92
!     LAST MODIFIED FOR SPECIAL CASES THAT INVOLVE B MATRIX IS ALSO
!     A DIAGONAL MATRIX OR A ROW-VECOTR,  2/93                 ----
 
 
 REAL, INTENT(IN OUT)                     :: a(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: da(1)
 REAL, INTENT(IN OUT)                     :: b(1)
 DOUBLE PRECISION, INTENT(IN)             :: db(1)
 REAL, INTENT(OUT)                        :: c(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dc(1)
 IMPLICIT INTEGER (a-z)
 INTEGER :: NAME(2),ad(7),sd(7)
 
 
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /mpyadx/  filea(7),fileb(7),filec(7),filed(7),nzz,t,signab,  &
     signc,prec,scr
 COMMON /system/  sysbuf,nout
 COMMON /TYPE  /  prc(2),words(4)
 COMMON /names /  rd,rdrew,wrt,wrtrew,clsrew
 COMMON /packx /  typep,typout,ip,jp,incrp
 COMMON /unpakx/  typeu,iu,ju,incru
 COMMON /trnspx/  namea(7),nameat(7),lcore,nscr,iscr
 EQUIVALENCE      (filea(1),fa   ),(filea(4),forma),  &
     (filea(5),typea),(fileb(1),fb   ), (fileb(4),formb),(fileb(5),typeb),  &
     (filec(1),fc   ),(filec(4),formc), (filec(5),typec),(filed(1),fd   ),  &
     (filed(2),cold ),(filed(4),formd), (filed(5),typed)
 DATA    NAME   / 4HMPYA , 4HDRI /, diagnl , rowvec , ident / 3, 7, 8  /
 
!     MOVE TRUE ROWS AND COLUMNS INTO ROWA/B/C AND COLA/B/C
 
 cola = filea(2)
 rowa = filea(3)
 colb = fileb(2)
 rowb = fileb(3)
 colc = filec(2)
 rowc = filec(3)
 IF (forma == diagnl .OR. forma == rowvec) cola = rowa
 IF (forma == rowvec) rowa = 1
 IF (formb == diagnl .OR. formb == rowvec) colb = rowb
 IF (formb == rowvec) rowb = 1
 IF (formc == diagnl .OR. formc == rowvec) colc = rowc
 IF (formc == rowvec) rowc = 1
 
 IF (signab == 0 .AND. fc == 0) GO TO 1100
 IF (signab == 0 .AND. fc /= 0) GO TO 780
 buf1  = nzz  - sysbuf
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 cold  = 0
 rowd  = rowa
 IF (t == 1) rowd = cola
 IF (prec == 1 .AND. (typed == 2 .OR. typed == 4)) typed = typed -1
 typout= typed
 nwds  = words(typed)
 rowa2 = rowa*nwds
 rowb2 = rowb*nwds
 rowd2 = rowd*nwds
 colb2 = colb*2
 nz    = buf3 - 1
 sd(1) = scr
 IF (fc /= 0) GO TO 10
 sd(1) = fd
 nz    = buf2 - 1
 10 CALL makmcb (sd,sd,rowd,formd,typed)
 
!     REMEMBER, ONLY FILEA CAN HAVE TRANSPOSE, NOT FILEB.
!     IF FILEA IS DIAGONAL, ROW VECTOR, OR IDENTITY MATRIX, THE ACTUAL
!     TRANSPOSE OF FILEA HAS NEVER TAKEN PLACE.
 
!     FA, FB, FC AND FD ARE FILEA, FILEB, FILEC AND FILED RESPECTIVELY.
!     AD(1) IS EITHER FILEA OR FILED, AND
!     SD(1) IS EITHER SCRATCH FILE OR FILED
 
 IF (t == 1) GO TO 30
 DO  i = 1,7
   ad(i) = filea(i)
 END DO
 GO TO 50
 30 DO  i = 1,7
   ad(i) = filed(i)
 END DO
 50 ip    = 1
 jp    = rowd
 incrp = 1
 iu    = 1
 incru = 1
 IF (fa <= 0) GO TO 60
 FILE = fa
 CALL OPEN (*1010,fa,a(buf1),rdrew)
 CALL fwdrec (*1020,fa)
 60 IF (fb <= 0) GO TO 70
 FILE = fb
 CALL OPEN (*1010,fb,a(buf2),rdrew)
 CALL fwdrec (*1020,fb)
 
 70 IF (fa <= 0) GO TO 80
 IF (forma == diagnl) GO TO  90
 IF (forma == rowvec) GO TO 200
 IF (forma == ident ) GO TO 400
 80 IF (t == 1)  GO TO 990
 IF (formb == diagnl) GO TO 490
 IF (formb == rowvec) GO TO 600
 IF (formb == ident ) GO TO 750
 FILE = 0
 GO TO 1070
 
!                                         D   G   J   M
!     FILEA IS                            E   H   K   N
!     DIAGONAL -                          F   I   L   O
!                      a      a  0  0    aD  aG  aJ  aM
!                      b      0  b  0    bE  bH  bK  bN
!                      c ==>  0  0  c    cF  cI  cL  cO
 
!     SPECIAL CASE NEEDS TO BE CONSIDERED -
!     FILEB IS ALSO A DIAGONAL MATRIX. (FILEB CANNOT BE A ROW VECTOR)
 
 90 FILE  = fa
 ju    = rowa
 typeu = typed*signab
 CALL unpack (*1050,fa,a)
 CALL CLOSE (fa,clsrew)
 CALL gopen (sd,a(buf1),wrtrew)
 FILE  = fb
 ju    = rowb
 typeu = typed
 IF (formb == diagnl) GO TO 150
 DO  i = 1,colb
   CALL unpack (*1050,fb,b)
   SELECT CASE ( typeb )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO 110
     CASE (    3)
       GO TO 120
     CASE (    4)
       GO TO 130
   END SELECT
   100 DO  j = 1,rowb
     c(j) = a(j)*b(j)
   END DO
   GO TO 140
   110 DO  j = 1,rowb
     dc(j) = da(j)*db(j)
   END DO
   GO TO 140
   120 DO  j = 1,rowb2,2
     c(j  ) = a(j)*b(j  ) - a(j+1)*b(j+1)
     c(j+1) = a(j)*b(j+1) + a(j+1)*b(j  )
   END DO
   GO TO 140
   130 DO  j = 1,rowb2,2
     dc(j  ) = da(j)*db(j  ) - da(j+1)*db(j+1)
     dc(j+1) = da(j)*db(j+1) + da(j+1)*db(j  )
   END DO
   140 CALL pack (c,sd,sd)
 END DO
 GO TO 190
 
!     SPECIAL CASE - FILEB IS ALSO A DIAGONAL MATRIX
 
 150 CALL unpack (*1050,fb,b)
 IF (typeb >= 3) GO TO 165
 DO  j = 1,rowb
   c(j) = 0.0
 END DO
 DO  j = 1,rowb
   IF (typeb == 1)  c(j) =  a(j)* b(j)
   IF (typeb == 2) dc(j) = da(j)*db(j)
   CALL pack (c,sd,sd)
   c(j) = 0.0
   IF (typeb == 2) dc(j) = 0.0D+0
 END DO
 GO TO 190
 
 165 DO  j = 1,rowb2
   c(j) = 0.0
 END DO
 DO  j = 1,rowb2,2
   IF (typeb == 4) GO TO 175
   c(j  ) = a(j)*b(j  ) - a(j+1)*b(j+1)
   c(j+1) = a(j)*b(j+1) + a(j+1)*b(j  )
   CALL pack (c,sd,sd)
   c(j  ) = 0.0
   c(j+1) = 0.0
   CYCLE
   175 dc(j  ) = da(j)*db(j  ) - da(j+1)*db(j+1)
   dc(j+1) = da(j)*db(j+1) + da(j+1)*db(j  )
   CALL pack (c,sd,sd)
   dc(j  ) = 0.0D+0
   dc(j+1) = 0.0D+0
 END DO
 
 190 CALL CLOSE (fb,clsrew)
 CALL CLOSE (sd,clsrew)
 GO TO 800
!                                         E       I      M
!     FILEA IS A ROW     a                F       J      N
!     VECTOR -           b                G       K      O
!     RESULT IN FILED,   c                H       L      P
!     A (Nx1) RECT.      d ==> a b c d  aE+bF+  aI+bJ+  aM+bN+
!     MATRIX or A ROW-                  cG+dH   cK+dL   cO+dP
!     VECTOR
 
!     SPECIAL CASE NEEDS TO BE CONSIDERED -
!     FILEB IS A DIAGONAL MATRIX. (FILEB CANNOT BE A ROW VECTOR)
 
 
!     TRANSPOSE OF FILEA,                 E       F       G
!     A ROW VECTOR -               a     aE      aF      aG
!                                  b     bE      bF      bG
!                                  c     cE      cF      cG
!                                  d     dE      dF      dG
 
!     SPECIAL CASES NEED TO BE CONSIDERED -
!     FILEB MUST BE A (Nx1) RECTANGULAR MATRIX, OR A ROW VECTOR
 
 200 FILE  = fa
 ju    = rowa
 typeu = typed*signab
 CALL unpack (*1050,fa,a)
 CALL CLOSE (fa,clsrew)
 CALL gopen (sd,a(buf1),wrtrew)
 FILE  = fb
 typeu = typed
 IF (t == 1) GO TO 350
 
!     FILEA IS A ROW-VECTOR, RESULT IS ALSO A ROW-VECTOR, OR A
!     (Nx1) RECTANGULAR MATRIX
 
 ju    = rowb
 IF (formb == diagnl) GO TO 260
 IF (rowb  /=   rowa) GO TO 1030
 colb4 = colb*4
 DO  j = 1,colb4
   c(j) = 0.0
 END DO
 DO  j = 1,colb
   CALL unpack (*290,fb,b)
   SELECT CASE ( typeb )
     CASE (    1)
       GO TO 210
     CASE (    2)
       GO TO 220
     CASE (    3)
       GO TO 230
     CASE (    4)
       GO TO 240
   END SELECT
   210 DO  k = 1,rowb
     c(j) = c(j) + a(k)*b(k)
   END DO
   CYCLE
   220 DO  k = 1,rowb
     dc(j) = dc(j) + da(k)*db(k)
   END DO
   CYCLE
   230 DO  k = 1,rowb2,2
     c(j  ) = c(j  ) + a(k)*b(k  ) - a(k+1)*b(k+1)
     c(j+1) = c(j+1) + a(k)*b(k+1) + a(k+1)*b(k  )
   END DO
   CYCLE
   240 DO  k = 1,rowb,2
     dc(j  ) = dc(j  ) + da(k)*db(k  ) - da(k+1)*db(k+1)
     dc(j+1) = dc(j+1) + da(k)*db(k+1) + da(k+1)*db(k  )
   END DO
 END DO
 GO TO 300
 
!     SPECIAL CASE - FILEB IS A DIAGONAL MATRIX.
 
 260 CALL unpack (*1050,fb,b)
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 270
   CASE (    2)
     GO TO 280
   CASE (    3)
     GO TO 290
   CASE (    4)
     GO TO 300
 END SELECT
 270 DO  j = 1,colb
   c(j) = a(j)*b(j)
 END DO
 GO TO 310
 280 DO  j = 1,colb
   dc(j) = da(j)*db(j)
 END DO
 GO TO 310
 290 DO  j = 1,colb2,2
   c(j  ) = a(j)*b(j  ) - a(j+1)*b(j+1)
   c(j+1) = a(j)*b(j+1) + a(j+1)*b(j  )
 END DO
 GO TO 310
 300 DO  j = 1,colb2,2
   dc(j  ) = da(j)*db(j  ) - da(j+1)*db(j+1)
   dc(j+1) = da(j)*db(j+1) + da(j+1)*db(j  )
 END DO
 
 310 CALL CLOSE (fb,clsrew)
 IF (fc == 0) GO TO 340
 FILE  = fc
 typeu = typec*signc
 CALL gopen (fc,a(buf2),rdrew)
 IF (formc /= rowvec) GO TO 311
 CALL unpack (*1050,fc,a(1))
 GO TO 314
 311 ip = 1
 jp = 1
 DO  j = 1,colc
   CALL unpack (*312,fc,a(j*nwds-1))
   CYCLE
   312 a(j*nwds-1) = 0.
   a(j*nwds  ) = 0.
 END DO
 
 314 CALL CLOSE (fc,clsrew)
 SELECT CASE ( typed )
   CASE (    1)
     GO TO 315
   CASE (    2)
     GO TO 325
 END SELECT
 315 DO  j = 1,rowd2
   c(j) = c(j) + a(j)
 END DO
 GO TO 340
 325 DO  j = 1,rowd2
   dc(j) = dc(j) + da(j)
 END DO
 
 340 CALL pack (c,sd,sd)
 formd = rowvec
 GO TO 970
 
!     FILEA (A ROW VECTOR) TRANSFPOSE
 
 350 IF (formb == rowvec) GO TO 390
 IF (rowb /= 1) GO TO 1030
 iu = 0
 j  = 1
 DO  i = 1,rowb
   CALL unpack (*360,fb,b(j))
   IF (iu /= i) GO TO 1030
   GO TO 380
   360 je = j + nwds
   DO  k = j,je
     b(k) = 0.0
   END DO
   380 j  = j + nwds
 END DO
 CALL CLOSE (fb,clsrew)
 iu = 1
 GO TO 610
 
!     SPECAIL CASE - FILE B IS A ROW VECTOR
 
 390 IF (rowb /= 1) GO TO 1030
 ju = colb
 CALL unpack (*1030,fb,b(1))
 CALL CLOSE (fb,clsrew)
 GO TO 610
 
!     FILEA IS IDENTITY -
 
!     SPECIAL CASEs NEED TO BE CONSIDERED -
!     SIGNAB IS NEGATIVE, OR FILEB IS A DIAGONAL MATRIX
!     (FILEB CANNOT BE A ROW-VECTOR)
 
 400 CALL CLOSE (fa,clsrew)
 IF (formb == diagnl .OR. signab < 0) GO TO 420
 FILE = sd(1)
 CALL OPEN (*1010,fa,a(buf1),wrtrew)
 CALL REWIND (fb)
 CALL cpyfil (fb,sd,a(1),nz,k)
 CALL CLOSE (fb,clsrew)
 CALL CLOSE (sd,clsrew)
 IF (fc == 0) GO TO 410
 DO  i = 2,7
   sd(i) = fileb(i)
 END DO
 GO TO 800
 410 DO  i = 2,7
   filed(i) = fileb(i)
 END DO
 GO TO 1100
 
!     SPECIAL CASE - FILEB IS A DIAGONAL MATRIX
!                    OR SIGNAB IS NEGATIVE
 
 420 CALL gopen (sd,a(buf1),wrtrew)
 ju    = rowb
 FILE  = fb
 typeu = typed*signab
 IF (formb /= diagnl) GO TO 430
 CALL unpack (*1050,fb,b)
 CALL CLOSE (fb,clsrew)
 j  = 1
 DO  i = 1,rowa
   ip = i
   jp = i
   CALL pack (b(j),sd,sd)
   j = j + nwds
 END DO
 CALL CLOSE (sd,clsrew)
 IF (fc == 0.0) THEN
   GO TO   950
 ELSE
   GO TO   800
 END IF
 
!     SPECIAL CASE - SIGNAB IS NEGATIVE
 
 430 FILE = fb
 DO  i = 1,colb
   CALL unpack (*1050,fb,b)
   CALL pack (b,sd,sd)
 END DO
 CALL CLOSE (sd,clsrew)
 CALL CLOSE (fb,clsrew)
 IF (fc == 0.0) THEN
   GO TO   950
 ELSE
   GO TO   800
 END IF
 
!     FILEA IS A COLUMN MATRIX -
!     i.e. A (1,N) RECTANGULAR MATRIX OR A (Nx1) TRANSPOSE
 
!     FILEB MUST BE A (Nx1) RECTANGULAR MATRIX
 
!     CURRENTLY THIS CASE IS HANDLED IN MPYAD SUBROUINTE
 
!     HOWEVER, IF FILEB IS A ROW VECTOR,  IT IS HANDLED IN 600
!     IF FILEA IS A ROW VECTOR TRANSPOSE, IT IS HANDLED IN 200/350
 
! 440 CONTINUE
 
!                                         X   0   0      X
!     FILEB IS DIAGONAL -                 0   Y   0      Y
!                                         0   0   Z <==  Z
!                             a  e  i    aX  eY  iZ
!                             b  f  j    bX  fY  jZ
!                             c  g  k    cX  gY  kZ
!                             d  h  l    dX  hY  lZ
 
 490 FILE  = fb
 ju    = colb
 typeu = typed*signab
 CALL unpack (*1050,fb,b)
 CALL CLOSE (fb,clsrew)
 CALL gopen (sd,a(buf2),wrtrew)
 FILE  = fa
 ju    = rowa
 typeu = typed
 DO  i = 1,cola
   CALL unpack (*1050,fa,a)
   SELECT CASE ( typeb )
     CASE (    1)
       GO TO 500
     CASE (    2)
       GO TO 520
     CASE (    3)
       GO TO 540
     CASE (    4)
       GO TO 560
   END SELECT
   500 DO  j = 1,rowa
     c(j)  = a(j)*b(i)
   END DO
   GO TO 580
   520 DO  j = 1,rowa
     dc(j) = da(j)*db(i)
   END DO
   GO TO 580
   540 DO  j = 1,rowa2,2
     c(j  ) = a(j)*b(j  ) - a(j+1)*b(j+1)
     c(j+1) = a(j)*b(j+1) + a(j+1)*b(j  )
   END DO
   GO TO 580
   560 DO  j = 1,rowa2,2
     dc(j  ) = da(j)*db(j  ) - da(j+1)*db(j+1)
     dc(j+1) = da(j)*db(j+1) + da(j+1)*db(j  )
   END DO
   580 CALL pack (c,sd,sd)
 END DO
 CALL CLOSE (ad,clsrew)
 CALL CLOSE (sd,clsrew)
 GO TO 800
 
!     FILEB IS A ROW VECTOR -                            E
!                                                        F
!     NOTE - FILEA MUST BE A               E   F   G <== G
!     ONE-COLUMN MATRIX.             a    aE  aF  aG
!     i.e. A(1xN) OR                 b    bE  bF  bG
!          A(Nx1) TRNASPOSE          c    cE  cF  cG
!                                    d    dE  dF  dG
!     WE ALREADY HANDLED THE CASE
!     WHERE FILEA IS A ROW-VECTOR TRANSPOSE IN 200
 
 600 FILE  = fb
 ju    = colb
 typeu = typed*signab
 IF (t == 1) GO TO 602
 IF (cola /= 1) GO TO 1030
 CALL unpack (*1050,fb,b)
 GO TO 608
 602 IF (rowa /= 1) GO TO 1030
 j = cola*nwds
 DO  i = 1,j
   b(i) = 0.0
 END DO
 j = 1
 DO  i = 1,cola
   CALL unpack (*606,fb,b(j))
   j = j + nwds
 END DO
 608 CALL CLOSE (fb,clsrew)
 FILE  = fa
 ju    = rowa
 typeu = typed
 CALL unpack (*1050,fa,a)
 CALL CLOSE (ad,clsrew)
 CALL gopen (fd,a(buf1),wrtrew)
 610 DO  j = 1,colb
   SELECT CASE ( typea )
     CASE (    1)
       GO TO 620
     CASE (    2)
       GO TO 640
     CASE (    3)
       GO TO 660
     CASE (    4)
       GO TO 680
   END SELECT
   620 DO  i = 1,rowa
     c(i) = a(i)*b(j)
   END DO
   GO TO 700
   640 DO  i = 1,rowa
     da(i) = da(i)*db(j)
   END DO
   GO TO 700
   660 DO  i = 1,rowa2,2
     c(i  ) = a(i)*b(j  ) - a(i+1)*b(j+1)
     c(i+1) = a(i)*b(j+1) + a(i+1)*b(j  )
   END DO
   GO TO 700
   680 DO  i = 1,rowa2,2
     dc(i  ) = da(i)*db(j  ) - da(i+1)*db(j+1)
     dc(i+1) = da(i)*db(j+1) + da(i+1)*db(j  )
     kx = kx + nwds
   END DO
   700 CALL pack (c,fd,filed)
 END DO
 CALL CLOSE (fd,clsrew)
 GO TO 800
 
!     FILEB IS IDENTITY -
 
!     SPECIAL CASE NEEDS TO BE CONSIDERED -
!     NEGATIVE SIGNAB
 
 750 CALL CLOSE (fb,clsrew)
 FILE = sd(1)
 CALL OPEN (*1010,sd,a(buf2),wrtrew)
 IF (signab < 0) GO TO 760
 CALL REWIND (fa)
 CALL cpyfil (fa,sd,a(1),nz,k)
 GO TO 770
 
 760 typeu = typed*signab
 ju    = rowa
 FILE  = fa
 DO  i = 1,cola
   CALL unpack (*1050,fa,a)
   CALL pack (a,sd,sd)
 END DO
 770 CALL CLOSE (fa,clsrew)
 CALL CLOSE (sd,clsrew)
 IF (fc == 0.0) THEN
   GO TO   950
 ELSE
   GO TO   800
 END IF
 
!     NULL MATRIX PRODUCT A*B, COPY FILEC TO FILED
 
 780 FILE = fc
 CALL OPEN (*1010,fc,a(buf1),rdrew)
 FILE = fd
 CALL OPEN (*1010,fd,a(buf2),wrtrew)
 CALL cpyfil (fc,fd,a(1),nz,k)
 CALL CLOSE (fc,clsrew)
 CALL CLOSE (fd,clsrew)
 DO  i = 2,7
   filed(i) = filec(i)
 END DO
 GO TO 1100
 
!     ADD PRODUCT OF A*B TO C
 
 800 IF (fc == 0) GO TO 950
 CALL gopen (fd,a(buf3),wrtrew)
 FILE = fc
 CALL OPEN (*1010,fc,a(buf2),rdrew)
 CALL fwdrec (*1020,fc)
 FILE = sd(1)
 CALL OPEN (*1010,sd,a(buf1),rdrew)
 CALL fwdrec (*1020,sd)
 ju   = rowc
 typep = typed
 DO  i = 1,colc
   typeu = typed*signc
   CALL unpack (*810,fc,c)
   GO TO 830
   810 DO  j = 1,rowd2
     c(j) = 0.0
   END DO
   830 typeu = typed
   CALL unpack (*840,sd,b)
   GO TO 860
   840 DO  j = 1,rowd2
     b(j) = 0.0
   END DO
   860 SELECT CASE ( typed )
     CASE (    1)
       GO TO 870
     CASE (    2)
       GO TO 890
     CASE (    3)
       GO TO 870
     CASE (    4)
       GO TO 890
   END SELECT
   870 DO  j = 1,rowd2
     a(j) = b(j) + c(j)
   END DO
   GO TO 910
   890 DO  j = 1,rowd2
     da(j) = db(j) + dc(j)
   END DO
   910 CALL pack (a,fd,filed)
 END DO
 CALL CLOSE (fc,clsrew)
 CALL CLOSE (sd,clsrew)
 
 950 IF (cold /= 0) GO TO 970
 DO  i = 2,7
   filed(i) = sd(i)
 END DO
 970 CALL CLOSE  (fd,clsrew)
 CALL wrttrl (filed)
 GO TO 1100
 
!     ERROR
 
 990 WRITE  (nout,1000) sfm
 1000 FORMAT (a25,'. MPYDRI DOES NOT HANDLE A-TRANSPOSE. SHOULD NOT BE',  &
     ' CALLED BY MPYAD')
 GO TO 1070
 1010 j = -1
 GO TO 1080
 1020 j = -2
 GO TO 1080
 1030 WRITE  (nout,1040) ufm
 1040 FORMAT (a23,' FROM MPYAD/MPYDRI.  FILES NOT COMPATIBLE')
 GO TO 1070
 1050 WRITE  (nout,1060) ufm
 1060 FORMAT (a23,' FROM MPYAD/MPYDRI.  NULL COLUMN ENCOUNTERED DURING',  &
     ' MATRIX UNPACK')
 1070 j = -37
 1080 CALL mesage (j,FILE,NAME)
 
 1100 RETURN
END SUBROUTINE mpydri
