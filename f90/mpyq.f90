SUBROUTINE mpyq (z      )
     
!     MPYQ IS CALLED ONCE PER EXECUTION OF MPYAD. IT PERFORMS GENERAL
!     INITIALIZATION FOR EACH OF THE ENTRY POINTS
!     I.E.  SETTING UP MPYAD GO-TO-BRANCHES, ARITH AND BPICK FOR METHOD
!           1, AND ARITH2, APICK2 AND BPICK2 FOR METHOD 2
 
!     ENTRY POINTS -
!           MPY1V   (PERFORMS THE INNER LOOP FOR MPYAD, METHOD 1,
!                    TRANSPOSE AND NON-TRANSPOSE. IT IS CALLED ONCE FOR
!                    EACH COLUMNN OF THE A MATRIX)
!           MPY2NV  (PERFORMS THE INNER LOOP FOR THE NON-TRANSPOSE CASE
!                    OF METHOD 2. IT IS CALLED ONCE FOR EACH COUMN OF
!                    THE B MATRIX)
!           MPY2TV  (SAME AS MPY2NV, EXECPT IT IS FOR THE TRANSPOSE
!                    CASE)
!           MPY3T   (PERFORMS THE INNER LOOP FOR THE TRANSPOSE CASE OF
!                    METHOD 3. IT IS CALLED ONCE FOR EACH COLUMN OF THE
!                    B MATRIX)
!     (WHERE V STANDS FOR VAX VERSION, AND T IS THE TRANSPOSE FLAG)
 
!     THE MPYi ROUTINES PERFORM THE MATRIX MULTIPLICATION AND ADDITION
!     FOR THE MPYAD INNER LOOPS
 
!           (+/-)A * B (+/-)C = D   OR (+/-)A(T) * B (+/-)C = D
 
 
!     LAST REVISED BY G.CHAN/UNISYS  1/91
!     (1) MPY3T WAS PREVIOUSLY A .MDS SUBROUTINE. IT IS NOW AN ENTRY
!         POINT IN THIS MPYQ ROUTINE
!         (MPY3T IS AN ENTRY POINT IN MPYQ, IN IBM AND CDC VERSIONS)
!     (2) TO IMPROVE MPYAD INNER LOOP LOGIC FOR THE COMMON CASES
 
 
 
 REAL, INTENT(IN OUT)                     :: z(1)
 REAL :: a(4)    ,b(4)    ,d(4)    , aa(4) ,aaa  ,  &
     bsr     ,aas(1)  ,dds(1)  ,bbb     ,bsi   ,bbs
 DOUBLE PRECISION :: ad(2)   ,bd(2)   ,dd(2)   ,zd(1)   ,add(2),bdr  ,  &
     bdi     ,aad(1)  ,ddd(1)  ,bbd(1)
 INTEGER :: all , arith, arith2, apick2,  bpick2, bpick,  typeb,  typea, typed
 DIMENSION        zz(1)
 COMMON /mpyadx/  filea(7),fileb(7),filec(7),filed(7),nz    ,t    ,  &
     signab  ,signc   ,prec1   ,scrtch  ,time /system/  ksystm(65)  &
     /TYPE  /  prc(2)  ,nwds(4) ,rc(4)  &
     /names /  rd      ,rdrew   ,wrt     ,wrtrew  ,clsrew,cls  &
     /zblpkx/  d       ,drow /zntpkx/  a       ,i       ,eol     ,eor  &
     /packx /  typed   ,typd1   ,one1    ,pp1     ,incr1  &
     /unpakx/  typebd  ,one2    ,pp2     ,incr2
 COMMON /mpyadz/  rcb     ,rcd     ,ll      ,lll     ,jb    ,nbx  ,  &
     ndx     ,jmax1x  ,acol    ,acol1   ,acoln ,acore,  &
     apoint  ,bcol    ,crow    ,firstl  ,na    ,nb   ,  &
     nd      ,nwda    ,nwdb    ,nwdd    ,prec  ,jmax , incra   ,BLOCK(20)
 COMMON /mpyqt4/  rca     ,prca    ,all4    ,jump4   ,prec4
 COMMON /zzzzzz/  bbs(17000)
 EQUIVALENCE      (ksystm( 1),sysbuf)  ,(ksystm( 2),mout )  ,  &
     (ksystm(58),ipass )
 EQUIVALENCE      (a(1)    ,ad(1)   )  ,(b(1)      ,bd(1))  ,  &
     (d(1)    ,dd(1)   )  ,(filea(2),m      )  ,  &
     (filea(3),n       )  ,(filea(5),typea  )  ,  &
     (fileb(2),q       )  ,(fileb(3),r      )  ,  &
     (fileb(5),typeb   )  ,(filec(5),typec  )  ,  &
     (filed(5),typd    )  ,(nzz     ,buf1   )  ,  &
     (acoln   ,arown   )  ,(aa(1)   ,add(1) )  ,  &
     (acol1   ,arow1   )  ,(acol    ,arow   )  , (bbs(1)  ,bbd(1)  )
 EQUIVALENCE      (BLOCK(2),TYPE    )  ,(BLOCK(3),FORM   )  ,  &
     (BLOCK(4),row, j3 )  ,(BLOCK(5),point  )  ,  &
     (BLOCK(6),nbrstr  )  ,(BLOCK(8),flag   )
 EQUIVALENCE      (BLOCK(5),jb31    )  ,(BLOCK(6),nterm3 )  ,  &
     (BLOCK(7),jb3n    )  ,(BLOCK(8),b3flag )
!     DATA    MASK6F / X'00FFFFFF'  /
 DATA    mask6f / 16777215 /
 
!     MASK6F= '00FFFFFF'X (OR X'00FFFFFF') = 16777215
!     RCB   = 1 IF B IS REAL,2 IF B IS COMPLEX
!     NWDB  = NUMBER OF WORDS PER ELEMENT OF B
!     NBX   = NUMBER OF ELEMENTS PER COLUMN OF B
!     NB    = NUMBER OF WORDS PER COLUMN OF B
!     NDX   = NUMBER OF ELEMENTS PER COLUMN OF C AND D
!     ND    = NUMBER OF WORDS PER COLUMN OF C AND D
!     NZZ   = BUF1 = POINTER TO FIRST GINO BUFFER
!     BUF2  = POINTER TO SECOND GINO BUFFER
!     BUF3  = POINTER TO THIRD GINO BUFFER
!     JJ    = MAX. NO. OF COLUMNS OF B AND D THAT MAY BE HELD IN CORE
!     MPASS1= NUMBER OF PASSES REQUIRED FOR METHOD ONE
!     JZB   = POINTER TO FIRST ELEMENT OF B FOR SP REFERENCE
!     JZDB  = POINTER TO FIRST ELEMENT OF B FOR DP REFERENCE
!     JB    = POINTER TO FIRST ELEMENT OF B FOR PRECISION OF PROBLEM
!     NWDA  = NUMBER OF WORDS PER ELEMENT OF A
!     NWDD  = NUMBER OF WORDS PER ELEMENT OF D
!     ACORE = POINTER TO FIRST WORD FOR STORAGE OF PACKED COLUMNS
!             OF A MATRIX FOR METHOD TWO
!*****
 rca  = rc(typea)
 modb = MOD(typeb,2)
 moda = MOD(typea,2)
 prca = prc(typea)
 fa3  = filea(3)
 
!     IF DIAG 43 IS ON, SKIP ALL (1991) SPEED IMPROVEMENT LOGIC
!     (THIS IS ONLY TEMPORARY)
 
 all  = 0
 IF (typea == typeb .AND. typea == typed) all = typea
 all4 = all
 IF (all4 == 0) all4 = 5
 IF (typed >= 3 .AND. typea <= 2 .AND. typeb <= 2) all4 = 6
 CALL sswtch (43,j)
 IF (j == 1) all = 0
 jump4 = typeb + (typea-1)*4
 prec4 = prec1 - 1
 
!     RCA, PRCA, ALL4 AND JUMP4 ARE USED IN MPY4T
 
 SELECT CASE ( typed )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 40
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 300
 END SELECT
 
!     REAL SINGLE PRECISION
 
 20 ASSIGN 840 TO arith
 ASSIGN 750 TO bpick
 IF (t /= 0)  GO TO 30
 
 ASSIGN 1210 TO arith2
 ASSIGN 1010 TO bpick2
 ASSIGN 1100 TO apick2
 GO TO 600
 
 30 ASSIGN 1600 TO arith2
 ASSIGN 1400 TO bpick2
 ASSIGN 1500 TO apick2
 GO TO 600
 
!     REAL DOUBLE PRECISION
 
 40 ASSIGN 850 TO arith
 ASSIGN 770 TO bpick
 IF (modb == 1) ASSIGN 790 TO bpick
 IF (t /= 0) GO TO 60
 
 ASSIGN 1220 TO arith2
 ASSIGN 1030 TO bpick2
 IF (modb == 1) ASSIGN 1050 TO bpick2
 ASSIGN 1120 TO apick2
 50 IF (moda == 1) ASSIGN 1140 TO apick2
 GO TO 600
 
 60 ASSIGN 1610 TO arith2
 ASSIGN 1420 TO bpick2
 IF (modb == 1) ASSIGN 1440 TO bpick2
 ASSIGN 1520 TO apick2
 70 IF (moda == 1) ASSIGN 1540 TO apick2
 GO TO 600
 
!     COMPLEX SINGLE PRECISION
 
 100 ASSIGN 860 TO arith
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 120
   CASE (    4)
     GO TO 130
 END SELECT
 110 ASSIGN 750 TO bpick
 GO TO 140
 120 ASSIGN 760 TO bpick
 GO TO 140
 130 ASSIGN 810 TO bpick
 140 IF (t /= 0)  GO TO 220
 
 ASSIGN 1230 TO arith2
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 150
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 160
   CASE (    4)
     GO TO 170
 END SELECT
 150 ASSIGN 1010 TO bpick2
 GO TO 180
 160 ASSIGN 1020 TO bpick2
 GO TO 180
 170 ASSIGN 1070 TO bpick2
 180 SELECT CASE ( typea )
   CASE (    1)
     GO TO 190
   CASE (    2)
     GO TO 190
   CASE (    3)
     GO TO 200
   CASE (    4)
     GO TO 210
 END SELECT
 190 ASSIGN 1100 TO apick2
 GO TO 600
 200 ASSIGN 1110 TO apick2
 GO TO 600
 210 ASSIGN 1160 TO apick2
 GO TO 600
 220 ASSIGN 1620 TO arith2
 
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 230
   CASE (    2)
     GO TO 230
   CASE (    3)
     GO TO 240
   CASE (    4)
     GO TO 250
 END SELECT
 230 ASSIGN 1400 TO bpick2
 GO TO 260
 240 ASSIGN 1410 TO bpick2
 GO TO 260
 250 ASSIGN 1460 TO bpick2
 260 SELECT CASE ( typea )
   CASE (    1)
     GO TO 270
   CASE (    2)
     GO TO 270
   CASE (    3)
     GO TO 280
   CASE (    4)
     GO TO 290
 END SELECT
 270 ASSIGN 1500 TO apick2
 GO TO 600
 280 ASSIGN 1510 TO apick2
 GO TO 600
 290 ASSIGN 1560 TO apick2
 GO TO 600
 
!     COMPLEX DOUBLE PRECISION
 
 300 ASSIGN 870 TO arith
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 310
   CASE (    2)
     GO TO 320
   CASE (    3)
     GO TO 330
   CASE (    4)
     GO TO 340
 END SELECT
 310 ASSIGN 790 TO bpick
 GO TO 350
 320 ASSIGN 770 TO bpick
 GO TO 350
 330 ASSIGN 800 TO bpick
 GO TO 350
 340 ASSIGN 780 TO bpick
 350 IF (t /= 0) GO TO 440
 
 ASSIGN 1240 TO arith2
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 360
   CASE (    2)
     GO TO 370
   CASE (    3)
     GO TO 380
   CASE (    4)
     GO TO 390
 END SELECT
 360 ASSIGN 1050 TO bpick2
 GO TO 400
 370 ASSIGN 1030 TO bpick2
 GO TO 400
 380 ASSIGN 1060 TO bpick2
 GO TO 400
 390 ASSIGN 1040 TO bpick2
 400 SELECT CASE ( typea )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 420
   CASE (    4)
     GO TO 430
 END SELECT
 410 ASSIGN 1120 TO apick2
 GO TO 600
 420 ASSIGN 1150 TO apick2
 GO TO 600
 430 ASSIGN 1130 TO apick2
 GO TO 600
 
 440 ASSIGN 1630 TO arith2
 SELECT CASE ( typeb )
   CASE (    1)
     GO TO 450
   CASE (    2)
     GO TO 460
   CASE (    3)
     GO TO 470
   CASE (    4)
     GO TO 480
 END SELECT
 450 ASSIGN 1440 TO bpick2
 GO TO 490
 460 ASSIGN 1420 TO bpick2
 GO TO 490
 470 ASSIGN 1450 TO bpick2
 GO TO 490
 480 ASSIGN 1430 TO bpick2
 490 SELECT CASE ( typea )
   CASE (    1)
     GO TO 70
   CASE (    2)
     GO TO 500
   CASE (    3)
     GO TO 510
   CASE (    4)
     GO TO 520
 END SELECT
 500 ASSIGN 1520 TO apick2
 GO TO 600
 510 ASSIGN 1550 TO apick2
 GO TO 600
 520 ASSIGN 1530 TO apick2
 
!     MPYQ INITIALIZATION DONE
 
 600 RETURN
 
 
 ENTRY mpy1v (zz      ,z      ,zd      )
!     =====================
 
!     METHOD 1  (TRANSPOSE AND NON-TRANSPOSE)
 
 700 b (2) = 0.
 bd(2) = 0.d0
 710 CALL zntpki
 i1 = i - 1
 IF (t == 0.0) THEN
   GO TO   720
 ELSE
   GO TO   730
 END IF
 720 k1 = ll
 k2 = i1*rcd + 1
 GO TO 740
 730 k1 = i1*rcb + jb
 k2 = lll
 740 k3 = k2 + jmax1x
 IF (all /= 0) THEN
    SELECT CASE ( all )
     CASE (    1)
       GO TO 900
     CASE (    2)
       GO TO 920
     CASE (    3)
       GO TO 940
     CASE (    4)
       GO TO 960
   END SELECT
 END IF
 DO  k = k2,k3,ndx
   j  = k1
   GO TO bpick, (750,760,770,780,790,800,810)
   750 IF (z(j) == 0.0) GO TO 880
   b(1) = z(j)
   GO TO 830
   760 IF (z(j) == 0.0 .AND. z(j+1) == 0.0) GO TO 880
   b(1) = z(j  )
   b(2) = z(j+1)
   GO TO 830
   770 IF (zd(j) == 0.d0) GO TO 880
   bd(1) = zd(j)
   GO TO 830
   780 IF (zd(j) == 0.d0 .AND. zd(j+1) == 0.d0) GO TO 880
   bd(1) = zd(j  )
   bd(2) = zd(j+1)
   GO TO 830
   790 IF (z(j) == 0.0) GO TO 880
   bd(1) = z(j)
   GO TO 830
   800 IF (z(j) == 0.0 .AND. z(j+1) == 0.0) GO TO 880
   bd(1) = z(j  )
   bd(2) = z(j+1)
   GO TO 830
   810 IF (zd(j) == 0.d0 .AND. zd(j+1) == 0.d0) GO TO 880
   b(1) = zd(j  )
   b(2) = zd(j+1)
   
   830 GO TO arith, (840,850,860,870)
   840 z(k) = z(k) + a(1)*b(1)
   GO TO 880
   850 zd(k) = zd(k) + ad(1)*bd(1)
   GO TO 880
   860 z(k  ) = z(k  ) + a(1)*b(1) - a(2)*b(2)
   z(k+1) = z(k+1) + a(1)*b(2) + a(2)*b(1)
   GO TO 880
   870 zd(k  ) = zd(k  ) + ad(1)*bd(1) - ad(2)*bd(2)
   zd(k+1) = zd(k+1) + ad(1)*bd(2) + ad(2)*bd(1)
   880 k1 = k1 + nbx
 END DO
 IF (eol == 0.0) THEN
   GO TO   710
 ELSE
   GO TO   980
 END IF
 
!     COMMON CASES (TYPEA=TYPEB=TYPED=PREC)
 
!     PREC=1, ARITH(840) AND BPICK(750)
!     PREC=2, ARITH(850) AND BPICK(770)
!     PREC=3, ARITH(860) AND BPICK(760)
!     PREC=4, ARITH(870) AND BPICK(780)
 
 900 DO  k = k2,k3,ndx
   z(k) = z(k) + a(1)*z(k1)
   k1 = k1 + nbx
 END DO
 IF (eol == 0.0) THEN
   GO TO   710
 ELSE
   GO TO   980
 END IF
 920 DO  k = k2,k3,ndx
   zd(k) = zd(k) + ad(1)*zd(k1)
   k1 = k1 + nbx
 END DO
 IF (eol == 0.0) THEN
   GO TO   710
 ELSE
   GO TO   980
 END IF
 940 DO  k = k2,k3,ndx
   z(k  ) = z(k  ) + a(1)*z(k1  ) - a(2)*z(k1+1)
   z(k+1) = z(k+1) + a(1)*z(k1+1) + a(2)*z(k1  )
   k1 = k1 + nbx
 END DO
 IF (eol == 0.0) THEN
   GO TO   710
 ELSE
   GO TO   980
 END IF
 960 DO  k = k2,k3,ndx
   zd(k  ) = zd(k  ) + ad(1)*zd(k1  ) - ad(2)*zd(k1+1)
   zd(k+1) = zd(k+1) + ad(1)*zd(k1+1) + ad(2)*zd(k1  )
   k1 = k1 + nbx
 END DO
 IF (eol == 0.0) THEN
   GO TO   710
 END IF
 980 RETURN
 
 
 ENTRY mpy2nv (zz      ,z      ,zd      )
!     ======================
 
!     METHOD 2 NON-TRANSPOSE CASE
 
 b(2)  = 0.
 bd(2) = 0.d0
 aa(2) = 0.
 add(2)= 0.d0
 l     = firstl
 acol  = acol1
 1000 CALL zntpki
 IF (i < acol1 .OR. i > acoln .OR. i < acol) GO TO 1290
 l     = l - 2*(i-acol)
 acol  = i
 apoint= zz(l)
 IF (apoint == 0) GO TO 1280
 nbr   = zz(l-1)
 IF (all /= 0) THEN
    SELECT CASE ( all )
     CASE (    1)
       GO TO 1260
     CASE (    2)
       GO TO 1265
     CASE (    3)
       GO TO 1270
     CASE (    4)
       GO TO 1275
   END SELECT
 END IF
 GO TO bpick2, (1010,1020,1030,1040,1050,1060,1070)
 1010 b(1)  = a(1)
 GO TO 1090
 1020 b(1)  = a(1)
 b(2)  = a(2)
 GO TO 1090
 1030 bd(1) = ad(1)
 GO TO 1090
 1040 bd(1) = ad(1)
 bd(2) = ad(2)
 GO TO 1090
 1050 bd(1) = a(1)
 GO TO 1090
 1060 bd(1) = a(1)
 bd(2) = a(2)
 GO TO 1090
 1070 b(1)  = ad(1)
 b(2)  = ad(2)
 
 1090 nbrstr = zz( apoint+1 )
 init   = zz( apoint )
 apoint = apoint + 2
 j      = apoint
 IF ( prca == 2 ) j = j/2 + 1
 apoint = apoint+ nbrstr*nwda
 irow   = init*rcd - rcd + 1
 nrow   = irow + nbrstr*rcd - 1
 DO  k = irow,nrow,rcd
   GO TO apick2, (1100,1110,1120,1130,1140,1150,1160)
   1100 aa(1) = z(j)
   GO TO 1200
   1110 aa(1) = z(j  )
   aa(2) = z(j+1)
   GO TO 1200
   1120 add(1)= zd(j)
   GO TO 1200
   1130 add(1)= zd(j  )
   add(2)= zd(j+1)
   GO TO 1200
   1140 add(1)= z(j)
   GO TO 1200
   1150 add(1)= z(j  )
   add(2)= z(j+1)
   GO TO 1200
   1160 aa(1) = zd(j  )
   aa(2) = zd(j+1)
   
   1200 GO TO arith2, (1210,1220,1230,1240)
   1210 z(k)  = z(k) + aa(1)*b(1)
   GO TO 1250
   1220 zd(k) = zd(k) + add(1)*bd(1)
   GO TO 1250
   1230 z(k  )= z(k  ) + aa(1)*b(1) - aa(2)*b(2)
   z(k+1)= z(k+1) + aa(1)*b(2) + aa(2)*b(1)
   GO TO 1250
   1240 zd(k  ) = zd(k  ) + add(1)*bd(1) - add(2)*bd(2)
   zd(k+1) = zd(k+1) + add(1)*bd(2) + add(2)*bd(1)
   1250 j = j + rca
 END DO
 nbr  = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1090
 ELSE
   GO TO  1280
 END IF
 
!     COMMON CASES (TYPEA=TYPEB=TYPED=PREC)
 
!     PREC=1, ARITH2(1210), APICK2(1100) AND BPICK2(1010)
!     PREC=2, ARITH2(1220), APICK2(1120) AND BPICK2(1030)
!     PREC=3, ARITH2(1230), APICK2(1110) AND BPICK2(1020)
!     PREC=4, ARITH2(1620), APICK2(1510) AND BPICK2(1410)
 
 1260 nbrstr = zz( apoint+1 )
 init   = zz( apoint )
 apoint = apoint + 2
 j      = apoint
 IF ( prca == 2 ) j = j/2 + 1
 apoint = apoint+ nbrstr*nwda
 irow  = init*rcd - rcd + 1
 nrow  = irow + nbrstr*rcd - 1
 DO  k = irow,nrow,rcd
   z(k)  = z(k) + z(j)*a(1)
   j     = j + rca
 END DO
 nbr   = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1260
 ELSE
   GO TO  1280
 END IF
 
 1265 nbrstr = zz( apoint+1 )
 init   = zz( apoint )
 apoint = apoint + 2
 j      = apoint
 IF ( prca == 2 ) j = j/2 + 1
 apoint = apoint+ nbrstr*nwda
 irow  = init*rcd - rcd + 1
 nrow  = irow + nbrstr*rcd - 1
 DO  k = irow,nrow,rcd
   zd(k) = zd(k) + zd(j)*ad(1)
   j     = j + rca
 END DO
 nbr   = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1265
 ELSE
   GO TO  1280
 END IF
 
 1270 nbrstr = zz( apoint+1 )
 init   = zz( apoint )
 apoint = apoint + 2
 j      = apoint
 IF ( prca == 2 ) j = j/2 + 1
 apoint = apoint+ nbrstr*nwda
 irow  = init*rcd - rcd + 1
 nrow  = irow + nbrstr*rcd - 1
 DO  k = irow,nrow,rcd
   z(k  )  = z(k  ) + z(j)*a(1) - z(j+1)*a(2)
   z(k+1)  = z(k+1) + z(j)*a(2) + z(j+1)*a(1)
   j     = j + rca
 END DO
 nbr   = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1270
 ELSE
   GO TO  1280
 END IF
 
 1275 nbrstr = zz( apoint+1 )
 init   = zz( apoint )
 apoint = apoint + 2
 j      = apoint
 IF ( prca == 2 ) j = j/2 + 1
 apoint = apoint+ nbrstr*nwda
 irow  = init*rcd - rcd + 1
 nrow  = irow + nbrstr*rcd - 1
 DO  k = irow,nrow,rcd
   zd(k  ) = zd(k  ) + zd(j)*ad(1) - zd(j+1)*ad(2)
   zd(k+1) = zd(k+1) + zd(j)*ad(2) + zd(j+1)*ad(1)
   j     = j + rca
 END DO
 nbr   = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1275
 END IF
 
 1280 l = l - 2
 acol = acol + 1
 1290 IF (eol == 0) GO TO 1000
 RETURN
 
 
 ENTRY mpy2tv (zz      ,z      ,zd      )
!     ======================
 
!     METHOD 2 - TRANSPOSE CASE
 
!     COMMENTS FROM G.CHAN/UNISYS      1/91
!     OBSERVE THAT THERE IS NO DO-LOOP IN THIS MPY2TV LOGIC. IT IS
!     THEREFORE CONCLUDED THAT THE TRANSPOSE CASE WOULD TAKE MUCH MORE
!     TIME THAN THE NON-TRANSPOSE CASE
 
 b(2)  = 0.
 bd(2) = 0.d0
 aa(2) = 0.
 add(2)= 0.d0
 dd(1) = 0.d0
 dd(2) = 0.d0
 l     = firstl
 apoint= zz(l)
 arow  = arow1
 IF (crow == mask6f) GO TO 1350
 GO TO 1330
 1300 apoint = zz(l)
 IF (crow-arow < 0.0) THEN
   GO TO  1320
 ELSE IF (crow-arow == 0.0) THEN
   GO TO  1340
 ELSE
   GO TO  1350
 END IF
 1310 drow = crow
 CALL zblpki
 1320 IF (eol /= 0) GO TO 1350
 1330 CALL zntpki
 crow  = i
 1340 dd(1) = ad(1)
 dd(2) = ad(2)
 IF (crow-arow < 0.0) THEN
   GO TO  1310
 ELSE IF (crow-arow == 0.0) THEN
   GO TO  1360
 END IF
 1350 dd(1) = 0.d0
 dd(2) = 0.d0
 IF (apoint == 0) GO TO 1690
 1360 drow  = arow
 IF (apoint == 0) GO TO 1680
 nbrstr= zz(l-1)
 1370 nbr   = zz( apoint+1 )
 nbr1  = nbr
 init  = zz( apoint )
 apoint = apoint + 2
 j = apoint
 IF ( prca > 1 ) j = j/2 + 1
 apoint = apoint + nbr*nwda
 k     = (init-1)*rcb + 1
 1380 SELECT CASE ( all )
   CASE (    1)
     GO TO 1640
   CASE (    2)
     GO TO 1645
   CASE (    3)
     GO TO 1650
   CASE (    4)
     GO TO 1655
 END SELECT
 GO TO bpick2, (1400,1410,1420,1430,1440,1450,1460)
 1400 b(1)  = z(k)
 GO TO 1470
 1410 b(1)  = z(k  )
 b(2)  = z(k+1)
 GO TO 1470
 1420 bd(1) = zd(k)
 GO TO 1470
 1430 bd(1) = zd(k  )
 bd(2) = zd(k+1)
 GO TO 1470
 1440 bd(1) = z(k)
 GO TO 1470
 1450 bd(1) = z(k  )
 bd(2) = z(k+1)
 GO TO 1470
 1460 b(1)  = zd(k  )
 b(2)  = zd(k+1)
 
 1470 GO TO apick2, (1500,1510,1520,1530,1540,1550,1560)
 1500 aa(1) = z(j)
 GO TO 1570
 1510 aa(1) = z(j  )
 aa(2) = z(j+1)
 GO TO 1570
 1520 add(1)= zd(j)
 GO TO 1570
 1530 add(1)= zd(j  )
 add(2)= zd(j+1)
 GO TO 1570
 1540 add(1)= z(j)
 GO TO 1570
 1550 add(1)= z(j  )
 add(2)= z(j+1)
 GO TO 1570
 1560 aa(1) = z(j  )
 aa(2) = z(j+2)
 
 1570 GO TO arith2, (1600,1610,1620,1630)
 1600 d(1)  = d(1) + aa(1)*b(1)
 GO TO 1660
 1610 dd(1) = dd(1) + add(1)*bd(1)
 GO TO 1660
 1620 d(1)  = d(1) + aa(1)*b(1) - aa(2)*b(2)
 d(2)  = d(2) + aa(1)*b(2) + aa(2)*b(1)
 GO TO 1660
 1630 dd(1) = dd(1) + add(1)*bd(1) - add(2)*bd(2)
 dd(2) = dd(2) + add(1)*bd(2) + add(2)*bd(1)
 GO TO 1660
 
!     COMMON CASES (TYPEA=TYPEB=TYPED=PREC)
 
!     PREC=1, ARITH2(1600), APICK2(1500) AND BPICK2(1400)
!     PREC=2, ARITH2(1610), APICK2(1520) AND BPICK2(1420)
!     PREC=3, ARITH2(1620), APICK2(1510) AND BPICK2(1410)
!     PREC=4, ARITH2(1630), APICK2(1530) AND BPICK2(1430)
 
 1640 d(1) = d(1) + z(j)*z(k)
 j = j + rca
 k = k + rcb
 nbr = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1640
 ELSE
   GO TO  1670
 END IF
 
 1645 dd(1) = dd(1) + zd(j)*zd(k)
 j = j + rca
 k = k + rcb
 nbr = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1645
 ELSE
   GO TO  1670
 END IF
 
 1650 d(1) = d(1) + z(j)*z(k  ) - z(j+1)*z(k+1)
 d(2) = d(2) + z(j)*z(k+1) + z(j+1)*z(k  )
 j = j + rca
 k = k + rcb
 nbr = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1650
 ELSE
   GO TO  1670
 END IF
 
 1655 dd(1) = dd(1) + zd(j)*zd(k  ) - zd(j+1)*zd(k+1)
 dd(2) = dd(2) + zd(j)*zd(k+1) + zd(j+1)*zd(k  )
 j = j + rca
 k = k + rcb
 nbr = nbr - 1
 IF (nbr > 0) THEN
   GO TO  1655
 ELSE
   GO TO  1670
 END IF
 
 1660 j = j + rca
 k = k + rcb
 nbr = nbr - 1
 IF (nbr > 0) GO TO 1380
 1670 nbrstr = nbrstr - 1
 IF (nbrstr > 0) GO TO 1370
 1680 CALL zblpki
 1690 l = l - 2
 arow = arow + 1
 IF (arow <= arown) GO TO 1300
 RETURN
 
 
 ENTRY mpy3t (*,aas      ,aad      ,dds      ,ddd      )
!     ===============================
 
!     METHOD 3 (TRANSPOSE ONLY)
 
 b3flag = -1
 CALL getstr (*2400,BLOCK)
!IBMNB 6/93
 IF ( BLOCK( 2 ) == typeb ) GO TO 1699
 typeb = BLOCK( 2 )
 rcb   = rc( typeb )
 all   = 0
 1699 CONTINUE
!IBMNE 6/93
 IF (all /= 0) THEN
    SELECT CASE ( all )
     CASE (    1)
       GO TO 1760
     CASE (    2)
       GO TO 1920
     CASE (    3)
       GO TO 2060
     CASE (    4)
       GO TO 2270
   END SELECT
 END IF
 SELECT CASE ( typed )
   CASE (    1)
     GO TO 1700
   CASE (    2)
     GO TO 1800
   CASE (    3)
     GO TO 2000
   CASE (    4)
     GO TO 2100
 END SELECT
 
!     PERFORM ARITHMETIC IN REAL SINGLE PRECISION
 
 1700 jb3n = jb31 + nterm3 - 1
 DO  jb3 = jb31,jb3n
   k = j3
   bbb = bbs(jb3)
   IF (BLOCK(2) == 2) bbb = bbd(jb3)
   IF (typea    /= 2) GO TO 1720
   DO  i = arow1,arown
     aaa = aad(k)
     dds(i) = dds(i) + aaa*bbb
     k  = k + na
   END DO
   GO TO 1740
   1720 DO  i = arow1,arown
     dds(i) = dds(i) + aas(k)*bbb
     k  = k + na
   END DO
   1740 j3 = j3 + 1
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 1700
 
!     COMMON CASE (TYPEA=TYPEB=TYPED=PREC=1)
 
 1760 jb3n = jb31 + nterm3 - 1
 DO  jb3 = jb31,jb3n
   k  = j3
   DO  i = arow1,arown
     dds(i) = dds(i) + aas(k)*bbs(jb3)
     k  = k  + na
   END DO
   j3 = j3 + 1
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 1760
 
!     PERFORM ARITHMETIC IN REAL DOUBLE PRECISION
 
 1800 k1 = 2*(prca-1) + prc(typeb)
 1810 jb3n = jb31 + nterm3 - 1
 DO  jb3 = jb31,jb3n
   k = j3
   SELECT CASE ( k1 )
     CASE (    1)
       GO TO 1820
     CASE (    2)
       GO TO 1840
     CASE (    3)
       GO TO 1860
     CASE (    4)
       GO TO 1880
   END SELECT
   1820 DO  i = arow1,arown
     ddd(i) = ddd(i) + aas(k)*bbs(jb3)
     k = k + fa3
   END DO
   GO TO 1900
   1840 DO  i = arow1,arown
     ddd(i) = ddd(i) + aas(k)*bbd(jb3)
     k = k + fa3
   END DO
   GO TO 1900
   1860 DO  i = arow1,arown
     ddd(i) = ddd(i) + aad(k)*bbs(jb3)
     k = k + fa3
   END DO
   GO TO 1900
   1880 DO  i = arow1,arown
     ddd(i) = ddd(i) + aad(k)*bbd(jb3)
     k = k + fa3
   END DO
   1900 j3 = j3 + 1
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 1810
 
!     COMMON CASE (TYPEA=TYPEB=TYPED=PREC=2)
 
 1920 jb3n = jb31 + nterm3 - 1
 DO  jb3 = jb31,jb3n
   k = j3
   DO  i = arow1,arown
     ddd(i) = ddd(i) + aad(k)*bbd(jb3)
     k = k  + fa3
   END DO
   j3 = j3 + 1
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 1920
 
!     PERFORM ARITHMETIC IN COMPLEX SINGLE PRECISION
 
 2000 bsi  = 0.
 i1   = 2*arow1 - 1
 in   = 2*arown - 1
 2010 IF (rca == 2) j3 = 2*j3 - 1
 jb3n = jb31 + rcb*nterm3 - rcb
 DO  jb3 = jb31,jb3n,rcb
   bsr = bbs(jb3)
   IF (rcb == 2) bsi = bbs(jb3+1)
   k = j3
   IF (rca == 2) GO TO 2030
   DO  i = i1,in,2
     dds(i  ) = dds(i  ) + aas(k)*bsr
     dds(i+1) = dds(i+1) + aas(k)*bsi
     k = k + na
   END DO
   GO TO 2050
   2030 DO  i = i1,in,2
     dds(i  ) = dds(i  ) + aas(k)*bsr - aas(k+1)*bsi
     dds(i+1) = dds(i+1) + aas(k)*bsi + aas(k+1)*bsr
     k = k + na
   END DO
   2050 j3 = j3 + rca
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 2010
 
!     COMMON CASE (TYPEA=TYPEB=TYPED=PREC=3)
 
 2060 i1 = 2*arow1 - 1
 in = 2*arown - 1
 2070 j3 = 2*j3 - 1
 jb3n = jb31 + rcb*nterm3 - rcb
 DO  jb3 = jb31,jb3n,rcb
   k = j3
   DO  i = i1,in,2
     dds(i  ) = dds(i  ) + aas(k)*bbs(jb3  ) - aas(k+1)*bbs(jb3+1)
     dds(i+1) = dds(i+1) + aas(k)*bbs(jb3+1) + aas(k+1)*bbs(jb3  )
     k = k + na
   END DO
   j3 = j3 + rca
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 2070
 
!     PERFORM ARITHMETIC IN COMPLEX DOUBLE PRECISION
 
 2100 bdi = 0.
 inca= rca*fa3
 i1  = 2*arow1 - 1
 in  = 2*arown - 1
 2110 IF (rca == 2) j3 = 2*j3 - 1
 jb3n = jb31 + rcb*nterm3 - rcb
 DO  jb3 = jb31,jb3n,rcb
   k = j3
   SELECT CASE ( typeb )
     CASE (    1)
       GO TO 2120
     CASE (    2)
       GO TO 2130
     CASE (    3)
       GO TO 2140
     CASE (    4)
       GO TO 2150
   END SELECT
   2120 bdr = bbs(jb3)
   GO TO 2160
   2130 bdr = bbd(jb3)
   GO TO 2160
   2140 bdr = bbs(jb3  )
   bdi = bbs(jb3+1)
   GO TO 2160
   2150 bdr = bbd(jb3  )
   bdi = bbd(jb3+1)
   2160 SELECT CASE ( typea )
     CASE (    1)
       GO TO 2170
     CASE (    2)
       GO TO 2190
     CASE (    3)
       GO TO 2210
     CASE (    4)
       GO TO 2230
   END SELECT
   2170 DO  i = i1,in,2
     ddd(i  ) = ddd(i  ) + aas(k)*bdr
     ddd(i+1) = ddd(i+1) + aas(k)*bdi
     k = k + inca
   END DO
   GO TO 2250
   2190 DO  i = i1,in,2
     ddd(i  ) = ddd(i  ) + aad(k)*bdr
     ddd(i+1) = ddd(i+1) + aad(k)*bdi
     k = k + inca
   END DO
   GO TO 2250
   2210 DO  i = i1,in,2
     ddd(i  ) = ddd(i  ) + aas(k)*bdr - aas(k+1)*bdi
     ddd(i+1) = ddd(i+1) + aas(k)*bdi + aas(k+1)*bdr
     k = k + inca
   END DO
   GO TO 2250
   2230 DO  i = i1,in,2
     ddd(i  ) = ddd(i  ) + aad(k)*bdr - aad(k+1)*bdi
     ddd(i+1) = ddd(i+1) + aad(k)*bdi + aad(k+1)*bdr
     k = k + inca
   END DO
   2250 j3  = j3 + rca
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 2110
 
!     COMMON CASE (TYPEA=TYPEB=TYPED=PREC=4)
 
 2270 inca= rca*fa3
 i1  = 2*arow1 - 1
 in  = 2*arown - 1
 2280 j3 = 2*j3 - 1
 jb3n = jb31 + rcb*nterm3 - rcb
 DO  jb3 = jb31,jb3n,rcb
   k = j3
   DO  i = i1,in,2
     ddd(i  ) = ddd(i  ) + aad(k)*bbd(jb3  ) - aad(k+1)*bbd(jb3+1)
     ddd(i+1) = ddd(i+1) + aad(k)*bbd(jb3+1) + aad(k+1)*bbd(jb3  )
     k = k + inca
   END DO
   j3  = j3 + rca
 END DO
 CALL endget (BLOCK)
 CALL getstr (*2500,BLOCK)
 GO TO 2280
 
 2400 RETURN 1
 2500 RETURN
END SUBROUTINE mpyq
