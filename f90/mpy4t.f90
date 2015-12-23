SUBROUTINE mpy4t (iz,z,dz)
     
!     INNER LOOP FOR MPYAD, METHOD 4 WITH TRANSPOSE
 
!          T
!         A * B + C = D
 
!     THIS ROUTINE IS CALLED ONLY BY MPYAD WHEN METHOD 2 TRANSPOSE,
!     MPY2T, IS SELECTED, AND DIAG 41 IS NOT TURNED ON BY USER.
 
!     MPY4T IS ABOUT 5 TIMES FASTER THAN MPY2T AS TESTED ON VAX
 
!     THERE IS A PICTORIAL DISCRIPTION ABOUT MPY4T IN MPYAD SUBROUTINE
 
!     THIS MACHINE INDEPENDENT ROUTINE CAN ACTUALLY BE INCORPORATED
!     INTO MPYQ, WHICH IS PRESENTLY A .MDS ROUTINE
 
!     IF MATRIX A, OR B, OR BOTH,  IS COMPLEX, MATRIX D IS COMPLEX.
!     MATRIX D CAN NOT BE COMPLEX, IF BOTH MATRICES A AND B ARE REAL.
 
 
!     WRITTEN BY G.CHAN/UNISYS   1/92
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN)                      :: iz(1)
 REAL, INTENT(IN OUT)                     :: z(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dz(1)
 REAL :: sumr    ,sumi
 INTEGER :: all
 DOUBLE PRECISION :: dsumr   ,dsumi   ,dzero
 DIMENSION  nam(2)
 COMMON /mpyadx/  filea(7),fileb(7),filec(7),filed(7)
 COMMON /TYPE  /  prc(2)  ,nwds(4) ,rc(4)
 COMMON /machin/  mach    ,ihalf   ,jhalf
 COMMON /unpakx/  typ     ,ii      ,jj
 COMMON /mpyadz/  rcb     ,rcd     ,ll      ,lll     ,jbb   ,  &
     nbx(3)  ,arow    ,arow1   ,arown   ,acore ,  &
     apoint  ,bcol    ,crow    ,firstl  ,na(3) , nwda
 COMMON /mpyqt4/  rca     ,prca    ,all     ,jump    ,prcd
 EQUIVALENCE      (dsumr  ,sumr  ) ,(dsumi  ,sumi)
 DATA    nam   /  4HMPY4  ,1HT   / ,dzero   / 0.0D+0 /
 
!*****
!     ANDF(I,J)   = IAND(I,J)
!     RSHIFT(I,J) = ISHFT(I,-J)
!     WHERE         ISHFT(I,-J) IS RIGHT-SHIFT I BY J BITS, ZERO FILL
!     AND           ISHFT IS SYSTEM ROUTINE
 
! UNIX:
!     REMOVE ABOVE 2 ON-LINE FUNCTIONS IF IAND AND ISHFT SYSTEM
!     FUNCTIONS ARE NOT AVAILABLE. ANDF AND RSHIFT ARE ALREADY ENTRY
!     POINTS IN SUBROUTINE MAPFNS.
!*****
 
!     METHOD 4T TRANSPOSE CASE
 
!     ARRAY Z(JBB) THRU Z(ACORE-1) HOLDS THE CURRENT COLUMN OF MATRIX B
!     ARRAY Z(1) THRU Z(JBB-1) IS A WORKING COLUMN SPACE FOR MATRIX D
 
!     ON EACH ROW OF A, WE WANT TO MULTIPLY
 
!        A(ROW,J)*B(J,COL) + C(ROW,COL) = D(ROW,COL)
 
!     NOTICE B(J,COL) RUNS FROM B(II,COL) THRU B(JJ,COL) WITHOUT
!     SKIPPING,
!     WHILE A(ROW,J) RUNS IN MULTIPLE STRING SEGMENTS ALONG J.
!     ALSO THE BEGINING OF J IN A(ROW,J) AND THE BEGINING OF J IN
!     B(J,COL) MOST LIKELY START DIFFERNTLY
 
!     NOW, ON EACH ROW, WE START FROM FIRST STRING. SKIP THIS STRING
!     IF IT IS NOT WITHIN B(II,) AND B(JJ,) RANGE. (ALSO, WE HAVE
!     SAVED PREVIOUSLY THE LAST TERM OF THE LAST STRING, AND THEREFORE
!     IF THE WHOLE ROW OF A(,J) WITH ITS STRINGS IS NOT WITHIN II,JJ
!     RANGE OF COLUMN B, WE SKIP THE WHOLE ROW-AND-COLUMN COMPUTATION.)
!     IF IT IS WITHIN THE RANGE, WE NEED TO SYNCHRONIZE THE J INDEX FOR
!     BOTH A(ROW,J) AND B(J,COL), THEN MULTIPLY, AND SUM ON AN ELEMENT
!     OF MATRIX D. THEN MOVE ON TO THE NEXT STRING, AND DO THE SAME.
!     REPEAT THIS PROCESS UNTIL J IS EXHAUST EITHER ON A(ROW,J) OR ON
!     B(J,COL).
!     WHEN ALL ROWS OF MATRIX A CURRENTLY IN CORE HAVE PASSED THRU, WE
!     HAVE ONE COLUMN OF MATRIX D DONE, FROM AROW1 THRU AROWN.
 
!     SINCE TRANSPOSE OF MATRIX A IS WHAT WE WANT, THE TERM 'ROW' IS
!     ACTUALLY 'COLUMN' WHEN THE DATA WAS MOVED INTO Z SPACE IN MPYAD
 
!     RCA,RCB    = 1, MATRIX A,B  IS REAL, = 2 MATRIX A,B IS COMPLEX
!     PRCA       = 1, MATRIX A IS IN S.P., = 2 MATRIX A IS IN D.P.
!     PRCD       = 0, MATRIX D IS IN S.P., = 1 MATRIX A IS IN D.P.
!     NWDA       = NUMBER OF WORDS PER ELEMENT OF MATRIX A
!     JBB        = POINTER TO FIRST WORD OF COLUMN B
!     II,JJ      = FIRST TO LAST NON-ZERO TERMS IN CURRENT COLUMN OF B
!     ALL        = 1,2,3,4 ALL MATRICES ARE OF THE SAME TYPE - S.P.,
!                  D.P., C.S.P., OR C.D.P. RESPECTIVELY
!                = 5, MATRICES ARE OF MIXED TYPES
!     JUMP       = BRANCHING INDEX TO MIXED TYPE MATRICES COMPUTATION.
 
!     APOINT     = POINTER TO STRING CONTROL WORD
!                = 0, CURRENT ROW OF A IS EXHAULTED
!     IZ(APOINT) = LEFT HALF OF WORD IS NBR, RIGHT HALF IS NBRSTR
!     NBR        = NO. OF WORDS   IN THIS STRING
!     NBRSTR     = NO. OF STRINGS IN THIS ROW A
!     INIT       = COLUMN POSITION OF 1ST STRING WORD
!     IF (INIT .GT. JJ) = 1ST STRING WORD IS BEYOND LAST WORD IN COLN B
!     IF (INIT+NBR .LT. II) = LAST STRING WORD IS BEFORE 1ST WORD IN
!                  COLUMN OF B
!     JB,JE      = BEGINNING AND ENDING J-INDEX FOR COLUMN A AND ROW B
!     IPOINT     = THE JB WORD POSITION IN ROW A
!     JA         = POINTER TO ROW A ELEMENT
!     KB         = POINTER TO COLUMN B ELEMENT
!     LAST       = POSITION OF LAST NON-ZERO COLUMN TERM IN ROW OF A
 
 
!     WE START FROM FIRST ROW AROW1, AND WILL RUN THRU TO LAST ROW AROWN
 
 arow = arow1
 l    = firstl
 10 apoint = iz(l)
 IF (apoint == 0) GO TO 510
 last = rshift(iz(l-1),ihalf)
 init = andf(iz(apoint),jhalf)
 IF (init > jj .OR. last < ii) GO TO 510
 nbrstr = andf(iz(l-1),jhalf)
 GO TO 30
 20 init = andf(iz(apoint),jhalf)
 30 nbr  = rshift(iz(apoint),ihalf)
 IF (init > jj) GO TO 510
 IF (init+nbr < ii) GO TO 500
 jb   = MAX0(init,ii)
 je   = MIN0(init+nbr-1,jj)
 IF (jb > je) GO TO 500
 ipoint = apoint + (jb-init+1)*prca
 ja   = (ipoint-1)/prca + 1
 kb   = (jb-ii)*rcb + jbb
 dsumr= dzero
 SELECT CASE ( all )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 80
   CASE (    4)
     GO TO 100
   CASE (    5)
     GO TO 120
   CASE (    6)
     GO TO 520
 END SELECT
 
 40 DO  j = jb,je
   sumr = sumr + z(ja)*z(kb)
   
!     DON'T BE SUPRISED TO SEE SOME Z(JA) ARE ZEROS
!     (VAX PACKING ROUTINE ALLOWS UP TO 3 ZEROS BETWEEN STRINGS)
   
   ja   = ja + rca
   kb   = kb + rcb
 END DO
 z(arow) = z(arow) + sumr
 GO TO 500
 
 60 DO  j = jb,je
   dsumr = dsumr + dz(ja)*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 dz(arow) = dz(arow) + dsumr
 GO TO 500
 
 80 sumi = 0.0
 DO  j = jb,je
   sumr = sumr + z(ja)*z(kb  ) - z(ja+1)*z(kb+1)
   sumi = sumi + z(ja)*z(kb+1) + z(ja+1)*z(kb  )
   ja   = ja + rca
   kb   = kb + rcb
 END DO
 z(arow  ) = z(arow  ) + sumr
 z(arow+1) = z(arow+1) + sumi
 GO TO 500
 
 100 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja)*dz(kb  ) - dz(ja+1)*dz(kb+1)
   dsumi = dsumi + dz(ja)*dz(kb+1) + dz(ja+1)*dz(kb  )
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 dz(arow  ) = dz(arow  ) + dsumr
 dz(arow+1) = dz(arow+1) + dsumi
 GO TO 500
 
 
 120 SELECT CASE ( jump )
   CASE (    1)
     GO TO 130
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 170
   CASE (    4)
     GO TO 190
   CASE (    5)
     GO TO  210
   CASE (    6)
     GO TO 230
   CASE (    7)
     GO TO 250
   CASE (    8)
     GO TO 270
   CASE (    9)
     GO TO 290
   CASE (   10)
     GO TO 310
   CASE (   11)
     GO TO 330
   CASE (   12)
     GO TO 350
   CASE (   13)
     GO TO  370
   CASE (   14)
     GO TO 390
   CASE (   15)
     GO TO 410
   CASE (   16)
     GO TO 430
 END SELECT
 
!                      +--------------- MATRIX  B -----------------+
!        MATRIX          REAL        REAL       COMPLEX     COMPLEX
!          A            SINGLE      DOUBLE      SINGLE      DOUBLE
!     ---------------  ----------  ---------  ----------  ----------
!     REAL SINGLE         130         150         170         190
!     REAL DOUBLE         210         230         250         270
!     COMPLEX SINGLE      290         310         330         350
!     COMPLEX DOUBLE      370         390         410         430
 
 
 130 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja)*z(kb))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 GO TO 460
 
 150 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja))*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   460
 ELSE
   GO TO   470
 END IF
 
 170 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja)*z(kb  ))
   dsumi = dsumi + DBLE(z(ja)*z(kb+1))
   ja   = ja + rca
   kb   = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 190 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja))*dz(kb  )
   dsumi = dsumi + DBLE(z(ja))*dz(kb+1)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 210 DO  j = jb,je
   dsumr = dsumr + dz(ja)*DBLE(z(kb))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   460
 ELSE
   GO TO   470
 END IF
 
 230 DO  j = jb,je
   dsumr = dsumr + dz(ja)*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 GO TO 470
 
 250 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja)*DBLE(dz(kb  ))
   dsumi = dsumi + dz(ja)*DBLE(dz(kb+1))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 270 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja)*dz(kb  )
   dsumi = dsumi + dz(ja)*dz(kb+1)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 290 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja  )*z(kb))
   dsumi = dsumi + DBLE(z(ja+1)*z(kb))
   ja   = ja + rca
   kb   = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 310 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja  ))*dz(kb)
   dsumi = dsumi + DBLE(z(ja+1))*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 330 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja)*z(kb  )) - DBLE(z(ja+1)*z(kb+1))
   dsumi = dsumi + DBLE(z(ja)*z(kb+1)) + DBLE(z(ja+1)*z(kb  ))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 GO TO 480
 
 350 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + DBLE(z(ja  ))*dz(kb)
   dsumi = dsumi + DBLE(z(ja+1))*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 370 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja  )*DBLE(z(kb))
   dsumi = dsumi + dz(ja+1)*DBLE(z(kb))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 390 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja  )*dz(kb)
   dsumi = dsumi + dz(ja+1)*dz(kb)
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 410 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja)*DBLE(z(kb  )) - dz(ja+1)*DBLE(z(kb+1))
   dsumi = dsumi + dz(ja)*DBLE(z(kb+1)) + dz(ja+1)*DBLE(z(kb  ))
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 IF (prcd > 0.0) THEN
   GO TO   480
 ELSE
   GO TO   490
 END IF
 
 430 dsumi = dzero
 DO  j = jb,je
   dsumr = dsumr + dz(ja)*dz(kb  ) - dz(ja+1)*dz(kb+1)
   dsumi = dsumi + dz(ja)*dz(kb+1) + dz(ja+1)*dz(kb  )
   ja    = ja + rca
   kb    = kb + rcb
 END DO
 GO TO 490
 
 460 dz(arow) = dz(arow) + dsumr
 GO TO 500
 470 z(arow)  = z(arow)  + SNGL(dsumr)
 GO TO 500
 480 dz(arow  ) = dz(arow  ) + dsumr
 dz(arow+1) = dz(arow+1) + dsumi
 GO TO 500
 490 z(arow  ) = z(arow  ) + SNGL(dsumr)
 z(arow+1) = z(arow+1) + SNGL(dsumi)
 
 
!     END OF STRING DATA. IF THIS IS NOT THE LAST STRING OF CURRENT
!     ROW OF A, RETURN FOR NEXT STRING
 
 500 nbrstr = nbrstr - 1
 apoint = apoint + nbr*nwda + prca
 IF (nbrstr > 0) GO TO 20
 
!     END OF A ROW OF MATRIX A.
!     RETURN FOR NEXT ROW IF THIS IS NOT THE LAST ROW IN OPEN CORE.
!     IF IT IS THE LAST ROW, RETURN TO CALLER FOR PACKING OUT THE
!     CURRENT COLUMN OF MATRIX D (IN C ARRAY)
 
 510 l    = l - 2
 arow = arow + 1
 IF (arow <= arown) GO TO 10
 RETURN
 
 520 CALL mesage (-37,0,nam)
 
 RETURN
END SUBROUTINE mpy4t
