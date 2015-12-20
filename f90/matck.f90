SUBROUTINE matck (mfile,pfile,a,z)
     
!     THIS ROUTINE CHECKS THE UNIQUENESS OF MATERIAL ID'S FOR
!          1. MAT1 (1)     8. MATT1  (MB)  15. MATS1  (MC)
!          2. MAT2         9. MATT2        16. MATPZ1 (MD)
!          3. MAT3        10. MATT3        17. MTTPZ1
!          4. MAT4        11. MATT4        18. MATPZ2
!          5. MAT5        12. MATT5        19. MTTPZ2 (ME)
!          6. MAT6        13. MATT6        20. DUMC
!          7. MAT8 (MA)   14. DUMB         21. DUMD  (NMAT)
!     AND THE MATERIAL ID SPECIFIED ON THE PROPERTY CARDS.
 
!     THIS ROUTINE SHOULD BE CALLED ONLY ONCE BY IFP.
!     THIS ROUTINE DOES NOT OPEN OR CLOSE MATERIAL FILE (MFILE) OR
!     ELEMENT PROPERTY FILE (PFILE)
 
!     WRITTEN BY G.CHAN/UNISYS,  OCT. 1982
 
 
 INTEGER, INTENT(IN OUT)                  :: mfile
 INTEGER, INTENT(IN OUT)                  :: pfile
 INTEGER, INTENT(IN)                      :: a(1)
 INTEGER, INTENT(IN)                      :: z(1)
 LOGICAL :: abort
 INTEGER :: ih(3),    NAME(2), mati(2,22)  &
     ,               group, epti(2,40),         matj(2,22)
 COMMON /system/ n1,       nout,     abort,    skip(42), kdum(9)
 DATA    matj  / 103,-12,  203,-17,  1403,-16, 2103,-3,  2203,-8,  &
     2503,-31,  603,-18, 703,-11,  803,-16,  1503,-16, 2303,-2,  2403,-7,  &
     2603,-31,  -11,-00, 503,-11, 1603,-07,  1803,-07, 1703,-44, 1903,-44,  &
     -11,-00,  -11,-00,   -11,-00/
 DATA    epti  /  52,191, 2502,071,  7002,071, 0502,041, 2202,041,  &
     5302,041, 0602,082,  0702,103, 0802,041, 0902,061,  &
     1002,041, 2102,041,  7052,171, 1102,082, 1202,103,  &
     1302,041, 7032,171,  1402,041, 1502,082, 1602,051,  &
     1702,041, 2002,031,  0152,243, 5102,241, 5802,174,  &
     5502,-49, 5602,-06,  5702,-06, 6102,001, 6202,001,  &
     6302,001, 6402,001,  6502,001, 6602,001, 6702,001,  &
     6802,001, 6902,001,     0,  0,    0,  0,    0,  0/
 DATA    nmat  / 21/,      group/ 7/
 DATA    nept  / 37/
 DATA    NAME  / 4HMATC,   4HK     /
 
!     FIRST WORDS ON THE EPTI TABLE ARE PROPERTY CARDS THAT SPECIFY
!     MATERIAL.  THE FIRST 2 DIGITS OF THE SECOND WORD INDICATE THE
!     NUMBER OF WORDS IN EACH PROPERTY INPUT CARD. AND THE 3RD DIGIT
!     INDICATES NUMBER OF MATERIAL ID'S SPECIFIED.
!     IF THIS SECOND WORD IS NEGATIVE, IT MEANS THE PROPERTY CARD IS
!     OPEN-ENDED. THE 3RD DIGIT INDICATES WHERE MID1 BEGINS, AND
!     REPEATING (FOR MID2, MID3,...) EVERY N WORDS WHERE N IS THE
!     ABSOLUTE VALUE OF THE FIRST 2 DIGITS. (NO REPEAT OF N=0)
 
!     ARRAY A CONTAINS A LIST OF ACTIVE PROPERTY IDS - SET UP BY PIDCK
 
 IF (abort) GO TO 220
 nomat = z(1)
 IF (nomat == 0) GO TO 145
 
!     UPDATE EPTI ARRAY IF DUMMY ELEMENT IS PRESENT
 
 DO  j = 1,9
   IF (kdum(j) == 0) CYCLE
   k = MOD(kdum(j),1000)/10
   epti(2,28+j) = k*10 + 1
 END DO
 
!     SET UP POINTERS FOR THE MATI TABLE
 
 ma = group
 mb = ma + 1
 mc = mb + group
 md = mc + 1
 me = mc + 4
 
!     READ MATERIAL ID INTO Z SPACE, AND SAVE APPROP. COUNT IN MATI(2,K)
 
 DO  j = 1,nmat
   mati(1,j) = matj(1,j)
   mati(2,j) = matj(2,j)
 END DO
 j = 1
 20   CALL fwdrec (*50,mfile)
 25   CALL READ (*50,*50,mfile,ih(1),3,0,kk)
 DO  k = 1,nmat
   IF (ih(1) == mati(1,k)) GO TO 35
 END DO
 GO TO 20
 35   nwds =-mati(2,k)
 IF (nwds < 0) CALL mesage (-37,0,NAME)
 mati(2,k) = 0
 40   CALL READ (*50,*25,mfile,z(j),nwds,0,kk)
 j = j + 1
 mati(2,k) = mati(2,k) + 1
 GO TO 40
 
!     INSTALL INITIAL COUNTERS IN MATI(1,K)
 
 50   jx = j
 IF (jx <= 1) GO TO 140
 mati(1,1) = 0
 DO  j = 1,nmat
   k = j + 1
   IF (mati(2,j) < 0) mati(2,j) = 0
   mati(1,k) = mati(1,j) + mati(2,j)
 END DO
 
!     NOTE - ORIGINAL DATA IN MATI TABLE IS NOW DESTROYED
 
!     CHECK MATERIAL ID UNIQUENESS AMONG MAT1, MAT2,..., MAT8
!     (MAT4 AND MAT5 ARE UNIQUE ONLY AMONG THEMSELVES)
 
 j = 0
 DO  k = 1,ma
   IF (mati(2,k) > 0) j = j + 1
 END DO
 IF (j <= 1) GO TO 90
 kk = mati(1,mb)
 k1 = kk - 1
 k4 = mati(1,4)
 loop80:  DO  k = 1,k1
   j  = z(k)
   ib = k + 1
   DO  i = ib,k1
     IF (j /= z(i)) CYCLE
     IF (k < k4 .AND. i >= k4) CYCLE
     CALL mesage (30,213,j)
     abort =.true.
     CYCLE loop80
   END DO
 END DO loop80
 
!     CHECK MATT1, MATT2,..., MATT6 AND MATS1 MATERIAL ID
!     AND THEIR CROSS REFERENCE TO MATI CARDS
 
 90   DO  k = mb,mc
   IF (mati(2,k) <= 0) CYCLE
   kk = MOD(k,ma)
   ib = mati(1,kk) + 1
   ie = mati(2,kk) + ib - 1
   jb = mati(1,k ) + 1
   je = mati(2,k ) + jb - 1
   loop105:  DO  j = jb,je
     k1 = z(j)
     IF (ie < ib) GO TO 100
     DO  i = ib,ie
       IF (z(i) == k1) CYCLE loop105
     END DO
     100  ih(1) = k1
     ih(2) = kk
     k1 = 217
     IF (k == 15) k1 = 17
     CALL mesage (30,k1,ih)
     abort =.true.
   END DO loop105
 END DO
 
!     CHECK MATERIAL ID UNIQUENESS AMONG MATPZI AND MTTPZI
 
 j = 0
 DO  k = md,me
   IF (mati(2,k) > 0) j = j + 1
 END DO
 IF (j <= 1) GO TO 140
 kk = mati(1,me+1)
 k1 = kk - 1
 nn = mati(1,md)
 loop130:  DO  k = nn,k1
   j  = z(k)
   ib = k + 1
   DO  i = ib,kk
     IF (j /= z(i)) CYCLE
     CALL mesage (30,213,j)
     abort =.true.
     CYCLE loop130
   END DO
 END DO loop130
 
!     NOW, WE CONTINUE TO CHECK MATERIAL ID'S ON MOST PROPERTY CARDS.
!     (MATERIAL ID'S ARE ON THE 2ND, 4TH, AND 6TH POSITIONS OF THE
!     PROPERTY CARDS, EXECPT THE OPEN-ENDED PCOMPI GROUP)
 
 140  je = mati(1,nmat)
 ii = a(1)
 145  CALL fwdrec (*220,pfile)
 150  CALL READ (*220,*220,pfile,ih(1),3,0,kk)
 DO  k = 1,nept
   IF (ih(1) == epti(1,k)) GO TO 170
 END DO
 GO TO 145
 170  IF (nomat == 0) GO TO 230
 nwds= epti(2,k)/10
 nn  = epti(2,k) - nwds*10
 ib  = 1
 ie  = nn*2
 ic  = 2
 komp= 0
 
!     CHANGE NWDS, IB, IE, AND IC IF PROPERTY CARD IS OPEN-ENDED
!     WHERE (IB+JX) POINTS TO THE FIRST MID POSITION
 
 IF (epti(2,k) > 0) GO TO 180
 komp= 1
 ib  =-nn - 1
 ic  =-nwds
 IF (nwds == 0) ic = 9999
 nwds= 10
 180  IF (komp == 1) ie = jx + nwds - 1
 
!     READ IN PROPERTY CARD. IF ID IS NOT ON ACTIVE LIST, SKIP IT.
!     SKIP IT TOO IF IT HAS NO MATERIAL ID REQUESTED.
!     (NO CORE SIZE CHECK HERE. SHOULD HAVE PLENTY AVAILABLE)
 
 CALL READ (*220,*150,pfile,z(jx),nwds,0,kk)
 IF (komp == 0) GO TO 182
 181  ie = ie + 1
 CALL READ (*220,*150,pfile,z(ie),1,0,kk)
 IF (z(ie) /= -1) GO TO 181
 ie = ie - 1 - jx
 182  DO  i = 2,ii
   IF (z(jx) == a(i)) GO TO 185
 END DO
 GO TO 180
 185 loop210: DO  i = ib,ie,ic
   kk = z(jx+i)
   IF (ie == 8 .AND. i == 7) kk = z(jx+i+3)
   IF (kk == 0) CYCLE loop210
   IF (jx <= 1) GO TO 200
   DO  j = 1,je
     IF (kk == z(j)) CYCLE loop210
   END DO
   200  ih(1) = kk
   ih(2) = z(jx)
   CALL mesage (30,215,ih)
   abort =.true.
 END DO loop210
 GO TO 180
 220  RETURN
 
 230  CALL mesage (30,16,ih)
 abort =.true.
 RETURN
END SUBROUTINE matck
