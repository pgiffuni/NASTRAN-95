SUBROUTINE anisop
     
    !     COMPUTES DIRECTION COSINES FOR RECTANGULAR COORD. SYSTEMS
    !     (W.R.T. BASIC COORD. SYSTEM) DESCRIBING ORIENTATION OF ANIS.
    !     MATERIAL FOR ISOPARAMETRIC SOLIDS
 
    !     ANISOP  GEOM1,EPT,BGPDT,EQEXIN,MPT/MPTA/S,N,ISOP $
    !     EQUIV   MPTA,MPT/ISOP $
    !     ISOP=-1 MEANS SUCH MATERIALS EXIST
 
    INTEGER :: geom1,ept,bgpdt,FILE,buf1,buf2,eqexin
    DIMENSION       iz(1),nam(2),ic1(2),ic2(2),ipi(2),imat6(2),  &
        idum(31),a(3),b(3),c(3),xp(3),yp(3),zp(3), store(9),xd(9),itrl(7),mat1(2)
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm
    COMMON /BLANK / isop
    COMMON /system/ ibuf,nout
    COMMON /zzzzzz/ z(1)
    EQUIVALENCE     (z(1),iz(1))
    DATA    geom1 , ept,bgpdt,eqexin,mpt,mpta /  &
        101   , 102,103  ,104   ,105,201  /
    DATA    ipi   / 7002,70/ ,ic1 / 1801,18   /, ic2 / 2101,21 /,  &
        imat6 / 2503,25/ ,mat1/  103, 1   /
    DATA    nam   / 4HANIS, 4HOP  /
 
    isop  = 1
    lcore = korsz(z)
    buf1  = lcore - ibuf - 1
    buf2  = buf1  - ibuf
    lcore = buf2 - 1
    IF (lcore <= 0) GO TO 1008
 
    !     GET LIST OF MAT1 AND MAT6 ID-S
 
    nmat1 = 0
    nmat6 = 0
    i     = 0
    !WKBI SPR93033 5/94
    FILE  = mpt
    CALL preloc (*1001,z(buf1),mpt)
    CALL locate (*2,z(buf1),mat1,idx)
    !WKBD SPR93022 5/94      FILE  = MPT
1   CALL READ (*1002,*2,mpt,idum,12,0,m)
    nmat1 = nmat1 + 1
    i = i + 1
    IF (i > lcore) GO TO 1008
    iz(i) = idum(1)
    GO TO 1
 
    !     DONE WITH MAT1
 
2   CALL locate (*4,z(buf1),imat6,idx)
3   CALL READ (*1002,*4,mpt,idum,31,0,m)
    nmat6 = nmat6 + 1
    i = i + 1
    IF (i > lcore) GO TO 1008
    iz(i) = idum(1)
    GO TO 3
4   CALL CLOSE (mpt,1)
 
    !     LOCATE PIHEX CARDS ON EPT AND FORM A LIST OF MATERIAL COORD.
    !     SYSTEM ID
 
    !WKBI SPR93033 5/94
    FILE = ept
    CALL preloc (*310,z(buf1),ept)
    CALL locate (*310,z(buf1),ipi,idx)
    !WKBD SPR93033 5/94      FILE = EPT
10  CALL READ (*1002,*40,ept,idum,7,0,m)
 
    icid = idum(3)
    mid  = idum(2)
 
    !     IF CID = 0, MID MUST BE MAT1
    !     IF CID IS NOT 0, MID MUST BE MAT6
 
    IF (icid  > 0) GO TO 11
    IF (nmat1 == 0) GO TO 14
    ii = 1
    nn = nmat1
    GO TO 12
11  IF (nmat6 == 0) GO TO 14
    ii = nmat1 + 1
    nn = nmat1 + nmat6
    12 DO  iii = ii,nn
        IF (mid == iz(iii)) GO TO 1141
    END DO
14  WRITE  (nout,1140) ufm,mid
1140 FORMAT (a23,', MATERIAL',i8,', SPECIFIED ON A PIHEX CARD, DOES ',  &
        'NOT REFERENCE THE PROPER MATERIAL TYPE', /5X,  &
        'CID = 0 MEANS MAT1, CID NOT 0 MEANS MAT6')
    GO TO 115
1141 CONTINUE
 
     !     STORE ALL CID,MID PAIRS WHERE CID IS NOT 0
 
     IF (icid == 0) GO TO 10
     i = i + 1
     iz(i) = icid
     i = i + 1
     iz(i) = mid
     GO TO 10
 
     !     LIST IS MADE. MOVE IT UP TO IZ(1)
 
40   CALL CLOSE (ept,1)
     ncid = i - nmat1 - nmat6
     IF (ncid == 0) RETURN
 
     DO  ii = 1,ncid
         jj = nmat1 + nmat6 + ii
         iz(ii) = iz(jj)
     END DO
 
     !     NOW MAKE A UNIQUE LIST OF CID-S
 
     ijk = ncid + 1
     num = 1
     iz(ijk) = iz(1)
     IF (ncid == 2) GO TO 44
     loop43:  DO  ii = 3,ncid,2
         icid = iz(ii)
         DO  jj = 1,num
             ncjj = ncid + jj
             IF (icid == iz(ncjj)) CYCLE loop43
         END DO
   
         !     UNIQUE - LIST IT
   
         ijk = ijk + 1
         num = num + 1
         IF (ijk > lcore) GO TO 1008
         iz(ijk) = icid
     END DO loop43
 
     !     UNIQUE LIST IS MADE- CHECK AGAINST CORD1R AND CORD2R ID-S
 
44   icord = ncid + num + 1
 
     ncord1 = 0
     ncord2 = 0
     FILE   = geom1
     CALL preloc (*1001,z(buf1),geom1)
     CALL locate (*70,z(buf1),ic1,idx)
45   IF (icord+12 > lcore) GO TO 1008
     CALL READ (*1002,*70,geom1,z(icord),6,0,m)
 
     !     COMPARE AGAINST CIDS ON PIHEX-S
 
     DO  jj = 1,num
         j = ncid + jj
         IF (iz(icord) == iz(j)) GO TO 60
     END DO
     GO TO 45
 
     !     MATCH- RESERVE 13 WORDS SINCE THIS CORD1R WILL BE CONVERTED TO
     !     CORD2R TYPE ENTRY LATER
 
60   iz(j)  =-iz(j)
     ncord1 = ncord1 + 1
     icord  = icord  + 13
     IF (ncord1 == num) GO TO 120
     GO TO 45
 
     !     TRY CORD2R
 
70   CALL locate (*100,z(buf1),ic2,idx)
75   IF (icord+12 > lcore) GO TO 1008
     CALL READ (*1002,*100,geom1,z(icord),13,0,m)
 
     !     COMPARE
 
     DO  jj = 1,num
         j = ncid + jj
         IF (iz(icord) == iz(j)) GO TO 90
     END DO
     GO TO 75
 
     !     MATCH ON CORD2R. CHECK FOR RID. MUST BE 0
 
90   IF (iz(icord+3) /= 0) GO TO 330
 
     iz(j)  =-iz(j)
     ncord2 = ncord2 + 1
     icord  = icord  + 13
     IF (ncord1+ncord2 == num) GO TO 120
     GO TO 75
 
     !     EXHAUSTED CORD2R-S, BUT NOT ALL  CID-S ARE LOCATED
 
     100 DO  jj = 1,num
         j = ncid + jj
         IF (iz(j) < 0) CYCLE
         WRITE  (nout,105) ufm,iz(j)
105      FORMAT (a23,', CID',i8,' ON A PIHEX CARD IS NOT DEFINED TO BE ',  &
             'CORD1R OR CORD2R')
     END DO
115  CALL mesage (-61,0,nam)
 
 
     !     MATCHING IS COMPLETE
 
120  CALL CLOSE (geom1,1)
 
     !     CID,MATID PAIRS ARE IN Z(1)-Z(NCID). UNIQUE CID LIST IS IN
     !     Z(NCID+1)-Z(NCID+NUM). THERE ARE NCORD1 CORD1R-S AND NCORD2
     !     CORD2R-S AT 13 WORDS EACH STARTING AT Z(NCID+NUM+1).
     !     NEXT AVAILABLE OPEN CORE IS AT Z(ICORD)
 
     DO  jj = 1,num
         j = ncid + jj
         iz(j) =-iz(j)
     END DO
 
     !     FOR CID-S ON CORD1R WE MUST OBTAIN THE BASIC COORDINATES OF EACH
     !     POINT FROM BGPDT. FIRST, THE EXTERNAL POINT NUMBERS ON CORD1R MUST
     !     BE CONVERTED TO  INTERNAL.
 
     lcore = lcore - (icord-1)
     IF (lcore <= 0) GO TO 1008
     mcore  = lcore
     ibgpdt = icord
     IF (ncord1 == 0) GO TO 200
     CALL gopen (bgpdt,z(buf1),0)
     FILE = bgpdt
     CALL READ (*1002,*140,bgpdt,z(ibgpdt),lcore,0,m)
     GO TO 1008
140  CALL CLOSE (bgpdt,1)
     ieq   = ibgpdt + m
     lcore = lcore  - m
     CALL gopen (eqexin,z(buf1),0)
     FILE = eqexin
     CALL READ (*1002,*150,eqexin,z(ieq),lcore,0,m)
     GO TO 1008
150  CALL CLOSE (eqexin,1)
     lcore = lcore - m
 
     !     FOR EACH CORD1R ENTRY, FIND THE BASIC COORDINATES FOR EACH POINT
     !     AND FORM A CORD2R ENTRY BACK WHERE THE CORD1R IS STORED
 
     DO  j = 1,ncord1
         ipoint = 13*(j-1) + ncid + num
         icid   = iz(ipoint+1)
         DO  k = 1,3
             isubk = ipoint + 3 + k
             k3 = 3*(k-1)
             igrid = iz(isubk)
             CALL bisloc (*350,igrid,z(ieq),2,m/2,jp)
     
             !     IM IS POINTER TO INTERNAL NUMBER. NOW FIND BGPDT ENTRY
     
             im = ieq + jp
             ip = 4*(iz(im)-1)
             DO  l = 1,3
                 isubb = ibgpdt + ip + l
                 isubl = k3 + l
                 store(isubl) = z(isubb)
             END DO
         END DO
   
         !     WE HAVE THE BASIC COORDINATES OF THE 3 POINTS. STORE IT BACK INTO
         !     THE CORD1R ENTRY. THE ENTRY STARTS AT Z(IPOINT+1)
   
         ip4 = ipoint + 4
         iz(ip4) = 0
         DO  l = 1,9
             isubl = ip4 + l
             z(isubl) = store(l)
         END DO
   
     !     GO BACK FOR ANOTHER CORD1R
   
     END DO
 
     !     FOR EACH COORDINATE SYSTEM, COMPUTE THE 9 DIRECTION COSINES FROM
     !     THE BASIC COORDINATE SYSTEM. Z(ICORD) IS THE NEXT AVAILABLE
     !     LOCATION OF OPEN CORE SINCE WE NO LONGER NEED EQEXIN OR BGPDT
     !     INFO.
 
200  lcore = mcore
     CALL gopen (mpt,z(buf1),0)
     CALL gopen (mpta,z(buf2),1)
     IF (icord+30 > lcore) GO TO 1008
 
     !     COPY MPT TO MPTA UNTIL MAT6 IS REACHED
 
     FILE = mpt
210  CALL READ (*280,*1003,mpt,z(icord),3,0,m)
     CALL WRITE (mpta,z(icord),3,0)
     IF (iz(icord) == imat6(1)) GO TO 240
220  CALL READ (*1002,*230,mpt,z(icord),lcore,0,m)
     CALL WRITE (mpta,z(icord),lcore,0)
     GO TO 220
230  CALL WRITE (mpta,z(icord),m,1)
     GO TO 210
 
     !     MAT6 RECORD FOUND. EACH MAT6 CONTAINS 31 WORDS. INCREASE THAT
     !     TO 40
 
240  CALL READ (*1002,*270,mpt,z(icord),31,0,m)
 
     !     SEE IF THIS ID MATCHES A CID ON PIHEX) IT NEED NOT
 
     DO  j = 2,ncid,2
         IF (iz(j) == iz(icord)) GO TO 258
     END DO
 
     !     NO MATCH. MAT6 NOT REFERENCED BY PIHEX. COPY IT TO MAT6 AND FILL
     !     IN ALL 3 DIRECTION COSINES. THIS MAT6 IS NOT REFERENCED BY PIHEX
 
     DO  k = 1,9
         xd(k) = 0.
     END DO
     GO TO 265
 
     !     MATCH. NOW FIND IT IN CORD1R,CORD2R LIST
 
258  icid = iz(j-1)
     DO  ii = 1,num
         ipoint = ncid + num + 13*(ii-1)
         IF (icid == iz(ipoint+1)) GO TO 260
     END DO
 
     !     LOGIC ERROR
 
     GO TO 370
 
260  iz(ipoint+1) = -iz(ipoint+1)
     a(1) = z(ipoint+ 5)
     a(2) = z(ipoint+ 6)
     a(3) = z(ipoint+ 7)
     b(1) = z(ipoint+ 8)
     b(2) = z(ipoint+ 9)
     b(3) = z(ipoint+10)
     c(1) = z(ipoint+11)
     c(2) = z(ipoint+12)
     c(3) = z(ipoint+13)
 
     !     ZP AXIS IS B-A. YP IS ZP X (C-A). XP IS YP X ZP
 
     zp(1) = b(1) - a(1)
     zp(2) = b(2) - a(2)
     zp(3) = b(3) - a(3)
     store(1) = c(1) - a(1)
     store(2) = c(2) - a(2)
     store(3) = c(3) - a(3)
     yp(1) = zp(2)*store(3) - zp(3)*store(2)
     yp(2) = zp(3)*store(1) - zp(1)*store(3)
     yp(3) = zp(1)*store(2) - zp(2)*store(1)
     xp(1) = yp(2)*zp(3) - yp(3)*zp(2)
     xp(2) = yp(3)*zp(1) - yp(1)*zp(3)
     xp(3) = yp(1)*zp(2) - yp(2)*zp(1)
 
     !     NOW COMPUTE DIRECTION COSINES BETWEEN XP,YP, ZP AND BASIC X,Y,Z
     !     X=(1,0,0),Y=(0,1,0),Z=(0,0,1)
     !     COS(THETA)=(DP.D)/(LENGTH OF DP)*(LENGTH OF D) WHERE DP=XP,YP,OR
     !     ZP   AND D=X,Y,OR Z. LENGTH OF D=1
 
     dl    = SQRT(xp(1)**2 + xp(2)**2 + xp(3)**2)
     xd(1) = xp(1)/dl
     xd(2) = xp(2)/dl
     xd(3) = xp(3)/dl
     dl    = SQRT(yp(1)**2 + yp(2)**2 + yp(3)**2)
     xd(4) = yp(1)/dl
     xd(5) = yp(2)/dl
     xd(6) = yp(3)/dl
     dl    = SQRT(zp(1)**2 + zp(2)**2 + zp(3)**2)
     xd(7) = zp(1)/dl
     xd(8) = zp(2)/dl
     xd(9) = zp(3)/dl
 
     !     WRITE OUT NEW MAT6 RECORD WITH DIRECTION COSINES APPENDED
 
265  CALL WRITE (mpta,z(icord),31,0)
     CALL WRITE (mpta,xd,9,0)
 
     !     GET ANOTHER MAT6
 
     GO TO 240
 
     !     MAT6 RECORD FINISHED. WRITE EOR, COPY REMAINDER OF MPT, AND CHECK
     !     TO SEE THAT ALL PIHEX CID-S HAVE BEEN ACCOUNTED FOR.
 
270  CALL WRITE (mpta,0,0,1)
     GO TO 210
 
     !     MPT EXHAUSTED
 
280  CALL CLOSE (mpt,1)
     CALL CLOSE (mpta,1)
     itrl(1) = mpt
     CALL rdtrl (itrl)
     itrl(1) = mpta
     CALL wrttrl (itrl)
     isop = -1
285  RETURN
 
     !WKBDB 5/94 SPR 93033
     !  310 WRITE  (NOUT,320) UWM
     !  320 FORMAT (A25,', EITHER EPT IS PURGED OR NO PIHEX CARDS FOUND ON ',
     !     1       'EPT IN ANISOP')
     !WKBDE 5/94 SPR 93033
     !WKBI  SPR 93033 5/94
310  CALL CLOSE ( ept, 1 )
     GO TO 285
330  WRITE  (nout,340) ufm,iz(j)
340  FORMAT (a23,', CORD2R',i8,' DEFINES A PIHEX CID BUT HAS NONZERO', ' RID')
     GO TO 115
350  WRITE  (nout,360) ufm,igrid
360  FORMAT (a23,', EXTERNAL GRID',i8,' CANNOT BE FOUND ON EQEXIN IN ',  &
         'ANISOP')
     GO TO 115
370  WRITE  (nout,380) ufm
380  FORMAT (a23,', NON-UNIQUE COORDINATE SYSTEMS ON PIHEX CARDS', /5X,  &
         '(SEE USER MANUAL P.2.4-233(05/30/86))')
     GO TO 115
 
     !WKBDB SPR93033 5/94
     ! 1001 N = -1
     !      GO TO 1010
     !WKBDE SPR93033 5/94
     !WKBIB SPR93033 5/94
1001 CALL CLOSE ( mpt, 1 )
     CALL CLOSE ( geom1, 1 )
     GO TO 285
     !WKBIE SPR93033 5/94
1002 n = -2
     GO TO 1010
1003 n = -3
     GO TO 1010
1008 FILE = 0
     n = -8
1010 CALL mesage (n,FILE,nam)

     RETURN
END SUBROUTINE anisop
