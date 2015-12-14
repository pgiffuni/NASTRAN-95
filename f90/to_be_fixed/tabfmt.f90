SUBROUTINE tabfmt
     
!     MODULE MAIN PROGRAM FOR DMAP MODULE TABPRT
 
!     THE CALL TO THIS MODULE IS
 
!     TABPRT TDB // C,N,KEY / C,N,OPT1 / C,N,OPT2 $
 
!            TDB IS THE TABLE DATA BLOCK TO BE PRINTED.
 
!            KEY IS THE BCD VALUE WHICH DETERMINES THE FORMAT BY
!                 WHICH THE TABLE IS PRINTED.
!                 THERE IS NO DEFAULT VALUE FOR KEY.
 
!            OPT1 IS A SKIP FACTOR BETWEEN DATA LINES.
!                 OPT1.EQ.0 MEANS  NO SPACE BETWEEN DATA LINES.
!                 OPT1.NE.0 MEANS ONE SPACE BETWEEN DATA LINES.
!                 THE DEFAULT VALUE FOR OPT1 IS 0
 
!            OPT2 IS ZERO BY DEFAULT.
!                 SKIP FILE-NAME AND KEY CHECKING IF OPT2 IS NON-ZERO.
 
 INTEGER :: p,p2,p3,na,r,x(14),subnam(2),NAME(2),NONE(2),  &
     nam(2),re,f,t(7),wd,rl,y,z(2),eid,a,b, h1,h2,h3,hx,zero,one,two
 REAL :: rx(14)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / p(2),p2,p3
 COMMON /system/ nb,no,junk1(6),nlpp,junk2(2),line
 COMMON /output/ t1(32),t2(32),t3(32),h1(32),h2(32),h3(32)
 COMMON /tabftx/ la,na(2,21)  ,  hx(32,40)  , re(21)
 COMMON /zzzzzz/ ix(1)
 EQUIVALENCE     (rx(1),x(1),ix(1))
 DATA    NONE  / 4H (no , 4HNE)  /
 DATA    subnam/ 4HTABF , 4HMT   /
 DATA    f     / 101   /
 DATA    zero  / 4H  0 /, one    / 4H  1 / , two / 4H  2 /
 
 
 1 FORMAT (1H )
 
 lc = korsz(x) - nb
 ib = lc + 1
 IF (lc <= 0) CALL mesage (-8,lc,subnam)
 ls = 1
 IF (p2 /= 0) ls = 2
 
 DO  i = 1,la
   IF (p(1) == na(1,i) .AND. p(2) == na(2,i)) GO TO 200
 END DO
 GO TO 9901
 
 200 CONTINUE
 CALL fname (f,NAME)
 IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) GO TO 9902
 t(1) = f
 CALL rdtrl (t)
 IF (t(1) <= 0) GO TO 9902
 CALL OPEN (*9902,f,x(ib),0)
 CALL READ (*9903,*9904,f,nam,2,re(i),kf)
 IF (nam(1) == p(1) .AND. nam(2) == p(2)) GO TO 250
 IF (p3 == 0) GO TO 9901
 250 CONTINUE
 
 GO TO (1100,1200,1300,1400,1500,1600,1700,1800,1900,2000  &
     ,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000 ,3100), i
 
!     PRINT CONTENTS OF TABLE DATA BLOCK BGPDT.
 
 1100 CONTINUE
 m1 = 2
 m2 = 3
 m3 = 4
 ASSIGN 1110 TO r
 GO TO 8000
 1110 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = one
 IF (lc < 4) GO TO 9905
 j = 0
 1120 CALL READ (*9903,*1180,f,x,4,0,kf)
 j = j + 1
 line = line + ls
 IF (line <= nlpp) GO TO 1140
 CALL page
 WRITE (no,1)
 line = ls
 1140 IF (p2 /= 0) WRITE (no,1)
 WRITE  (no,1150) j,x(1),(rx(l),l=2,4)
 1150 FORMAT (20X,i10,i13,1X,1P,3E20.5)
 GO TO 1120
 1180 CONTINUE
 IF (kf == 0) GO TO 9000
 GO TO 9906
 
!     PRINT CONTENTS OF TABLE DATA BLOCK GPL.
 
 1200 CONTINUE
 
!     RECORD 1
 
 m1 = 2
 m2 = 5
 m3 = 6
 ASSIGN 1205 TO r
 GO TO 8000
 1205 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = one
 IF (lc < 5) GO TO 9905
 j = -4
 1210 CALL READ (*9903,*1230,f,x,5,0,kf)
 j = j + 5
 line = line + ls
 IF (line <= nlpp) GO TO 1220
 CALL page
 WRITE (no,1)
 line = ls
 1220 IF (p2 /= 0) WRITE (no,1)
 WRITE  (no,1225) j,(x(l),l=1,5)
 1225 FORMAT (11X,i8,5(9X,i8,3X))
 GO TO 1210
 1230 CONTINUE
 IF (kf == 0) GO TO 1250
 j = j + 5
 line = line + ls
 IF (line <= nlpp) GO TO 1240
 CALL page
 WRITE (no,1)
 line = ls
 1240 IF (p2 /= 0) WRITE (no,1)
 WRITE (no,1225) j,(x(l),l=1,kf)
 
!     RECORD 2
 
 1250 IF (i == 4) GO TO 9000
 m1 = 2
 m2 = 7
 m3 = 8
 ASSIGN 1255 TO r
 GO TO 8000
 1255 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = two
 IF (lc < 6) GO TO 9905
 j = -2
 1260 CALL READ (*9903,*1280,f,x,6,0,kf)
 j = j + 3
 line = line + ls
 IF (line <= nlpp) GO TO 1270
 CALL page
 WRITE (no,1)
 line = ls
 1270 IF (p2 /= 0) WRITE (no,1)
 WRITE  (no,1275) j,(x(l),l=1,6)
 1275 FORMAT (11X,i8,2X,3(9X,i8,4X,i12))
 GO TO 1260
 1280 CONTINUE
 IF (kf == 0) GO TO 1292
 j = j + 3
 line = line + ls
 IF (line <= nlpp) GO TO 1290
 CALL page
 WRITE (no,1)
 line = ls
 1290 IF (p2 /= 0) WRITE (no,1)
 WRITE (no,1275) j,(x(l),l=1,kf)
 
 1292 IF (MOD(kf,2) == 0) GO TO 9000
 GO TO 9906
 
!     PRINT CONTENTS OF TABLE DATA BLOCK CSTM
 
 1300 CONTINUE
 m1 = 2
 m2 = 9
 m3 = 10
 ASSIGN 1310 TO r
 GO TO 8000
 1310 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = one
 IF (lc < 14) GO TO 9905
 j = 0
 1320 CALL READ (*9903,*1380,f,x,14,0,kf)
 j = j + 1
 line = line + ls + 2
 IF (line <= nlpp) GO TO 1340
 CALL page
 WRITE (no,1)
 line = ls + 2
 1340 IF (p2 /= 0) WRITE (no,1)
 WRITE (no,1350) j,x(1),x(2),rx( 6),rx( 7),rx( 8),rx( 3)  &
     ,rx( 9),rx(10),rx(11),rx( 4) ,rx(12),rx(13),rx(14),rx( 5)
 1350 FORMAT( 10X,i10,i10,i10, 1P,3E20.8,10X,1P,e20.8  &
     /40X,             1P,3E20.8,10X,1P,e20.8  &
     /40X,             1P,3E20.8,10X,1P,e20.8)
 GO TO 1320
 1380 CONTINUE
 IF (kf == 0) GO TO 9000
 GO TO 9906
 
!     PRINT CONTENTS OF TABLE DATA BLOCK GPLD
 
 1400 CONTINUE
 GO TO 1200
 
!     PRINT CONTENTS OF TABLE DATA BLOCK EQEXIN
 
 1500 CONTINUE
 m1 = 2
 m2 = 11
 m3 = 12
 ASSIGN 1501 TO r
 GO TO 8000
 1501 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = one
 1502 IF (lc < 8) GO TO 9905
 j = -3
 1510 CALL READ (*9903,*1530,f,x,8,0,kf)
 j = j + 4
 line = line + ls
 IF (line <= nlpp) GO TO 1520
 CALL page
 WRITE (no,1)
 line = ls
 1520 IF (p2 /= 0) WRITE (no,1)
 WRITE  (no,1525) j,(x(l),l=1,8)
 1525 FORMAT (7X,i8,4(7X,i8,6X,i8))
 GO TO 1510
 1530 CONTINUE
 IF (kf == 0) GO TO 1550
 j = j + 4
 line = line + ls
 IF (line <= nlpp) GO TO 1540
 CALL page
 WRITE (no,1)
 line = ls
 1540 IF (p2 /= 0) WRITE (no,1)
 WRITE (no,1525) j,(x(l),l=1,kf)
 
!     RECORD 2
 
 1550 IF (h1(24) == two) GO TO 9000
 m1 = 2
 m2 = 11
 m3 = 22
 ASSIGN 1555 TO r
 GO TO 8000
 1555 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = two
 GO TO 1502
 
!     PRINT CONTENTS OF TABLE DATA BLOCK EQDYN
 
 1600 CONTINUE
 GO TO 1500
 
!     PRINT CONTENTS OF TABLE DATA BLOCK GPDT
 
 1700 CONTINUE
 m1 = 2
 m2 = 13
 m3 = 14
 ASSIGN 1710 TO r
 GO TO 8000
 1710 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = one
 IF (lc < 7) GO TO 9905
 j = 0
 1720 CALL READ (*9903,*1750,f,x,7,0,kf)
 j = j + 1
 line = line + ls
 IF (line <= nlpp) GO TO 1730
 CALL page
 WRITE (no,1)
 line = ls
 1730 IF (p2 /= 0) WRITE (no,1)
 WRITE  (no,1740) x(1),x(2),rx(3),rx(4),rx(5),x(6),x(7)
 1740 FORMAT (7X,i8,10X,i8,10X,3(1P,e12.5,5X),5X,i8,10X,i8)
 GO TO 1720
 1750 CONTINUE
 IF (kf == 0) GO TO 9000
 GO TO 9906
 
!     PRINT CONTENTS OF TABLE DATA BLOCK GPTT
 
 1800 CONTINUE
 
!     RECORD 0
 
 m1 = 2
 m2 = 15
 m3 = 16
 ASSIGN 1805 TO r
 GO TO 8000
 1805 CONTINUE
 h1(19) = p(1)
 h1(20) = p(2)
 h1(24) = zero
 IF ((lc/3)*3 == 0) GO TO 9905
 ival = (lc/3)*3
 CALL READ (*9903,*1812,f,x,ival,0,kf)
 GO TO 9905
 1812 wd = ((kf-1)/3) + 1
 IF (kf == 0) wd = ((lc-1)/3) + 1
 DO  j = 1,wd
   line = line + ls
   IF (line <= nlpp) GO TO 1815
   CALL page
   WRITE (no,1)
   line = ls
   1815 IF (p2 /= 0) WRITE (no,1)
   IF (x(3*j-1) == -1) GO TO 1822
   WRITE  (no,1821) j,x(3*j-2),rx(3*j-1),x(3*j)
   1821 FORMAT (7X,i8,10X,i8,14X,1P,e12.5,19X,i8)
   CYCLE
   1822 WRITE  (no,1823) j,x(3*j-2),x(3*j-1),x(3*j)
   1823 FORMAT (7X,i8,10X,i8,14X,6X,i3,22X,i8)
 END DO
 
!     RECORD 1 AND ALL OTHERS
 
 m1 = 17
 m2 = 1
 m3 = 18
 ASSIGN 1835 TO r
 GO TO 8000
 1835 CONTINUE
 DO  rl = 1,j
   IF (x(3*rl)   == 0) CYCLE
   IF ((lc-3*wd) < 4) GO TO 9905
   CALL READ (*9903,*9904,f,y,1,0,kf)
   1840 CALL READ (*9903,*1895,f,z,2,0,kf)
   1845 CALL READ (*9903,*9904,f,eid,1,0,kf)
   
!     ELEMENT ID EQUALS ZERO INDICATES THE END OF DATA FOR CURRENT TYPE
   
   IF ((lc-3*wd) < z(2)) GO TO 9905
   
!     ELEMENT ID LESS THAN ZERO INDICATES NONEXISTENT TEMPERATURE VALUES
   
   IF (eid < 0.0) THEN
     GO TO  1848
   ELSE IF (eid == 0.0) THEN
     GO TO  1840
   END IF
   1847 ival = 3*wd + 1
   CALL READ (*9903,*9904,f,x(ival),z(2),0,kf)
   1848 a = 3*wd + 1
   b = a + 7
   IF (b >= (a+z(2))) b = a + z(2) - 1
   line = line + ls
   IF (line <= nlpp) GO TO 1865
   CALL page
   WRITE  (no,1850) x(3*rl),y
   1850 FORMAT (9X,i8,11X,i8,11X,i8,13X,i8)
   WRITE  (no,1860)
   1860 FORMAT (14H0   element id,8X,5H( 1 ),9X,5H( 2 ),9X,5H( 3 ),9X,  &
       5H( 4 ),9X,5H( 5 ),9X,5H( 6 ),9X, 5H( 7 ),9X,5H( 8 ) )
   WRITE (no,1)
   line = ls + 3
   1865 IF (eid < 0.0) THEN
     GO TO  1866
   ELSE IF (eid == 0.0) THEN
     GO TO  1840
   ELSE
     GO TO  1867
   END IF
   1866 WRITE (no,1870) eid
   GO TO 1845
   1867 WRITE  (no,1870) eid,(rx(l),l=a,b)
   1870 FORMAT (4X,i8,3X,8(2X,1P,e12.5) )
   IF (p2 /= 0) WRITE (no,1)
   IF (b == (a+z(2)-1)) GO TO 1845
   a = a + 8
   1875 b = a + 7
   IF (b >= (a+z(2))) b = a + z(2) - 1
   line = line + ls
   IF (line <= nlpp) GO TO 1880
   CALL page
   WRITE (no,1850) x(3*rl),y
   WRITE (no,1)
   line = ls + 3
   1880 WRITE  (no,1885) (rx(l),l=a,b)
   1885 FORMAT (17X,1P,e12.5,2X,1P,e12.5,2X,1P,e12.5,2X,1P,e12.5,  &
       2X,1P,e12.5,2X,1P,e12.5,2X,1P,e12.5,2X,1P,e12.5)
   IF (p2 /= 0) WRITE (no,1)
   IF (b == (a+z(2)-1)) GO TO 1845
   a = a + 8
   GO TO 1875
 END DO
 IF (kf == 0) GO TO 9000
 GO TO 9906
 
!     PRINT CONTENTS OF TABLE DATA BLOCK GPCT
 
 1900 CONTINUE
 m1 = 19
 m2 = 20
 m3 = 21
 ASSIGN 1910 TO r
 GO TO 8000
 1910 CONTINUE
 IF (lc < 12) GO TO 9905
 j = 0
 1920 CALL READ (*1990,*9904,f,pi,1,0,kf)
 j  = j + 1
 wd = 10
 line = line + ls
 IF (line <= nlpp) GO TO 1930
 CALL page
 WRITE (no,1)
 line = ls
 1930 IF (p2 /= 0) WRITE (no,1)
 CALL READ (*9903,*1980,f,m,1,0,kf)
 CALL READ (*9903,*1940,f,x,10,0,kf)
 GO TO 1945
 1940 wd = kf
 1945 WRITE  (no,1950) j,pi,m,(x(l),l=1,wd)
 1950 FORMAT (3X,i8,12(2X,i8))
 IF (m <= 10) GO TO 1920
 1960 line = line + ls
 IF (line <= nlpp) GO TO 1965
 CALL page
 WRITE (no,1)
 line = ls
 1965 IF (p2 /= 0) WRITE (no,1)
 CALL READ (*9903,*1975,f,x,10,0,kf)
 WRITE  (no,1970) (x(l),l=1,wd)
 1970 FORMAT (31X,10(2X,i8))
 GO TO 1960
 1975 wd = kf
 WRITE (no,1970) (x(l),l=1,wd)
 GO TO 1920
 1980 m = 0
 WRITE (no,1950) pi,m
 GO TO 1920
 1990 IF (j == 0) GO TO 9903
 GO TO 9000
 
!     PRINT CONTENTS OF
 
 2000 CONTINUE
 GO TO 9000
 
 2100 CONTINUE
 GO TO 9000
 
 2200 CONTINUE
 GO TO 9000
 
 2300 CONTINUE
 GO TO 9000
 
 2400 CONTINUE
 GO TO 9000
 
 2500 CONTINUE
 GO TO 9000
 
 2600 CONTINUE
 GO TO 9000
 
 2700 CONTINUE
 GO TO 9000
 
 2800 CONTINUE
 GO TO 9000
 
 2900 CONTINUE
 GO TO 9000
 
 3000 CONTINUE
 GO TO 9000
 
 3100 CONTINUE
 GO TO 9000
 
 
!     INTERNAL ROUTINE TO SET HEADINGS AND INITIALIZE LINE COUNTER.
!     -------------------------------------------------------------
 
 8000 CONTINUE
 DO  m = 1,32
   h1(m) = hx(m,m1)
   h2(m) = hx(m,m2)
   h3(m) = hx(m,m3)
 END DO
 line = nlpp
 GO TO r, (1110,1205,1255,1310,1501,1555,1710,1805,1835,1910)
 
 
!     PRINT TRAILER OF TABLE DATA BLOCK
!     ---------------------------------
 
 9000 CONTINUE
 WRITE  (no,9010) (t(l),l=2,7)
 9010 FORMAT (15H0*** trailer = ,6I18)
 GO TO 9999
 
 
 9901 WRITE  (no,9951) uwm,p
 9951 FORMAT (a25,' 2094, SUBROUTINE TABFMT, KEYNAME ',2A4,  &
     ' NOT IN LIST OF AVAILABLE KEYNAMES.')
 GO TO 9993
 
 9902 WRITE  (no,9952) uwm
 9952 FORMAT (a25,' 2095, SUBROUTINE TABFMT, PURGED INPUT.')
 GO TO 9995
 
 9903 WRITE  (no,9953) uwm
 9953 FORMAT (a25,' 2096, SUBROUTINE TABFMT, EOF ENCOUNTERED.')
 GO TO 9995
 
 9904 WRITE  (no,9954) uwm
 9954 FORMAT (a25,' 2097, SUBROUTINE TABFMT, EOR ENCOUNTERED.')
 GO TO 9995
 
 9905 WRITE  (no,9955) uwm
 9955 FORMAT (a25,' 2098, SUBROUTINE TABFMT, INSUFFICIENT CORE.')
 GO TO 9995
 
 9906 WRITE  (no,9956) uwm,kf
 9956 FORMAT (a25,' 2099, SUBROUTINE TABFMT, KF =',i10)
 GO TO 9995
 
 9993 WRITE  (no,9994) (na(1,l),na(2,l),l=1,la)
 9994 FORMAT ('0*** LIST OF RECOGNIZED KEYNAMES FOLLOWS...', /(20X,2A4))
 
 9995 CONTINUE
 
!     DO NOT CALL PEXIT SINCE THIS IS AN OUTPUT PROCESSOR.
 
 9999 CALL CLOSE (f,1)
 RETURN
 
END SUBROUTINE tabfmt
