SUBROUTINE giggks
     
!     THIS SUBROUTINE READS THE SPLINE CARDS AND DETERMINES THE
!     POINTS IN THE G AND K S
 
 INTEGER :: buff,buff1,sysbuf,buff2,eqt(7),  &
     scard(5),ccard(5),ss1(3),ls2(3),caero(3),set1(3),  &
     st2(3),ns(2),gkset,trl(7),TYPE,out,spl3(3),  &
     atab(2),pcstm,pbgpt,prcp,ptcp,pte,pre,ctyp,  &
     spline,useta,cstm,bagpdt,sila,ecta,gm,GO,scr1,  &
     scr2,scr3,scr4,scr5,ns1,ns2,ksize,gsize,gtka
 DIMENSION       c(18),x1b(3),x4b(3),temp(3),temp1(6),x1e(3),  &
     x4e(3),cb(18),b(6),z(28)
 DIMENSION       set2(8),crard(16)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ sysbuf,out
 COMMON /gicom / spline,useta,cstm,bagpdt,sila,ecta ,gm,GO,gtka,  &
     ksize,gsize,scr1,scr2,scr3,scr4,scr5
 COMMON /condas/ dum(3),degra
 COMMON /zzzzzz/ iz(1)
 
!     CHANGE IN EQUIV FOR SIZE OF SCARD OR CCARD
 
 EQUIVALENCE    (iz(1),z(1),scard(1),set2(1)),(z(28),nkset)
 EQUIVALENCE    (z(11),ccard(1),crard(1))    ,(z(27),ngset),  &
     (set2(3),sp1), (set2(4),sp2) ,(set2(5),ch1),  &
     (set2(6),ch2), (set2(7),z1 ) ,(set2(8),z2 )
 DATA       c / 18*0.0      /,  set1   / 3502,35,999/
 DATA     ss1 / 3302,33, 6  /,  ls2    / 3402,34,10 /,  &
     st2 / 3602,36, 8  /,  caero  / 3002,30,16 /
 DATA    spl3 / 4901,49, 1  /,  atab   / 200 ,2     /
 DATA      ns / 4HGIGG,  4HKS          /,iz2 /2     /
!     DATA    IECT / 3002,46     /
 
!     INITILIZE
 
 CALL sswtch (18,i18)
 nwds = korsz(iz)
 nogo = 0
 ns1  = 0
 ns2  = 0
 ns3  = 0
 
!     BUFF  HAS SPLINE
!     BUFF1 HAS CSTM,BGPT,EQAERO,SILA,SCR1
!     BUFF2 HAS SCR2
 
 buff  = nwds - sysbuf - 1
 buff1 = buff - sysbuf - 1
 buff2 = buff1- sysbuf
 
!     PROCESS SET CARDS AND WRITE G LISTS ON SCR2
 
 ifil = scr2
 CALL OPEN (*999,scr2,iz(buff2+1),1)
 ifil = spline
 CALL preloc (*999,iz(buff+1),spline)
 
!     SET1 CARDS
 
 CALL locate (*340,iz(buff+1),set1,idum)
 n   = 1
 nco = buff2 - n
 CALL READ (*998,*310,spline,iz(n),nco,1,nwr)
 GO TO 993
 310 i  = n - 1
 n1 = 0
 ASSIGN 335 TO TYPE
 320 i = i + 1
 IF (iz(i) == -1) GO TO 330
 IF (i    == nwr) GO TO 990
 n1 = n1 + 1
 GO TO 320
 330 IF (n1 < 2) GO TO 9971
 CALL WRITE (scr2,iz(n),n1,1)
 335 IF (i == nwr) GO TO 340
 n  = i + 1
 n1 = 0
 GO TO 320
 
!     SET 2 CARDS
 
 340 CALL locate (*490,iz(buff+1),st2,idum)
 
!     READ IN BAGPDT AND CSTM
 
 n = ls2(3) + caero(3) + 1
 trl(1) = cstm
 CALL rdtrl (trl)
 IF (trl(1) < 0) trl(3) = 0
 ncstm = (trl(3)+1)*14
 pcstm = buff2 - ncstm
 trl(1)= bagpdt
 CALL rdtrl (trl)
 nbg   = (trl(2)-trl(3))*4
 pbgpt = pcstm - nbg
 IF (pbgpt < n+150) GO TO 993
 
!     READ IN CSTM AT PCSTM + 14 ADD BASIC COORD SYSTEM
 
 iz(pcstm  ) = 0
 iz(pcstm+1) = 1
 DO  i = 2,13
   z(pcstm+i) = 0.0
 END DO
 z(pcstm+5 ) = 1.0
 z(pcstm+9 ) = 1.0
 z(pcstm+13) = 1.0
 IF (ncstm == 14) GO TO 7
 ifil = cstm
 CALL gopen (cstm,iz(buff1+1),0)
 CALL READ  (*998,*998,cstm,iz(pcstm+14),ncstm-14,1,nwr)
 CALL CLOSE (cstm,1)
 7 CONTINUE
 
!     READ IN BAGPDT AT PBGPT
 
 ifil = bagpdt
 CALL gopen (bagpdt,iz(buff1+1),0)
 CALL READ  (*998,*998,bagpdt,iz(pbgpt),nbg,1,nwr)
 CALL CLOSE (bagpdt,1)
 
!     READ IN SET2 CARDS WITH CAERO1 APPENDED
 
 ifil = spline
 lca  = 0
 ASSIGN 350 TO TYPE
 350 CALL READ (*998,*490,spline,iz(1),n-1,0,nwr)
 n1 = 1
 IF (ccard(1) == lca) GO TO 4001
 lca= ccard(1)
 k  = pcstm
 j  = pcstm + ncstm - 1
 IF (ccard(3) == 0) GO TO 371
 DO  i = k,j,14
   IF (ccard(3) == iz(i)) GO TO 370
 END DO
 GO TO 990
 370 prcp = i + 2
 ptcp = i + 5
 ctyp = iz(i+1)
 
!     LOCATE POINTS 1 AND 4 AS INPUT
 
 SELECT CASE ( ctyp )
   CASE (    1)
     GO TO 371
   CASE (    2)
     GO TO 372
   CASE (    3)
     GO TO 373
 END SELECT
 371 x1b(1) = crard(9)
 x1b(2) = crard(10)
 x1b(3) = crard(11)
 x4b(1) = crard(13)
 x4b(2) = crard(14)
 x4b(3) = crard(15)
 IF (ccard(3) == 0) GO TO 390
 GO TO 374
 372 x1b(1) = crard( 9)*COS(crard(10)*degra)
 x1b(2) = crard( 9)*SIN(crard(10)*degra)
 x1b(3) = crard(11)
 x4b(1) = crard(13)*COS(crard(14)*degra)
 x4b(2) = crard(13)*SIN(crard(14)*degra)
 x4b(3) = crard(15)
 GO TO 374
 373 x1b(1) = crard( 9)*SIN(crard(10)*degra)*COS(crard(11)*degra)
 x1b(2) = crard( 9)*SIN(crard(10)*degra)*SIN(crard(11)*degra)
 x1b(3) = crard( 9)*COS(crard(10)*degra)
 x4b(1) = crard(13)*SIN(crard(14)*degra)*COS(crard(15)*degra)
 x4b(2) = crard(13)*SIN(crard(14)*degra)*SIN(crard(15)*degra)
 x4b(3) = crard(13)*COS(crard(14)*degra)
 374 CALL gmmats (z(ptcp),3,3,0, x1b,3,1,0, temp)
 x1b(1) = temp(1) + z(prcp  )
 x1b(2) = temp(2) + z(prcp+1)
 x1b(3) = temp(3) + z(prcp+2)
 CALL gmmats (z(ptcp),3,3,0, x4b,3,1,0, temp)
 x4b(1) = temp(1) + z(prcp  )
 x4b(2) = temp(2) + z(prcp+1)
 x4b(3) = temp(3) + z(prcp+2)
 390 IF (ccard(2) == 0) GO TO 399
 
!     FIND ELEMENT COORDINATE SYSTEM
 
 DO  i = k,j,14
   IF (ccard(2) == iz(i)) GO TO 392
 END DO
 GO TO 990
 392 pre = i + 2
 pte = i + 5
 x1b(1) = x1b(1) - z(pre  )
 x1b(2) = x1b(2) - z(pre+1)
 x1b(3) = x1b(3) - z(pre+2)
 x4b(1) = x4b(1) - z(pre  )
 x4b(2) = x4b(2) - z(pre+1)
 x4b(3) = x4b(3) - z(pre+2)
 CALL gmmats (z(pte),3,3,1, x1b(1),3,1,0, x1e)
 CALL gmmats (z(pte),3,3,1, x4b(1),3,1,0, x4e)
 GO TO 400
 399 x1e(1) = x1b(1)
 x1e(2) = x1b(2)
 x4e(1) = x4b(1)
 x4e(2) = x4b(2)
 400 x2e = x1e(1) + crard(12)
 y2e = x1e(2)
 x3e = x4e(1) + crard(16)
 y3e = x4e(2)
 
!     FIND PRISM POINTS
 
 4001 CONTINUE
 px1 = (1.0-sp1)*(1.0-ch1)*x1e(1) + (1.0-sp1)*ch1*x2e +  &
     sp1*ch1*x3e + sp1*(1.0-ch1)*x4e(1)
 px2 = (1.0-sp1)*(1.0-ch2)*x1e(1) + (1.0-sp1)*ch2*x2e +  &
     sp1*ch2*x3e + sp1*(1.0-ch2)*x4e(1)
 px3 = (1.0-sp2)*(1.0-ch2)*x1e(1) + (1.0-sp2)*ch2*x2e +  &
     sp2*ch2*x3e + sp2*(1.0-ch2)*x4e(1)
 px4 = (1.0-sp2)*(1.0-ch1)*x1e(1) + (1.0-sp2)*ch1*x2e +  &
     sp2*ch1*x3e + sp2*(1.0-ch1)*x4e(1)
 
!     CHECK FOR BAD GEOMETRY
 
 IF (px1 > px2 .OR. px4 > px3) GO TO 997
 py1 = (1.0-sp1)*(1.0-ch1)*x1e(2) + (1.0-sp1)*ch1*y2e +  &
     sp1*ch1*y3e + sp1*(1.0-ch1)*x4e(2)
 py2 = (1.0-sp1)*(1.0-ch2)*x1e(2) + (1.0-sp1)*ch2*y2e +  &
     sp1*ch2*y3e + sp1*(1.0-ch2)*x4e(2)
 py3 = (1.0-sp2)*(1.0-ch2)*x1e(2) + (1.0-sp2)*ch2*y2e +  &
     sp2*ch2*y3e + sp2*(1.0-ch2)*x4e(2)
 py4 = (1.0-sp2)*(1.0-ch1)*x1e(2) + (1.0-sp2)*ch1*y2e +  &
     sp2*ch1*y3e + sp2*(1.0-ch1)*x4e(2)
 
!     BUILD PRISM INEQUALITY MATRICES
 
 c(1) = py1 - py2
 c(2) = px2 - px1
 c(4) = py2 - py3
 c(5) = px3 - px2
 c(7) = py3 - py4
 c(8) = px4 - px3
 c(10)= py4 - py1
 c(11)= px1 - px4
 c(15)= 0.0
 c(18)= 0.0
 b(1) = px2*py1 - px1*py2
 b(2) = px3*py2 - px2*py3
 b(3) = px4*py3 - px3*py4
 b(4) = px1*py4 - px4*py1
 nr   = 4
 IF (z1 == 0.0) GO TO 401
 c(15)=-1.0
 b(5) =-z1
 nr   = 5
 401 IF (z2 == 0.0) GO TO 404
 IF (z1 == 0.0) GO TO 403
 c(18)= 1.0
 b(6) = z2
 nr   = 6
 GO TO 404
 403 c(15)= 1.0
 b(5) = z2
 nr   = 5
 
!     CONVERT TO BASIC
 
 404 IF (ccard(2) == 0) GO TO 406
 CALL gmmats (c,nr,3,0, z(pte),3,3,1, cb)
 CALL gmmats (z(pte),3,3,1, z(pre),3,1,0, temp)
 CALL gmmats (c,nr,3,0, temp,3,1,0, temp1)
 b(1) = b(1) + temp1(1)
 b(2) = b(2) + temp1(2)
 b(3) = b(3) + temp1(3)
 b(4) = b(4) + temp1(4)
 IF (nr == 4) GO TO 405
 b(5) = b(5) + temp1(5)
 IF (nr == 5) GO TO 405
 b(6) = b(6) + temp1(6)
 GO TO 405
 406 DO  i = 1,18
   cb(i) = c(i)
 END DO
 405 CONTINUE
 
!     FINALLY TEST ALL GRID POINTS TO SEE IF THEY ARE IN PRISM
 
 kk = pbgpt
 kkk= kk + nbg - 1
 loop440:  DO  k = kk,kkk,4
   IF (iz(k) == -1) CYCLE loop440
   jj = 0
   DO  i = 1,nr
     sum = 0.0
     DO  j = 1,3
       jj  = jj + 1
       sum = sum + cb(jj)*z(k+j)
     END DO
     IF (sum < b(i)) CYCLE loop440
   END DO
   
!     FOUND ONE
   
   n1 = n1 + 1
   iz(n1) = (k-pbgpt)/4 + 1
 END DO loop440
 IF (n1  < 2) GO TO 997
 IF (i18 == 0) GO TO 446
 WRITE  (out,445) (iz(ii),ii=1,n1)
 445 FORMAT (5H0SET2 ,i8,2X,(/,10I9))
 446 CONTINUE
 CALL WRITE (scr2,iz(1),n1,1)
 GO TO 350
 490 CALL CLOSE (scr2,1)
 CALL OPEN  (*999,scr2,iz(buff2+1),0)
 neq = ksize*3
 eqt(1) = sila
 CALL rdtrl (eqt)
 nsil = eqt(2)
 ieq  = buff2 - neq - nsil
 
!     INITIAL CORE CHECK  PLUS FUDGE FACTOR
 
 IF (ieq-150 < 0) GO TO 993
 
!     READ SPLINE FOR K POINT POINTERS
 
!     READ SILA
 
 CALL locate (*990,iz(buff+1),atab,idum)
 CALL READ (*998,*11,spline,iz(ieq),neq+1,0,nwr)
 GO TO 990
 11 neq  = nwr
 ifil = sila
 CALL gopen (sila,iz(buff1+1),0)
 CALL READ  (*998,*998,sila,iz(ieq+neq),nsil,1,nwr)
 CALL CLOSE (sila,1)
 ifil  = spline
 trl(1)= scr1
 MAX   = 0
 CALL gopen (scr1,iz(buff1+1),1)
 
!     N = LENGTH OF LONGEST SPLINE CARD + CAERO1 CARD + 3
!     N  POINTS TO 1 ST LOCATION OF CORE AVAILABLE SEE EQIV
 
 n   = ls2(3) + caero(3) + 3
 nco = ieq - n
 
!     READ SPLINE1 CARDS
 
 CALL locate (*100,iz(buff+1),ss1,idum)
 ASSIGN 10 TO TYPE
 nr = ls2(3) + caero(3)
 10 CALL READ (*998,*100,spline,iz(1),nr,0,nwr)
 ns1 = ns1 + 1
 ASSIGN 30 TO gkset
 GO TO 300
 
!     G AND K SET ARE IN CORE SORTED  BY INTERNAL NUMBERS
!     A SECOND SET OF G   ARE SORTED  BY SIL NUMBERS
!     A SECOND SET OF K   ARE IN CORE BY K NUMBER
!     NK POINTS TO K SET
!     N1 IS FIRST LOCATION OF OPEN CORE
!     NGSET IS THE NUMBER OF G  NKSET FOR K
 
 30 IF (nogo == 1) GO TO 10
 
!     WRITE ALL SPLINE1 DATA ON SCR1 AS PROCESSED
!     ID OF SPLINE1 = 1
 
 iz(iz2) = 1
 nw  = n1 - 1
 MAX = MAX0(MAX,nw)
 CALL WRITE (scr1,iz(1),nw,1)
 GO TO 10
 
!     END OF SPLINE1 CARDS
 
!     READ SPLINE2 CARDS
 
 100 CALL locate (*190,iz(buff+1),ls2,idum)
 ASSIGN 110 TO TYPE
 nr = ls2(3) + caero(3)
 110 CALL READ (*998,*190,spline,iz(1),nr,0,nwr)
 ns2 = ns2 + 1
 ASSIGN 120 TO gkset
 GO TO 300
 
!     ID OF SPLINE2 = 2
 
 120 IF (nogo == 1) GO TO 110
 iz(iz2) = 2
 nw  = n1 - 1
 MAX = MAX0(MAX,nw)
 CALL WRITE (scr1,iz(1),nw,1)
 GO TO 110
 
!     END OF SPLINE2 CARDS
 
 190 CALL CLOSE (scr1,1)
 CALL CLOSE (scr2,1)
 CALL gopen (scr3,iz(buff1+1),1)
 
!     SPLINE 3 CARDS TO SCR3
 
 CALL locate (*290,iz(buff+1),spl3,idum)
 CALL READ (*998,*200,spline,iz,ieq,0,ns3)
 GO TO 993
 200 n = ns3 + 1
 
!     CONVERT AERO IDS TO K COLUMN NUMBERS, BUILD A LIST OF SPLINE CARD
!     POINTERS, SORT ON K COLUMNS, PROCESS CARDS IN SORTED ORDER GET
!     G POINTS TO SILS
 
 n1 = 1
 nw = ieq - 1
 ASSIGN 240 TO TYPE
 i = n
 210 k = iz(n1+3)
 DO  j = 1,neq,3
   IF (k == iz(nw+j)) GO TO 230
 END DO
 GO TO 992
 230 iz(n1+3) = iz(nw+j+2)
 iz(i   ) = n1
 iz(i +1) = iz(n1+3)
 i  = i+2
 240 n1 = n1 + iz(n1) + 1
 IF (n1 >= ns3) GO TO 250
 GO TO 210
 250 nw  = i - n
 ns3 = nw/2
 IF (ns3 == 0) GO TO 1001
 IF (ns3 == 1) GO TO 255
 CALL sort (0,0,2,2,iz(n),nw)
 
!     PROCESS BY SORTED ORDER
 
 255 n  = n - 1
 j  = ieq + neq - 1
 jj = 5
 DO  i = 1,nw,2
   n1 = iz(n+i)
   jjj= iz(n1) - caero(3)
   DO  k = jj,jjj,3
     l  = iz(n1+k)
     iz(n1+k) = iz(j+l)
   END DO
   CALL WRITE (scr3,iz(n1+1),iz(n1),1)
 END DO
 290 CALL CLOSE (spline,1)
 CALL CLOSE (scr3,1)
 CALL dmpfil (scr1,z,nwds)
 CALL dmpfil (scr3,z,nwds)
 trl(2) = MAX
 trl(3) = ns1 + ns2
 CALL wrttrl (trl)
 IF (nogo == 1) GO TO 1001
 IF (ns1 == 0 .AND. ns2 == 0 .AND. ns3 == 0) GO TO 990
 GO TO 1000
 
!     SET 1 CARDS
!     SET 2 CARDS
 
 300 ngset = 0
 ifil  = scr2
 301 CALL READ (*996,*996,scr2,iz(n),1,0,nwr)
 IF (scard(5) == iz(n)) GO TO 305
 CALL fwdrec (*998,scr2)
 GO TO 301
 305 CALL READ (*998,*306,scr2,iz(n),nco,1,nwr)
 GO TO 993
 306 CALL REWIND (scr2)
 ifil  = spline
 ngset = nwr
 n1 = n+ngset
 CALL sort (0,0,1,1,iz(n),ngset)
 
!     GET K SET
 
 nk   = n1 -1
 nkset= 0
 nmin = scard(3)
 nmax = scard(4)
 ncord= ccard(5)
 ifrst= ccard(1)
 IF (nmin > nmax) GO TO 990
 j1 = ncord*ccard(4) + ifrst - 1
 IF (nmin < ifrst .OR. nmax > j1) GO TO 990
 j1 = (nmin-ifrst)/ncord + 1
 i1 = (nmin-ifrst) - ncord*(j1-1) + 1
 jl = (nmax-ifrst)/ncord + 1
 il = (nmax-ifrst) - ncord*(jl-1) + 1
 DO  j = j1,jl
   DO  i = i1,il
     iz(n1) = ifrst + (i-1) + ncord*(j-1)
     n1 = n1 + 1
     nkset = nkset + 1
   END DO
 END DO
 
!     MAKE A LIST OF SIL NUMBERS   FOR G SET
 
 nw = ngset
 j = ieq + neq - 1
 DO  i = 1,nw
   k = iz(n+i-1)
   iz(n1) = iz(k+j)
   n1 = n1 + 1
 END DO
 
!     FIND INTERNAL K POINT NUMBER  FOR BGPT PLUS K NUMBER
 
 jj = 1
 nw = ieq - 1
 DO  i = 1,nkset
   DO  j = jj,neq,3
     IF (iz(nk+i) == iz(nw+j)) GO TO 550
   END DO
   GO TO 991
   550 jj = j
   iz(nk+i) = iz(nw+j+1)
   iz(n1  ) = iz(nw+j+2)
   n1 = n1 + 1
 END DO
 GO TO gkset, (30,120)
 
!     ERROR MESSAGES
 
 999 CALL mesage (-1,ifil,ns)
 998 CALL mesage (-3,ifil,ns)
 993 CALL mesage (-8,0,ns)
 990 CALL mesage (-7,0,ns)
 9971 scard(1) = iz(n)
 997 WRITE  (out,9970) uwm,scard(5),scard(1)
 9970 FORMAT (a25,' 2257, SET',i9,' REFERENCED ON SPLINE CARD',i9,  &
     ' IS EMPTY.')
 GO TO 901
 996 WRITE  (out,9960) ufm,scard(5),scard(1)
 9960 FORMAT (a23,' 2258, SET',i9,' REFERENCED ON SPLINE CARD',i9,  &
     ' NOT FOUND OR IT IS EMPTY.')
 CALL REWIND (scr2)
 GO TO 900
 991 WRITE  (out,9910) sfm,iz(nk+i-1),ccard(1)
 9910 FORMAT (a25,' 2259, POINT ASSIGNED TO BOX',i9,' FOR CAER01',i9,  &
     ' NOT IN EQAERO.')
 GO TO 900
 992 WRITE (out,9910) k,iz(n1+2)
 GO TO 900
 1001 CALL mesage (-61,0,ns)
 900 nogo = 1
 901 GO TO TYPE, (10,100,240,335,350)
 1000 RETURN
END SUBROUTINE giggks
