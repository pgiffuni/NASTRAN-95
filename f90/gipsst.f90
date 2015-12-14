SUBROUTINE gipsst
     
!     THIS SUBROUTINE LOCATES ALL THE G AND K SET POINTS IN THE SPLINE
!     COORDINATE SYSTEM AND FORMS G FOR EACH SET THEN
!     INSERTS THE G INTO THE FULL SIZED G MATRIX
 
 LOGICAL :: oxr,oyr,zap,kcol
 INTEGER :: sysbuf,out,buff,buff1,trl(7),tgkg(7),oldid,iz(28),  &
     pcstm,pbgpt,nwr,TYPE,ns(2),proe,pte,isng,slope,  &
     pg,pk,prol,ptl,buff2,ctype,psil
 INTEGER :: scard(10),ccard(16)
 INTEGER :: spline,useta,cstm,bagpdt,sila,eqaero,gm,GO,scr1,  &
     scr2,scr3,scr4,scr5,ksize,gsize,gtka
 DIMENSION       tl(9),rol(3),an(6),BLOCK(20),tgs(18),t(3),tg(9)
 DIMENSION       tt(9),z(1),srard(10)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zntpkx/ a,dum(3),nr,ieol,IEOR
 COMMON /system/ sysbuf,out
 COMMON /gicom / spline,useta,cstm,bagpdt,sila,eqaero,gm,GO,gtka,  &
     ksize,gsize,scr1,scr2,scr3,scr4,scr5
 COMMON /unpakx/ itc,ii,j1,incr
 COMMON /zzzzzz/ izx(1)
 
!     CHANGE IN EQUIV FOR SIZE OF SCARD OR CCARD
!     NEED TO CHANGE PENDC
 
 EQUIVALENCE     (izx(1),iz(1),z(1),scard(1),srard(1)),  &
     (iz(11),ccard(1)) ,(iz(27),ngset) ,(iz(28),nkset)
 DATA     ns   / 4HGIPS,4HST  /
 DATA     tgs  / 18*0.0       /
 
 pendc = 28
 nogo  = 0
 oldid =-1
 lc    =-1
 nwds  = korsz(iz)
 buff  = nwds - sysbuf
 buff1 = buff - sysbuf
 buff2 = buff1- sysbuf
 trl(1)= cstm
 CALL rdtrl (trl)
 IF (trl(1) < 0) trl(3) = 0
 ncstm = (trl(3)+1)*14
 pcstm = buff2 - ncstm
 trl(1)= bagpdt
 CALL rdtrl (trl)
 nbg   = trl(2)*4
 pbgpt = pcstm - nbg
 trl(1)= scr1
 CALL rdtrl (trl)
 MAX   = trl(2)
 IF (trl(3) == 0) GO TO 1000
 i     = scr2
 scr2  = scr3
 scr3  = i
 ipass = 0
 
!     INITIAL CORE CHECK
 
 IF (pbgpt-2*MAX < 0) GO TO 993
 
!     OPEN SCR1 TO LOOP ON G AND K SET RECORDS
 
 CALL gopen (scr1,iz(buff+1),0)
 
!     READ IN CSTM AT PCSTM + 14 ADD BASIC COORD SYSTEM
 
 1 iz(pcstm  ) = 0
 iz(pcstm+1) = 1
 DO  i = 2,13
   z(pcstm+i) = 0.0
 END DO
 z(pcstm+5) = 1.0
 z(pcstm+9) = 1.0
 z(pcstm+13)= 1.0
 IF (ncstm == 14) GO TO 7
 ifil = cstm
 CALL gopen (cstm,iz(buff1+1),0)
 CALL READ  (*999,*999,cstm,iz(pcstm+14),ncstm-14,1,nwr)
 CALL CLOSE (cstm,1)
 7 CALL pretrs (iz(pcstm),ncstm)
 
!     READ IN BAGPDT AT PBGPT
 
 ifil = bagpdt
 CALL gopen (bagpdt,iz(buff1+1),0)
 CALL READ  (*999,*999,bagpdt,iz(pbgpt),nbg,1,nwr)
 CALL CLOSE (bagpdt,1)
 
!     READ SCR1 AND PROCESS A SPLINE DEPENDING ON TYPE
 
 10 n1  = MAX + 1
 ifil= scr1
 CALL READ (*500,*20,scr1,iz(1),n1,1,nwr)
 20 j   = 2
 TYPE= iz(j)
 pg  = pendc
 pk  = pg + ngset
 psil= pk + nkset
 ipk = psil + ngset
 np  = ngset + nkset
 
!     USE A K POINT TO PICK UP POINTER TO BAGPDT FOR
!     COORDINATE SYSTEM ID OF SPLINE
 
 newid = ccard(2)
 ctype = ccard(8)
 k = pcstm
 j = pcstm + ncstm - 1
 IF (newid == oldid) GO TO 40
 DO  i = k,j,14
   IF (iz(i) /= newid) CYCLE
   proe = i + 2
   pte  = i + 5
   oldid = newid
   GO TO 40
 END DO
 ic = newid
 GO TO 997
 40 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 100
 END SELECT
 
!     SURFACE SPLINE
 
 50 SELECT CASE ( ctype )
   CASE (    1)
     GO TO 51
   CASE (    2)
     GO TO 998
   CASE (    3)
     GO TO 51
   CASE (    4)
     GO TO 51
   CASE (    5)
     GO TO 51
 END SELECT
 51 CONTINUE
 is  = 1
 ptl = pte
 DO  i = 1,9
   tl(i) = z(pte+i-1)
 END DO
 DO  i = 1,np
   k = (iz(pg+i)-1)*4
   
!     BASIC COORDINATES
   
   bx =  z(pbgpt+k+1)
   by =  z(pbgpt+k+2)
   bz =  z(pbgpt+k+3)
   IF (newid == 0) GO TO 55
   
!     X AND Y OF SPLINE
   
   t1 = bx -  z(proe  )
   t2 = by -  z(proe+1)
   t3 = bz -  z(proe+2)
   z(n1  ) = z(pte)*t1   + z(pte+3)*t2 + z(pte+6)*t3
   z(n1+1) = z(pte+1)*t1 + z(pte+4)*t2 + z(pte+7)*t3
   GO TO 59
   55 z(n1  ) = bx
   z(n1+1) = by
   59 n1 = n1 + 2
 END DO
 k = MAX + 1
 j = k + 2*ngset
 ncore = pbgpt - n1
 
!     CORE CHECK
 
 n  = ngset + 3
 nd = nkset*2
 nn = n*n + 3*n + n*nd + nd*ngset
 IF (nn < ncore) GO TO 70
 ncore = buff2 - n1
 IF (nn > ncore) GO TO 992
 zap =.true.
 
!     GET G FOR A SURFACE SPLINE
 
 70 CALL ssplin (ngset,iz(k),nkset,iz(j),0,0,1,1,scard(6),iz(n1), ncore,isng)
 IF (isng == 2) GO TO 998
 IF (nogo == 1) GO TO 10
 
!     REVERSE SIGN OF SLOPE COLUMN
 
 k = n1
 DO  i = 1,nkset
   k = k + ngset
   DO  j = 1,ngset
     z(k) = -z(k)
     k = k + 1
   END DO
 END DO
 GO TO 300
 
!     LINEAR SPLINE
 
 100 SELECT CASE ( ctype )
   CASE (    1)
     GO TO 101
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 101
   CASE (    4)
     GO TO 101
   CASE (    5)
     GO TO 101
 END SELECT
 
!     CAERO2 PROSESSING   BODIES
 
 130 scard( 8) = newid
 scard( 9) = scard(10)
 scard(10) = -1.0
 ibtyp = ccard(16)
 kd    = 1
 DO  i = 2,8
   tl(i) = 0.0
 END DO
 tl(1) = 1.0
 tl(5) = 1.0
 tl(9) = 1.0
 GO TO 102
 101 kd = 2
 102 CONTINUE
 
!     FIND CORD SYSTEM OF LINEAR SPLINE
 
 IF (scard(8) == lc) GO TO 120
 DO  i = k,j,14
   IF (scard(8) /= iz(i)) CYCLE
   lc   = scard(8)
   prol = i + 2
   ptl  = i + 5
   GO TO 120
 END DO
 ic = scard(8)
 GO TO 997
 120 IF (newid == 0 .AND. scard(8) == 0) GO TO 145
 t1 = z(prol  ) - z(proe)
 t2 = z(prol+1) - z(proe+1)
 t3 = z(prol+2) - z(proe+2)
 t1 = z(pte+2)*t1 + z(pte+5)*t2 + z(pte+8)*t3
 t2 = z(pte+5)*t1
 t3 = z(pte+8)*t1
 t1 = z(pte+2)*t1
 rol(1) = z(prol  ) - t1
 rol(2) = z(prol+1) - t2
 rol(3) = z(prol+2) - t3
 t1 = z(ptl+4)*z(pte+8) - z(ptl+7)*z(pte+5)
 t2 = z(ptl+7)*z(pte+2) - z(ptl+1)*z(pte+8)
 t3 = z(ptl+1)*z(pte+5) - z(ptl+4)*z(pte+2)
 t4 = SQRT(t1*t1 + t2*t2 + t3*t3)
 IF (t4 == 0.0) GO TO 996
 tl(1) = t1/t4
 tl(4) = t2/t4
 tl(7) = t3/t4
 tl(2) = z(pte+5)*tl(7) - z(pte+8)*tl(4)
 tl(5) = z(pte+8)*tl(1) - z(pte+2)*tl(7)
 tl(8) = z(pte+2)*tl(4) - z(pte+5)*tl(1)
 tl(3) = z(pte+2)
 tl(6) = z(pte+5)
 tl(9) = z(pte+8)
 145 DO  i = 1,np
   
!     BASIC CORD
   
   k  = (iz(pg+i)-1)*4
   bx =  z(pbgpt+k+1)
   by =  z(pbgpt+k+2)
   bz =  z(pbgpt+k+3)
   IF (newid == 0 .AND. scard(8) == 0) GO TO 150
   t1 = bx - rol(1)
   t2 = by - rol(2)
   t3 = bz - rol(3)
   z(n1  ) = tl(1)*t1 + tl(4)*t2 + tl(7)*t3
   z(n1+1) = tl(2)*t1 + tl(5)*t2 + tl(8)*t3
   GO TO 155
   150 z(n1  ) = bx
   z(n1+1) = by
   155 n1 = n1 + 2
 END DO
 IF (ctype /= 2) GO TO 169
 n1 = MAX + 1
 DO  i = 1,np
   z(n1+1) = z(n1)
   z(n1  ) = 0.0
   n1 = n1 + 2
 END DO
 
!     CHECK CORE
 
 169 k = MAX + 1
 j = k + 2*ngset
 ncore = pbgpt - n1
 oyr = .false.
 oxr = .false.
 IF (srard( 9) < 0.0) oxr = .true.
 IF (srard(10) < 0.0) oyr = .true.
 is = 3
 IF (oxr) is = is - 1
 IF (oyr) is = is - 1
 n  = is*ngset + 3
 nd = nkset*(1+kd)
 nn = n*n + 3*n + n*nd + nd*ngset*is
 IF (nn < ncore) GO TO 170
 ncore = buff2 - n1
 IF (nn > ncore) GO TO 992
 zap =.true.
 
!     GET G FOR A LINEAR SPLINE
 
 170 CALL lsplin (ngset,iz(k),nkset,iz(j),0,kd,1,scard(6),scard(9),  &
     scard(10),scard(7),iz(n1),ncore,isng)
 IF (isng == 2) GO TO 998
 IF (nogo == 1) GO TO 10
 IF (ctype == 2) GO TO 300
 
!     TRANSFORM G TO SPLINE COORDINATES
 
 tyl = 1.0
 txl = 0.0
 IF (newid == 0 .AND. scard(8) == 0) GO TO 190
 tyl = z(pte+1)*tl(2) + z(pte+4)*tl(5) + z(pte+7)*tl(8)
 txl = z(pte+1)*tl(1) + z(pte+4)*tl(4) + z(pte+7)*tl(7)
 
!     MOVE COLUMNS UP
 
 190 nrgs = ngset*is
 k2   = nrgs + nrgs
 k3   = k2 + nrgs
 ncore= n1
 n1   = n1 + nrgs - 1
 n2   = n1
 DO  i = 1,nkset
   DO  k = 1,nrgs
     z(n2+k) = z(n1+k)*txl + z(n1+nrgs+k)*tyl
     z(n2+nrgs+k) = z(n1+k2+k)
   END DO
   n1 = n1 + k3
   n2 = n2 + k2
 END DO
 n1 = ncore
 
!     TRANSFORM G INTO GLOBAL
 
 300 CONTINUE
 
!                         T
!     OPEN SCR2 TO WRITE G   MATRIX
!                         KG
 
 CALL gopen (scr2,iz(buff1+1),1)
 CALL gopen (scr3,iz(buff2+1),0)
 tgkg(3) = gsize
 tgkg(4) = 2
 tgkg(5) = 1
 tgkg(1) = scr2
 tgkg(2) = 0
 tgkg(6) = 0
 tgkg(7) = 0
 ibcc    = 1
 SIGN    = 1.0
 slope   = 1
 kcol    = .false.
 kn      = 1
 kcoln   = iz(ipk+kn)
 
!     KCOLN PICKS UP COLUMN NUMBER TO INSERT
!     KN POINT TO COLUMN OF G MATRIX
!     SLOPE IS FLIP FLOP SWITCH FOR SLOPE COLUMN (KEEPS KCOL TRUE)
 
 
!     LOOP THROUGH COLUMNS OF GKT
 
 DO  i = 1,ksize
   CALL bldpk (1,1,scr2,BLOCK,1)
   IF (kcoln == i) kcol = .true.
   
!     COPY A COLUMN OR OUTPUT A NULL COLUMN
   
   CALL intpk (*340,scr3,0,1,0)
   IF (kcol) GO TO 995
   330 CALL zntpki
   CALL bldpki (a,nr,scr2,BLOCK)
   IF (ieol == 0) GO TO 330
   GO TO 390
   340 IF (.NOT.kcol) GO TO 390
   
!     LOOP THROUGH COLUMN OF G BUILDING COLUMN OF GKT
   
   DO  j = 1,ngset
     nr = iz(psil+j)
     k  = (iz(pg+j)-1)*4
     CALL transs (iz(pbgpt+k),tt)
     CALL gmmats (tt,3,3,1,tl,3,3,0,tg)
     SELECT CASE ( TYPE )
       CASE (    1)
         GO TO 350
       CASE (    2)
         GO TO 360
     END SELECT
     
!     TERMS OF SURFACE SPLINE
     
     350 CONTINUE
     DO  jj = 3,9,3
       a = tg (jj)*z(n1)
       CALL bldpki (a,nr,scr2,BLOCK)
       nr = nr + 1
     END DO
     n1 = n1 + 1
     CYCLE
     
!     TERMS OF LINEAR SPLINE
     
     360 IF (ctype == 2) GO TO 370
     IF (is == 1) GO TO 350
     tgs( 1) = tg(3)
     tgs( 4) = tg(6)
     tgs( 7) = tg(9)
     tgs(11) = tg(1)
     tgs(12) = tg(2)
     tgs(14) = tg(4)
     tgs(15) = tg(5)
     tgs(17) = tg(7)
     tgs(18) = tg(8)
     GO TO 365
     
!     BODIES
     
     370 SELECT CASE ( ibtyp )
       CASE (    1)
         GO TO 372
       CASE (    2)
         GO TO 371
       CASE (    3)
         GO TO 373
     END SELECT
     371 SELECT CASE ( ibcc )
       CASE (    1)
         GO TO 373
       CASE (    2)
         GO TO 372
       CASE (    3)
         GO TO 372
       CASE (    4)
         GO TO 373
     END SELECT
     372 tgs( 1) = tg(3)*SIGN
     tgs( 4) = tg(6)*SIGN
     tgs( 7) = tg(9)*SIGN
     tgs(11) =-tg(2)*SIGN
     tgs(12) = tg(1)*SIGN
     tgs(14) =-tg(5)*SIGN
     tgs(15) = tg(4)*SIGN
     tgs(17) =-tg(8)*SIGN
     tgs(18) = tg(7)*SIGN
     GO TO 365
     373 tgs( 1) = tg(2)
     tgs( 4) = tg(5)
     tgs( 7) = tg(8)
     tgs(11) = tg(3)
     tgs(12) = tg(1)
     tgs(14) = tg(6)
     tgs(15) = tg(4)
     tgs(17) = tg(9)
     tgs(18) = tg(7)
     365 t(1) = z(n1)
     n1   = n1 + 1
     t(2) = 0.0
     t(3) = 0.0
     IF (oxr) GO TO 361
     t(2) = z(n1)
     n1   = n1 + 1
     361 IF (oyr) GO TO 362
     t(3) = z(n1)
     n1   = n1 + 1
     362 CALL gmmats (tgs,6,3,0,t,3,1,0,an)
     DO  jj = 1,6
       CALL bldpki (an(jj),nr,scr2,BLOCK)
       nr = nr + 1
     END DO
   END DO
   
!     COLUMN FINISHED CHECKSLOPE COLUMN NEXT OR END OF G
   
   IF (ctype /= 3) GO TO 382
   n1 = n1 + ngset*is
   GO TO 384
   382 IF (ctype /= 2) GO TO 383
   IF (ibtyp == 1) SIGN = -SIGN
   IF (ibtyp /= 2) GO TO 383
   ibcc = ibcc + 1
   IF (ibcc == 3) SIGN = -SIGN
   IF (ibcc == 5) SIGN = -SIGN
   IF (ibcc == 5) ibcc = 1
   
!     KEEP SLOPE NEG FOR ZY BODIES AND REPROCESS SAME COLUMN TWICE
   
   IF (ibcc == 2 .OR. ibcc == 4) n1 = n1 - ngset*is
   IF (ibcc > 2) GO TO 390
   383 slope = -slope
   IF (slope /= 1) GO TO 390
   384 kn = kn + 1
   IF (kn > nkset) GO TO 385
   kcoln = iz(ipk+kn)
   385 kcol  = .false.
   390 CALL bldpkn (scr2,BLOCK,tgkg)
 END DO
 
!     SWITCH FILES FOR ANOTHER SPLINE
 
 CALL CLOSE (scr2,1)
 CALL wrttrl (tgkg)
 CALL CLOSE (scr3,1)
 i    = scr2
 scr2 = scr3
 scr3 = i
 ipass= ipass + 1
 IF (zap) GO TO 1
 GO TO 10
 
!     FINISHED SWITCH FILES SO OUTPUT IS SCR2
 
 
!     IF ALL DONE BE SURE SCR2 IS GTKA
 
 500 i    = scr2
 scr2 = scr3
 scr3 = i
 IF (scr3 /= 201) GO TO 520
 CALL gopen (scr2,z(buff1),0)
 CALL gopen (scr3,z(buff2),1)
 tgkg(1) = scr2
 CALL rdtrl (tgkg)
 n    = tgkg(2)
 tgkg(1) = scr3
 tgkg(2) = 0
 tgkg(6) = 0
 tgkg(7) = 0
 incr = 1
 itc  = 1
 CALL cyct2b (scr2,scr3,n,z,tgkg)
 CALL CLOSE (scr2,1)
 CALL CLOSE (scr3,1)
 CALL wrttrl (tgkg)
 520 CONTINUE
 CALL CLOSE (scr1,1)
 IF (nogo == 0) GO TO 1000
 
!     ERROR MESSAGES
 
 CALL mesage (-61,0,ns)
 999 CALL mesage (-3,ifil,ns)
 993 CALL mesage (-8,0,ns)
 998 WRITE  (out,9980) ufm,scard(1)
 9980 FORMAT (a23,' 2260, SINGULAR MATRIX DEVELOPED WHILE PROCESSING ',  &
     'SPLINE',i9)
 GO TO  1001
 997 CALL mesage (30,25,ic)
 GO TO 1001
 996 WRITE  (out,9960) ufm,scard(1),ccard(1)
 9960 FORMAT (a23,' 2261, PLANE OF LINEAR SPLINE',i9,  &
     ' PERPENDICULAR TO PLANE OF AERO ELEMENT',i9)
 GO TO 1001
 995 WRITE  (out,9950) ufm,scard(1)
 9950 FORMAT (a23,' 2262, SPLINE',i9,' INCLUDES AERO BOX INCLUDED ON A',  &
     ' EARLIER SPLINE')
 GO TO 1001
 992 WRITE  (out,9920) ufm,scard(1)
 9920 FORMAT (a23,' 2263, INSUFFICIENT CORE TO PROCESS SPLINE',i9)
 1001 nogo = 1
 GO TO 10
 1000 RETURN
END SUBROUTINE gipsst
