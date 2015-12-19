SUBROUTINE em3d (eltype,istart,itype,ncount,ido,iwords,nbdys,all,  &
        nelout)
     
!     E  AND  M LOADS FOR 3-D ELEMENTS
!     TETRA  39   WEDGE  40   HEXA1 41  HEXA2  42
!     IHEX1  65   IHEX2  66   IHEX3 67
 
 
 INTEGER, INTENT(IN)                      :: eltype
 INTEGER, INTENT(IN)                      :: istart
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER, INTENT(IN)                      :: ncount
 INTEGER, INTENT(IN)                      :: ido
 INTEGER, INTENT(IN)                      :: iwords
 INTEGER, INTENT(IN)                      :: nbdys
 INTEGER, INTENT(IN OUT)                  :: all
 INTEGER, INTENT(IN)                      :: nelout
 LOGICAL :: onlyc
 INTEGER :: scr6, typold,elid,  outpt,sysbuf,pointr(7,7),frstgd,tmap(88)
 REAL :: ll(4,5),w(5)
 DIMENSION       isc(5),sc(5),xlacc(3),buf(50),ibuf(50),hcx3(60),  &
     hcx(4),hcy(4),hcz(4), g(9),necpt(1),dr(24),iz(1),ip(4),r(3,8),xload(8),  &
     gpt(32),bxyz(3,32),s(4),h(4),gauss(8),f(32),  &
     shp(32),dshp(3,32),xjacob(3,3),dshpb(3,32),hc(96), hcxyz(3),gh(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ksystm(2)
 COMMON /zzzzzz/ z(1)
 COMMON /emecpt/ ecpt(200)
 COMMON /matin / matid,inflag,eltemp
 COMMON /hmtout/ xmat(6)
 EQUIVALENCE     (ksystm(1),sysbuf), (ksystm(2),outpt),  &
     (ecpt(1),necpt(1)), (z(1),iz(1)), (i1,ip(1)),  &
     (i2,ip(2)),(i3,ip(3)), (i4,ip(4))
 EQUIVALENCE     (buf(1),ibuf(1))  , (isc(1),sc(1))
 
!     GRID POINT NO FOR EACH ELEMENT
 
 DATA    tmap  / 1, 2, 3, 4,    1, 2, 3, 5,    1, 2, 3, 6,  &
     1, 4, 5, 6,    2, 4, 5, 6,    3, 4, 5, 6,  &
     1, 2, 4, 6,    2, 3, 4, 6,    1, 3, 4, 5,  &
     2, 3, 4, 5,    1, 3, 5, 6,    1, 2, 5, 6,  &
     1, 2, 3, 6,    1, 3, 4, 8,    1, 3, 8, 6,  &
     1, 5, 6, 8,    3, 6, 7, 8,    2, 3, 4, 7,  &
     1, 2, 4, 5,    2, 4, 5, 7,    2, 5, 6, 7, 4, 5, 7, 8/
!     DATA    NAM   / 4HEM3D,4H       /
 DATA    typold/ 0               /
 DATA    scr6  / 306             /
 
!     SET UP GAUSSIAN INTEGRATION POINTS
 
 DATA    gauss / 0.57735027,     0.55555556, 0.77459667,     0.88888889,  &
     0.34785484,     0.86113631, 0.65214515,     0.33998104/
 
!     SET UP POINTR ARRAY ONTO EST
 
!                    TYPE  MID   FRSTGD  ISYS1   NIP     ITEMP     NELS
 DATA    pointr/ 39,   2,    3    ,  7,      0,      23,       1,  &
     40,   2,    3,      9,      0,      33,      12,  &
     41,   2,    3,      11,     0,      43,       5,  &
     42,   2,    3,      11,     0,      43,       10,  &
     65,   10,   2,      16,     12,     48,       1,  &
     66,   22,   2,      28,     24,     108,      1,  &
     67,   34,   2,      40,     36,     168,      1 /
 
 
 onlyc  =.false.
 nopts  = 6
 IF (eltype == typold) GO TO 30
 typold = eltype
 DO  l = 1,7
   ilis   = l
   IF (eltype-pointr(1,l) < 0.0) THEN
     GO TO  1200
   ELSE IF (eltype-pointr(1,l) == 0.0) THEN
     GO TO    20
   ELSE
     GO TO    10
   END IF
   10 CONTINUE
 END DO
 GO TO 1200
 
!     SET UPPOINTERS INTO EST(ECPT) DATA
 
 20 mid = pointr(2,ilis)
 
!     MATERIAL ID
 
 frstgd = pointr(3,ilis)
 
!     FIRST SIL
 
 isys1 = pointr(4,ilis)
 
!     FIRST CSIL
 
 nip = pointr(5,ilis)
 
!     NO OF INTEGRATION POINTS (ISOPARAMETRICS ONLY)
 
 itemp = pointr(6,ilis)
 
!     TEMPERATURE DATA
 
 nels = pointr(7,ilis)
 
!     NO. OF ELEMENTS
 
!     GO TO SECTION 190 FOR ISOPARAMETRICS
 
!     CHECK FOR ZERO LOAD
 
 30 ngrid = isys1 - frstgd
 IF (eltype >= 65) ngrid = ngrid - 6
!                     65 TO 67 ??
 isc(1) = necpt(1)
 isc(2) = 1
 IF (eltype == 65) isc(2) = 9
 IF (eltype == 66 .OR. eltype == 67) isc(2) = 21
 
!     CHECK TO SEE IF THIS ELEMENT CONTAINS A GRID POINT ON A PERMBDY
!     CARD. IF SO, OR IF NO PERMBDY CARD EXISTS, COMPUTE LOADS FOR THE
!     ELEMENT. IF NOT, COMPUTE HC CENTROIDAL VALUE ONLY. (ONLYC=.TRUE.)
!     THE PERMBDY SILS START AT Z(ISTART-NBDYS-1)
 
 IF (nbdys == 0) GO TO 50
 
 DO  i = 1,ngrid
   ng = necpt(frstgd+i-1)
   DO  j = 1,nbdys
     IF (ng == iz(istart-nbdys-nelout+j-1)) GO TO 50
   END DO
 END DO
 
!     ELEMENT HAS NO GRIDS ON PERMBDY
 
 onlyc =.true.
 nopts = 1
 50 IF (onlyc .AND. itype == 24) RETURN
 
!     IF ONLYC=TRUE, CHECK TO SEE IF THE ELEMENT HAD AN ELFORCE REQUEST.
!     IF SO, CONTINUE. IF NOT, JUST WRITE ZEROS TO HCCEN,SCR6) AND
!     RETURN.
 
 IF (.NOT.onlyc) GO TO 70
 IF (all == 1) GO TO 70
 IF (nelout == 0) GO TO 100
 
 DO  i = 1,nelout
   IF (necpt(1) == iz(istart-nelout+i-1)) GO TO 70
 END DO
 GO TO 100
 70 IF (itype /= 20 .AND. itype /= 24) GO TO 130
 g1 = 0.
 g2 = 0.
 g3 = 0.
 h1 = 0.
 h2 = 0.
 h3 = 0.
 DO  i = 1,ngrid
   isub = istart + 3*necpt(frstgd+i-1) - 3
   IF (itype == 24) isub = istart + 3*ncount - 3
   h1 = h1 + ABS(z(isub  ))
   h2 = h2 + ABS(z(isub+1))
   h3 = h3 + ABS(z(isub+2))
   g1 = g1 + z(isub  )
   g2 = g2 + z(isub+1)
   g3 = g3 + z(isub+2)
   IF (itype == 24) EXIT
 END DO
 90 hl = h1 + h2 + h3
 IF (hl /= 0.) GO TO 120
 IF (itype == 24) RETURN
 
 100 sc(3) = 0.
 sc(4) = 0.
 sc(5) = 0.
 CALL WRITE (scr6,sc,2,0)
 isc2  = isc(2)
 DO  i = 1,isc2
   CALL WRITE (scr6,sc(3),3,0)
 END DO
 RETURN
 
 120 IF (itype == 24) GO TO 130
 
!     AVERGAGE SPCFLD
 
 ahcx = g1/FLOAT(ngrid)
 ahcy = g2/FLOAT(ngrid)
 ahcz = g3/FLOAT(ngrid)
 
 130 IF (eltype >= 65) GO TO 500
 IF (onlyc) GO TO 140
 
!     GET MATERIAL INFO
!     INFLAG = 3  RETURNS A 3X3 MATRIX
 
 ll(1,1) = .25
 ll(2,1) = .25
 ll(3,1) = .25
 ll(4,1) = .25
 ll(1,2) = .5
 ll(2,2) = 1./6.
 ll(3,2) = ll(2,2)
 ll(4,2) = ll(2,2)
 ll(1,3) = 1./6.
 ll(2,3) = .5
 ll(3,3) = ll(1,3)
 ll(4,3) = ll(1,3)
 ll(1,4) = 1./6.
 ll(2,4) = ll(1,4)
 ll(3,4) = .5
 ll(4,4) = ll(1,4)
 ll(1,5) = 1./6.
 ll(2,5) = ll(1,5)
 ll(3,5) = ll(1,5)
 ll(4,5) = .5
 w(1)    =-.8
 w(2)    = 9./20.
 w(3)    = w(2)
 w(4)    = w(2)
 w(5)    = w(2)
 inflag  = 3
 matid   = necpt(mid)
 eltemp  = ecpt(itemp)
 CALL hmat (necpt(1))
 
!     G STORED BY ROW
 
 g(1) = xmat(1)
 g(2) = xmat(2)
 g(3) = xmat(3)
 g(4) = xmat(2)
 g(5) = xmat(4)
 g(6) = xmat(5)
 g(7) = xmat(3)
 g(8) = xmat(5)
 g(9) = xmat(6)
 
!     PUT COORDINATES OF GRID POINTS INTO ARRAY
!     FOR HEXA2  DIVIDE VOLUME BY 2.
 
 xm = 1.
 IF (eltype == 42) xm = 2.
 
!     TETRA   4 GRID PTS    1 ELEMENT
!     WEDGE   6 GRID PTS   18 ELEMENTS(6 ARE DUPLICATES-4 POINTS AT A
!     HEXA1   8 GRID PTS    5 ELEMENT (4 PTS AT A TIME)
!     HEXA2   8 GRID PTS    10ELEMENT (4 PTS AT A TIME)
!     SET UP PROPER POINTERS VIA TMAP
!     R ARRAY CONTAINS COORDINATE INFO
 
 140 DO  i = 1,ngrid
   itt    = isys1 + 4*i - 4
   r(1,i) = ecpt(itt+1)
   r(2,i) = ecpt(itt+2)
   r(3,i) = ecpt(itt+3)
 END DO
 
!     SET UP POINTER TO GRID PT NO
 
 irow = 0
 IF (eltype == 41 .OR. eltype == 42) irow = 12
 DO  i = 1,8
   xload(i) = 0.0
 END DO
 
!     SET UP POINTS FOR AVERAGE COORDINATES
 
 xxc = 0.
 yyc = 0.
 zzc = 0.
 DO  i = 1,ngrid
   xxc = xxc + r(1,i)
   yyc = yyc + r(2,i)
   zzc = zzc + r(3,i)
 END DO
 xxc = xxc/FLOAT(ngrid)
 yyc = yyc/FLOAT(ngrid)
 zzc = zzc/FLOAT(ngrid)
 
!     PRINCIPAL LOOP OVER ELEMENT OF THE GIVEN TYPE
 
 DO  iel = 1,nels
   IF (onlyc) GO TO 200
   
!     RESET XM FOR WEDGES. 1ST 12 CONFIGURATIONS ARE MULTIPLIED BY 2.
!     ALL 18 ARE DIVIDED BY 6.(SINCE XM IS A DIVISOR, USE RECIPROCALS)
   
   IF (eltype == 40 .AND. iel <= 6) xm = 6./2.
   IF (eltype == 40 .AND. iel > 6) xm = 6.
   isub = (irow+iel-1)*4
   DO  i = 1,4
     f(i)  = 0.
     ip(i) = i
     IF (eltype >= 40) ip(i) = tmap(isub+i)
   END DO
   
!     NEED DET TO COMPUTE VOL
   
   term1 =  r(3,i4)*((r(1,i2)-r(1,i1))*r(2,i3) +  &
       (r(1,i1)-r(1,i3))*r(2,i2) + (r(1,i3)-r(1,i2))*r(2,i1))
   term2 =  r(3,i3)*((r(1,i1)-r(1,i2))*r(2,i4) +  &
       (r(1,i4)-r(1,i1))*r(2,i2) + (r(1,i2)-r(1,i4))*r(2,i1))
   term3 =  r(3,i2)*((r(1,i3)-r(1,i1))*r(2,i4) + (r(1,i1)-r(1,i4))*  &
       r(2,i3) + (r(1,i4)-r(1,i3))*r(2,i1))
   term4 =  r(3,i1)*((r(1,i2)-r(1,i3))*r(2,i4) + (r(1,i4)-r(1,i2))*  &
       r(2,i3) + (r(1,i3)-r(1,i4))*r(2,i2))
   det   =  term1 + term2 + term3 + term4
   vol   =  ABS(det)/6.
   
!     GRADIENTS OF SHAPE FUNCTIONS
   
   dr( 1) = r(3,i3)*r(2,i4) - r(3,i4)*r(2,i3) + r(2,i2)*(r(3,i4)-  &
       r(3,i3)) - r(3,i2)*(r(2,i4)-r(2,i3))
   dr( 2) = r(1,i3)*r(3,i4) - r(1,i4)*r(3,i3) - r(1,i2)*(r(3,i4)-  &
       r(3,i3)) + r(3,i2)*(r(1,i4)-r(1,i3))
   dr( 3) = r(2,i3)*r(1,i4) - r(1,i3)*r(2,i4) + r(1,i2)*(r(2,i4)-  &
       r(2,i3)) - r(2,i2)*(r(1,i4)-r(1,i3))
   dr( 4) = r(2,i3)*r(3,i4) - r(2,i4)*r(3,i3) - r(2,i1)*(r(3,i4)-  &
       r(3,i3)) + r(3,i1)*(r(2,i4)-r(2,i3))
   dr( 5) = r(1,i4)*r(3,i3) - r(1,i3)*r(3,i4) + r(1,i1)*(r(3,i4)-  &
       r(3,i3)) - r(3,i1)*(r(1,i4)-r(1,i3))
   dr( 6) = r(1,i3)*r(2,i4) - r(2,i3)*r(1,i4) - r(1,i1)*(r(2,i4)-  &
       r(2,i3)) + r(2,i1)*(r(1,i4)-r(1,i3))
   dr( 7) = r(3,i2)*r(2,i4) - r(2,i2)*r(3,i4) + r(2,i1)*(r(3,i4)-  &
       r(3,i2)) - r(3,i1)*(r(2,i4)-r(2,i2))
   dr( 8) = r(1,i2)*r(3,i4) - r(1,i4)*r(3,i2) - r(1,i1)*(r(3,i4)-  &
       r(3,i2)) + r(3,i1)*(r(1,i4)-r(1,i2))
   dr( 9) = r(2,i2)*r(1,i4) - r(1,i2)*r(2,i4) + r(1,i1)*(r(2,i4)-  &
       r(2,i2)) - r(2,i1)*(r(1,i4)-r(1,i2))
   dr(10) = r(2,i2)*r(3,i3) - r(3,i2)*r(2,i3) - r(2,i1)*(r(3,i3)-  &
       r(3,i2)) + r(3,i1)*(r(2,i3)-r(2,i2))
   dr(11) = r(3,i2)*r(1,i3) - r(1,i2)*r(3,i3) + r(1,i1)*(r(3,i3)-  &
       r(3,i2)) - r(3,i1)*(r(1,i3)-r(1,i2))
   dr(12) = r(1,i2)*r(2,i3) - r(2,i2)*r(1,i3) - r(1,i1)*(r(2,i3)-  &
       r(2,i2)) + r(2,i1)*(r(1,i3)-r(1,i2))
   
   DO   k = 1,12
     dr(k) = dr(k)/det
   END DO
   
!     MULTIPLY SHAPE FUNCTION  BY G
   
   IF (itype /= 24) CALL gmmats (dr(1),4,3,0,g,3,3,0,dr(13))
   
!     COMPUTE HC
   
   IF (itype /= 24) GO TO 200
   
!     REMFLUX
   
   nsubx = istart + 3*ncount - 3
   hc(1) = z(nsubx)
   hc(2) = z(nsubx+1)
   hc(3) = z(nsubx+2)
   GO TO 360
   
!     INTEGRATE TO GET HC
   
   200 ktype = itype - 19
   xlacc(1) = 0.
   xlacc(2) = 0.
   xlacc(3) = 0.
   
!     START INTEGRATION PROCEDURE-NEED 5 POINTS FOR CUBIC + CENTROID
   
   DO  npts = 1,nopts
     
!     DO CENTROID FOR ONLY 1ST TETRA
     
     IF (npts == nopts .AND. iel > 1) CYCLE
     
!     COMPUTE BASIC COORDS OF INTEGRATION POINT
     
     IF (npts /= nopts) GO TO 210
     
!     CENTROID
     
     xx = xxc
     yy = yyc
     zz = zzc
     IF (itype /= 20) GO TO 220
     
!     AVERAGE SPCFLD
     
     hc(1) = ahcx
     hc(2) = ahcy
     hc(3) = ahcz
     GO TO 310
     210 xx = ll(1,npts)*r(1,i1) + ll(2,npts)*r(1,i2) + ll(3,npts)*r(1,i3)  &
         + ll(4,npts)*r(1,i4)
     yy = ll(1,npts)*r(2,i1) + ll(2,npts)*r(2,i2) + ll(3,npts)*r(2,i3)  &
         + ll(4,npts)*r(2,i4)
     zz = ll(1,npts)*r(3,i1) + ll(2,npts)*r(3,i2) + ll(3,npts)*r(3,i3)  &
         + ll(4,npts)*r(3,i4)
     220 hc(1) = 0.
     hc(2) = 0.
     hc(3) = 0.
     
!     COMPUTE HC AT THIS PPOINT FOR ALL LOADS OF THIS TYPE IN THIS
!     SUBCASE
     
     DO  ijk = 1,ido
       IF (itype == 20) GO TO 240
       isub = istart + (ijk-1)*iwords - 1
       DO  i = 1,iwords
         buf(i) = z(isub+i)
       END DO
       
       SELECT CASE ( ktype )
         CASE (    1)
           GO TO 240
         CASE (    2)
           GO TO 260
         CASE (    3)
           GO TO 270
         CASE (    4)
           GO TO 280
       END SELECT
       
!     SPCFLD
       
       240 DO  i = 1,4
         isil = frstgd - 1 + ip(i)
         ist  = istart + 3*necpt(isil) - 3
         hcx(i) = z(ist  )
         hcy(i) = z(ist+1)
         hcz(i) = z(ist+2)
       END DO
       hc1 = ll(1,npts)*hcx(1) + ll(2,npts)*hcx(2) + ll(3,npts)*hcx(3) +  &
           ll(4,npts)*hcx(4)
       hc2 = ll(1,npts)*hcy(1) + ll(2,npts)*hcy(2) + ll(3,npts)*hcy(3) +  &
           ll(4,npts)*hcy(4)
       hc3 = ll(1,npts)*hcz(1) + ll(2,npts)*hcz(2) + ll(3,npts)*hcz(3) +  &
           ll(4,npts)*hcz(4)
       GO TO 290
       
!     CEMLOOP,GEMLOOP,MDIPOLE
       
       260 CALL axloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 290
       270 CALL geloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 290
       280 CALL dipole (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
       290 hc(1) = hc(1) + hc1
       hc(2) = hc(2) + hc2
       hc(3) = hc(3) + hc3
     END DO
     310 IF (npts /= nopts) GO TO 320
     sc(3) = hc(1)
     sc(4) = hc(2)
     sc(5) = hc(3)
     CALL WRITE (scr6,sc,5,0)
     CYCLE
     
!     WE HAVE HC AT THIS POINT. MULT. BY  WEIGHT AND ACCUMULATE
     
     320 DO  i = 1,3
       xlacc(i) = xlacc(i) + hc(i)*w(npts)
     END DO
     
!     GET ANOTHER INTEGRATTION POINT
     
   END DO
   
   IF (onlyc) RETURN
   DO  i = 1,3
     hc(i) = xlacc(i)
   END DO
   360 CONTINUE
   
!     MULTIPLY HC BY GRADIENTS AND MATERIALS
   
   isubx = 13
   IF (itype == 24) isubx = 1
   CALL gmmats (dr(isubx),4,3,0,hc,3,1,0,f(1))
   
   DO  k = 1,4
     kk = ip(k)
     xload(kk) = xload(kk) + f(k)*vol/xm
   END DO
   
!     XLOAD   IS SUM OF ALL LOADS FOR ALL THE ELEMENTS
!     F COMPUTED FOR A GIVEN TETRA OF THE TOTAL SHAPE
!     SO MULTIPLY BY VOL
   
 END DO
 
 DO   i = 1,ngrid
   is   = frstgd - 1 + i
   isil = necpt(is)
   
!     IF PERMBDY EXISTS AND IF GRID IS NOT ON IT, IGNORE ITS LOAD
   
   IF (nbdys == 0) GO TO 420
   DO  j = 1,nbdys
     IF (isil /= iz(istart-nbdys-nelout+j-1)) CYCLE
     GO TO 420
   END DO
   CYCLE
   420 z(isil) = z(isil) - xload(i)
 END DO
 RETURN
 
!     ISOPARAMETRIC SOLIDS
 
 500 jtype = itype
 itype = eltype - 64
 inip  = necpt(nip)
 IF (inip == 0) inip = itype/2 + 2
 np    = 12*itype - 4
 elid  = necpt(1)
 
!     SET UP FOR FETCHING SHAPE FUNCTIONS
 
 DO  i = 1,np
   gpt(i) = ecpt(itemp-1+i)
   DO  j = 1,3
     bxyz(j,i) = ecpt(np+4+4*i+j)
   END DO
 END DO
 IF (onlyc) GO TO 570
 i = inip - 1
 SELECT CASE ( i )
   CASE (    1)
     GO TO 520
   CASE (    2)
     GO TO 530
   CASE (    3)
     GO TO 540
 END SELECT
 520 h(1)  = 1.
 s(1)  = gauss(1)
 h(2 ) = h(1)
 s(2)  = -s(1)
 GO TO 550
 530 h(1)  = gauss(2)
 s(1)  = gauss(3)
 h(2 ) = gauss(4)
 s(2)  = 0.
 h(3 ) = h(1)
 s(3 ) = -s(1)
 GO TO 550
 540 h(1 ) = gauss(5)
 s(1)  = gauss(6)
 h(2)  = gauss(7)
 s(2)  = gauss(8)
 h(3)  = h(2)
 s(3)  = -s(2)
 h(4) = h(1)
 s(4) = -s(1)
 550 DO  i = 1,32
   f(i) = 0.0
 END DO
 
!     SET UP HC ARRAY GIVING HC AT EACH GRID
 
 IF (jtype /= 24) GO TO 570
 
!     REMFLUX
 
 isub  = istart + 3*ncount - 3
 gh(1) = z(isub)
 gh(2) = z(isub+1)
 gh(3) = z(isub+2)
 GO TO 610
 
!     IF SPCFLD,PICK UP GRID VALUES HERE. IF NOT, PICK UP INTEGRATION
!     POINT VALUES LATER.(THERE IS ONLY ONE SPCFLD CARD AT THIS POINT)
 
 570 IF (jtype /= 20) GO TO 590
 DO  i = 1,np
   isil = 3*necpt(i+1)
   hc(3*i-2) = z(istart+isil-3)
   hc(3*i-1) = z(istart+isil-2)
   hc(3*i  ) = z(istart+isil-1)
 END DO
 590 inflag = 3
 matid  = necpt(mid)
 610 ktype  = jtype - 20
 IF (onlyc) GO TO 850
 
!     START INTEGRATION
 
 DO  i = 1,inip
   DO  j = 1,inip
     DO  k = 1,inip
       
!     FETCH SHAPE FUNCTIONS FOR THIS INTEGRATION POINT
       
       CALL ihexss(itype,shp,dshp,xjacob,detj,elid,s(i),s(j),s(k),bxyz)
       
!     COMPUTE NI W.R.T. X,Y,Z(REVERVSE CALLING SEQUENCE,SINCE COL STOR)
       
       CALL gmmats(dshp,np,3,0,xjacob,3,3,0,dshpb)
       
!     COMPUTE TEMPERATURES AND HC  AT THIS INTEGRSTION POINT
       
       eltemp = 0
       DO  l = 1,np
         eltemp = eltemp + shp(l)*gpt(l)
       END DO
       IF (jtype == 24) GO TO 730
       hcxyz(1) = 0.
       hcxyz(2) = 0.
       hcxyz(3) = 0.
       IF (jtype == 20) GO TO 700
       
!     FOR LOOPS AND DIPOLES, COMPUTE BASIC COORDS FOR THIS INTEGRATION
!     POINT
       
       xx = 0.
       yy = 0.
       zz = 0.
       DO  l = 1,np
         xx = xx + shp(l)*bxyz(1,l)
         yy = yy + shp(l)*bxyz(2,l)
         zz = zz + shp(l)*bxyz(3,l)
       END DO
       DO  ijk = 1,ido
         isub = istart + (ijk-1)*iwords - 1
         DO  l = 1,iwords
           buf(l) = z(isub+l)
         END DO
         
!     COMPUTE HC AT THIS POINT DUE TO ALL LOADS OF PRESENT TYPE
         
         SELECT CASE ( ktype )
           CASE (    1)
             GO TO 650
           CASE (    2)
             GO TO 660
           CASE (    3)
             GO TO 670
         END SELECT
         650 CALL axloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
         GO TO 680
         660 CALL geloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
         GO TO 680
         670 CALL dipole (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
         680 hcxyz(1) = hcxyz(1) + hc1
         hcxyz(2) = hcxyz(2) + hc2
         hcxyz(3) = hcxyz(3) + hc3
       END DO
       GO TO 720
       
!     SPCFLD
       
       700 DO   l = 1,np
         hcxyz(1) = hcxyz(1) + shp(l)*hc(3*l-2)
         hcxyz(2) = hcxyz(2) + shp(l)*hc(3*l-1)
         hcxyz(3) = hcxyz(3) + shp(l)*hc(3*l  )
       END DO
       
       CALL hmat(elid)
       
       720 g(1) = xmat(1)
       g(2) = xmat(2)
       g(3) = xmat(3)
       g(4) = xmat(2)
       g(5) = xmat(4)
       g(6) = xmat(5)
       g(7) = xmat(3)
       g(8) = xmat(5)
       g(9) = xmat(6)
       
       CALL gmmats (g,3,3,0,hcxyz,3,1,0,gh)
       
       730 sfact = h(i)*h(j)*h(k)*detj
       DO  l = 1,np
         f(l) = f(l) + (dshpb(1,l)*gh(1) + dshpb(2,l)*gh(2) +  &
             dshpb(3,l)*gh(3))*sfact
       END DO
       
!     GET ANOTHER INTEGRATIONPOINT
       
     END DO
   END DO
 END DO
 
!     ADD LOADS INTO LOAD ARRAY
 
 DO  l = 1,np
   isil = necpt(frstgd+l-1)
   
!     IF PERMBDY EXISTS AND IF GRID IS NOT ON IT, IGNORE ITS LOAD
   
   IF (nbdys == 0) GO TO 830
   DO  i = 1,nbdys
     IF (isil /= iz(istart-nbdys-nelout+i-1)) CYCLE
     GO TO 830
   END DO
   CYCLE
   830 z(isil) = z(isil) - f(l)
 END DO
 850 itype = jtype
 
!     BEFORE LEAVING, WE MUST COMPUTE HC VALUES AT GRIDS OF ISOPARA-
!     METRICS AND WRITE TO SCR6
 
 IF (jtype == 24) GO TO 1150
 CALL WRITE (scr6,isc,2,0)
 IF (jtype /= 20) GO TO 1010
 
!     FOR SPCFLD THE VALUES ARE IN CORE(EXCEPT FOR MIFPOINTS OF IHEX3)
 
 IF (eltype == 67) GO TO 880
 CALL WRITE (scr6,hc,3*np,0)
 
!     CENTROID/ XI = ETA = ZETA = 0
 
 860 CALL ihexss (eltype-64,shp,dshp,xjacob,detj,elid,0.,0.,0.,bxyz)
 hcx3(1) = 0.
 hcx3(2) = 0.
 hcx3(3) = 0.
 DO  l = 1,np
   hcx3(1) = hcx3(1) + shp(l)*hc(3*l-2)
   hcx3(2) = hcx3(2) + shp(l)*hc(3*l-1)
   hcx3(3) = hcx3(3) + shp(l)*hc(3*l  )
 END DO
 CALL WRITE (scr6,hcx3,3,0)
 GO TO 1150
 880 isub1 = 1
 isub2 = 10
 j = -5
 890 DO  i = isub1,isub2,3
   j = j + 6
   k = 3*i - 2
   hcx3(j  ) = hc(k  )
   hcx3(j+1) = hc(k+1)
   hcx3(j+2) = hc(k+2)
   hcx3(j+3) = .5*(hc(k+3) + hc(k+6))
   hcx3(j+4) = .5*(hc(k+4) + hc(k+7))
   hcx3(j+5) = .5*(hc(k+5) + hc(k+8))
 END DO
 IF (isub1 == 21) GO TO 1000
 j = 22
 DO  i = 13,16
   j = j + 3
   k = 3*i - 2
   hcx3(j  ) = .5*(hc(k  ) + hc(k+12))
   hcx3(j+1) = .5*(hc(k+1) + hc(k+13))
   hcx3(j+2) = .5*(hc(k+2) + hc(k+14))
 END DO
 isub1 = 21
 isub2 = 30
 j = 31
 GO TO 890
 
!     DONE - WRITE RESULTS
 
 1000 CALL WRITE (scr6,hcx3,60,0)
 GO TO 860
 
!     CEMLOOP, GEMLOOP, MDIPOLE
 
 1010 nx = np + 1
 IF (eltype == 67) nx = 21
 DO  j = 1,nx
   IF (j /= nx) GO TO 1030
   
!     CENTROID
   
   CALL ihexss (eltype-64,shp,dshp,xjacob,detj,elid,0.,0.,0.,bxyz)
   xx = 0.
   yy = 0.
   zz = 0.
   DO  l = 1,np
     xx = xx + shp(l)*bxyz(1,l)
     yy = yy + shp(l)*bxyz(2,l)
     zz = zz + shp(l)*bxyz(3,l)
   END DO
   GO TO 1070
   
   1030 IF (eltype /= 67) GO TO 1060
   
!     IHEX3
   
   IF (j == 1) k1 =-1
   IF (j == 13) k1 = 7
   IF (j == 1) k2 =-1
   IF (j == 13) k2 = 7
   IF (j < 9 .OR. j > 12) GO TO 1040
   xx = .5*(bxyz(1,j+4) + bxyz(1,j+8))
   yy = .5*(bxyz(2,j+4) + bxyz(2,j+8))
   zz = .5*(bxyz(3,j+4) + bxyz(3,j+8))
   GO TO 1070
   1040 IF ((j/2)*2 /= j)GO TO 1050
   k1 = k1 + 1
   xx = .5*(bxyz(1,j+k1) + bxyz(1,j+k1+1))
   yy = .5*(bxyz(2,j+k1) + bxyz(2,j+k1+1))
   zz = .5*(bxyz(3,j+k1) + bxyz(3,j+k1+1))
   GO TO 1070
   1050 k2 = k2 + 1
   xx = bxyz(1,j+k2)
   yy = bxyz(2,j+k2)
   zz = bxyz(3,j+k2)
   GO TO 1070
   
   1060 xx = bxyz(1,j)
   yy = bxyz(2,j)
   zz = bxyz(3,j)
   1070 hc(1) = 0.
   hc(2) = 0.
   hc(3) = 0.
   DO  ijk = 1,ido
     isub = istart + (ijk-1)*iwords - 1
     DO  i = 1,iwords
       buf(i) = z(isub+i)
     END DO
     
!     COMPUTE HC AT THIS POINT
     
     SELECT CASE ( ktype )
       CASE (    1)
         GO TO 1090
       CASE (    2)
         GO TO 1100
       CASE (    3)
         GO TO 1110
     END SELECT
     1090 CALL axloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 1120
     1100 CALL geloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 1120
     1110 CALL dipole (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     1120 hc(1) = hc(1) + hc1
     hc(2) = hc(2) + hc2
     hc(3) = hc(3) + hc3
   END DO
   
   CALL WRITE (scr6,hc,3,0)
 END DO
 
 1150 RETURN
 
 1200 WRITE  (outpt,1210) ufm
 1210 FORMAT (a23,' - WRONG ELEMENT TYPE IN EM3D PROBLEM.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE em3d
