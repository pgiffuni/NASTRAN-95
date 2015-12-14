SUBROUTINE em2d (itype,istart,jtype,ncount,ido,iwords,nbdys,all,  &
        nelout)
     
!     COMPUTES ADDITIONAL E AND M LOADS FOR TWO DIMENSIONAL ELEMENTS
 
!     THIS ROUTINE HANDLES THE FOLLOWING 2-D ELEMENTS
 
!     TRIA1 -6-   TRMEM -9-   QDMEM-16-  TRIA2-17-  QUAD2-18-  QUAD1-19-
!     TRIARG-36-  TRAPRG-37   IS2D8-80-
 
 
 INTEGER, INTENT(IN)                      :: itype
 INTEGER, INTENT(IN)                      :: istart
 INTEGER, INTENT(IN)                      :: jtype
 INTEGER, INTENT(IN)                      :: ncount
 INTEGER, INTENT(IN)                      :: ido
 INTEGER, INTENT(IN)                      :: iwords
 INTEGER, INTENT(IN)                      :: nbdys
 INTEGER, INTENT(IN OUT)                  :: all
 INTEGER, INTENT(IN)                      :: nelout
 LOGICAL :: onlyc
 INTEGER :: otpe, pointr(9,9),typold,scr6
 REAL :: l(3,4),w(4)
 DIMENSION       buf(50),jbuf(50),xlacc(3),iz(1),nam(2),necpt(10),  &
     r(3,8),ip(3),hc(3),xload(3),d12(3),d13(3),xn(18),  &
     g(9),dxx(3),zi(3),zj(3),zk(3),et(9),xng(9),hcx(3),  &
     hcy(3),hcz(3),isc(5),sc(5),pt(3),h(3),z14(3),  &
     xz(16),vec(3),vvec(3),hci(24),f(8),gh(3),dn(8),  &
     dnxi(1),dneta(1),dnc(16),dnl(16),dnx(1),dny(1),  &
     xi(8),eta(8),xjb(4),xxjb(2,2),iws(2,3),hcxyz(3), ddnl(24),ddnlb(24)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,otpe,idum(78)
 COMMON /BLANK / nrowsp
 COMMON /emecpt/ ecpt(200)
 COMMON /zzzzzz/ z(1)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /hmtout/ xmat(6)
 EQUIVALENCE     (buf(1),jbuf(1)),(sc(1),isc(1)),(z(1),iz(1)),  &
     (ecpt(1),necpt(1)),(i1,ip(1)),(i2,ip(2)),  &
     (i3,ip(3)),(dnc(1),dnxi(1)),(dnc(9),dneta(1)),  &
     (dnl(1),dnx(1)),(dnl(9),dny(1))
 DATA    xi    / -1., 1., 1.,-1., 0., 1., 0.,-1./
 DATA    eta   / -1.,-1., 1., 1.,-1., 0., 1., 0./
 DATA    twopi3/ 2.094395103  /
 DATA    nam   / 4HEM2D,4H    /
 DATA    typold/ 0 /,   scr6  / 306/
 
!     EST STARTING POINTERS
 
!     ISIL   = 1ST SIL NUMBER
!     ITH    = MATERIAL ANGLE
!     MID    = MATERIAL ID
!     IA     = AREA FACTOR (TO COMPUTE VOLUME)
!     ISYS   = 1ST OUTPUT CORRDINATE SYSYTEM NUMBER
!     NGRIDS = NUMBER OF GRID POINTS
!     ITEMP  = ELEMENT TEMPERATURE
!     NEL    = NUMBER OF TRIANGLES USED TO FORM ELEMENT
 
!              ITYPE ISIL ITH MID IA ISYS NGRIDS ITEMP NEL
 
 DATA    pointr/  6,    2,  5,  6,  7, 15,    3,    27,  1,  &
     9,    2,  5,  6,  7,  9,    3,    21,  1,  &
     16,    2,  6,  7,  8, 10,    4,    26,  4,  &
     17,    2,  5,  6,  7,  9,    3,    21,  1,  &
     18,    2,  6,  7,  8, 10,    4,    26,  4,  &
     19,    2,  6,  7,  8, 16,    4,    32,  4,  &
     36,    2,  5,  6,  0,  7,    3,    19,  1,  &
     37,    2,  6,  7,  0,  8,    4,    24,  4,  &
     80,    2, 11, 12, 13, 14,    8,    46,  1/
 
 onlyc  = .false.
 IF (itype == 80) GO TO 10
 l(1,1) = 1./3.
 l(2,1) = l(1,1)
 l(3,1) = l(1,1)
 l(1,2) = .6
 l(2,2) = .2
 l(3,2) = .2
 l(1,3) = .2
 l(2,3) = .6
 l(3,3) = .2
 l(1,4) = .2
 l(2,4) = .2
 l(3,4) = .6
 w(1)   =-27./48.
 w(2)   = 25./48.
 w(3)   = w(2)
 w(4)   = w(2)
 nopts  = 4
 10 CONTINUE
 isc(1) = necpt(1)
 isc(2) = 1
 IF (itype == 80) isc(2) = 9
 
!     FIND ELEMENT TYPE TO PICK UP POINTERS
 
 IF (itype == typold) GO TO 40
 typold = itype
 DO  i = 1,9
   jel = i
   IF (itype-pointr(1,i) < 0) THEN
     GO TO  1600
   ELSE IF (itype-pointr(1,i) == 0) THEN
     GO TO    30
   ELSE
     GO TO    20
   END IF
 END DO
 GO TO 1600
 
 30 isil  = pointr(2,jel)
 ith   = pointr(3,jel)
 mid   = pointr(4,jel)
 ia    = pointr(5,jel)
 isys  = pointr(6,jel)
 ngrids= pointr(7,jel)
 itemp = pointr(8,jel)
 nel   = pointr(9,jel)
 
!     CHECK TO SEE IF THIS ELEMENT CONTAINS A GRID POINT ON A PERMBDY
!     CARD. IF SO, OR IF NO PERMBDY CARD EXISTS, COMPUTE LOADS FOR THE
!     ELEMENT. IF NOT, COMPUTE HC CENTROIDAL VALUE ONLY. (ONLYC=.TRUE.)
!     THE PERMBDY SILS START AT Z(ISTART-NBDYS-1)
 
 40 IF (nbdys == 0) GO TO 60
 
 DO  i = 1,ngrids
   ng = necpt(isil+i-1)
   DO  j = 1,nbdys
     IF (ng == iz(istart-nbdys-nelout+j-1)) GO TO 60
   END DO
 END DO
 
!     ELEMENT HAS NO GRIDS ON PERMBDY
 
 onlyc = .true.
 nopts = 0
 60 IF (onlyc .AND. jtype == 24) RETURN
 
!     IF ONLYC=TRUE, CHECK TO SEE IF THE ELEMENT HAD AN ELFORCE REQUEST.
!     IF SO, CONTINUE. IF NOT, JUST WRITE ZEROS TO HCCEN,SCR6) AND
!     RETURN.
 
 IF(.NOT.onlyc) GO TO 80
 IF(all == 1) GO TO 80
 IF(nelout == 0) GO TO 110
 
 DO  i = 1,nelout
   IF (necpt(1) == iz(istart-nelout+i-1)) GO TO 80
 END DO
 GO TO 110
 
!     CHECK FOR ZERO LOAD
 
 80 IF (jtype /= 20 .AND. jtype /= 24) GO TO 210
 h1 = 0.
 h2 = 0.
 h3 = 0.
 g1 = 0.
 g2 = 0.
 g3 = 0.
 DO  i = 1,ngrids
   isub = istart + 3*necpt(isil+i-1) - 3
   IF (jtype == 24) isub = istart + 3*ncount - 3
   h1 = h1 + ABS(z(isub  ))
   h2 = h2 + ABS(z(isub+1))
   h3 = h3 + ABS(z(isub+2))
   g1 = g1 + z(isub  )
   g2 = g2 + z(isub+1)
   g3 = g3 + z(isub+2)
   IF (jtype == 24) EXIT
 END DO
 100 hl = h1 + h2 + h3
 IF (hl /= 0.) GO TO 200
 IF (jtype == 24) RETURN
 
!     ALL ZEROS - WRITE TO SCR6
 
 110 sc(3) = 0.
 sc(4) = 0.
 sc(5) = 0.
 CALL WRITE (scr6,sc,2,0)
 isc2  = isc(2)
 DO  i = 1,isc2
   CALL WRITE (scr6,sc(3),3,0)
 END DO
 RETURN
 
 200 IF (jtype == 24) GO TO 210
 
!     AVERAGE SPCFLD
 
 ahcx = g1/FLOAT(ngrids)
 ahcy = g2/FLOAT(ngrids)
 ahcz = g3/FLOAT(ngrids)
 
 210 IF (onlyc) GO TO 310
 
!     PICK UP MATERIAL INFO
!     INFLAG = 3 MEANS A 3 X 3 MATERIAL MATRIX WILL BE RETURNED. THE
!     REASON FOR DOING THIS FOR A 2-D ELEMENT IS THAT HC CAN HAVE A
!     COMPONENT NORMAL TO THE PLANE OF THE ELEMENT. PARTIAL DERIVATIVE
!     W.R.T Z IS 0.  BUT IF THE MATERIAL IS ANISOTROPIC, THEN A
!     CONTRIBUTION TO THE SCALAR LOAD IS POSSIBLE IF MATERIAL CONTAINS
!     A NON-ZERO X-Z TERM. FOR ISOTROPIC MATERIALS, THE NORMAL COMPONENT
!     OF HC WILL BE IGNORED W.R.T ITS CONTRIBUTION TO THE LOAD. IF ALL
!     TERMS OF MATERIAL MATRIX W.R.T.Z ARE 0, AND IF ANISOTROPIC ANGLE
!     IS NOT 0, THEN WE MUST TRANSFORM MATERIALS TO ELEMENT SYSTEM HERE.
 
 inflag = 3
 IF (jtype == 24) GO TO 260
 matid  = necpt(mid)
 eltemp = ecpt(itemp)
 angle  = ecpt(ith)*0.017453293
 sinth  = SIN(angle)
 costh  = COS(angle)
 CALL hmat (necpt(1))
 
!     CHECK FOR 3-D ANISOTROPY
 
 IF (xmat(3) == 0. .AND. xmat(5) == 0.) GO TO 230
 
 220 g(1) = xmat(1)
 g(2) = xmat(2)
 g(3) = xmat(3)
 g(5) = xmat(4)
 g(6) = xmat(5)
 g(9) = xmat(6)
 GO TO 240
 
!     CHECK FOR 2-D ANISOTROPY
 
 230 IF (ABS(angle) <= .0001) GO TO 220
 
!     2-D ANISOTROPY
 
 csq  = costh*costh
 ssq  = sinth*sinth
 cs   = costh*sinth
 g(1) = csq*xmat(1) - 2.*cs*xmat(2) + ssq*xmat(4)
 g(2) = cs*(xmat(1) - xmat(4)) + (csq-ssq)*xmat(2)
 g(3) = 0.
 g(5) = ssq*xmat(1) + 2.*cs*xmat(2) + csq*xmat(4)
 g(6) = 0.
 g(9) = xmat(6)
 
 240 IF (itype /= 36 .AND. itype /= 37) GO TO 250
 
!     SWITCH Y-Z MATERIALS FOR TRAPRG AND TRIARG
 
 temp = g(5)
 g(5) = g(9)
 g(9) = temp
 temp = g(2)
 g(2) = g(3)
 g(3) = temp
 
!     FILL IN SYMMETRIC PART
 
 250 g(4) = g(2)
 g(7) = g(3)
 g(8) = g(6)
 
!     SINCE QUADRILATERALS ARE COVERED BY 4 OVERLAPPING TRIANGLES,
!     MUST DIVIDE QUAD RESULTS BY 2
 
 260 xmul = 1.
 IF (ngrids == 4) xmul = .5
 
!     PICK UP COORDINATES OF GRID POINTS
 
 DO  i = 1,ngrids
   isubi = isys + 4*i - 4
   DO  j = 1,3
     isub  = isubi + j
     r(j,i)= ecpt(isub)
   END DO
 END DO
 310 IF (itype == 80) GO TO 900
 
!     COMPUTE COORDINATES OF CENTROID (OR, AT LEAST, AVERAGE ELEMENT
!     COORDS)
 
 xxc = 0.
 yyc = 0.
 zzc = 0.
 DO  i = 1,ngrids
   xxc = xxc + r(1,i)
   yyc = yyc + r(2,i)
   zzc = zzc + r(3,i)
 END DO
 xxc = xxc/FLOAT(ngrids)
 yyc = yyc/FLOAT(ngrids)
 zzc = zzc/FLOAT(ngrids)
 
!     NOW COMPUTE PROPER LOADS FOR EACH TRIANGLE
 
 DO  iel = 1,nel
   IF (onlyc) GO TO 500
   
!     1ST SET UP AN ARRAY TO PICK UP GRID POINTS IN A PARTICULAR ORDER.
!     FOR TRIANGLES, IT IS 1,2,3. FOR QUADRILATERALS, FORM 4 TRIANGLES
!     BY TAKING GRIDS 1,2,3, 2,3,4, 3,4,1, AND 4,1,2
   
   DO  i = 1,3
     ip(i) = i + iel - 1
     IF (ip(i) > 4) ip(i) = ip(i) - 4
   END DO
   
!     COMPUTE VECTORS FROM 1ST GRID TO 2ND AND FROM 1ST TO 3RD
   
   DO  i = 1,3
     d12(i) = r(i,i2) - r(i,i1)
     d13(i) = r(i,i3) - r(i,i1)
   END DO
   
!     SET UP GRADIENTS FOR AXISYMMETRIC ELEMENTS SEPARATELY
   
   IF (itype /= 36 .AND. itype /= 37) GO TO 360
   
!     THE LENGTH OF THE CROSS PRODUCT VECTOR IS TWICE THE AREA OF THE
!     TRIANG
   
   CALL saxb (d12(1),d13(1),d12(1))
   area = .5*SQRT(d12(1)**2 + d12(2)**2 + d12(3)**2)
   vol  = area*twopi3*(r(1,i1) + r(1,i2) + r(1,i3))
   
!     NOW SET UP GRADIENT OF THE SHAPE FUNCTION AT EACH GRID POINT.
!     SET UP A 3 X3 MATRIX ROW-STORED FOR GMMATS
   
   d     = (r(1,i2)-r(1,i1))*r(3,i3) + (r(1,i1)-r(1,i3))*r(3,i2) +  &
       (r(1,i3)-r(1,i2))*r(3,i1)
   xn(1) = r(3,i2) - r(3,i3)
   xn(2) = 0.
   xn(3) = r(1,i3) - r(1,i2)
   xn(4) = r(3,i3) - r(3,i1)
   xn(5) = 0.
   xn(6) = r(1,i1) - r(1,i3)
   xn(7) = r(3,i1) - r(3,i2)
   xn(8) = 0.
   xn(9) = r(1,i2) - r(1,i1)
   
   DO  i = 1,9
     xn(i) = xn(i)/d
   END DO
   
!     FOR ALL EXCEPT REMFLUX, MULT. GRADIENTS INTO MATERIALS
   
   IF (jtype /= 24) CALL gmmats (xn(1),3,3,0,g,3,3,0,xn(10))
   GO TO 420
   
!     FIRST, CONVERT COORDINATES TO ELEMNT COORDINATE SYSTEM
   
   360 zlen = SQRT(d12(1)**2 + d12(2)**2 + d12(3)**2)
   DO  i = 1,3
     zi(i) = d12(i)/zlen
   END DO
   
   CALL saxb (zi(1),d13(1),dxx(1))
   
   x2 = zlen
   x3 = d13(1)*zi(1) + d13(2)*zi(2) + d13(3)*zi(3)
   y3 = SQRT(dxx(1)**2 + dxx(2)**2 + dxx(3)**2)
   
   area = .5*x2*y3
   vol  = area*ecpt(ia)
   
!     GET J AND K VECTORS FOR LATER USE
   
   DO  i = 1,3
     zk(i) = dxx(i)/y3
   END DO
   
   CALL saxb (zk(1),zi(1),zj(1))
   zlen = SQRT(zj(1)**2 + zj(2)**2 + zj(3)**2)
   DO  i = 1,3
     zj(i) = zj(i)/zlen
   END DO
   DO  i = 1,3
     et(i  ) = zi(i)
     et(i+3) = zj(i)
     et(i+6) = zk(i)
   END DO
   
!     SHAPE FUNCTION GRADIENTS
   
   xn(1) = -1./x2
   xn(2) = (x3-x2)/(x2*y3)
   xn(3) = 0.
   xn(4) = -xn(1)
   xn(5) = -x3/(x2*y3)
   xn(6) = 0.
   xn(7) = 0.
   xn(8) = 1./y3
   xn(9) = 0.
   
!     TRANSFORM SHAPE FN GRADIENTS FROM LOCAL TO BASIC
   
   CALL gmmats (et,3,3,1,xn(1),3,3,1,xng(1))
   
!     FOR ALL EXCEPT REMFLUX, MULT. GRADIENTS OF SHAPE FNS INTO
!     MATERIALS
   
   IF (jtype == 24) GO TO 410
   CALL gmmats (xng(1),3,3,1,g,3,3,0,xn(10))
   GO TO 420
   410 xn(1) = xng(1)
   xn(2) = xng(4)
   xn(3) = xng(7)
   xn(4) = xng(2)
   xn(5) = xng(5)
   xn(6) = xng(8)
   xn(7) = xng(3)
   xn(8) = xng(6)
   xn(9) = xng(9)
   420 IF (jtype == 24) GO TO 740
   
!     START INTEGRATION PROCEDURE- 4 POINTS FOR CUBIC PLUS ONE AT
!     CENTROID
   
   500 ktype    = jtype - 19
   xlacc(1) = 0.
   xlacc(2) = 0.
   xlacc(3) = 0.
   noptsp   = nopts + 1
   DO  npts = 1,noptsp
     
!     DO CENTROID FOR ONLY 1ST TRIANGLE
     
     IF (npts == noptsp .AND. iel > 1) CYCLE
     
!     COMPUTE BASIC COORDS OF INTEGRATION POINT
     
     IF (npts /= noptsp) GO TO 510
     
!     CENTROID
     
     xx = xxc
     yy = yyc
     zz = zzc
     IF (jtype /= 20) GO TO 520
     
!     AVERAGE SPCFLD
     
     hc(1) = ahcx
     hc(2) = ahcy
     hc(3) = ahcz
     GO TO 610
     510 xx    = l(1,npts)*r(1,i1) + l(2,npts)*r(1,i2) + l(3,npts)*r(1,i3)
     yy    = l(1,npts)*r(2,i1) + l(2,npts)*r(2,i2) + l(3,npts)*r(2,i3)
     zz    = l(1,npts)*r(3,i1) + l(2,npts)*r(3,i2) + l(3,npts)*r(3,i3)
     520 hc(1) = 0.
     hc(2) = 0.
     hc(3) = 0.
     
!     COMPUTE HC AT THIS POINT FOR ALL LOADS OF THIS TYPE
     
     DO  ijk = 1,ido
       IF (jtype == 20) GO TO 540
       isub = istart + (ijk-1)*iwords - 1
       DO  i = 1,iwords
         buf(i) = z(isub+i)
       END DO
       SELECT CASE ( ktype )
         CASE (    1)
           GO TO 540
         CASE (    2)
           GO TO 560
         CASE (    3)
           GO TO 570
         CASE (    4)
           GO TO 580
       END SELECT
       540 DO  i = 1,3
         ipi  = ip(i)
         nsil = necpt(isil+ipi-1)
         ipt  = istart + 3*nsil - 3
         hcx(i) = z(ipt  )
         hcy(i) = z(ipt+1)
         hcz(i) = z(ipt+2)
       END DO
       hc1 = l(1,npts)*hcx(1) + l(2,npts)*hcx(2) + l(3,npts)*hcx(3)
       hc2 = l(1,npts)*hcy(1) + l(2,npts)*hcy(2) + l(3,npts)*hcy(3)
       hc3 = l(1,npts)*hcz(1) + l(2,npts)*hcz(2) + l(3,npts)*hcz(3)
       GO TO 590
       
!     CEMLOOP, GEMLOOP, MDIPOLE
       
       560 CALL axloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 590
       570 CALL geloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 590
       580 CALL dipole (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       590 hc(1) = hc(1) + hc1
       hc(2) = hc(2) + hc2
       hc(3) = hc(3) + hc3
     END DO
     
     610 IF (npts /= noptsp) GO TO 700
     sc(3) = hc(1)
     sc(4) = hc(2)
     sc(5) = hc(3)
     CALL WRITE (scr6,sc,5,0)
     CYCLE
     
!     WE HAVE HC AT THIS INTEG. PT. MULT. BY WEIGHT AND ACCUMULATE
     
     700 DO  i = 1,3
       xlacc(i) = xlacc(i) + hc(i)*w(npts)
     END DO
     
!     GET ANOTHER INTEGRATION POINT
     
   END DO
   
   IF (onlyc) RETURN
   DO  i = 1,3
     hc(i) = xlacc(i)
   END DO
   GO TO 750
   
!     REMFLUX
   
   740 ipt   = istart + 3*ncount - 3
   hc(1) = z(ipt  )
   hc(2) = z(ipt+1)
   hc(3) = z(ipt+2)
   
!    TAKE XMUL MULTIPLIER INTO ACCOUNT
   
   750 DO  i = 1,3
     hc(i) = hc(i)*xmul
   END DO
   
!     MAKE FINAL COMPUTATION. MULTIPLY PRODUCT OF SHAPE FUNCTION
!     GRADIENTS AND MATERIAL MATRIX INTO HC AND MULTIPLY BY VOLUME
   
   isub = 10
   IF (jtype == 24) isub = 1
   CALL gmmats (xn(isub),3,3,0,hc,3,1,0,xload(1))
   
!     ADD THIS ELEMENT LOAD VECTOR IN OVERALL VECTOR. USE NSIL AND IP TO
!     POI
   
   DO  j = 1,3
     ipi  = ip(j)
     nsil = necpt(isil+ipi-1)
     
!     IF PERMBDY EXISTS AND IF GRID IS NOT ON IT, IGNORE ITS LOAD
     
     IF (nbdys == 0) GO TO 780
     DO  i = 1,nbdys
       IF (nsil /= iz(istart-nbdys-nelout+i-1)) CYCLE
       GO TO 780
     END DO
     CYCLE
     780 z(nsil) = z(nsil) - xload(j)*vol
   END DO
   
!     DONE FOR THIS TRIANGLE. GO BACK FOR ANOTHER
   
 END DO
 RETURN
 
!     IS2D8
 
!     SET UP QUADRATURE POINTS AND WEIGHTS
 
 900 IF (onlyc) GO TO 1000
 pt(1) = -0.57735027
 pt(2) = -pt(1)
 h(1)  = 1.
 h(2)  = 1.
 IF (necpt(10) == 2) GO TO 910
 pt(1) = -0.77459667
 pt(2) = 0.
 pt(3) = -pt(1)
 h(1)  = 5./9.
 h(2)  = 8./9.
 h(3)  = h(1)
 
!     COMPUTE I,J,K VECTORS- I IS 1 TO 2
 
 910 DO  i = 1,3
   zi(i) = r(i,2) - r(i,1)
   z14(i)= r(i,4) - r(i,1)
 END DO
 zlen  = SQRT(zi(1)**2 + zi(2)**2 + zi(3)**2)
 DO  i = 1,3
   zi(i) = zi(i)/zlen
 END DO
 
!     GET K BY CROSSING I INTO VECTOR FROM 1 TO 4
 
 zk(1) = zi(2)*z14(3) - zi(3)*z14(2)
 zk(2) = zi(3)*z14(1) - zi(1)*z14(3)
 zk(3) = zi(1)*z14(2) - zi(2)*z14(1)
 zklen = SQRT(zk(1)**2 + zk(2)**2 + zk(3)**2)
 DO  i = 1,3
   zk(i) = zk(i)/zklen
 END DO
 
!     GET J BY CROSSING K INTO I AND STORE INTO TRANSFORMATION MATRIX
 
 zj(1) = zk(2)*zi(3) - zk(3)*zi(2)
 zj(2) = zk(3)*zi(1) - zk(1)*zi(3)
 zj(3) = zk(1)*zi(2) - zk(2)*zi(1)
 zjlen = SQRT(zj(1)**2 + zj(2)**2 + zj(3)**2)
 DO  i = 1,3
   zj(i) = zj(i)/zjlen
 END DO
 
 DO  i = 1,3
   et(i  ) = zi(i)
   et(i+3) = zj(i)
   et(i+6) = zk(i)
 END DO
 
!     COMPUTE ELMENT COORDS FOR 1 AND 2
 
 xz(1) = 0.
 xz(2) = 0.
 xz(3) = zlen
 xz(4) = 0.
 
!     FOR 3-8, X IS DOT PRODUCT OF VECTOR FROM 1 TO GRID WITH I.
!     Y IS THE LENFTH OF THE VECTOR RESULTING FROM CROSSING I INTO
!     VECTOR FROM 1 TO GRID
 
 DO  i = 3,8
   ixx = 2*i - 1
   DO  j = 1,3
     vec(j)  = r(j,i) - r(j,1)
   END DO
   xz(ixx) = vec(1)*zi(1) + vec(2)*zi(2) + vec(3)*zi(3)
   vvec(1) = zi(2)*vec(3) - zi(3)*vec(2)
   vvec(2) = zi(3)*vec(1) - zi(1)*vec(3)
   vvec(3) = zi(1)*vec(2) - zi(2)*vec(1)
   xz(ixx+1) = SQRT(vvec(1)**2 + vvec(2)**2 + vvec(3)**2)
 END DO
 
 DO  i = 1,8
   f(i) = 0.
 END DO
 
!     GET HC AT EACH GRID
 
 IF (jtype /= 24) GO TO 1000
 
!     REMFLUX
 
 isub  = istart + 3*ncount - 3
 gh(1) = z(isub  )
 gh(2) = z(isub+1)
 gh(3) = z(isub+2)
 GO TO 1020
 
!     IF SPCFLD, PICK UP GRID VALUES HERE. IF NOT, PICK UP INTEGRATION
!     POINT VALUES LATER
 
 1000 IF (jtype /= 20) GO TO 1020
 DO  i = 1,ngrids
   isil = 3*necpt(i+1)
   hci(3*i-2) = z(istart+isil-3)
   hci(3*i-1) = z(istart+isil-2)
   hci(3*i  ) = z(istart+isil-1)
 END DO
 1020 inip  = necpt(10)
 ktype = jtype - 20
 IF (onlyc) GO TO 1340
 
!     START INTEGRATION
 
 DO  iii = 1,inip
   DO  jjj = 1,inip
     
!     COMPUTE DERIVATIVES WITH RESPECT TO XI AND ETA
!     EACH GRID POINT
     
     DO  n = 1,4
       dn(n)   = .25*(1.+pt(iii)*xi(n))*(1.+pt(jjj)*eta(n))*  &
           (pt(iii)*xi(n)+pt(jjj)*eta(n)-1.)
       dnxi(n) = .25*xi(n)*(1.+pt(jjj)*eta(n))*  &
           (2.*pt(iii)*xi(n)+pt(jjj)*eta(n))
       dneta(n)= .25*eta(n)*(1.+pt(iii)*xi(n))*  &
           (pt(iii)*xi(n)+2.*pt(jjj)*eta(n))
     END DO
     DO  n = 5,7,2
       
       dn(n)   = .5*(1.-pt(iii)*pt(iii))*(1.+pt(jjj)*eta(n))
       dnxi(n) = -pt(iii)*(1.+pt(jjj)*eta(n))
       dneta(n)= .5*(1.-pt(iii)*pt(iii))*eta(n)
     END DO
     
     DO  n = 6,8,2
       dn(n)   = .5*(1.+pt(iii)*xi(n))*(1.-pt(jjj)*pt(jjj))
       dnxi(n) = .5*xi(n)*(1.-pt(jjj)*pt(jjj))
       dneta(n)= -pt(jjj)*(1.+pt(iii)*xi(n))
     END DO
     
!     COMPUTE JACOBEAN
     
!           N1XI   N2XI   N3XI   N4XI   N5XI   N6XI   N7XI   N8XI
!     DNC = N1ETA  N2ETA  N3ETA  N4ETA  N5ETA  N6ETA  N7ETA  N8ETA
     
!          X1  Y1
!          X2  Y2
!          X3  Y3
!     XX = X4  Y4
!          X5  Y5
!          X6  Y6
!          X7  Y7
!          X8  Y8
     
     CALL gmmats (dnc,2,8,0,xz,8,2,0,xjb)
     
!     XJB IS ROW-STORED-IT MUST BE COLUMN-STORED AND DOUBLY DIMENSIONED
!     FOR INVERSION
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xxjb(i,j) = xjb(k)
       END DO
     END DO
     
!     COMPUTE INVERSE AND DETERMINANT OF JACOBEAN
     
     CALL invers (2,xxjb,2,dumarg,0,determ,ising,iws)
     
!     COMPUTE DERIVATIVES WITH RESPECT TO X AND Y
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xjb(k) = xxjb(i,j)
       END DO
     END DO
     CALL gmmats (xjb,2,2,0,dnc,2,8,0,dnl)
     
!           N1X N2X N3X N4X N5X N6X N7X N8X
!     DNL = N1Y N2Y N3Y N4Y N5Y N6Y N7Y N8Y
     
     IF (jtype == 24) GO TO 1190
     
!     INITIALIZE HC AT PRESENT UNTEGRATION POINT
     
     DO  i = 1,3
       hcxyz(i) = 0.
     END DO
     IF (jtype == 20) GO TO 1160
     
!     FOR LOOPS AND DIPOLES, COMPITE BASIC COORDS. FOR THIS INTEGRATION
!     PT
     
     xx = 0.
     yy = 0.
     zz = 0.
     DO  m = 1,ngrids
       xx = xx + dn(m)*r(1,m)
       yy = yy + dn(m)*r(2,m)
       zz = zz + dn(m)*r(3,m)
     END DO
     
     DO  ijk = 1,ido
       isub = istart + (ijk-1)*iwords - 1
       DO  m = 1,iwords
         buf(m) = z(isub+m)
       END DO
       SELECT CASE ( ktype )
         CASE (    1)
           GO TO 1110
         CASE (    2)
           GO TO 1120
         CASE (    3)
           GO TO 1130
       END SELECT
       1110 CALL axloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 1140
       1120 CALL geloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       GO TO 1140
       1130 CALL dipole (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
       1140 hcxyz(1) = hcxyz(1) + hc1
       hcxyz(2) = hcxyz(2) + hc2
       hcxyz(3) = hcxyz(3) + hc3
     END DO
     GO TO 1180
     
!     SPCFLD
     
     1160 DO  m = 1,ngrids
       hcxyz(1) = hcxyz(1) + dn(m)*hci(3*m-2)
       hcxyz(2) = hcxyz(2) + dn(m)*hci(3*m-1)
       hcxyz(3) = hcxyz(3) + dn(m)*hci(3*m)
     END DO
     
!     MULTIPLY MATERIAL INTO HC AT THIS INTEGRATION POINT
     
     1180 CALL gmmats (g,3,3,0,hcxyz,3,1,0,gh)
     1190 sfact = h(iii)*h(jjj)*determ
     
!     TRANSFORM DNL FROM LOCAL TO BASIC
!     1 ST EXPAND TO ADD IN ZEROS CORRESPONDING TO Z DIRECTION
     
     DO  i = 1,16
       ddnl(i) = dnl(i)
     END DO
     DO  i = 17,24
       ddnl(i) = 0.
     END DO
     
     CALL gmmats (et,3,3,1,ddnl,3,8,0,ddnlb)
     
     DO  m = 1,ngrids
       f(m) = f(m) + (ddnlb(m)*gh(1) + ddnlb(m+8)*gh(2) +  &
           ddnlb(m+16)*gh(3))*sfact
     END DO
     
!     GET ANOTHER INTEGRATION POINT
     
   END DO
 END DO
 
!     ADD LOAD INTO LOAD ARRAY
 
 DO  m = 1,ngrids
   isil = necpt(m+1)
   
!     IF PERMBDY EXISTS AND IF GRID IS NOT ON IT, IGNORE ITS LOAD
   
   IF (nbdys == 0) GO TO 1320
   DO  i = 1,nbdys
     IF (isil /= iz(istart-nbdys-nelout+i-1)) CYCLE
     GO TO 1320
   END DO
   CYCLE
   1320 z(isil) = z(isil)-f(m)*ecpt(ia)
 END DO
 
!     BEFORE LEAVING COMPUTE HC AT GRIDS AND CENTROID AND WRITE TO SCR6
 
 1340 IF (jtype == 24) RETURN
 CALL WRITE (scr6,isc,2,0)
 
!     SET UP SHAPE FUNCTIONS AT CENTROID
 
 DO  i = 1,4
   dn(i) = -.25
 END DO
 DO  i = 5,8
   dn(i) = .5
 END DO
 
 IF (jtype /= 20) GO TO 1400
 
!     FOR SPCFLD HC VALUES AT GRIDS ARE IN CORE
 
 CALL WRITE (scr6,hci,24,0)
 
 DO  i = 1,3
   hcxyz(i) = 0.
 END DO
 DO  m = 1,ngrids
   hcxyz(1) = hcxyz(1) + dn(m)*hci(3*m-2)
   hcxyz(2) = hcxyz(2) + dn(m)*hci(3*m-1)
   hcxyz(3) = hcxyz(3) + dn(m)*hci(3*m  )
 END DO
 
 CALL WRITE (scr6,hcxyz,3,0)
 RETURN
 
!     NOT SPCFLD
 
 1400 DO  j = 1,9
   IF (j /= 9) GO TO 1420
   
!     CENTROID
   
   xx = 0.
   yy = 0.
   zz = 0.
   DO  m = 1,8
     xx = xx + dn(m)*r(1,m)
     yy = yy + dn(m)*r(2,m)
     zz = zz + dn(m)*r(3,m)
   END DO
   GO TO 1430
   1420 xx = r(1,j)
   yy = r(2,j)
   zz = r(3,j)
   1430 hc(1) = 0.
   hc(2) = 0.
   hc(3) = 0.
   DO  ijk = 1,ido
     isub  = istart + (ijk-1)*iwords - 1
     DO  i = 1,iwords
       buf(i) = z(isub+i)
     END DO
     SELECT CASE ( ktype )
       CASE (    1)
         GO TO 1450
       CASE (    2)
         GO TO 1460
       CASE (    3)
         GO TO 1470
     END SELECT
     1450 CALL axloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 1480
     1460 CALL geloop (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 1480
     1470 CALL dipole (buf,jbuf,xx,yy,zz,hc1,hc2,hc3)
     1480 hc(1) = hc(1) + hc1
     hc(2) = hc(2) + hc2
     hc(3) = hc(3) + hc3
   END DO
   
   CALL WRITE (scr6,hc,3,0)
 END DO
 
 RETURN
 
 1600 WRITE  (otpe,1610) ufm,nam,itype
 1610 FORMAT (a23,', IN SUBROUTINE',2A4,' ELEMENT TYPE',i8,' IS NOT ',  &
     'LEGAL')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE em2d
