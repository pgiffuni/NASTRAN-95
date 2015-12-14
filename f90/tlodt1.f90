SUBROUTINE tlodt1 (treal,tint)
     
!     THERMAL LOAD VECTOR FOR TRPLT1 (HIGHER ORDER PLATE BENDING ELEMENT
 
!     ECPT ENTRIES
!     AS IN STIFFNESS ROUTINE KTRPL1
 
 
 REAL, INTENT(IN)                         :: treal(6)
 INTEGER, INTENT(IN OUT)                  :: tint(6)
 LOGICAL :: nogo,nots,uniben,unitem
 INTEGER :: xpower(20),ypower(20),xthk(10),ythk(10),pt(3), qt(3),sil(6),sil1
 REAL :: ivect,jvect,kvect
 DIMENSION       f(10,10),xc(6),yc(6),zc(6),qqq(20,20),qqinv(360),  &
     ts1(60),ts2(60),iest(42), trand(9),dd(3), ics(6),ge1(9),nam(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /ssgwrk/ x,y,z,dista,distb,distc,a1,a2,a3,b1,b2,b3,g1(3),  &
     d(3),e(18),ivect(3),jvect(3),kvect(3),cc(10),g(9),  &
     ptem(20),ptele(18),ptglb(36),psub(3),psubt(6),  &
     psubt1(6),ts6(40),NAME(2),INDEX(20,3),nl(6),tl(3), balotr(36)
 COMMON /system/ sysbuf,iout
 COMMON /trimex/ est(100)
 COMMON /zzzzzz/ pg(1)
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin / matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/ em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy,  &
     rj11,rj12,rj22
 
!     EQUIVALENCE IECPT WITH ECPT IN COMMON BLOCK /SMA1ET/ SINCE ECPT IS
!     A MIXED INTEGER AND REAL ARRAY
 
 EQUIVALENCE     (thk1,tmem1),(thk2,tmem3),(thk3,tmem5),  &
     (a,dista),(b,distb),(c,distc),(iest(1),est(1)),  &
     (c1,cc(1)),(c2,cc(2)),(c3,cc(3)),(c4,cc(4)),  &
     (c5,cc(5)),(c6,cc(6)),(c7,cc(7)),(c8,cc(8)),  &
     (c9,cc(9)),(c10,cc(10)),(d(1),d1),(d(2),d2), (d(3),d3),(dd(1),d(1))
 DATA    BLANK , nam  / 4H    , 4HTRPL, 4HT1    /
 DATA    xpower/ 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/
 DATA    ypower/ 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/
 DATA    xthk  / 0,1,0,2,1,0,3,2,1,0 /
 DATA    ythk  / 0,0,1,0,1,2,0,1,2,3 /
 DATA    pt    / 0,1,0 /,  qt / 0,0,1/
 DATA    degra / 0.0174532925 /
 
 
 nots  = .false.
 idele = iest(1)
 DO  i = 1,6
   nl(i) = iest(i+1)
 END DO
 thetam = est(8)
 matid1 = iest(9)
 tmem1  = (est(10)*12.0)**0.333333333333
 tmem3  = (est(11)*12.0)**0.333333333333
 tmem5  = (est(12)*12.0)**0.333333333333
 tshr1  = est(14)
 tshr3  = est(15)
 tshr5  = est(16)
 j      = 0
 DO  i = 24,44,4
   j      = j + 1
   ics(j) = iest(i)
   xc(j)  = est(i+1)
   yc(j)  = est(i+2)
   zc(j)  = est(i+3)
 END DO
 temp1  = treal(1)
 temp3  = treal(1)
 temp5  = treal(1)
 t1prim =-treal(2)
 t3prim =-treal(2)
 t5prim =-treal(2)
 
!     IF TMEM3 OR TMEM5 EQUAL TO ZERO OR BLANK,THEY WILL BE SET EQUAL TO
!     SO ALSO FOR TEMP3 AND TEMP5
 
 IF (tmem3 == 0.0 .OR. tmem3 == BLANK) tmem3  = tmem1
 IF (tmem5 == 0.0 .OR. tmem5 == BLANK) tmem5  = tmem1
 IF (temp3 == 0.0 .OR. temp3 == BLANK) temp3  = temp1
 IF (temp5 == 0.0 .OR. temp5 == BLANK) temp5  = temp1
 IF (t3prim == .0 .OR. t3prim == BLANK) t3prim = t1prim
 IF (t5prim == .0 .OR. t5prim == BLANK) t5prim = t1prim
 IF (tshr3 == 0.0 .OR. tshr3 == BLANK) tshr3  = tshr1
 IF (tshr5 == 0.0 .OR. tshr5 == BLANK) tshr5  = tshr1
 eltemp = est(48)
 avthk  = (tmem1+tmem3+tmem5)/3.0
 aviner = avthk**3/12.0
 IF (tshr1 == 0.0) nots = .true.
 theta1 = thetam*degra
 sinth  = SIN(theta1)
 costh  = COS(theta1)
 IF (ABS(sinth) <= 1.0E-06) sinth = 0.0
 
!     EVALUATE MATERIAL PROPERTIES
 
 matflg = 2
 matid = matid1
 CALL mat (idele)
 
 g(1) = em(1)
 g(2) = em(2)
 g(3) = em(3)
 g(4) = em(2)
 g(5) = em(4)
 g(6) = em(5)
 g(7) = em(3)
 g(8) = em(5)
 g(9) = em(6)
 
!     IF TINT(6).NE.1,G1 IS G AND T1PRIME IS ALPHA TIMES T1PRIME
!     IF TINT(6).EQ.1,G1 IS G TIMES ALPHA AND T1PRIME IS T1PRIME
 
 IF (tint(6) /= 1) GO TO 147
 
!     G1 IS G TIMES ALPHA
 
 CALL gmmats (g,3,3,0, alf,3,1,0, g1)
 GO TO 149
 147 CONTINUE
 DO  i = 1,9
   ge1(i) = g(i)*aviner
 END DO
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (3,ge1(1),3,ts1(1),0,determ,ising,INDEX)
 IF (ising == 2) GO TO 901
 CALL gmmats (ge1,3,3,0, treal(2),3,1,0, tl(1))
 
!     CALCULATIONS FOR THE TRIANGLE
 
 149 CALL trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,iest(1),nam )
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e( i) = 0.0
 END DO
 e( 1) = kvect(1)
 e( 4) = kvect(2)
 e( 7) = kvect(3)
 e(11) = ivect(1)
 e(14) = ivect(2)
 e(17) = ivect(3)
 e(12) = jvect(1)
 e(15) = jvect(2)
 e(18) = jvect(3)
 
!     EVALUATE CONSTANTS D1,D2,D3 IN THE LINEAR EQUATION FOR TEMPERATURE
!     GRADIENT VARIATION OVER THE ELEMENT
 
 CALL af (f,10,a,b,c,d1,d2,d3,thk1,thk2,thk3,1)
 unitem = .false.
 IF (ABS(d2) <= 1.0E-06 .AND. ABS(d3) <= 1.0E-06) unitem =.true.
 
 distab = dista + distb
 a1   = (thk1*dista+thk2*distb)/distab
 a2   = (thk2-thk1)/distab
 a3   = (thk3-a1)/distc
 a1sq = a1*a1
 a2sq = a2*a2
 a3sq = a3*a3
 c1   = a1sq*a1
 c2   = 3.0*a1sq*a2
 c3   = 3.0*a1sq*a3
 c4   = 3.0*a1*a2sq
 c5   = 6.0*a1*a2*a3
 c6   = 3.0*a3sq*a1
 c7   = a2sq*a2
 c8   = 3.0*a2sq*a3
 c9   = 3.0*a2*a3sq
 c10  = a3*a3sq
 CALL af (f,10,a,b,c,b1,b2,b3,tshr1,tshr3,tshr5,1)
 uniben =.false.
 IF (ABS(a2) <= 1.0E-06 .AND. ABS(a3) <= 1.0E-06) uniben =.true.
 
!     COMPUTE THE AREA INTEGRATION FUNCTION F
 
 CALL af (f,10,a,b,c,0,0,0,0,0,0,-1)
 
!     CALCULATIONS FOR QMATRIX (QQQ) AND ITS INVERSE
 
 DO  i = 1,400
   qqq(i,1) = 0.0
 END DO
 DO  i = 1,6
   i1 = (i-1)*3 + 1
   i2 = (i-1)*3 + 2
   i3 = (i-1)*3 + 3
   qqq(i1, 1) = 1.0
   qqq(i1, 2) = xc(i)
   qqq(i1, 3) = yc(i)
   qqq(i1, 4) = xc(i)*xc(i)
   qqq(i1, 5) = xc(i)*yc(i)
   qqq(i1, 6) = yc(i)*yc(i)
   qqq(i1, 7) = qqq(i1, 4)*xc(i)
   qqq(i1, 8) = qqq(i1, 4)*yc(i)
   qqq(i1, 9) = qqq(i1, 5)*yc(i)
   qqq(i1,10) = qqq(i1, 6)*yc(i)
   qqq(i1,11) = qqq(i1, 7)*xc(i)
   qqq(i1,12) = qqq(i1, 7)*yc(i)
   qqq(i1,13) = qqq(i1, 8)*yc(i)
   qqq(i1,14) = qqq(i1, 9)*yc(i)
   qqq(i1,15) = qqq(i1,10)*yc(i)
   qqq(i1,16) = qqq(i1,11)*xc(i)
   qqq(i1,17) = qqq(i1,12)*yc(i)
   qqq(i1,18) = qqq(i1,13)*yc(i)
   qqq(i1,19) = qqq(i1,14)*yc(i)
   qqq(i1,20) = qqq(i1,15)*yc(i)
   qqq(i2, 3) = 1.0
   qqq(i2, 5) = xc(i)
   qqq(i2, 6) = yc(i)*2.0
   qqq(i2, 8) = qqq(i1, 4)
   qqq(i2, 9) = qqq(i1, 5)*2.0
   qqq(i2,10) = qqq(i1, 6)*3.0
   qqq(i2,12) = qqq(i1, 7)
   qqq(i2,13) = qqq(i1, 8)*2.0
   qqq(i2,14) = qqq(i1, 9)*3.0
   qqq(i2,15) = qqq(i1,10)*4.0
   qqq(i2,17) = qqq(i1,12)*2.0
   qqq(i2,18) = qqq(i1,13)*3.0
   qqq(i2,19) = qqq(i1,14)*4.0
   qqq(i2,20) = qqq(i1,15)*5.0
   qqq(i3, 2) =-1.0
   qqq(i3, 4) =-2.0*xc(i)
   qqq(i3, 5) =-yc(i)
   qqq(i3, 7) =-qqq(i1, 4)*3.0
   qqq(i3, 8) =-qqq(i1, 5)*2.0
   qqq(i3, 9) =-qqq(i1, 6)
   qqq(i3,11) =-qqq(i1, 7)*4.0
   qqq(i3,12) =-qqq(i1, 8)*3.0
   qqq(i3,13) =-qqq(i1, 9)*2.0
   qqq(i3,14) =-qqq(i1,10)
   qqq(i3,16) =-qqq(i1,11)*5.0
   qqq(i3,17) =-qqq(i1,13)*3.0
   qqq(i3,18) =-qqq(i1,14)*2.0
   qqq(i3,19) =-qqq(i1,15)
   
!     IF NO TRANSVERSE SHEAR GO TO 113
   
   IF (nots) CYCLE
   x = xc(i)
   y = yc(i)
   CALL tlodt3 (ts6,nots)
   DO  jj = 1,20
     qqq(i2,jj) = qqq(i2,jj) - ts6(20+jj)
     qqq(i3,jj) = qqq(i3,jj) + ts6(   jj)
   END DO
 END DO
 
 qqq(19,16) = 5.0*a**4*c
 qqq(19,17) = 3.0*a**2*c**3 - 2.0*a**4*c
 qqq(19,18) =-2.0*a*c**4 + 3.0*a**3*c**2
 qqq(19,19) = c**5 - 4.0*a**2*c**3
 qqq(19,20) = 5.0*a*c**4
 qqq(20,16) = 5.0*b**4*c
 qqq(20,17) = 3.0*b**2*c**3 - 2.0*b**4*c
 qqq(20,18) = 2.0*b*c**4 - 3.0*b**3*c**2
 qqq(20,19) = c**5 - 4.0*b**2*c**3
 qqq(20,20) =-5.0*b*c**4
 
!     FOURTH ARGUMENT IS A DUMMY LOCATION FOR INVERSE AND HENCE TS1(1)
!     IS U
 
!     AGAIN SET ISING = -1
 
 ising = -1
 CALL  invers (20,qqq,20,ts1(1),0,determ,ising,INDEX)
 
!     ISING EQUAL TO 2 IMPLIES THAT QQQ IS SINGULAR
 
!     FIRST 18 COLUMNS OF QQQ INVERSE IS THE QQQINV FOR USE IN STIFFNESS
!     MATRIX CALCULATIONS
 
 DO  i = 1,20
   DO  j = 1,18
     ij = (i-1)*18 + j
     qqinv (ij) = qqq(i,j)
   END DO
 END DO
 
 DO  i = 1,20
   mx   = xpower(i)
   rmx  = mx
   nx   = ypower(i)
   rnx  = nx
   rmnx = rmx*rnx
   rmx1 = rmx*(rmx-1.0D0)
   rnx1 = rnx*(rnx-1.0D0)
   ptemp= 0.0
   mx01 = mx - 1
   mx1  = mx + 1
   nx01 = nx - 1
   nx1  = nx + 1
   DO  k = 1,10
     mx01x= mx01+ xthk(k)
     nx1y = nx1 + ythk(k)
     mx1x = mx1 + xthk(k)
     nx01y= nx01+ ythk(k)
     mxx  = mx  + xthk(k)
     nxy  = nx  + ythk(k)
     IF (tint(6) /= 1) GO TO 213
     DO  l = 1,3
       mx01xp= mx01x+ pt(l)
       nx1yq = nx1y + qt(l)
       mx1xp = mx1x + pt(l)
       nx01yq= nx01y+ qt(l)
       mxxp  = mxx  + pt(l)
       nxyq  = nxy  + qt(l)
       IF (mx01xp > 0 .AND. nx1yq > 0)  &
           ptemp = ptemp+cc(k)*dd(l)*g1(1)*rmx1*f(mx01xp,nx1yq)
       IF (mx1xp > 0 .AND. nx01yq > 0)  &
           ptemp = ptemp+cc(k)*dd(l)*g1(2)*rnx1*f(mx1xp,nx01yq)
       IF (mxxp > 0 .AND. nxyq > 0)  &
           ptemp = ptemp+cc(k)*dd(l)*g1(3)*rmnx*f(mxxp,nxyq)
       IF (unitem) GO TO 213
     END DO
     
     213 IF (tint(6) == 1) GO TO 214
     IF (mx01x > 0) ptemp = ptemp + cc(k)*rmx1*(tl(1)*g(1)  &
         + tl(2)*g(2) + tl(3)*g(3))*f(mx01x,nx1y)
     IF (nx01y > 0) ptemp = ptemp + cc(k)*rnx1*(tl(1)*g(4)  &
         + tl(2)*g(5) + tl(3)*g(6))*f(mx1x,nx01y)
     IF (mxx > 0 .AND. nxy > 0) ptemp = ptemp + cc(k)*rmnx*  &
         (tl(1)*g(7)+tl(2)*g(8)+tl(3)*g(9))*f(mxx,nxy)
     214 IF (uniben) GO TO 216
   END DO
   
   216 ptem(i) = ptemp/12.0
 END DO
 
!     IF NO TRANSVERSE SHEAR GO TO 230
 
!     IF TSHR EQUAL TO ZERO OR MATID3 EQUAL TO ZERO, SKIP THESE
!     CALCULATIONS
 
 IF (nots) GO TO 230
 
 CALL tlodt2 (ts1,ts2)
 DO  i = 1,20
   ptem(i) = ptem(i) + ts2(i)
 END DO
 
!     (QQQINV) TRANSPOSE (KTR3)  (QQQINV)
 
 230 CALL  gmmats (qqinv,20,18,+1, ptem,20,1,0, ptele)
 
!     LOCATE THE TRANSFORMATION MATRICES FROM BASIC TO LOCAL (THAT IS
!     COORDINATE AT ANY GRID POINT IN WHICH DISPLACEMENT AND STRESSES
!     ARE R - NOT NEEDED IF FIELD 7 IN GRID CARD IS ZERO)
 
 DO  i = 1,36
   ptglb(i) = 0.0
 END DO
 DO  i = 1,6
   sil(i) = i
 END DO
 DO  i = 1,6
   DO  ii = 1,36
     balotr(ii) = 0.0D0
   END DO
   sil1 = sil(i)
   DO  k = 1,3
     k1 = (sil1-1)*3 + k
     psub(k) = ptele(k1)
   END DO
   CALL gmmats (e,6,3,0, psub,3,1,0, psubt)
   
!     TRANSFORM THE PSUBT(6) FROM BASIC TO DISPLACEMENT COORDINATES
   
   IF (nl(i) == 0 .OR. ics(i) == 0) GO TO 330
   jj = 4*i + 20
   CALL transs (iest(jj),trand)
   DO  jj = 1,3
     l = 6*(jj-1) + 1
     m = 3*(jj-1) + 1
     balotr(l   ) = trand(m  )
     balotr(l+1 ) = trand(m+1)
     balotr(l+2 ) = trand(m+2)
     balotr(l+21) = trand(m  )
     balotr(l+22) = trand(m+1)
     balotr(l+23) = trand(m+2)
   END DO
   CALL gmmats (balotr(1),6,6,1, psubt,6,1,0, psubt1)
   DO  k = 1,6
     psubt(k) = psubt1(k)
   END DO
   
!     INSERT PTGLB IN PG
   
   330 DO  ii = 1,6
     i1 = (i-1)*6 + ii
     i2 = iest(i+1) + ii - 1
     ptglb(i1) = psubt(ii)
     pg(i2) = pg(i2) + psubt(ii)
   END DO
 END DO
 GO TO 999
 
 901 WRITE  (iout,905) ufm,iest(1)
 905 FORMAT (a23,' 2412, A SINGULAR MATERIAL MATRIX FOR ELEMENT ID =',  &
     i9,' HAS BEEN DETECTED BY SUBROUTINE TLODT1', /26X,'WHILE',  &
     ' TRYING TO COMPUTE THERMAL LOADS WITH TEMPP2 CARD DATA.')
 nogo=.true.
 999 RETURN
END SUBROUTINE tlodt1
