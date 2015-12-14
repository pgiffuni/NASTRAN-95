SUBROUTINE hdlin (x,y,z,np,nc,  &
        xcc,icount,irct,x21,y21,z21,iia,xe,ye,xu,yu,xi,yi,zi,  &
        di,ibeg,iend,ict,icct,  &
        ind,nind,xxx,ccc,in,in1,in2,tgm,tgmt,tgi,zm,zmi,rv,  &
        rvi,nno,noct,ymin,zmin,coord,sndt,neh,KEEP)
     
 
!     THIS SUBROUTINE IS THE EXECUTIVE.
 
 
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(IN)                         :: y(1)
 REAL, INTENT(IN)                         :: z(1)
 INTEGER, INTENT(IN)                      :: np
 INTEGER, INTENT(IN OUT)                  :: nc
 REAL, INTENT(IN OUT)                     :: xcc(1)
 INTEGER, INTENT(OUT)                     :: icount(1)
 INTEGER, INTENT(OUT)                     :: irct(1)
 REAL, INTENT(OUT)                        :: x21(1)
 REAL, INTENT(OUT)                        :: y21(1)
 REAL, INTENT(OUT)                        :: z21(1)
 INTEGER, INTENT(OUT)                     :: iia(1)
 REAL, INTENT(OUT)                        :: xe(1)
 REAL, INTENT(OUT)                        :: ye(1)
 REAL, INTENT(OUT)                        :: xu(1)
 REAL, INTENT(OUT)                        :: yu(1)
 REAL, INTENT(OUT)                        :: xi(1)
 REAL, INTENT(OUT)                        :: yi(1)
 REAL, INTENT(OUT)                        :: zi(1)
 REAL, INTENT(OUT)                        :: di(1)
 INTEGER, INTENT(OUT)                     :: ibeg(1)
 INTEGER, INTENT(OUT)                     :: iend(1)
 INTEGER, INTENT(OUT)                     :: ict(1)
 INTEGER, INTENT(OUT)                     :: icct(1)
 INTEGER, INTENT(OUT)                     :: ind(1)
 INTEGER, INTENT(OUT)                     :: nind(1)
 REAL, INTENT(IN OUT)                     :: xxx(1)
 REAL, INTENT(IN)                         :: ccc(1)
 INTEGER, INTENT(OUT)                     :: in(1)
 INTEGER, INTENT(OUT)                     :: in1(1)
 INTEGER, INTENT(OUT)                     :: in2(1)
 REAL, INTENT(OUT)                        :: tgm(1)
 REAL, INTENT(OUT)                        :: tgmt(1)
 REAL, INTENT(OUT)                        :: tgi(1)
 REAL, INTENT(OUT)                        :: zm(1)
 REAL, INTENT(OUT)                        :: zmi(1)
 REAL, INTENT(OUT)                        :: rv(1)
 REAL, INTENT(OUT)                        :: rvi(1)
 INTEGER, INTENT(OUT)                     :: nno(1)
 INTEGER, INTENT(OUT)                     :: noct(1)
 REAL, INTENT(OUT)                        :: ymin(1)
 REAL, INTENT(OUT)                        :: zmin(1)
 REAL, INTENT(OUT)                        :: coord(1)
 REAL, INTENT(OUT)                        :: sndt(1)
 INTEGER, INTENT(OUT)                     :: neh(1)
 INTEGER, INTENT(OUT)                     :: KEEP(1)
 DIMENSION  i2(2),i3(2),rrx(20),ngx(15),h(15), u(6),v(6),w(6),x1(10),y1(10)
 
 
 COMMON /go3 / l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12, l13
 COMMON /hdsc/ scx,yaw,roll,pit,lz,vp,jjj,icore
 COMMON /hedg/ jat,me,jt,vx,vx1,vx2,vx3,nn
 
 IF (vp < 0.) GO TO 20
 hxx = .015
 ava = .0
 hx1 = .001
 lc  = 10**6
 ixxx= 0
 IF (scx < 0.) ixxx = 1
 scx = ABS(scx)
 
!     INITIALIZE VARIABLES.
 
 lz  = lz*5
 sw1 = 0
 sw  = 0
 idav= 0
 
!     CALCULATE MAXIMUM ALLOWABLE ELEMENTS.
 
 iabc = icore/(25+lz+4*jjj)
 sct  = 1.
 vp   = vp/sct
 vpx  = ABS(vp)
 isave= nc
 nc   = iabc
 l5   = 0
 l6   = nc
 l7   = 2*nc
 l8   = 3*nc
 l2   = 4*nc
 l3   = 5*nc
 l4   = 6*nc
 l00  = 7*nc
 l01  = 8*nc
 l1   = 9*nc
 l0   = 10*nc
 l9   = 11*nc
 l10  = 12*nc
 l11  = 13*nc
 l15  = 14*nc
 l16  = 15*nc
 l17  = 16*nc
 l18  = 19*nc
 l12  = 20*nc
 l13  = 25*nc
 l14  = l13 + lz*nc
 DO  j = 1,nc
   rvi(l8+j) = 10**6
   tgm(l5+j) = 10**6
   rv (l7+j) =-rvi(l8+j)
   tgi(l6+j) =-tgm(l5+j)
   noct(l9+j)= 0
   zm  (l2+j)= rv(l7+j)
   zmi(l3+j) = rvi(l8+j)
   nind(l16+j) = 0
   ind (l15+j) = j
   KEEP(l18+j) = 0
 END DO
 nc  = isave
 ik  = 0
 ikt = 0
 kr  = jjj
 pi  = 3.1416/180.
 u(6)= scx
 v(6)= scx
 vp  =-vp
 
!     STORE EULERIAN ANGLES.
 
 xx   = yaw*pi
 yy   = roll*pi
 zz   = pit*pi
 cosy = COS(yy)
 siny = SIN(yy)
 cosz = COS(zz)
 sinz = SIN(zz)
 cosx = COS(xx)
 sinx = SIN(xx)
 20 CONTINUE
 nt  = np-1
 ikk = ik+1
 ik  = ik+1
 
!     SET ERROR CODES, IF NECESSARY.
 
 IF (ikk <= iabc) GO TO 30
 sw = 1
 30 CONTINUE
 IF (nc == 0) GO TO 40
 idav = 1
 nc   =-sw1
 IF (sw == 0.) GO TO 50
 icore = (25+lz+4*jjj)*ikk
 nc    =-(sw+sw1)
 40 CONTINUE
 50 CONTINUE
 DO  j = 1,np
   x21(j) = x(j)
   y21(j) = y(j)
   z21(j) = z(j)
 END DO
 
!     STORE COORDINATES AND SET PEN POSITION WHENEVER ABS(Z)=9999.
 
 DO  j = 1,nt
   iia(j) = 0
   IF (z21(j) /= 9999.) CYCLE
   iia(j) = 1
   ixu = j - 2
   ibb = j - ISIGN(1,ixu)
   x21(j) = x21(ibb)
   y21(j) = y21(ibb)
   z21(j) = z21(ibb)
 END DO
 iia(np) = 1
 z21(np) = z21(nt)
 y21(np) = y21(nt)
 x21(np) = x21(nt)
 jxx = ikk
 i   = 1
 vl  = ABS(vp)
 
!     LOOP THAT DOES THE THREE DIMENSIONAL TRANSFORMATION ON THE
!     COORDINATES.
 
 jv = l14 + (ikk-1)*4*jjj
 jt = 1
 DO  j = 1,np
   xj = x21(j)/sct
   yj = y21(j)/sct
   zj = z21(j)/sct
   u(i) = zj*(cosy*sinx) + xj*(cosy*cosx) - yj*siny
   tw = yj*cosy*cosz
   tz = xj*( sinz*sinx+siny*cosz*cosx)
   ty = zj*(-sinz*cosx+siny*cosz*sinx)
   v(i) = tz + tw + ty
   pt = yj*cosy*sinz
   pk = zj*( cosz*cosx+siny*sinz*sinx)
   ps = xj*(-cosz*sinx+siny*sinz*cosx)
   zj = pk + ps + pt
   IF (zj < vl) GO TO 80
   sw1 = 2
   vpx = AMAX1(zj,vpx)
   vpx = vpx + (.5/sct)
   80 CONTINUE
   t = sw + sw1
   IF (t /= 0.) CYCLE
   
!     CALCULATES PERSPECTIVE BASED ON VALUE VP(DV) FROM CALLING PROGRAM.
   
   hh = vl/(vl-zj)
   x21(j) = u(i)*hh
   y21(j) = v(i)*hh
   z21(j) = zj*hh
   
!     CALCULATES MAX/MIN VALUES OF EACH ELEMENT ON THE X,Y,Z DIMENSION
   
   rv (l7+jxx) = AMAX1(rv (l7+jxx),y21(j))
   rvi(l8+jxx) = AMIN1(rvi(l8+jxx),y21(j))
   tgi(l6+jxx) = AMAX1(tgi(l6+jxx),x21(j))
   tgm(l5+jxx) = AMIN1(tgm(l5+jxx),x21(j))
   zm (l2+jxx) = AMAX1(zm (l2+jxx),z21(j))
   zmi(l3+jxx) = AMIN1(zmi(l3+jxx),z21(j))
   coord(jt+jv  ) = x21(j)
   coord(jt+jv+1) = y21(j)
   coord(jt+jv+2) = z21(j)
   coord(jt+3+jv) = iia(j)
   jt = jt + 4
 END DO
 IF (idav == 1) vp = vpx*sct
 IF (t   /= 0.) GO TO 400
 noct(l9+ikk) = noct(l9+ikk) + np
 ns  = np
 ava = ava + (tgi(l6+jxx)-tgm(l5+jxx))*(rv(l7+jxx)-rvi(l8+jxx))
 IF (ixxx == 1) GO TO 95
 
!     CALL SUBROUTINE WHICH CALCULATES BOTH THE EQUATIONS OF THE LINE
!     SEGMENTS AND POLYGONS.
 
 CALL hdcoef (x21,y21,z21,xxx,jxx,ns,ccc,lz)
 
!     CHECKS TO SEE IF ALL ELEMENTS(SETS) HAVE BEEN PASSED.
 
 95 CONTINUE
 IF (idav == 1) GO TO 100
 GO TO 400
 100 CONTINUE
 ava = ava/ikk
 DO  j = 1,100
   icct(j) = 0
   ict (j) = 0
   irct(j) = j - 1
   ibeg(j) = 1
   iend(j) = 0
 END DO
 iaug  = 50 + (ikk/10000)*2
 amaxx =-999999.
 amaxy =-999999.
 aminx = 999999.
 aminy = 999999.
 DO  j = 1,ikk
   amaxx = AMAX1(amaxx,tgi(l6+j))
   amaxy = AMAX1(amaxy,rv (l7+j))
   aminx = AMIN1(aminx,tgm(l5+j))
   aminy = AMIN1(aminy,rvi(l8+j))
 END DO
 tmax = (amaxx-aminx)*(amaxy-aminy)
 ibl  = tmax/ava
 ibl  = ibl/4
 
!     DETERMINES THE NUMBER OF GRID POINTS IN THE GRID.
 
 
 en = ikk
 k  = (ALOG(en)/ALOG(2.)) + .01
 k  = k + iaug
 k  = MIN0(k,ibl)
 IF (k <= 1) k = 1
 t  = k
 r  = t**.5
 ks = r + .5
 s  = t/ks
 ms = s + .5
 n  = ks*ms
 mnd= n + 1
 xmd= mnd
 t  = 3./(mnd-1)
 igy= t*ikk
 k  = ks
 k1 = ms
 crx= (amaxx-aminx)/k
 cry= (amaxy-aminy)/k1
 
 
!     DETERMINES THE RELEVANT ELEMENTS VIA THE GRID BLOCKS.
 
 
 loop3:  DO  j = 1,ikk
   ia = 0
   xmat = tgi(l6+j)
   xmit = tgm(l5+j)
   ymat =  rv(l7+j)
   ymit = rvi(l8+j)
   m = 0
   DO  i = 1,k1
     DO  l = 1,k
       m  = m + 1
       s  = xmat - ((l-1)*crx+aminx)
       s1 = xmat - (l*crx+aminx)
       r  = xmit - ((l-1)*crx+aminx)
       r1 = xmit - (l*crx+aminx)
       a  = ymat - ((i-1)*cry+aminy)
       a1 = ymat - (i*cry+aminy)
       b  = ymit - ((i-1)*cry+aminy)
       b1 = ymit - (i*cry+aminy)
       IF (s <= 0. .OR. r1 >= 0.) CYCLE
       IF (a <= 0. .OR. b1 >= 0.) CYCLE
       IF (s*s1 > 0. .OR. r*r1 > 0.) GO TO 4
       IF (a*a1 > 0. .OR. b*b1 > 0.) GO TO 4
       nind(l16+j) = m
       CYCLE loop3
       4 CONTINUE
       ia = ia + 1
       IF (ia <= 4) GO TO 8000
       nind(j+l16) = 0
       GO TO 8001
       8000 CONTINUE
       nind(l16+j) = nind(l16+j) + m*(mnd**(ia-1))
       8001 CONTINUE
       IF (icct(m) < 0) CYCLE
       icct(m) = icct(m) + 1
       jk      = (m-1)*igy+icct(m) + l17
       neh(jk) = j
       IF (icct(m) >= igy) icct(m) = -1
     END DO
   END DO
 END DO loop3
 CALL hdvs1 (nind(l16+1),ik,ind(l15+1))
 sw = 0
 l  = 1
 DO  i = 1,ikk
   11 CONTINUE
   IF (nind(l16+i) /= irct(l)) GO TO 6
   sw = sw + 1
   IF (sw == 1.) LT = i
   ict(l) = ict(l) + 1
   CYCLE
   6 CONTINUE
   IF (sw /= 0.) GO TO 8
   l = l + 1
   GO TO 11
   8 CONTINUE
   ibeg(l) = LT
   iend(l) = LT + ict(l) - 1
   sw = 0
   IF (nind(l16+i) >= mnd) GO TO 2110
   l = l + 1
   GO TO 11
 END DO
 ibeg(l) = LT
 iend(l) = LT + ict(l) - 1
 2110 CONTINUE
 DO  j = 1,ikk
   sndt(l4+j) = ind(l15+j)
 END DO
 CALL hdvsr (sndt(l4+1),ik,nind(l16+1))
 en  = ikk
 igx = (ALOG(en)/ALOG(2.)) + 1.
 DO  j = 1,igx
   rrx(j) = 2**(igx-j)
 END DO
 u(6) = scx
 v(6) = scx
 w(6) = scx
 ikt  = nc
 t    = aminy
 t1   = aminx
 v(5) = t
 u(5) = t1
 ij   = 0
 x1(3)= u(5)
 y1(3)= v(5)
 x1(4)= u(6)
 y1(4)= v(6)
 x1(4)= x1(4)/sct
 y1(4)= y1(4)/sct
 DO  j = 1,ikk
   in(l11 +j) = j
   in1(l0 +j) = j
   in2(l00+j) = j
   tgmt(l10+j) = tgm(l5+j)
   ymin(l1 +j) = rvi(l8+j)
   zmin(l01+j) = zm(l2+j)
 END DO
 
!     CALL SUBROUTINE WHICH WILL SORT ON X,Y AND Z.
 
 CALL hdvsr (tgmt(l10+1),ik,in(l11+1))
 CALL hdvsr (ymin(l1+1),ik,in1(l0+1))
 CALL hdvsr (zmin(l01+1),ik,in2(l00+1))
 h(8) = 0
 DO  j = 1,ikk
   ks = ikk
   jj = l14 + (j-1)*4*jjj
   jh = 1
   ii = 0
   ixr= noct(l9+j)
   nit= 0
   jt = l12 + 5*(j-1)
   jo = l13 + lz*(j-1)
   IF (ixxx == 1) GO TO 200
   ns = xxx(5+jt)
   ng = ns*5
   a3 = xxx(1+jt)
   b3 = xxx(2+jt)
   c3 = xxx(3+jt)
   d3 = xxx(4+jt)
   i  = 0
   DO  ix = 1,ng,5
     IF (ixr <= 3) CYCLE
     i = i + 1
     xe(i) = ccc(ix+3+jo)
     IF (ccc(ix+jo) /= 0.) GO TO 120
     xe(i) =-ccc(ix+2+jo)
     ye(i) = ccc(ix+3+jo)
     CYCLE
     120 CONTINUE
     ye(i) =-ccc(ix+2+jo) - ccc(ix+1+jo)*xe(i)
   END DO
   
!     THIS LOOP DETERMINES THE RELEVANT ELEMENTS AS THEY RELATE TO A
!     PARTICULAR ELEMENT.  THAT IS, EACH ELEMENT HAS ASSOCIATED WITH IT
!     THOSE OTHER ELEMENTS WHICH COULD POSSIBLY HIDE SOME PORTION
!     OF THE GIVEN ELEMENT.
   
   k  = 2**igx
   k1 = k
   k2 = k
   
!     DO LOGARITHMIC SEARCH TO DETERMINE RELEVANT ELEMENTS.
   
   s = -1
   DO  i = 1,igx
     k = k + SIGN(rrx(i),s)
     IF (k > ikk) k = ikk
     s  = tgi(l6+j) - tgmt(l10+k  )
     s1 = tgi(l6+j) - tgmt(l10+k-1)
     IF (s*s1 <= 0.) GO TO 132
   END DO
   k = ikk
   132 CONTINUE
   s = -1
   DO  i = 1,igx
     k1 = k1 + SIGN(rrx(i),s)
     IF (k1 > ikk) k1 = ikk
     s  = rv(l7+j) - ymin(l1+k1  )
     s1 = rv(l7+j) - ymin(l1+k1-1)
     IF (s*s1 <= 0.) GO TO 134
   END DO
   k1 = ikk
   134 CONTINUE
   s = -1
   DO  i = 1,igx
     k2 = k2 + SIGN(rrx(i),s)
     IF (k2 <=   1) k2 = 2
     IF (k2 > ikk) k2 = ikk
     s  = zmi(l3+j) - zmin(l01+k2  )
     s1 = zmi(l3+j) - zmin(l01+k2-1)
     IF (s*s1 <= 0.) GO TO 136
   END DO
   k2 = 1
   136 CONTINUE
   i1 = ikk - k2 + 1
   
!     RETRIEVE THE RELEVANT ELEMENTS DETERMINED FROM SCHEME 1.
   
   IF  (nind(l16+j) == 0) GO TO 1270
   ir = nind(l16+j)
   vx = nind(l16+j)
   t  = ALOG(vx)
   IF (nind(l16+j) <= lc) GO TO 1800
   e  = lc
   lg = nind(l16+j)/lc
   mu = MOD(ir,lc)
   ux = lg + (mu/e)
   t  = ALOG(ux) + ALOG(e)
   1800 CONTINUE
   ixt = 0
   iexp= (t/ALOG(xmd)) + 1
   DO  l = 1,iexp
     iv = ir/(mnd**(iexp-l))
     ir = ir - iv*(mnd**(iexp-l))
     iv = iv + 1
     IF (icct(iv-1) == 0) EXIT
     IF (icct(iv-1) > 0) EXIT
     GO TO 1270
     4001 CONTINUE
     ke = icct(iv-1)
     il = 0
     jtt= (iv-2)*igy + l17
     DO  i = 1,ke
       kv = neh(i+jtt)
       IF (KEEP(l18+kv) == j) EXIT
       il = il + 1
       nno(l4+ixt+il) = kv
       KEEP(l18+kv) = j
     END DO
     ixt = ixt + il
     4000 CONTINUE
     ix  = ibeg(iv)
     ix1 = iend(iv)
     DO  i = ix,ix1
       nno(l4+ixt+i-ix+1) = ind(l15+i)
     END DO
     ixt = ixt + ix1 - ix + 1
   END DO
   ks = ixt
   1270 CONTINUE
   im = MIN0(i1,k,k1)
   
!     PICK MINIMUM COUNT FROM BOTH SCHEMES.
   
   IF (ks < im) GO TO 129
   IF (im == i1) GO TO 1000
   IF (im ==  k) GO TO 1001
   IF (im == k1) GO TO 1002
   1000 CONTINUE
   ks = i1
   DO  i = 1,ks
     nno(l4+i) = in2(l00+ikk-i+1)
   END DO
   GO TO 129
   1001 CONTINUE
   ks = k
   DO  i = 1,ks
     nno(l4+i) = in(l11+i)
   END DO
   GO TO 129
   1002 CONTINUE
   ks = k1
   DO  i = 1,ks
     nno(l4+i) = in1(l0+i)
   END DO
   129 CONTINUE
   loop170:  DO  i = 1,ks
     it = 0
     jb = nno(l4+i)
     IF (j == jb) CYCLE loop170
     jk = l13 + lz*(jb-1)
     js = l12 +  5*(jb-1)
     IF (tgm(l5+j) >= tgi(l6+jb) .OR. tgi(l6+j) <= tgm(l5+jb)) CYCLE loop170
     IF (rv(l7+j) <= rvi(l8+jb) .OR. rvi(l8+j) >= rv(l7+jb)) CYCLE loop170
     IF (zmi(l3+j) >= zm(l2+jb)) CYCLE loop170
     nv = xxx(5+js)
     IF (xxx(js+3) == 0.) CYCLE loop170
     IF (xxx(3+jt) == 0.) GO TO 165
     nb = 5*nv
     
     
!     TEST TO SEE IF ALL VERTICES LIE EITHER BEHIND OR IN FRONT OF
!     THE GIVEN POLYGON.
     
     
     m = 0
     DO  ix = 1,nb,5
       m = m + 1
       a = ccc(ix+3+jk)
       IF (ccc(ix+jk) /= 0.) GO TO 130
       a =-ccc(ix+2+jk)
       b = ccc(ix+3+jk)
       GO TO 140
       130 CONTINUE
       b =-ccc(ix+2+jk) - ccc(ix+1+jk)*a
       140 CONTINUE
       xu(m) = a
       yu(m) = b
       vx  = xxx(4+js)
       vx1 = xxx(2+js)*b
       vx2 = xxx(1+js)*a
       zs  =-(vx+vx1+vx2)/xxx(3+js)
       vx  = xxx(4+jt)
       vx1 = xxx(2+jt)*b
       vx2 = xxx(1+jt)*a
       zs1 =-(vx+vx1+vx2)/xxx(3+jt)
       IF (ABS(zs-zs1) < hxx) CYCLE
       it = it + 1
       icount(it) = 0
       IF (zs > zs1) icount(it) = 1
     END DO
     
     
!     TESTS FOR SEMI-RELEVANT PLANES.  THAT IS,NEGATIVE INDEXES
!     INDICATE ELEMENT IS TO BE USED FOR VISIBILITY TEST, BUT NOT FOR
!     INTERSECTION LINE DETERMINATION.
     
     
     IF (it == 0) CYCLE loop170
     l = 0
     DO  m = 1,it
       l = l + icount(m)
     END DO
     IF (l  ==  0) CYCLE loop170
     IF (l  == it) jb = -jb
     IF (ii /=  0) GO TO 165
     
     
!     INTERROGATE THE RELATIONSHIP OF THE CANDIDATE POLYGON TO THE
!     GIVEN POLYGON BY DETERMINING IF THE PROJECTION OF ONE POLYGON
!     CAN BE SEPARATED BY AN EDGE FROM THE OTHER'S PROJECTION
     
     
     c3 = xxx(3+jt)
     c4 = xxx(3+js)
     sd = 0
     i3(1) = jk
     i3(2) = jo
     i2(1) = nv*5
     i2(2) = ns*5
     DO  ku = 1,2
       is = i3(ku)
       ib = i2(ku)
       DO  l = 1,ib,5
         151 CONTINUE
         IF (sd == 1.) GO TO 152
         a = xxx(2+jt)*c4 - xxx(2+js)*c3
         b = xxx(1+jt)*c4 - xxx(1+js)*c3
         c = xxx(4+jt)*c4 - xxx(4+js)*c3
         GO TO 153
         152 CONTINUE
         a = ccc(l+is  )
         b = ccc(l+is+1)
         c = ccc(l+is+2)
         153 CONTINUE
         IF (a == 0. .AND. b == 0.) GO TO 162
         IF (a /= 0.) GO TO 154
         a = 0
         c = c/b
         b = 1
         GO TO 155
         154 CONTINUE
         b = b/a
         c = c/a
         a = 1
         155 CONTINUE
         m = 0
         r1= 0
         DO  ix = 1,nv
           m = m + 1
           yg= yu(m)
           IF (a /= 0.) GO TO 156
           dy = -c/b
           yg = xu(m)
           GO TO 157
           156 CONTINUE
           dy = -c - b*xu(m)
           157 IF (ABS(dy-yg) < hxx) CYCLE
           r = yg - dy
           IF (r*r1 < 0.) GO TO 162
           r1 = r
         END DO
         m  = 0
         r2 = 0
         DO  ix = 1,ns
           m  = m + 1
           yg = ye(m)
           IF (a /= 0.) GO TO 159
           dy = -c/b
           yg = xe(m)
           GO TO 160
           159 CONTINUE
           dy = -c - b*xe(m)
           160 IF (ABS(dy-yg) < hxx) CYCLE
           r  = yg - dy
           IF (r*r2 < 0.) GO TO 162
           r2 = r
         END DO
         IF (r1*r2 < 0.) CYCLE loop170
         162 CONTINUE
         IF (sd /= 0.) CYCLE
         sd = 1
         GO TO 151
       END DO
     END DO
     165 CONTINUE
     ii = ii + 1
     nno(l4+ii) = jb
   END DO loop170
   js  = 1
   jat =-4
   jt  = l12 + (j-1)*5
   nn  = xxx(jt+5)
   vx  = xxx(jt+4)
   vx1 = xxx(2+jt)
   vx2 = xxx(1+jt)
   vx3 = xxx(3+jt)
   IF (ixr <= 2) GO TO 200
   IF (ii  == 0) GO TO 190
   
!     CALL SUBROUTINE WHICH SOLVES FOR THE LINES OF INTERSECTION,IF ANY,
!     OF THE JTH ELEMENT WITH OTHER ELEMENTS.
   
   CALL hdsolv(ixr,j,xxx,ccc,ii,nno,nit,x21,y21,z21,iia,nc,zm,zmi,lz)
   190 CONTINUE
   200 CONTINUE
   DO  jm = 1,ixr
     x21(jm) = coord(jh  +jj)
     y21(jm) = coord(jh+1+jj)
     z21(jm) = coord(jh+2+jj)
     iia(jm) = coord(jh+3+jj)
     jh = jh + 4
   END DO
   ixr = ixr + 3*nit
   IF (ii   == 0) GO TO 220
   IF (ixxx /= 1) GO TO 240
   220 CONTINUE
   DO  jm = 1,ixr
     x1(2) = x21(jm)
     y1(2) = y21(jm)
     im = iia(jm)
     CALL hdplt (x1,y1,ij,im)
   END DO
   GO TO 390
   240 CONTINUE
   jx = 1
   250 CONTINUE
   
!     PLOTS IF IIA(JX+1) IS EQUAL TO 1.
   
   IF (iia(jx) == 0 .AND. iia(jx+1) == 0) GO TO 260
   im    = iia(jx+1)
   x1(2) = x21(jx+1)
   y1(2) = y21(jx+1)
   CALL hdplt (x1,y1,ij,im)
   jx = jx + 2
   IF (jx >= ixr) GO TO 390
   GO TO 250
   260 CONTINUE
   jat = jat + 5
   me  = 0
   
!     CALL SUBROUTINE WHICH DETERMINES THE POINTS OF INTERSECTIONS
!     OF THE LINES OF THE JTH SET WITH THE RELEVANT LINES AND PLANES
!     OF OTHER ELEMENTS.
   
   CALL hdchk (xxx,ccc,nno,ii,xi,yi,ngx,zm,zmi,rv,rvi,tgm,tgi,zi,lz, xcc)
   IF (js /= 1) STOP 'my GOSH. JS IS NOT 1 /HDLIN'
   ng = ngx(js) + 2
   xi(1)  = x21(jx)
   yi(1)  = y21(jx)
   zi(1)  = z21(jx)
   xi(ng) = x21(jx+1)
   yi(ng) = y21(jx+1)
   zi(ng) = z21(jx+1)
   IF (ng <= 3) GO TO 340
   
!     THE FOLLOWING CODE SORTS THE INTERSECTION POINTS IN ASCENDING
!     ORDER OF OCCURENCE AND THEN SHRINKS THE LIST IF REDUNDANCY EXIST.
   
   ni  = ng - 2
   nii = ni
   DO  m = 1,ng
     di(m) = (xi(m)-xi(1))**2
     ppppp = (yi(m)-yi(1))**2
     di(m) = di(m) + ppppp
   END DO
   DO  m = 2,ni
     DO  mx= 2,nii
       IF (di(mx) <= di(mx+1)) CYCLE
       hold   = di(mx)
       hold1  = xi(mx)
       hold2  = yi(mx)
       hold3  = zi(mx)
       xi(mx) = xi(mx+1)
       yi(mx) = yi(mx+1)
       zi(mx) = zi(mx+1)
       di(mx) = di(mx+1)
       di(mx+1) = hold
       xi(mx+1) = hold1
       yi(mx+1) = hold2
       zi(mx+1) = hold3
     END DO
     nii = nii - 1
   END DO
   lx  = 1
   npx = ng
   300 npx = npx - 1
   i   = lx
   DO  m = i,npx
     rx  = 0
     t   = xi(m) - xi(m+1)
     t1  = yi(m) - yi(m+1)
     t   = (t**2+t1**2)**.5
     IF (t > hx1) CYCLE
     ix  = m
     ix1 = npx
     DO  mx = ix,ix1
       xi(mx) = xi(mx+1)
       yi(mx) = yi(mx+1)
       zi(mx) = zi(mx+1)
     END DO
     rx = 1
     lx = m
     IF (lx == npx) EXIT
     GO TO 300
   END DO
   330 CONTINUE
   IF (rx == 1.) npx = npx - 1
   ng = npx + 1
   340 CONTINUE
   
!     THIS CODE DETERMINES THE HDSTUS(VISIBILITY) OF EVERY OTHER POINT
!     AS SUGGESTED BY THE THEOREM IN THE TECHNICAL REPORT.
   
   DO  l = 1,ng,2
     
     oj  = xi(l)
     tmj = yi(l)
     zj  = zi(l)
     CALL hdstus (oj,tmj,xxx,tgm,rv,rvi,tgi,zm,nno,ii,h,im,jxt,zj,nc,  &
         zmi,ccc,lz)
     di(l) = im
   END DO
   DO  l = 1,ng,2
     IF (l == ng  ) CYCLE
     IF (l == ng-1) GO TO 360
     c = di(l) + di(l+2)
     IF (c /= 2.) GO TO 360
     di(l+1) = di(l)
     CYCLE
     360 oj  = xi(l+1)
     tmj = yi(l+1)
     zj  = zi(l+1)
     CALL hdstus (oj,tmj,xxx,tgm,rv,rvi,tgi,zm,nno,ii,h,im,jxt,zj,nc,  &
         zmi,ccc,lz)
     di(l+1) = im
   END DO
   
!     THE FOLLOWING CODE ACTUALLY PLOTS THE POINTS ON A GIVEN LINE
!     GOVERNED BY THE VALUE(IM) RETURNED BY HDSTUS SUBROUTINE.
!     1 MEANS INVISIBLE,...0 MEANS VISIBLE.
   
   DO  l = 1,ng
     x1(2) = xi(l)
     y1(2) = yi(l)
     im = di(l)
     CALL hdplt (x1,y1,ij,im)
     IF (l == ng) CYCLE
     c = di(l) + di(l+1)
     IF (c > 0.) CYCLE
     h(8) = 1
     oj   = (xi(l)+xi(l+1))/2
     tmj  = (yi(l)+yi(l+1))/2
     zj   = (zi(l)+zi(l+1))/2
     CALL hdstus (oj,tmj,xxx,tgm,rv,rvi,tgi,zm,nno,ii,h,im,jxt,zj,nc,  &
         zmi,ccc,lz)
     h(8) = 0
     x1(2)= oj
     y1(2)= tmj
     CALL hdplt (x1,y1,ij,im)
   END DO
   jx = jx + 1
   GO TO 250
   390 CONTINUE
   
!     DECREMENTS THE COUNT OF THE NUMBER OF LINES IN THE JTH SET
!     SINCE THE LINES OF INTERSECTIONS WERE ADDED TO THIS ELEMENT
!     BY THE SUBROUTINE SOLVE.
   
   xxx(5+jt) = xxx(5+jt) - nit
 END DO
 400 CONTINUE
 RETURN
END SUBROUTINE hdlin
