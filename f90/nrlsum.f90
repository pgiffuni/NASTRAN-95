SUBROUTINE nrlsum
     
!     NRLSUM   OES2,OEF2/NRLSTR,NRLFOR/V,N,NMODES/V,N,NSHOCK(NDIR)/
!              C,Y,DIRECT=123/C,Y,SQRSS=0 $
 
!     NRLSUM COMPUTES NRL SUM STRESSES AND FORCES FOR DDAM. IT IS
!     ASSUMED THAT THE USER HAS REQUESTED STRESSES AND FORCES IN SORT2
!     FORMAT (BUT RESULTS WILL BE SORT1). NRLSUM READS ITEMS FOR AN
!     ELEMENT (FOR ALL SUBCASES) AND COMPUTES THE NRL SUM.  UP TO 3
!     SCRATCH FILES ARE USED TO STORE THE SUMS FOR EACH SHOCK DIRECTION.
!     PRINCIPAL STRESSES WILL BE COMPUTED BASED ON THE SUMS. THE NUMBER
!     OF SUBCASES IS NMODES*NSHOCK WITH THE ORDER 1-NMODES,
!     NMODES+1 - 2*NMODES, ... NSHOCK*NMODES.
 
!     (IF (SQRSS.EQ.1), SQUARE ROOT OF THE SUM OF THE SQUARES IS USED
!     INSTEAD OF NRL SUM
 
 INTEGER :: FILE,buf1,buf2,buf3,buf4,oes2,sysbuf,eltype,scr1,  &
     scr2,scr3,elid,oldtyp,scr(3), oef2,ofil,idir(2),inum(3),nsub(3),DIRECT,sqrss
 DIMENSION       sig(6),sigp(3),smat(3,3),DCOS(3,3)
 DIMENSION       z(20),nam(2),stress(146),istres(146),mcb(7)
 COMMON /system/ sysbuf
 COMMON /BLANK / nmodes,nshock,DIRECT,sqrss
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE    (z(1),iz(1)), (stress(1),istres(1))
 EQUIVALENCE    (sigp(1),sa), (sigp(2),sb), (sigp(3),sc)
 EQUIVALENCE    (sig(1) ,sx), (sig(2) ,sy), (sig(3) ,sz),  &
     (sig(4),sxy), (sig(5),syz), (sig(6),szx)
 DATA    oes2  , nrlstr,scr1,scr2,scr3 / 101,201,301,302,303/
 DATA    oef2  , nrlfor /102,202 /
 DATA    scr   / 301,302,303     /
 DATA    idir  / 4HDIRE,4HCTIO   /
 DATA    inum  / 4HN 1 ,4HN 2 ,4HN 3 /
 DATA    dtor  / 0.0174532925E0  /
 DATA    nam   / 4HNRLS,4HUM     /, i0 / 0 /
 
 lcore= korsz(z)
 buf1 = lcore- sysbuf + 1
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 IF (nshock == 3) GO TO 20
 IF (nshock == 2) GO TO 10
 buf4 = buf2
 buf3 = buf2
 GO TO 20
 
 10 buf4 = buf3
 
 20 lcore = buf4 - 1
 IF (lcore <= 0) GO TO 1008
 ndir = nshock
 IF (ndir  > 1) GO TO 11
 nsub(1) = DIRECT
 GO TO 14
 11 IF (ndir   >  2) GO TO 13
 IF (DIRECT == 23) GO TO 12
 nsub(1) = 1
 nsub(2) = 2
 IF (DIRECT == 13) nsub(2) = 3
 GO TO 14
 12 nsub(1) = 2
 nsub(2) = 3
 GO TO 14
 13 nsub(1) = 1
 nsub(2) = 2
 nsub(3) = 3
 14 CONTINUE
 ifil = oes2
 ofil = nrlstr
 
 15 FILE = ifil
 oldtyp = 0
 CALL OPEN (*710,ifil,z(buf1),0)
 CALL fwdrec (*1002,ifil)
 CALL gopen (scr1,z(buf2),1)
 IF (nshock > 1) CALL gopen (scr2,z(buf3),1)
 IF (nshock == 3) CALL gopen (scr3,z(buf4),1)
 
 30 CALL READ (*410,*1003,ifil,stress,146,1,iwords)
 eltype= istres(3)
 nwds  = istres(10)
 i5    = istres(5)
 elid  = istres(5)/10
 
!     REFORMULATE TO SORT 1 FORMAT
 
 istres(2) = 5
 IF (ifil == oef2) istres(2) = 4
 istres(141) = idir(1)
 istres(142) = idir(2)
 
!     WRITE ONTO SCRATCH ONLY FOR NEW ELEMENT TYPE
 
 IF (eltype == oldtyp)  GO TO 45
 DO  i = 1,nshock
   istres(4) = i
   istres(5) = i
   istres(8) = i
   isub = nsub(i)
   istres(143) = inum(isub)
   IF (oldtyp /= 0) CALL WRITE (scr(i),0,0,1)
   CALL WRITE (scr(i),stress,146,1)
 END DO
 
 oldtyp = eltype
 
!     READ STRESS INFO FOR NUMBER OF MODES AND SHOCK DIRECTIONS
 
 IF (nmodes*nwds > lcore)  GO TO 1008
 45 DO  ns = 1,nshock
   iscr = 300 + ns
   
   CALL fread (ifil,z(1),nwds*nmodes,0)
   
!     GO TO PROPER SECTION FOR EACH ELEMENT TYPE
   
!     FOR FORCES, COMPUTATIONS ARE EASIER. SO LETS NOT HAVE A COMPUTED
!     GO TO
   
   IF (ifil == oes2) GO TO 46
   
   IF (eltype >= 20 .AND. eltype <= 33) CYCLE
   IF (eltype >= 39 .AND. eltype <= 52) CYCLE
   IF (eltype == 62 .OR.  eltype == 68 .OR. eltype == 69 .OR.  &
       eltype == 72) CYCLE
   IF (eltype >= 65 .AND. eltype <= 67) CYCLE
   IF (eltype == 9 .OR.  eltype == 16 .OR. eltype == 73 .OR.  &
       eltype == 76) CYCLE
   i3 = 1
   i2 = nwds
   i1 = 2
   IF (eltype == 35 .OR.  eltype == 70 .OR. eltype == 71) i1 = 3
   GO TO 105
   
   46 CONTINUE
   
   GO TO ( 50, 60, 50, 70, 70, 80, 80, 80, 90, 50,  &
       100,100,100,400, 80, 90, 80, 80, 80,400,  &
       400,400,400,400,400,400,400,400,400,400,  &
       400,400,400,110,120,130,140,150,160,160,  &
       160,160,400,400,400,400,400,400,400,400,  &
       400,400,170,170,170,170,170,170,170,170,  &
       170, 90, 90, 80,220,220,220,400,400,180,  &
       190,400,200,200,200,210,400,400,400,400, 400,400, 80), eltype
   
!     ROD, TUBE, CONROD
   
   50 i1 = 2
   i2 = 4
   i3 = 2
   ASSIGN 55 TO iret
   GO TO 390
   
!     IGNORE MARGINS OF SAFETY
   
   55 iz(i0+3) = 1
   iz(i0+5) = 1
   GO TO 395
   
!     BEAM
   
   60 i1 = 2
   i2 = 5
   i3 = 1
   ASSIGN 65 TO iret
   GO TO 390
   65 z(6) = z(5) + AMAX1(z(2),z(3),z(4))
   z(7) = z(5) + AMIN1(z(2),z(3),z(4))
   iz(i0+8) = 1
   i1 = 9
   i2 = 11
   i3 = 1
   ASSIGN 66 TO iret
   GO TO 390
   66 z(12) = z(5) + AMAX1(z(9),z(10),z(11))
   z(13) = z(5) + AMIN1(z(9),z(10),z(11))
   iz(i0+14) = 1
   GO TO 395
   
!     SHEAR
   
   70 i1 = 2
   i2 = 3
   i3 = 1
   ASSIGN 75 TO iret
   GO TO 390
   75 iz(i0+4) = 1
   GO TO 395
   
!     TRBSC, TRPLT, QDPLT, TRIA1, TRIA2, TRIA3, QUAD1, QUAD2, QUAD4
   
   80 i1 = 3
   i2 = 5
   i3 = 1
   j3 = 3
   j4 = 4
   j5 = 5
   j6 = 6
   j7 = 7
   j8 = 8
   j9 = 9
   ASSIGN 85 TO iret
   GO TO 390
   85 ss = .5*(z(j3) + z(j4))
   st = z(j3) - z(j4)
   sq = SQRT(.25*st**2 + z(j5)**2)
   z(j7) = ss + sq
   z(j8) = ss - sq
   z(j9) = sq
   sd = 2.*z(j5)
   IF (ABS(sd) < 1.e-15 .AND. ABS(st) < 1.e-15) GO TO 87
   z(j6) = ATAN2(sd,st)*28.6478898
   GO TO 88
   87 z(j6) = 0.
   88 IF (j3 == 11) GO TO 395
   IF (eltype == 9 .OR. eltype == 16) GO TO 395
!                 TRMEM             QDMEM
   IF (eltype == 62 .OR. eltype == 63) GO TO 395
!                QDMEM1            QDMEM2
   IF (eltype == 35) GO TO 125
!                  CONEAX
   i1 = 11
   i2 = 13
   i3 = 1
   j3 = 11
   j4 = 12
   j5 = 13
   j6 = 14
   j7 = 15
   j8 = 16
   j9 = 17
   GO TO 390
   
!     TRMEM, QDMEM,  QDMEM1, QDMEM2
   
   90 i1 = 2
   i2 = 4
   i3 = 1
   j3 = 2
   j4 = 3
   j5 = 4
   j6 = 5
   j7 = 6
   j8 = 7
   j9 = 8
   ASSIGN 85 TO iret
   GO TO 390
   
!     CELAS1,2,3
   
   100 i1 = 2
   i2 = 2
   i3 = 1
   105 ASSIGN 395 TO iret
   GO TO 390
   
!     BAR - ADD AXIAL STRESS TO EXTENSIONAL STRESSES DUE TO BENDING
!           BEFORE COMPUTING NRL SUMS. THEN ZERO OUT AXIAL STRESS
!           AND MAX AND MIN STRESSES
   
   110 i1 = 2
   i2 = 5
   i3 = 1
   DO  j = 1,nmodes
     isub = nwds*(j-1)
     DO  i = 2,5
       z(isub+i) = z(isub+i)+z(isub+6)
     END DO
     DO  i = 10,13
       z(isub+i) = z(isub+i)+z(isub+6)
     END DO
   END DO
   ASSIGN 115 TO iret
   GO TO 390
   115 z(6) = 0.
   z(7) = 0.
   z(8) = 0.
   iz(i0+9) = 1
   i1   = 10
   i2   = 13
   i3   = 1
   ASSIGN 116 TO iret
   GO TO 390
   116 z(14) = 0.
   z(15) = 0.
   iz(i0+16) = 1
   GO TO 395
   
!     CONEAX
   
   120 i1 = 4
   i2 = 6
   i3 = 1
   j3 = 4
   j4 = 5
   j5 = 6
   j6 = 7
   j7 = 8
   j8 = 9
   j9 = 10
   ASSIGN 85 TO iret
   GO TO 390
   125 IF (j3 == 12) GO TO 395
   i1 = 12
   i2 = 14
   i3 = 1
   j3 = 12
   j4 = 13
   j5 = 14
   j6 = 15
   j7 = 16
   j8 = 17
   j9 = 18
   GO TO 390
   
!     TRIARG
   
   130 i1 = 2
   i2 = 5
   i3 = 1
   GO TO 105
   
!     TRAPRG
   
   140 i1 = 2
   i2 = 21
   i3 = 1
   GO TO 105
   
!     TORDRG
   
   150 i1 = 2
   i2 = 16
   i3 = 1
   GO TO 105
   
!     TETRA, WEDGE, HEXA1, HEXA2
   
   160 i1 = 2
   i2 = 7
   i3 = 1
   ASSIGN 165 TO iret
   GO TO 390
   165 z(8) = SQRT((z(2)-z(3))**2 + (z(3)-z(4))**2 + (z(4)-z(2))**2 +  &
       6.*(z(5)**2 + z(6)**2 + z(7)**2)) / 3.
   z(9) = -(z(2)+z(3)+z(4)) / 3.
   GO TO 395
   
!     DUM1 - DUM9
   
   170 i1 = 2
   i2 = 10
   i3 = 1
   GO TO 105
   
!     TRIAAX
   
   180 i1 = 3
   i2 = 11
   i3 = 1
   GO TO 105
   
!     TRAPAX
   
   190 i1 = 3
   i2 = 47
   i3 = 1
   GO TO 105
   
!     TRIM6, TRPLT1, TRSHL
   
   200 iend = 8
   iskip= 8
   IF (eltype /= 73) GO TO 201
   iend = 4
   iskip= 7
   201 j2 = -5
   ij = 0
   202 ij = ij + 1
   j2 = j2 + iskip
   j4 = j2 + 2
   i1 = j2
   i2 = j4
   i3 = 1
   ASSIGN 205 TO iret
   GO TO 390
   205 ss = .5*(z(j2)+z(j2+1))
   st = z(j2) - z(j2+1)
   sq = SQRT(.25*st**2 + z(j4)**2)
   z(j4+2) = ss + sq
   z(j4+3) = ss - sq
   z(j4+4) = sq
   sd = 2.*z(j4)
   IF (ABS(sd) < 1.e-15 .AND. ABS(st) < 1.e-15) GO TO 206
   z(j4+1) = ATAN2(sd,st) * 28.6478898
   GO TO 207
   206 z(j4+1) = 0.
   207 IF (ij < iend) GO TO 202
   GO TO 395
   
!     IS2D8
   
   210 ij = 0
   j2 = 1
   211 ij = ij + 1
   j2 = j2 + 5
   j4 = j2 + 2
   i1 = j2
   i2 = j4
   i3 = 1
   ASSIGN 215 TO iret
   GO TO 390
   215 IF (ij < 8) GO TO 211
   GO TO 395
   
!     IHEX1,2,3
   
   220 i1 = 3
   i2 = 4
   i3 = 1
   ASSIGN 221 TO iret
   GO TO 390
   221 i1 = 11
   IF (eltype == 67) i1 = 12
!                   IHEX3
   i2 = i1+1
   ASSIGN 222 TO iret
   GO TO 390
   222 i1 = i1 + 6
   i2 = i1 + 1
   ASSIGN 223 TO iret
   GO TO 390
   
!     COMPUTE PRINCIPAL STRESSES
   
   223 sig(1) = z( 3)
   sig(2) = z(11)
   sig(3) = z(17)
   sig(4) = z( 4)
   sig(5) = z(12)
   sig(6) = z(18)
   IF (eltype /= 67) GO TO 224
!                   IHEX3
   sig(2) = z(12)
   sig(3) = z(18)
   sig(5) = z(13)
   sig(6) = z(19)
   224 CONTINUE
!*****
!     SOLVE CUBIC EQUATION FOR PRINCIPAL STRESSES
!*****
   
!     S**3 + P*S**2 + Q*S + R = 0.0
   
!     REF. -- CRC STANDARD MATH TABLES 14TH ED., PP. 392,3
   
   rm = 0.0
   DO  i = 1,6
     IF (ABS(sig(i)) > rm) rm = ABS(sig(i))
   END DO
   IF (rm <= 0.0) GO TO 267
   thresh = 1.0E-5
   264 DO  i = 1,6
     IF (ABS(sig(i)/rm) < thresh) sig(i) = 0.0
   END DO
   rx = sx/rm
   ry = sy/rm
   rz = sz/rm
   rxy= sxy/rm
   ryz= syz/rm
   rzx= szx/rm
   p  =-rx - ry - rz
   q  = rx*ry + ry*rz + rz*rx - rxy**2 - ryz**2 - rzx**2
   r  =-(rx*ry*rz +2.0*rxy*ryz*rzx -rx*ryz**2 -ry*rzx**2 -rz*rxy**2)
   a  = (3.0*q - p**2)/3.0
   b  = (2.0*p**3 - 9.0*p*q + 27.0*r)/27.0
   x  =-a**3/27.0
   IF (x > 0.0) GO TO 270
   
!     CHECK FOR IMAGINARY ROOTS
   
   IF (ABS(x) > rm*1.0E-6) GO TO 265
   
!     CHECK FOR 3 EQUAL ROOTS
   
   IF (ABS(b) > 1.0E-6) GO TO 265
   x  = 0.0
   phi= 0.0
   GO TO 275
   265 thresh = 10.0*thresh
   IF (thresh < 1.1E-3) GO TO 264
   267 sa = 0.0
   sb = 0.0
   sc = 0.0
   GO TO 280
   270 cosphi =-(b/2.0)/SQRT(x)
   IF (ABS(cosphi) > 1.0) GO TO 265
   phi= ACOS(cosphi)
   x  = 2.0*SQRT(-a/3.0)
   275 sa = (x*COS(phi/3.0)-p/3.0)*rm
   sb = (x*COS(phi/3.0+120.0*dtor)-p/3.0)*rm
   sc = (x*COS(phi/3.0+240.0*dtor)-p/3.0)*rm
   rm = 0.0
   DO  i = 1,3
     IF (ABS(sigp(i)) > rm) rm = ABS(sigp(i))
   END DO
   DO  i = 1,3
     IF (ABS(sigp(i)/rm) < 1.0E-5) sigp(i) = 0.0
   END DO
!*****
!     COMPUTE MEAN STRESS OR PRESSURE
!*****
   280 sn =-(sa+sb+sc)/3.0
!*****
!     COMPUTE OCTAHEDRAL SHEAR STRESS
!*****
   so = SQRT(((sa+sn)**2 + (sb+sn)**2 + (sc+sn)**2)/3.0)
!*****
!     COMPUTE DIRECTION COSINES OF THE PRINCIPAL PLANES
!*****
   rm = 1.0E-6
   DO  i = 1,3
     IF (sigp(i) == 0.0) GO TO 580
     smat(1,1) = 1.0 - sx/sigp(i)
     smat(2,1) =-sxy/sigp(i)
     smat(3,1) =-szx/sigp(i)
     smat(1,2) = smat(2,1)
     smat(2,2) = 1.0 - sy/sigp(i)
     smat(3,2) =-syz/sigp(i)
     smat(1,3) = smat(3,1)
     smat(2,3) = smat(3,2)
     smat(3,3) = 1.0 - sz/sigp(i)
     CALL saxb (smat(1,1),smat(1,2),DCOS(1,i))
     rx = sadotb(DCOS(1,i),DCOS(1,i))
     j  = 1
     CALL saxb (smat(1,2),smat(1,3),DCOS(1,i))
     ry = sadotb(DCOS(1,i),DCOS(1,i))
     IF (ry > rx) j = 2
     CALL saxb (smat(1,3),smat(1,1),DCOS(1,i))
     rz = sadotb(DCOS(1,i),DCOS(1,i))
     IF (rz > ry .AND. rz > rx) j = 3
     p = smat(1,j)
     q = smat(2,j)
     r = smat(3,j)
     IF (j-2 < 0) THEN
       GO TO   450
     ELSE IF (j-2 == 0) THEN
       GO TO   460
     ELSE
       GO TO   470
     END IF
     450 j = 2
     GO TO 480
     460 j = 3
     GO TO 480
     470 j = 1
     480 s = smat(1,j)
     t = smat(2,j)
     v = smat(3,j)
     IF (ABS(q)  <= rm) GO TO 500
     rx = v - t*r/q
     IF (ABS(rx) <= rm) GO TO 490
     rz =-(s - t*p/q)/rx
     ry =-(p + r*rz)/q
     485 x  = 1.0 + rz*rz + ry*ry
     DCOS(1,i) = 1.0/SQRT(x)
     DCOS(2,i) = ry*DCOS(1,i)
     DCOS(3,i) = rz*DCOS(1,i)
     CYCLE
     490 rx = s - t*p/q
     IF (ABS(rx) <= rm) GO TO 580
     ry =-r/q
     x  = 1.0 + ry*ry
     DCOS(1,i) = 0.0
     DCOS(3,i) = 1.0/SQRT(x)
     DCOS(2,i) = ry*DCOS(3,i)
     CYCLE
     500 IF (ABS(r) <= rm) GO TO 520
     rz = -p/r
     IF (ABS(t) <= rm) GO TO 510
     ry =-(s - v*p/r)/t
     GO TO 485
     510 IF (ABS(s-v*p/r) <= rm) GO TO 580
     DCOS(1,i) = 0.0
     DCOS(2,i) = 1.0
     DCOS(3,i) = 0.0
     CYCLE
     520 IF (ABS(p) <= rm) GO TO 580
     IF (ABS(v) <= rm) GO TO 530
     rz =-t/v
     x  = 1.0 + rz*rz
     DCOS(1,i) = 0.0
     DCOS(2,i) = 1.0/SQRT(x)
     DCOS(3,i) = rz*DCOS(2,i)
     CYCLE
     530 IF (ABS(t) <= rm) GO TO 580
     DCOS(1,i) = 0.0
     DCOS(2,i) = 0.0
     DCOS(3,i) = 1.0
     CYCLE
     580 DCOS(1,i) = 0.0
     DCOS(2,i) = 0.0
     DCOS(3,i) = 0.0
   END DO
   ipts = 0
   IF (eltype == 67) ipts = 1
!                   IHEX3
   z(5) = sa
   z(9) = sn
   z(10)= so
   z(ipts+13) = sb
   z(ipts+19) = sc
   DO  i = 1,3
     z(      5+i) = DCOS(1,i)
     z(ipts+13+i) = DCOS(2,i)
     z(ipts+19+i) = DCOS(3,i)
   END DO
   GO TO 395
   
!     PERFORM NRL SUMS
   
   390 DO  i = i1,i2,i3
     sum  = 0.
     rmax = 0.
     DO  j = 1,nmodes
       isub = nwds*(j-1) + i
       sum  = sum + z(isub)**2
       IF (ABS(z(isub)) > rmax) rmax = ABS(z(isub))
     END DO
     IF (sqrss == 1) rmax = 0.
     sum = sum  - rmax**2
     sum = rmax + SQRT(sum)
     z(i)= sum
   END DO
   
   GO TO iret, (55,65,66,75,85,115,116,165,205,215,221,222,223,395)
   
!     WRITE NRL SUMS TO APPROPRIATE SCRATCH FILE
   
   395 iz(1) = i5
   CALL WRITE (iscr,z,nwds,0)
   400 CONTINUE
 END DO
 
!     DONE WITH THIS ELEMENT.  SINCE WE ARE WRITING IN SORT1, EOR IS
!     NEEDED ON SCRATCH FILE ONLY IF ELEMENT TYPE CHANGES.  THIS WILL BE
!     CHECKED ABOVE.  SKIP EOR ON OES2 AND GO BACK.
 
 FILE = ifil
 CALL fwdrec (*1002,ifil)
 GO TO 30
 
!     EOF ON OES2.  WRITE EOR ON SCRATCH FILE AND COPY THEM TO OUTPUT
!     DATA BLOCK.
 
 410 CALL CLOSE (ifil,1)
 
 DO  i = 2,7
   mcb(i) = 1
 END DO
 DO  i = 1,nshock
   CALL WRITE (scr(i),0,0,1)
   CALL CLOSE (scr(i),1)
   mcb(1) = scr(i)
   CALL wrttrl (mcb)
 END DO
 
 lcore = buf2 - 1
 CALL gopen (ofil,z(buf1),1)
 DO  i = 1,nshock
   CALL gopen (scr(i),z(buf2),0)
   
   430 CALL READ (*690,*440,scr(i),z,lcore,0,iwords)
   CALL WRITE (ofil,z,lcore,0)
   GO TO 430
   
!     EOR
   
   440 CALL WRITE (ofil,z,iwords,1)
   GO TO 430
   
!     EOF
   
   690 CALL CLOSE (scr(i),1)
   
 END DO
 
 CALL CLOSE (ofil,1)
 mcb(1) = ofil
 CALL wrttrl (mcb)
 
!     GO BACK FOR FORCES
 
 710 IF (ifil == oef2) RETURN
 ifil = oef2
 ofil = nrlfor
 GO TO 15
 
 1002 n = -2
 GO TO 1010
 1003 n = -3
 GO TO 1010
 1008 n = -8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE nrlsum
