SUBROUTINE opt2c (pt,iel,ipr,pr,rr)
     
 
 INTEGER, INTENT(IN)                      :: pt(2,1)
 INTEGER, INTENT(IN)                      :: iel(1)
 INTEGER, INTENT(IN OUT)                  :: ipr(1)
 REAL, INTENT(IN)                         :: pr(1)
 REAL, INTENT(IN)                         :: rr(1)
 LOGICAL :: kpun
 INTEGER :: b1, count,eid,eject,est1,est2,etyp,headng,  &
     outtap, iz(100),NAME(2),neop(21), sysbuf,wdopt(42),ycor,zcor,mcb(7),iy(1),  &
     tube,quad4,trim6,tria3
 REAL :: y(1),blk,pcd(2,21),g(2,10),parm(8)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / skp1(2),count,ncard,skp2,ycor,b1,nelop,nwdse,  &
     nwdsp,skp3(2),est1,skp4,est2,nelw,nprw,nklw,ntotl, conv
 COMMON /optpw2/ zcor,z(100)
 COMMON /zzzzzz/ core(1)
 COMMON /names / nrd,noeor,nwrt,nweor
 COMMON /system/ sysbuf,outtap,skps1(6),nlpp,skps2(2),nlines, skps3(78),lpch
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 EQUIVALENCE     (iz(1),z(1)),   (eid,z(1)), (core(1),parm(1),MAX),  &
     (g(1,1),iz(100)), (g(1,10),ig10), (iprnt,parm(7)), (iy(1),y(1),parm(8))
!     EQUIVALENT ARE  (IPR,PR)
 
 
!     NOTE - CHANGE EQUIVALENCE IF AN ELEMENT TO BE OPTIMIZED HAS EST
!     (EPT ONLY) ENTRIES BEYOND 100 WORDS.
 
 DATA   NAME   / 4H opt, 4H2C   /
 DATA   nmes   , yes,plus,blk   /  0, 4HYES , 4H+aaa, 4H    /
 DATA   tube   , quad4,trim6,tria3  / 3, 64 , 73, 83        /
 DATA   pcd    /  &
     4HPBAR,4H    , 4HPELB,4HOW  , 4HPIS2,4HD8  , 4HPQDM,4HEM  ,  &
     4HPQDM,4HEM1 , 4HPQDM,4HEM2 , 4HPQDP,4HLT  , 4HPQUA,4HD1  ,  &
     4HPQUA,4HD2  , 4HPROD,4H    , 4HPSHE,4HAR  , 4HPTRB,4HSC  ,  &
     4HPTRI,4HA1  , 4HPTRI,4HA2  , 4HPTRI,4HM6  , 4HPTRM,4HEM  ,  &
     4HPTRP,4HLT  , 4HPTUB,4HE   , 4HPSHE,4HLL  , 4HPSHE,4HLL  , 4HYYYY,4H    /
 
!     POINTERS TO WORDS ON EST TO CONVERT.  NEOP(ITP) IS POINTER INTO
!     -WDOPT- ARRAY.  THE -WDOPT- FIRST ENTRY FOR THE ELEMENT IS THE
!     NUMBER OF ENTRIES ON -EST- TO CONVERT FOLLOWED BY THE WORD NUMBERS
!     TO OPTIMIZE.
 
 DATA   neop   / 21,30,39,15,15,  15,27,17,15,1,  6,12,8,6,35,  &
     6,12, 4,41,41,   0/
 DATA wdopt / & 
!     ROD (A,J)  &
 2,   5,6 & 
!     TUBE (O.D.)  &
 , 1,   5 &  
!     SHEAR(T), TRMEM(T), TRIA2(T)  &
 , 1,   7 & 
!     TRIA1(T1,T2,I)  &
 , 3,   7,9,11 & 
!     TRBSC(T2,I),TRPLT(T2,I)  &
 , 2,   7,9 & 
!     QDMEM(T), QDMEM1(T), QDMEM2(T), QUAD2(T)  &
 , 1,   8 & 
!     QUAD1(T1,T2,I)  &
 , 3,   8,10,12 & 
!     BAR(A,J,I1,I2,I12)  &
 , 5,   17,18,19,20,33 & 
!     QDPLT(T2,I)  &
 , 2,   8,10 & 
!     ELBOW(A,J,I1,I2)  &
 , 4,   9,10,11,12 & 
!     TRIM6(T1,T3,T5)  &
 , 3,   10,11,12 & 
!     IS2D8(T)  &
 , 1,   13 & 
!     QUAD4(T), TRIA3(T) PSHELL ONLY  &
 , 1,   14 /
 
!     DETERMINE IF PROPTETY CARDS ARE TO BE PUNCHED
 
 kpun = .false.
 kount  = 0
 headng = 0
 ch     = 1.0
 icp    = ntotl
 IF (count == MAX .OR. conv == 2.0) kpun =.true.
 IF (parm(5) /= yes) kpun = .false.
 IF (iprnt /= 0) nlines = nlpp
 ie2 = 1
 lel = 0
 
!     READ EST1 ELEMENT TYPE
 
 10 CALL READ (*400,*360,est1,etyp,1,noeor,i)
 CALL WRITE (est2,etyp,1,noeor)
 itp = iy(etyp)
 IF (itp == 0) GO TO 20
 ie1 = pt(1,itp)
 
!     CHECK IF CORE ELEMENTS SKIPPED BECAUSE TYPE NOT ON EST
 
 IF (ie1 > ie2) GO TO 60
 ie2 = pt(1,itp+1)
 lel = iel(ie1)
 ip1 = pt(2,itp) - 1
 IF (ie2 > ie1) GO TO 40
 
!     SKIP THIS ELEMENT TYPE.  COPY RECORD TO EST2
 
 20 j = 1
 n = zcor
 CALL READ (*30,*30,est1,z,zcor,noeor,n)
 j = 0
 30 CALL WRITE (est2,z(1),n,j)
 IF (j == 0) THEN
   GO TO    20
 ELSE
   GO TO    10
 END IF
 
!     ELEMENT TYPE HAS CORE ENTRIES
 
 40 CONTINUE
 nwds = incr*(etyp-1) + 12
 nwds = NE(nwds)
 npcard = 0
 IF (nwds > zcor) CALL mesage (-8,zcor,NAME)
 
!     READ ONE EST1 ELEMENT INTO CORE
 
 50 CALL READ (*350,*340,est1,z,nwds,noeor,i)
 IF (eid-lel < 0.0) THEN
   GO TO    55
 ELSE IF (eid-lel == 0.0) THEN
   GO TO    80
 ELSE
   GO TO    60
 END IF
 
!     ELEMENT ID NOT IN CORE
 
 55 CALL WRITE (est2,iz(1),nwds,noeor)
 GO TO 50
 
!     ELEMENT IN CORE NOT ON EST
 
 60 i = eject(2)
 IF (i == 0) GO TO 68
 IF (count == MAX .OR. conv == 2.0) GO TO 66
 WRITE  (outtap,65) count
 65 FORMAT (1H0,8X,45HPROPERTIES used during intermediate iteration,  &
     i5, 10H by optpr2/)
 GO TO 68
 66 WRITE  (outtap,67) count
 67 FORMAT (1H0,8X,38HPROPERTIES used during final iteration,  &
     i5, 10H by optpr2/)
 68 WRITE  (outtap,70) sfm,etyp,lel,NAME
 70 FORMAT (a25,' 2297, INCORRECT LOGIC FOR ELEMENT TYPE',i4,  &
     ', ELEMENT',i8,2H (,2A4,2H).)
 CALL mesage (-61,lel,NAME)
 
!     ELEMENT IN CORE - CONVERT THE ENTRIES
 
 80 ipl = iel(ie1+4) + ip1
 ie1 = ie1 + nwdse
 lel = iel(ie1)
 IF (ie1 > ie2) lel = 100000000
 a = pr(ipl+4)
 IF (a > 0.0) GO TO 100
 nmes = nmes + 1
 IF (iprnt == 0 .OR. nmes > 100) GO TO 55
 i = eject (2)
 IF (i == 0) GO TO 88
 IF (count == MAX .OR. conv == 2.0) GO TO 84
 WRITE (outtap,65) count
 GO TO 88
 84 WRITE  (outtap,65) count
 88 WRITE  (outtap,90) uim,eid
 90 FORMAT (a29,' 2305, OPTPR2 DETECTED NEGATIVE ALPHA FOR ELEMENT', i8)
 GO TO 55
 
 100 locf = neop(itp)
 j = locf
 k = wdopt(locf)
 irr = (ipl+nwdsp)/nwdsp
 IF (ABS(parm(3)-1.0) < 0.0001) ch = 0.25*rr(irr) + 0.75
 c = (a/(a+(1.0-a)*parm(3)))**ch
 IF (etyp /= trim6) GO TO 105
 
!     SPECIAL HANDLING FOR TRIM6
!     IF THICKNESS-3 OR THICKNESS-5 IS ZERO, SET EQUAL TO THICKNESS-1
 
 DO  jj = 1,k
   j = j +1
   l = wdopt(j)
   IF (jj /= k .AND. ABS(z(l+1)) < 1.e-7) z(l+1) = z(l)
   pc = y(icp+jj)
   z(l) = z(l)*(pc/(pc+(1.0-pc)*parm(3)))
 END DO
 icp = icp + 4
 GO TO 115
 
 105 DO  i = 1,k
   j = j + 1
   l = wdopt(j)
   z(l) = c*z(l)
 END DO
 IF (etyp /= quad4 .AND. etyp /= tria3) GO TO 112
 z(l+6) =  0.5*z(l)
 z(l+7) = -0.5*z(l)
 112 IF (etyp == tube .AND. z(l) < 2.*z(l+1)) z(l+1) = .5*z(l)
 115 CALL WRITE (est2,z(1),nwds,noeor)
 
!     PUNCH AND/OR PRINT PROPERTY CARDS
 
 IF (iprnt == 0 .OR. ipr(ipl) <= 0) GO TO 50
 SELECT CASE ( itp )
   CASE (    1)
     GO TO 120
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 140
   CASE (    4)
     GO TO 150
   CASE (    5)
     GO TO 150
   CASE (    6)
     GO TO 150
   CASE (    7)
     GO TO 160
   CASE (    8)
     GO TO 170
   CASE (    9)
     GO TO 150
   CASE (   10)
     GO TO 180
   CASE (   11)
     GO TO 150
   CASE (   12)
     GO TO 160
   CASE (   13)
     GO TO 170
   CASE (   14)
     GO TO 150
   CASE (   15)
     GO TO 180
   CASE (   16)
     GO TO 150
   CASE (   17)
     GO TO 160
   CASE (   18)
     GO TO 190
   CASE (   19)
     GO TO 170
   CASE (   20)
     GO TO 170
 END SELECT
 
!     PBAR
 
 120 k1 = 02222211
 k2 = 22222222
 k3 = 00000222
 GO TO 250
 
!     PELBOW
 
 130 k1 = 02222211
 k2 = 22222222
 k3 = 22222222
 GO TO 250
 
!     PIS2D8
 
 140 k1 = 00000211
 GO TO 230
 
!     PQDMEM, PQDMEM1, PQDMEM2, PQUAD2, PSHEAR, PTRIA2, PTRMEM
 
 150 k1 = 00002211
 GO TO 230
 
!     PQDPLT, PTRBSC, PTRPLT
 
 160 k1 = 22221211
 GO TO 230
 
!     PQUAD1, PTRIA1, PSHELL
 
 170 k1 = 22121211
 k2 = 00000022
 GO TO 240
 
!     PROD, PTRIM6
 
 180 k1 = 00222211
 GO TO 230
 
!     PTUBE
 
 190 k1 = 00022211
 
!     OUTPUT THE CARD(S)
 
 230 k2 = 0
 240 k3 = 0
 250 ii = wdopt(locf+1) - 4
 kk = k1
 g(1,1)   = pcd(1,itp)
 g(2,1)   = pcd(2,itp)
 iz(ii+2) = ipr(ipl)
 ipr(ipl) =-ipr(ipl)
 260 DO  i = 2,9
   g(1,i) = blk
   g(2,i) = blk
   j = MOD(kk,10)
   IF (j == 0) GO TO 270
   IF (j == 1) CALL int2a8 (*370,iz(i+ii),g(1,i))
   IF (j == 2) CALL  fp2a8 (*380, z(i+ii),g(1,i))
   270 kk = kk/10
 END DO
 g(1,10) = blk
 g(2,10) = blk
 IF (k2 == 0 .OR. (k2 == -1 .AND. k3 == 0) .OR. k3 == -1) GO TO 320
 kount = kount + 1
 CALL int2a8 (*375,kount,g(1,10))
 g(2,10) = g(1,10)
 ig10 = khrfn3(g(1,1),plus,-3,1)
 IF (headng == 0) GO TO 320
 280 WRITE  (outtap,290) g
 290 FORMAT (5X,10(2A4,1X))
 IF (.NOT.kpun) GO TO 300
 WRITE  (lpch,295) g
 295 FORMAT (20A4)
 ncard = ncard + 1
 
!     SET UP FOR CONTINUATION CARD(S)
 
 300 IF (k2 == 0 .OR. (k2 == -1 .AND. k3 == 0) .OR. k3 == -1) GO TO 50
 g(1,1) = g(1,10)
 g(2,1) = g(2,10)
 ii     = ii + 8
 IF (k2 < 0) THEN
   GO TO   315
 ELSE IF (k2 == 0) THEN
   GO TO    50
 END IF
 310 kk = k2
 k2 = -1
 GO TO 260
 315 kk = k3
 k3 = -1
 GO TO 260
 
!     PRINT HEADING
 
 320 headng = 1
 IF (eject(1) == 0) GO TO 280
 IF (count == MAX .OR. conv == 2.0) GO TO 330
 WRITE (outtap,65) count
 GO TO 280
 330 WRITE (outtap,67) count
 GO TO 280
 
!     EOR ON EST1
 
 340 CALL WRITE (est2,0,0,nweor)
 IF (ie1-ie2 < 0) THEN
   GO TO    60
 ELSE
   GO TO    10
 END IF
 
!     ERRORS
 
 350 CALL mesage (-2,est1,NAME)
 360 CALL mesage (-3,est1,NAME)
 370 j = 370
 GO TO 390
 375 j = 375
 i = kount
 GO TO 390
 380 j = 380
 390 WRITE  (outtap,395) j,g(1,1),g(2,1),i,ii,iz(i+ii),z(i+ii)
 395 FORMAT (16H0*** opt2c/error,i5,9X,5HELEM ,2A4,3I9,e10.4 )
 GO TO 50
 
 400 CALL eof (est2)
 mcb(1) = est1
 CALL rdtrl(mcb)
 mcb(1) = est2
 CALL wrttrl(mcb)
 RETURN
END SUBROUTINE opt2c
