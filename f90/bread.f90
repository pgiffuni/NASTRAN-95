SUBROUTINE bread (ig,inv,ii3,norig,kg)
     
!      THIS ROUTINE IS USED ONLY IN BANDIT MODULE
!      IT READS THE CONNECTING ELEMENTS AND GENEL ELEM. FROM GEOM2 FILE
!      AND PREPROCESS THE MPC CARDS AND THE RIGID ELEMENTS FROM GEOM4
 
!      REVISED BY G.CHAN/UNISYS
!      12/89, TO INCLUDE NEW RIGID ELEMENTS CRROD, CRBAR, CRTRPLT,
!      CRBE1, CREB2, CRBE3 AND CRSPLINE
!      03/92, TO INCLUDE DUMMY ELEMENTS, CDUM1,...,CDUM9
 
 IMPLICIT INTEGER (a-z) 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(OUT)                     :: inv(ii3,1)
 INTEGER, INTENT(IN OUT)                  :: ii3
 INTEGER, INTENT(IN OUT)                  :: norig(1)
 INTEGER, INTENT(OUT)                     :: kg(7)
 LOGICAL :: debug
 DIMENSION        sub(2),   xxx(3),   iz(3)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,      uwm
 COMMON /banda /  ibuf1,    nompc
 COMMON /bandb /  nbitin,   kore,     ifl,      ngrid,    ipnw(2), kdim
 COMMON /bandd /  dum6(6),  nel,      neq,      neqr
 COMMON /bands /  nn(10)
 COMMON /geomx /  geom1,    geom2,    geom4,    scr1
 COMMON /names /  rd,       rdrew,    wrt,      wrtrew,   rew
 COMMON /system/  ibuf,     nout,     dum43(43),kdum(9)
 COMMON /gpta1 /  NE,       last,     incr,     ke(1)
 COMMON /zzzzzz/  z(1)
 DATA             crigdr,   crigd1,   crigd2,   crigd3,   genel  /  &
     8210,     5310,     5410,     8310,     4301   /
 DATA             chbdy,    plotel,   crrod,    crbar,    crtrpt /  &
     4208,     5201,     6510,     6610,     6710   /
 DATA             crbe1,    crbe2,    crbe3,    crspln,   mset   /  &
     6810,     6910,     7010,     7110,     4HMSET /
 DATA             sub,                mpc,      maxmpc,   debug  /  &
     4HBREA,   4HD   ,   4901,     150,      .false./
 
 
!     CHECK THE PRESENCE OF GEOM2 FILE
 
 kg(1) = geom2
 CALL rdtrl (kg(1))
 j = kg(2) + kg(3) + kg(4) + kg(5) + kg(6) + kg(7)
 IF (kg(1) < 0 .OR. j == 0) GO TO 370
 DO  i = 1,7
   kg(i) = 0
 END DO
 
!     UPDATE /GPTA1/ IF DUMMY ELEMENTS ARE PRESENT
 
 DO  i = 1,9
   IF (kdum(i) == 0) CYCLE
   k = kdum(i)/10000000
   l = (kdum(i)-k*10000000)/10000
   j = (i+51)*incr
   ke(j+ 6) = 2 + k + l
   ke(j+10) = k
 END DO
 
!     CHECK THE PRESENCE OF MPC CARDS AND RIGID ELEMENTS.  SAVE THEIR
!     GRID DATA IN SCR1 FILE FOR TIGER AND UPDATE NEQ AND NEQR COUNTERS
 
 IF (nompc == 0) GO TO 200
 z(1) = geom4
 CALL rdtrl (z(1))
 j = 0
 DO  i = 2,7
   j = j + z(i)
 END DO
 IF (z(1) < 0 .OR. j == 0) GO TO 200
 
 ibuf2 = ibuf1 - ibuf
 CALL OPEN (*290,scr1,z(ibuf2),wrtrew)
 ifile  = geom4
 CALL preloc (*190,z(ibuf1),geom4)
 
 IF (nompc == 1) GO TO 40
 
 xxx(1) = mpc
 xxx(2) = xxx(1)/100
 CALL locate (*40,z(ibuf1),xxx,j)
 25   j = 1
 CALL READ (*300,*40,geom4,iz,1,0,m)
 30   j = j + 1
 CALL READ (*300,*40,geom4,kg(j),3,0,m)
 IF (kg(j) /= -1) IF (j+3-maxmpc) 30,30,320
 j = j - 1
 kg(1) = j - 1
 CALL WRITE (scr1,kg,j,1)
 neq = neq + 1
 GO TO 25
 
!     LOCATE ANY CRIGDR AND CRROD ELEMENTS, AND SAVE THE GRID DATA IN
!     SCR1. (DEPENDENT GRID FIRST, AND ONE INDEPENDENT GRID LAST)
 
!     FOR ALL RIGID ELEMENTS, THE FIRST WORD OF KG ARRAY CONTAINS
!     (NO. OF DEPENDENT + INDEP. GRIDS)*1000 + (NO. OF INDEP. GRIDS)
!     THE DATA IN SCR1 WILL BE PROCESSED BY TIGER
 
 40   IF (nompc == 3) GO TO 180
 xxx(1) = crigdr
 50   xxx(2) = xxx(1)/100
 CALL locate (*60,z(ibuf1),xxx,j)
 55   CALL READ (*300,*60,geom4,iz,1,0,m)
 CALL READ (*300,*60,geom4,kg(3),3,0,m)
 kg(1) = 2*1000 + 1
 kg(2) = kg(4)
 CALL WRITE (scr1,kg,3,1)
 neqr  = neqr + 1
 GO TO 55
 
 60   IF (xxx(1) == crrod) GO TO 70
 xxx(1) = crrod
 GO TO 50
 
!     LOCATE ANY CRIGD1, CRIGD2  AND CRBE2  ELEMENTS, AND SAVE GRID
!     DATA IN SCR1. PUT THE ONE INDEPENDENT GRID LAST
 
 70   xxx(1) = crigd1
 75   xxx(2) = xxx(1)/100
 CALL locate (*90,z(ibuf1),xxx,j)
 80   j = 1
 CALL READ (*300,*90,geom4,iz,2,0,m)
 iz2 = iz(2)
 85   j = j + 1
 CALL READ (*300,*90,geom4,kg(j),1,0,m)
 CALL READ (*300,*90,geom4,   0,-6,0,m)
 IF (kg(j) /= -1) IF (j-maxmpc) 85,85,320
 kg(j) = iz2
 kg(1) = (j-1)*1000 + 1
 CALL WRITE (scr1,kg,j,1)
 neqr  = neqr + 1
 GO TO 80
 90   IF (xxx(1) == crbe2) GO TO 110
 
!     LOCATE ANY CRIGD2 ELEMENT
 
 IF (xxx(1) == crigd2) GO TO 100
 xxx(1) = crigd2
 GO TO 75
 
!     LOCATE ANY CRBE2 ELEMENT
 
 100  xxx(1) = crbe2
 GO TO 75
 
!     LOCATE ANY CRIGD3, CRBE1, CRBAR AND CRTRPLT ELEMENTS, AND SAVE
!     GRID DATA IN SCR1 FILE. PUT THE INDEPENDENT GRID LAST
 
 110  xxx(1) = crbar
 ASSIGN 115 TO irtn
 GO TO 130
 115  xxx(1) = crtrpt
 ASSIGN 120 TO irtn
 GO TO 130
 120  xxx(1) = crbe1
 ASSIGN 125 TO irtn
 GO TO 130
 125  xxx(1) = crigd3
 ASSIGN 150 TO irtn
 130  xxx(2) = xxx(1)/100
 CALL locate (*145,z(ibuf1),xxx,j)
 133  j = 2
 k = 1
 CALL READ (*300,*145,geom4,iz,1,0,m)
 135  CALL READ (*300,*145,geom4,iz(k),1,0,m)
 IF (iz(k) == mset) GO TO 137
 CALL READ (*300,*145,geom4,0,-6,0,m)
 k = k + 1
 IF (k > 999) GO TO 340
 GO TO 135
 137  CALL READ (*300,*145,geom4,kg(j),1,0,m)
 CALL READ (*300,*145,geom4,0,-6,0,m)
 IF (kg(j) == -1) GO TO 140
 j = j + 1
 IF (j > maxmpc) GO TO 320
 GO TO 137
 140  k = k - 1
 DO  i = 1,k
   kg(j) = iz(i)
   j = j + 1
 END DO
 j = j - 1
 kg(1) = (j-1)*1000 + k
 CALL WRITE (scr1,kg,j,1)
 neqr = neqr + 1
 GO TO 133
 
!     LOCATE ANY CRSPLINE ELEMENTS, AND SAVE GRID DATA IN SCR1 FILE.
!     PUT THE INDEPENDENT GRIDS LAST
 145  GO TO irtn, (115,120,125,150)
 
!     LOCATE ANY CRBE3 ELEMENTS, AND SAVE GRID DATA IN SCR1 FILE. PUT
!     THE INDEPENDENT GRID LAST
 
 150  xxx(1) = crbe3
 xxx(2) = xxx(1)/100
 CALL locate (*165,z(ibuf1),xxx,j)
 151  CALL READ (*300,*165,geom4,iz,3,0,m)
 iz2 = iz(2)
 j = 2
 CALL READ (*300,*165,geom4,0,-2,0,m)
 153  CALL READ (*300,*165,geom4,kg(j),1,0,m)
 k = -kg(j)
 IF (k > 0) THEN
    SELECT CASE ( k )
     CASE (    1)
       GO TO 155
     CASE (    2)
       GO TO 157
     CASE (    3)
       GO TO 160
   END SELECT
 END IF
 j = j + 1
 IF (j-maxmpc > 0) THEN
   GO TO   320
 ELSE
   GO TO   153
 END IF
 155  CALL READ (*300,*165,geom4,i,1,0,m)
 IF (i == -2) GO TO 157
 CALL READ (*300,*165,geom4,0,-1,0,m)
 GO TO 153
 157  CALL READ (*300,*165,geom4,kg(j),1,0,m)
 IF (kg(j) < 0) GO TO 160
 CALL READ (*300,*165,geom4,0,-1,0,m)
 j = j + 1
 GO TO 157
 160  kg(j) = iz2
 kg(1) = (j-1)*1000 + 1
 CALL WRITE (scr1,kg,j,1)
 neqr = neqr + 1
 GO TO 151
 
!     LOCATE ANY CRSPLINE ELEMENTS, AND SAVE GRID DATA IN SCR1 FILE.
!     PUT THE INDEPENDENT GRIDS LAST
 
 165  xxx(1) = crspln
 xxx(2) = xxx(1)/100
 CALL locate (*180,z(ibuf1),xxx,j)
 167  CALL READ (*300,*180,geom4,iz,3,0,m)
 k = 1
 iz(k) = iz(3)
 j = 1
 170  j = j + 1
 173  CALL READ (*300,*175,geom4,kg(j),2,0,m)
 IF (kg(j) == -1) GO TO 175
 IF (j+2 > maxmpc) GO TO 320
 IF (kg(j+1) /=  0) GO TO 170
 k = k + 1
 IF (k > 999) GO TO 340
 iz(k) = kg(j)
 GO TO 173
 175  DO  i = 1,k
   kg(j) = iz(i)
   j = j + 1
 END DO
 j = j - 1
 kg(1) = (j-1)*1000 + k
 CALL WRITE (scr1,kg,j,1)
 neqr = neqr + 1
 GO TO 167
 
 180  DO  k = 1,maxmpc
   kg(k) = 0
 END DO
 190  CALL CLOSE (geom4,rew)
 CALL CLOSE (scr1,rew)
 
!     PROCESS ELEMENT CARDS AND FILL UP CONNECTION TABLE IG
 
 200  ifile = geom2
 CALL preloc (*300,z(ibuf1),geom2)
 ielem = 1 - incr
 205  ielem = ielem + incr
 IF (ielem > last) GO TO 250
 IF (ke(ielem+3) == chbdy ) GO TO 205
 IF (ke(ielem+3) == plotel) GO TO 205
 scalar = ke(ielem+10)
 IF (scalar == -1) GO TO 205
 CALL locate (*205,z(ibuf1),ke(ielem+3),j)
 nwds  = ke(ielem+ 5)
 ngpts = ke(ielem+ 9)
 ngpt1 = ke(ielem+12)
 ncon  = ngpts
 210  CALL READ (*300,*205,geom2,kg(1),nwds,0,m)
 IF (scalar == 0) GO TO 220
 IF (kg(5) == 0 .OR. kg(6) == 0) GO TO 210
!     THE ABOVE CONDITIONS HOLD TRUE FOR CDAMPI, CELASI, AND CMASSI
!     WHERE I = 1,2
 220  nel = nel + 1
 CALL scat (kg(ngpt1),ncon,inv,ii3,norig)
 IF (ngrid == -1) GO TO 270
 IF (ncon  <=  1) GO TO 240
 ngpt2 = ngpt1 + ncon - 1
 k = ngpt2 - 1
 DO  i = ngpt1,k
   l = i + 1
   DO  j = l,ngpt2
     CALL setig (kg(i),kg(j),ig,norig)
   END DO
 END DO
 240  IF (ielem-last > 0) THEN
   GO TO   255
 ELSE
   GO TO   210
 END IF
 
!     SPECIAL TREATMENT FOR GENERAL ELEM.
!     (LIMITED TO KDIM*4 GRID POINTS PER GENEL)
 
 250  xxx(1) = genel
 xxx(2) = xxx(1)/100
 CALL locate (*270,z(ibuf1),xxx,j)
 kdim4  = kdim*4
 255  ntot   = 0
 CALL READ (*300,*270,geom2,k,1,0,m)
 k    = 0
 kgpv = 0
 GO TO 263
 260  IF (kg(ncon) == kgpv) GO TO 265
 kgpv = kg(ncon)
 263  ntot = ntot + 1
 IF (ntot < kdim4) ncon = ntot
 265  CALL READ (*300,*270,geom2,kg(ncon),2,0,m)
 IF (kg(ncon) /= -1) IF (kg(ncon+1)) 260,  265,  260
!                                           GRD  SCALAR GRD
!                                           PT.   PT.   PT.
 k = k + 1
 xxx(k) = kg(ncon+1)
 IF (k < 2) GO TO 265
 ncon = ncon - 1
 m    = xxx(1)
 nwds = 1 + (m*m-m)/2 + m
 CALL READ (*300,*270,geom2,k,-nwds,0,m)
 CALL READ (*300,*270,geom2,k,   1,0,m)
 ngpt1 = 1
 IF (k == 0) GO TO 220
 nwds  = m*xxx(2)
 CALL READ (*300,*270,geom2,k,-nwds,0,m)
 GO TO 220
 270  CALL CLOSE (geom2,rew)
 IF (ntot > kdim4) GO TO 330
 IF (.NOT.debug) RETURN
 
 m = nn(1)
 WRITE  (nout,280) nn
 WRITE  (nout,285) ((inv(i,j),j=1,2),i=1,m)
 280  FORMAT (//21H /bands/ from bread =,10I8)
 285  FORMAT (/12H table inv =,(/10X,2I8))
 RETURN
 
 290  ifile = scr1
 300  CALL mesage (-1,ifile,sub)
 320  WRITE  (nout,325) uwm,iz(1),maxmpc
 325  FORMAT (a25,', MPC SET (OR CRIGID ID)',i9,  &
     ' IS TOO LONG,  ONLY THE FIRST',i4, /5X,  &
     ' GRID POINTS ARE USED IN THE BANDIT COMPUTATION')
 GO TO 180
 330  WRITE  (nout,335) ufm,ntot
 335  FORMAT (a23,', GENEL ELEMENT HAS TOO MANY GRID POINTS,',i7)
 j = ntot/400 + 1
 IF (j <= 9) WRITE (nout,336) j
 336  FORMAT (5X,'USER NEEDS TO ADD A ''NASTRAN BANDTDIM=',i1,  &
     ''' CARD AND RERUN JOB')
 GO TO 350
 340  WRITE  (nout,345)
 345  FORMAT ('0*** MORE THAN 1000 INDEPENDENT GRID POINTS USED IN A ',  &
     'RIGID ELEMENT')
 350  CALL mesage (-61,0,0)
 
 370  ngrid = 0
 RETURN
END SUBROUTINE bread
