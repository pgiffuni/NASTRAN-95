SUBROUTINE ta1cpd
     
!     G3 MATRIX CALCULATION WITH NEW FORMULATION
 
!     THIS ROUTINE IS CALLED IN TA1 IF PARAM COMPS IS SET TO -1
!     INDICATING PCOMP, PCOMP1 OR PCOMP2 BULK DATA ENTRIES ARE
!     PRESENT. IT'S PRIMARY FUNCTION IS TO -
!       1. CREATE FILE PCOMPS WHICH WILL CONTAIN THE ECHO OF THE
!          'PCOMPS' ENTRIES ALONG WITH INDIVIDUAL LAYER INTRINISIC
!          PROPERTY MATRICES.
!       2. CALCULATE OVERALL MATERIAL PROPERTIES IN THE FORM OF MAT2
!          ENTRIES AND WRITE TO FILE MPTX.
!       3. GENERATE EQUIVALENT PSHELL PROPERTY ENTRIES AND WRITE TO
!          FILE EPTX.
 
 EXTERNAL         andf,orf
 LOGICAL :: ok uai
 INTEGER :: pcomp(2),pcomp1(2),pcomp2(2),comps,pcbit(3),eptx,  &
     pshlpr,eptwds,pshbit,rd,rdrew,wrt,wrtrew,clsrew,  &
     cls,ipshel(17),pshnam(3),pcompr,typc,typc1,typc2,  &
     flag,ept,pcomps,eof,elid,pidloc,eoeloc,sym,z,  &
     symmem,eoe,sysbuf,pos,pos1,buf0,buf1,buf2,buf3,  &
     buf4,buf5,FILE,INDEX(6,3),indexx(3,3),andf,orf, BLANK
 DIMENSION        rz(1),npcmp(3),npcmp1(3),npcmp2(3),nam(2),  &
     nam1(2),nam2(2),matnam(3),ipcomp(7),imembr(17),  &
     ibendg(17),imembd(17),itrshr(17),imptx(7), ieptx(7)
 REAL :: glay(25),gmembr(17),gbendg(17),gmembd(17),  &
     gtrshr(17),exx,eyy,eixx,eiyy,zx,zy,rpshel(17), alfa1,alfa2,alfa12,tref,gsube
 DOUBLE PRECISION :: theta,thetar,c,c2,c4,s,s2,s4,pi,twopi,raddeg,  &
     degrad,t(9),GT(9),gbr(9),gbar(3,3),g(25),  &
     gd(9),gdt(9),gdbr(9),gdbar(3,3),gd2(3,3),u(9),  &
     g3i(9),g3iu(9),g3br(9),g3bar(3,3),g1(3,3),  &
     g2(3,3),g3(2,2),g4(3,3),tlam,zk,zk1,zref,  &
     zg1,zg2,zg4,zi,ti,rho,detrmn,const,zbarx,zbary,  &
     trflx(2,2),zbarxt,zbarxb,zbaryt,zbaryb,zbar(2),  &
     gtrflx(2,2),g3invd(2,2),ex,ey,e(2),fi(2),  &
     fii(2),ri(2),determ,dum(6),dummy(3),stiff(6,6), ei(2),gd1(3,3),gd4(3,3),epsi
 COMMON /BLANK /  luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /ta1com/  nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4
 COMMON /names /  rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /machin/  mach
 COMMON /system/  sysbuf,nout  ,nogo  ,dm(20),icfiat
 COMMON /matin /  matid ,inflag,eltemp
 COMMON /matout/  rmtout(25)
 COMMON /zzzzzz/  z(1)
 COMMON /condad/  pi    ,twopi ,raddeg,degrad
 COMMON /two   /  two(32)
 EQUIVALENCE      (z(1)     ,rz(1)    ), (ipshel(1),rpshel(1)),  &
     (imembr(1),gmembr(1)), (ibendg(1),gbendg(1)),  &
     (imembd(1),gmembd(1)), (itrshr(1),gtrshr(1))
!     DATA    MPT   /  107/
!     DATA    MPTX  /  206/
!     DATA    PCOMPS/  207/
!     DATA    EPTX  /  208/
 DATA    pcomp /  5502,55/
 DATA    pcomp1/  5602,56/
 DATA    pcomp2/  5702,57/
 DATA    npcmp /  5502,55,280/
 DATA    npcmp1/  5602,56,281/
 DATA    npcmp2/  5702,57,282/
 DATA    pshnam/  5802,58,283/
 DATA    matnam/  203, 2, 78 /
 DATA    pcbit /  55, 56, 57 /
 DATA    pshbit/  58/
 DATA    i1st  /  1 /
 DATA    sym   /  1 /
 DATA    mem   /  2 /
 DATA    symmem/  3 /
 DATA    mt2bit/  2 /
 DATA    eoe   /  -1/
 DATA    nam   /  4HTA1C, 4HPD   /
 DATA    nam2  /  4HPCOM, 4HPS   /
 DATA    BLANK /  4HBLNK         /
 DATA    ok uai/  .true.         /
 DATA    epsi  /  1.0D-15        /
 
 buf0 = korsz(z) - sysbuf - 2
 buf1 = buf0 - sysbuf - 2
 buf2 = buf1 - sysbuf - 2
 buf3 = buf2 - sysbuf - 2
 buf4 = buf3 - sysbuf - 2
 buf5 = buf4 - sysbuf - 2
 
!     PERFORM GENERAL INITILIZATION
 
 matwds = 0
 eof    = 0
 elid   = 0
 mat2pr = 0
 pshlpr = 0
 icount = 0
 rho    = 0.0D0
 IF (mach == 2) epsi = 1.0D-12
 
!     OPEN EPTX AND WRITE HEADER RECORD
 
 FILE = eptx
 CALL OPEN  (*1200,eptx,z(buf0),wrtrew)
 CALL fname (eptx,nam1)
 CALL WRITE (eptx,nam1,2,1)
 
!     OPEN MPTX AND WRITE HEADER RECORD
 
 FILE = mptx
 CALL OPEN  (*1200,mptx,z(buf1),wrtrew)
 CALL fname (mptx,nam1)
 CALL WRITE (mptx,nam1,2,1)
 
!     OPEN MPT AND POSITION FILE
 
 FILE = mpt
 CALL OPEN (*1200,mpt,z(buf2),rdrew)
 CALL fwdrec (*1200,mpt)
 
!     OPEN PCOMPS AND WRITE HEADER RECORD
!     WRITE TO IPCOMP(1), THE GINO FILE NAME OF PCOMPS
 
 FILE = pcomps
 CALL OPEN  (*1200,pcomps,z(buf3),wrtrew)
 CALL WRITE (pcomps,nam2,2,1)
 
 ipcomp(1) = pcomps
 DO  ll = 2,7
   ipcomp(ll) = 0
 END DO
 
!     COPY ALL EPT ENTRIES UP TO PSHELL TYPE TO FILE EPTX
!     IF NONE FOUND, MUST CREATE ONE BEFORE THE LAST RECORD IN FILE
 
!     SET AVAILABLE CORE
 
 n = buf5 - 1
 iept = i1st
 FILE = ept
 CALL OPEN (*1200,ept,z(buf4),rdrew)
 CALL fwdrec (*1200,ept)
 irec = 0
 20 CALL fwdrec (*30,ept)
 irec = irec + 1
 GO TO 20
 
 30 CALL REWIND (ept)
 CALL fwdrec (*1200,ept)
 ired = 0
 40 CALL READ (*1200,*50,ept,z(iept),n,1,eptwds)
 CALL mesage (-8,0,nam)
 50 IF (z(iept) == 4902) GO TO 60
 ired = ired + 1
 IF (ired == irec) GO TO 70
 CALL WRITE (eptx,z(iept),eptwds,1)
 eptwds = 0
 GO TO 40
 
 60 pshlpr = 1
 70 CALL bckrec (ept)
 CALL savpos (ept,pos1)
 CALL CLOSE  (ept,clsrew)
 
!     OPEN EPT
 
 FILE = ept
 CALL preloc (*1200,z(buf4),ept)
 
!     COPY ALL MAT ENTRIES UP TO MAT2 TYPE TO FILE MPTX
 
!     SET AVAILABLE CORE
 
 n = buf5 - 1
 imat = i1st
 80 CALL READ (*110,*90,mpt,z(imat),n,1,matwds)
 CALL mesage (-8,0,nam)
 90 IF (z(imat) >= 203) GO TO 100
 CALL WRITE (mptx,z(imat),matwds,1)
 matwds = 0
 GO TO 80
 100 CALL bckrec (mpt)
 CALL savpos (mpt,pos)
 IF (z(imat) == 203) mat2pr = 1
 GO TO 120
 
!     SET END OF FILE FLAG
 
 110 eof = 1
 
!     CLOSE MPT BEFORE CALLING PREMAT
 
 120 CALL CLOSE (mpt,1)
 
!     SET POINTERS AND PERFORM INITILIZATION
 
 ipc1  = 1
 npc   = 0
 npc1  = 0
 npc2  = 0
 typc  = 0
 typc1 = 0
 typc2 = 0
 
!     SET SIZE OF AVAILABLE CORE
 
 n   = buf5 - 1
 ipc = 1
 
!     LOCATE PCOMP DATA AND READ INTO CORE
 
 CALL locate (*140,z(buf4),pcomp,flag)
 
 CALL READ (*1200,*130,ept,z(ipc),n,0,npc)
 CALL mesage (-8,0,nam)
 130 IF (npc > 0) typc = 1
 ipc1 = ipc + npc
 IF (ipc1 >= buf5) CALL mesage (-8,0,nam)
 n = n - npc
 
!     LOCATE PCOMP1 DATA AND READ INTO CORE
 
 140 CALL locate (*160,z(buf4),pcomp1,flag)
 
 ipc1 = ipc + npc
 CALL READ (*180,*150,ept,z(ipc1),n,0,npc1)
 CALL mesage (-8,0,nam)
 150 IF (npc1 > 0) typc1 = 1
 ipc2 = ipc1 + npc1
 IF (ipc2 >= buf5) CALL mesage (-8,0,nam)
 n = n - npc1
 
!     LOCATE PCOMP2 DATA AND READ INTO CORE
 
 160 CALL locate (*180,z(buf4),pcomp2,flag)
 
 ipc2 = ipc1 + npc1
 CALL READ (*180,*170,ept,z(ipc2),n,0,npc2)
 CALL mesage (-8,0,nam)
 170 IF (npc2 > 0) typc2 = 1
 
!     SET SIZE OF LPCOMP. NUMBER OF WORDS READ INTO CORE
 
 180 lpcomp = ipc + npc + npc1 + npc2
 IF (lpcomp >= buf5) CALL mesage (-8,0,nam)
 
!     CLOSE EPT BEFORE PROCESSING PCOMPI
 
 CALL CLOSE (ept,1)
 
!     READ MATERIAL PROPERTY TABLE INTO CORE
 
 imat  = lpcomp + 1
 n1mat = buf5 - imat
 CALL premat (z(imat),z(imat),z(buf5),n1mat,n2mat,mpt,dit)
 IF (imat+n2mat >= buf5) CALL mesage (-8,0,nam)
 icore = imat + n2mat + 1
 
!     SET POINTERS
 
 itype  =-1
 istart = 0
 ifinis = 0
 
!     PROCESS ALL 'PCOMP' ENTRY TYPES SEQUENTIALLY
 
!     PCOMP ENTRIES
 
 IF (typc == 0) GO TO 190
 itype  = 0
 istart = ipc
 ifinis = ipc1 - 1
 nwdpc  = 8
 kpc    = 4
 pcompr = 1
 GO TO 220
 
!     PCOMP1 ENTRIES
 
 190 IF (typc1 == 0) GO TO 200
 itype  = 1
 istart = ipc1
 ifinis = ipc2 - 1
 nwdpc  = 8
 kpc    = 1
 pcompr = 1
 GO TO 220
 
!     PCOMP2 ENTRIES
 
 200 IF (typc2 == 0) GO TO 210
 itype  = 2
 istart = ipc2
 ifinis = lpcomp - 1
 nwdpc  = 8
 kpc    = 2
 
!     CHECK IF NO PCOMP DATA HAS BEEN READ INTO CORE
 
 210 IF (typc == 0 .AND. typc1 == 0 .AND. typc2 == 0) GO TO 1210
 
!     SET INFLAG = 12, SO THAT FOR LAMINA REFERENCING MAT1 OR MAT2
!     PROPERTY ENTRY WILL BE RETURNED IN MAT2 FORMAT. EXECPT FOR
!     THOSE REFERENCING MAT8 PROPERTY, IN WHICH CASE THE ENTRY
!     IS MERELY ECHOED.
 
 220 inflag = 12
 
!     SET POINTERS
 
!     WRITE 3-WORD IDENTITY FOR PCOMP DATA
 
!     PCOMP TYPE
 
 IF (itype /= 0) GO TO 230
 CALL WRITE (pcomps,npcmp,3,0)
 GO TO 250
 
!     PCOMP1 TYPE
 
 230 IF (itype /= 1) GO TO 240
 CALL WRITE (pcomps,npcmp1,3,0)
 GO TO 250
 
!     PCOMP2 TYPE
 
 240 CALL WRITE (pcomps,npcmp2,3,0)
 
!     PROCESS ALL 'PCOMP' ENTRIES
 
 250 LEN    = 0
 nlay   = 0
 eoeloc = 0
 pidloc = 1
 tlam   = 0.d0
 rho    = 0.d0
 zk     = 0.0D0
 zk1    = 0.0D0
 
! ... NEXT 5 TERMS ARE NEW IN 2/1990 UAI CODE
!     PICK THEM UP IF OK UAI FLAG IS .TRUE.
 
 IF (.NOT.ok uai) GO TO 255
 tref   = 0.0
 gsube  = 0.0
 alfa1  = 0.0
 alfa2  = 0.0
 alfa12 = 0.0
 
 255 DO  ii = istart,ifinis
   IF (z(ii) == -1) GO TO 270
 END DO
 
 270 eoeloc = ii
 pidloc = istart
 LEN  = eoeloc - pidloc
 nlay = (LEN - nwdpc)/kpc
 lamopt = z(pidloc+7)
 
!     DETERMINE LAMINATE THICKNESS
 
!     PCOMP DATA
 
 IF (itype > 0) GO TO 290
 DO  k = 1,nlay
   iik  = (pidloc+5) + 4*k
   tlam = tlam + rz(iik)
 END DO
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0D0*tlam
 GO TO 320
 
!     PCOMP1 DATA
 
 290 IF (itype > 1) GO TO 300
 iik  = pidloc + 6
 tlam = rz(iik)*nlay
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0D0*tlam
 GO TO 320
 
!     PCOMP2 DATA
 
 300 DO  k = 1,nlay
   iik  = (pidloc+6) + 2*k
   tlam = tlam + rz(iik)
 END DO
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0D0*tlam
 
!     WRITE TO PCOMPS
!      1. PID
!      2. NLAY - NUMBER OF LAYERS
!      3. REMAINDER OF PCOMP ENTRY
 
 320 CALL WRITE (pcomps,z(pidloc),1,0)
 CALL WRITE (pcomps,nlay,1,0)
 
!     SET LEN TO THE NO. WORDS TO BE WRITTEN TO PCOMPS
 
 LEN = LEN - 1
 CALL WRITE (pcomps,z(pidloc+1),LEN,0)
 
!     CALL MAT TO GET LAYER PROPERTIES AND WRITE TO PCOMPS
!     NOTE FOR PCOMP1 AND PCOMP2 ENTRIES THE PROPERTY MATRIX
!     IS ONLY WRITTEN TO PCOMPS ONCE. (ALL LAYER PER ENTRY HAVE
!     THE SAME MID.
!     SIMILARILY FOR PCOMP ENTRY, IF ALL LAYERS REFERENCE THE SAME
!     MID, THEN THE PROPERTY MATRIX IS ONLY WRITTEN ONCE TO PCOMPS.
 
!          ITYPE = 0 PCOMP  ENTRY
!          ITYPE = 1 PCOMP1 ENTRY
!          ITYPE = 2 PCOMP2 ENTRY
 
 mid = 0
 
!     INTILIZISE G1, G2, G3 AND G4 MATRICES
 
 DO  ll = 1,3
   DO  mm = 1,3
     g1 (ll,mm) = 0.0D0
     gd1(ll,mm) = 0.0D0
     g2 (ll,mm) = 0.0D0
     gd2(ll,mm) = 0.0D0
     g4 (ll,mm) = 0.0D0
     gd4(ll,mm) = 0.0D0
   END DO
 END DO
 
 DO  ll = 1,2
   fii(ll)  = 0.0D0
   fi(ll)   = 0.0D0
   ri(ll)   = 0.0D0
   zbar(ll) = 0.0D0
   DO  mm = 1,2
     g3(ll,mm)     = 0.0D0
     gtrflx(ll,mm) = 0.0D0
     trflx(ll,mm)  = 0.0D0
     g3invd(ll,mm) = 0.0D0
   END DO
 END DO
 
!     INTILIZISE ZBAR
 
 zbarx   = 0.0D0
 zbary   = 0.0D0
 zbarxt  = 0.0D0
 zbarxb  = 0.0D0
 zbaryt  = 0.0D0
 zbaryb  = 0.0D0
 zx      = 0.000
 zy      = 0.000
 
 eixx    = 0.000
 eiyy    = 0.000
 
!     LOOP OVER LAYERS
 
 DO  k = 1,nlay
   IF (itype == 0) matid = z(pidloc+4+4*k)
   IF (itype == 1 .OR. itype == 2) matid = z(pidloc+5)
   IF (k >= 2 .AND. (itype == 0 .AND. mid == matid)) GO TO 410
   IF (k >= 2 .AND. (itype == 1 .OR.  itype == 2)  ) GO TO 420
   
   mid = matid
   CALL mat (elid)
   
!     CALL LPROPD TO GET LAYER PROPERTY MATRICES
   
   CALL lpropd (g)
   
!     COPY G(25) TO GLAY(25), FOR WRITING TO PCOMPS
   
   DO  kk = 1,25
     glay(kk) = g(kk)
   END DO
   
! ... NEXT 20 LINES ARE NEW FROM 2/1990 UAI CODE
   
!     COPY ALFA1, ALFA2 AND ALFA12 FROM GLAY(14 THRU 16)
   
   IF (.NOT.ok uai) GO TO 410
   alfa1  = glay(14)
   alfa2  = glay(15)
   alfa12 = glay(16)
   
!     IF PCOMP, COPY TREF AND GE FROM THE MAIN CARD TO THE MATERIAL
!     PROPERTY DATA. THIS IS DONE FOR THE FIRST LAYER
   
   IF (k     > 1) GO TO 410
   IF (itype >= 1) GO TO 405
   tref  = rz(pidloc+5)
   gsube = rz(pidloc+6)
   glay(24) = tref
   glay(25) = gsube
   GO TO 410
   405 tref  = glay(24)
   gsube = glay(25)
   
!     WRITE THE LAYER PROPERTY MATRIX G TO FILE PCOMPS
   
   410 CALL WRITE (pcomps,glay(1),25,0)
   
   
!     CALCULATE CONTRIBUTION OF EACH LAYER TO OVERALL PROPERTY
!     MATRICES G1, G2, G4
   
!     BUILD TRANSFORMATION MATRIX T
   
   420 IF (itype == 0) theta = rz(pidloc+6+4*k)
   IF (itype == 1) theta = rz(pidloc+7+  k)
   IF (itype == 2) theta = rz(pidloc+7+2*k)
   c = DABS(theta)
   IF (c < 0.00002D0) c = 0.0D0
   IF (c > 89.9998D0 .AND. c < 90.0002D0) c =  90.0D0
   IF (c > 179.998D0 .AND. c < 180.002D0) c = 180.0D0
   IF (c > 269.998D0 .AND. c < 270.002D0) c = 270.0D0
   IF (c > 359.998D0 .AND. c < 360.002D0) c = 360.0D0
   IF (theta < 0.0D0) c = -c
   thetar = c*degrad
   
   c  = DCOS(thetar)
   IF (DABS(c) < epsi) c = 0.0D0
   c2 = c*c
   c4 = c2*c2
   s  = DSIN(thetar)
   IF (DABS(s) < epsi) s = 0.0D0
   s2 = s*s
   s4 = s2*s2
   
   t(1) = c2
   t(2) = s2
   t(3) = c*s
   t(4) = s2
   t(5) = c2
   t(6) =-c*s
   t(7) =-2.0*c*s
   t(8) = 2.0*c*s
   t(9) = c2 - s2
   
!                       T
!     CALCULATE GBAR = T  X G X T
   
!     MULTIPLY G X T AND WRITE TO GT
   
   CALL gmmatd (g(1),3,3,0, t(1),3,3,0, GT(1))
   
!               T
!     MULTIPLY T  X GT AND WRITE TO GBR
   
   CALL gmmatd (t(1),3,3,1, GT(1),3,3,0, gbr(1))
   
!     WRITE GBR IN TWO DIMENSIONED ARRAY GBAR
   
   DO  ll = 1,3
     DO  mm = 1,3
       nn = mm + 3*(ll-1)
       gbar(ll,mm) = gbr(nn)
     END DO
   END DO
   
!     PROCESSING FOR G3 MATRIX
   
!                        T
!     CALCULATE GDBAR = T  X GD X T
   
!     DETERMINE GD MATRIX, WHICH IS EQUAL TO G MATRIX WITH POISSONS
!     RATIO=0.0
!        GD(1) ---- YOUNGS MODULUS IN X-DIRN
!        GD(5) ---- YOUNGS MODULUS IN Y-DIRN
!        GD(9) ---- INPLANE SHEAR MODULUS
   
   DO  ll = 1,9
     gd(ll) = 0.0D0
   END DO
   const = 1.0D0 - (g(2)*g(4))/(g(5)*g(1))
   gd(1) = g(1)*const
   gd(5) = g(5)*const
   gd(9) = g(9)
   
!     MULTIPLY GD X T AND WRITE TO GDT
   
   CALL gmmatd (gd(1),3,3,0, t(1),3,3,0, gdt(1))
   
!               T
!     MULTIPLY T  X GDT AND WRITE TO GDBR
   
   CALL gmmatd (t(1),3,3,1, gdt(1),3,3,0, gdbr(1))
   
!     WRITE GDBR IN TWO DIMENSIONED ARRAY GDBAR
   
   DO  ll = 1,3
     DO  mm = 1,3
       nn = mm + 3*(ll-1)
       gdbar(ll,mm) = gdbr(nn)
     END DO
   END DO
   
!     *********************************************************
!     *   NOTE TO APPROXIMATE BEAM BEHAVIOUR THE CROSS AND    *
!     *   COUPLING TERMS IN THE GDBAR MATRIX NEED TO BE       *
!     *   DEGRADED I.E SET TO ZERO.                           *
!     *********************************************************
   
   gdbar(1,2) = 0.0D0
   gdbar(2,1) = 0.0D0
   gdbar(1,3) = 0.0D0
   gdbar(3,1) = 0.0D0
   gdbar(2,3) = 0.0D0
   gdbar(3,2) = 0.0D0
   
!     PERFORM INITIALIZATION
   
   zref = -tlam/2.0D0
   zk1  = zk
   IF (k == 1) zk1 = zref
   IF (itype == 0) zk = zk1 + rz(pidloc+5+4*k)
   IF (itype == 1) zk = zk1 + rz(pidloc+6    )
   IF (itype == 2) zk = zk1 + rz(pidloc+6+2*k)
   zg1 = zk - zk1
   zg4 =-(zk**2 - zk1**2)*0.5D0
   zg2 = (zk**3 - zk1**3)*0.33333333D0
   
!     CALCULATE LAYER CONTRIBUTION TO G1, G2, GD2 ,G4
   
   DO  ir = 1,3
     DO  ic = 1,3
       g1 (ir,ic) =  g1(ir,ic) +  gbar(ir,ic)*zg1
       gd1(ir,ic) = gd1(ir,ic) + gdbar(ir,ic)*zg1
       IF (lamopt == mem .OR. lamopt == symmem) CYCLE
       g2 (ir,ic) =  g2(ir,ic) +  gbar(ir,ic)*zg2
       gd2(ir,ic) = gd2(ir,ic) + gdbar(ir,ic)*zg2
       IF (lamopt == sym) CYCLE
       g4 (ir,ic) =  g4(ir,ic) +  gbar(ir,ic)*zg4
       gd4(ir,ic) = gd4(ir,ic) + gdbar(ir,ic)*zg4
     END DO
   END DO
   
!     CHECK LAMINATION OPTION AND IF SYMM OR SYMM.MEMB CALCULATE
!     LAYER CONTRIBUTION TO THE MEMBRANE, BENDING AND THE
!     MEMEBRANE-BENDING MATRICES
   
   IF (lamopt /= sym .AND. lamopt /= symmem) GO TO 480
   
   DO  ir = 1,3
     DO  ic = 1,3
       g1 (ir,ic) =  g1(ir,ic) +  gbar(ir,ic)*zg1
       gd1(ir,ic) = gd1(ir,ic) + gdbar(ir,ic)*zg1
       IF (lamopt == symmem) CYCLE
       g2 (ir,ic) =  g2(ir,ic) +  gbar(ir,ic)*zg2
       gd2(ir,ic) = gd2(ir,ic) + gdbar(ir,ic)*zg2
     END DO
   END DO
   
   480 CONTINUE
   
!     ************************************************************
!     CALCULATION OF ZBARX AND ZBARY
!            NEUTRAL SURFACE LOCATION IN X- AND Y- DIRECTION
   
!          TI  -  THICKNESS OF LAYER K
!          ZI  -  DISTANCE FROM REFERENCE SURFACE TO MID OF LAMINA K
!       EX,EY  -  APPARENT ENGINEERING PROPERTY. I.E YOUNGS MODULUS
!                 IN THE LONGITUDINAL AND TRANSVERSE DIRECTIONS IN
!                 THE MATERIAL COORDINATE SYSTEM.
!     ************************************************************
   
!     INVERT GDBAR TO DETERMINE EX AND EY
   
   ising = -1
   CALL inverd (3,gdbar,3,dummy,0,determ,ising,indexx)
   
!     THE YOUNGS MODULI EX AND EY IN THE MATERIAL COORD SYSTEM
   
   ex = 1.0D0/gdbar(1,1)
   ey = 1.0D0/gdbar(2,2)
   
   exx = ex
   eyy = ey
   
!     WRITE EXX AND EYY TO PCOMPS
   
   CALL WRITE (pcomps,exx,1,0)
   CALL WRITE (pcomps,eyy,1,0)
   
   IF (lamopt == sym) GO TO 490
   
   ti = zk - zk1
   zi = (zk + zk1)/2.0D0
   
   zbarxt = zbarxt + ex*ti*zi
   zbarxb = zbarxb + ex*ti
   zbaryt = zbaryt + ey*ti*zi
   zbaryb = zbaryb + ey*ti
   
!     CALCULATE CONTRIBUTION TO OVERALL DENSITY RHO
   
   490 IF (g(23) == 0.) CYCLE
   rho = rho + g(23)*zg1
   
!     PROCESS NEXT LAYER
   
 END DO
 
!     JUMP IF LAMOPT IS MEMBRANE OR SYMM.MEMBRANE
 
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 520
 
!     WRITE GD1, GD2 AND GD4 TO STIFF MATRIX AND INVERT
!     TO DETERMINE THE OVERALL BENDING PROPERTY FOR THE
!     LAMINATE.
 
 DO  ll= 1,3
   DO  mm = 1,3
     stiff(ll  ,mm  ) = gd1(ll,mm)
     stiff(ll  ,mm+3) = gd4(ll,mm)
     stiff(ll+3,mm  ) = gd4(ll,mm)
     stiff(ll+3,mm+3) = gd2(ll,mm)
   END DO
 END DO
 
!     INVERT STIFF
 
 ising = -1
 CALL inverd (6,stiff,6,dum,0,determ,ising,INDEX)
 
 ei(1) = 1.0D0/stiff(4,4)
 ei(2) = 1.0D0/stiff(5,5)
 
 eixx = ei(1)
 eiyy = ei(2)
 
!     WRITE EIXX AND EIYY TO PCOMPS
 
 520 CALL WRITE (pcomps,eixx,1,0)
 CALL WRITE (pcomps,eiyy,1,0)
 
!     ***************************************************************
!     *   THE MEMBRANE, BENDING, AND MEMEBRANE-BENDING MATRICES     *
!     *   G1, G2, AND G4 ARE GIVEN BY THE FOLLOWING                 *
!     ***************************************************************
 
 DO  ir = 1,3
   DO  ic = 1,3
     g1(ir,ic) = (1.0D0/tlam)*g1(ir,ic)
     IF (lamopt == mem .OR. lamopt == symmem) CYCLE
     g2(ir,ic) = (12.0D0/tlam**3)*g2(ir,ic)
     IF (lamopt == sym) CYCLE
     g4(ir,ic) = (1.0D0/tlam**2)*g4(ir,ic)
   END DO
 END DO
 
!     CALCULATE LOCATION OF NEUTRAL SURFACE ZBARX AND ZBARY
!     FOR LAMINATE
 
 IF (lamopt == sym .OR. lamopt == mem .OR. lamopt == symmem) GO TO 540
 zbarx   = zbarxt/zbarxb
 zbary   = zbaryt/zbaryb
 zbar(1) = zbarx
 zbar(2) = zbary
 
 zx = zbarx
 zy = zbary
 
!     WRITE ZX AND ZY TO PCOMPS
 
 540 CALL WRITE (pcomps,zx,1,0)
 CALL WRITE (pcomps,zy,1,0)
 
!     CALCULATE OVERALL DENSITY RHO
 
 IF (rho == 0.) GO TO 550
 IF (lamopt == sym .OR. lamopt == symmem) rho = 2.0D0*rho
 rho = rho/tlam
 
!     *****************************************************************
!     *    CHECK IF TRANSVERSE FLEXIBILITY MATRIX NEEDS TO CALCULATED *
!     *    OTHERWISE JUMP TO PROCEED AS PER NORMAL.                   *
!     *****************************************************************
 
 550 IF (lamopt == mem .OR. lamopt == symmem) GO TO 830
 IF (g(10) == 0.0D0) GO TO 830
 
!     LOOP OVER ALL THE LAYERS
 
 DO  k = 1,nlay
   IF (itype == 0) matid = z(pidloc+4+4*k)
   IF (itype == 1 .OR. itype == 2) matid = z(pidloc+5)
   IF (k >= 2 .AND. (itype == 0 .AND. mid == matid)) GO TO 560
   IF (k >= 2 .AND. (itype == 1 .OR.  itype == 2)  ) GO TO 560
   
   mid = matid
   CALL mat (elid)
   
!     CALL LPROPD TO GET LAYER PROPERTY MATRICES
   
   CALL lpropd (g)
   
!     BUILD TRANSFORMATION MATRIX T
   
   560 IF (itype == 0) theta = rz(pidloc+6+4*k)
   IF (itype == 1) theta = rz(pidloc+7+  k)
   IF (itype == 2) theta = rz(pidloc+7+2*k)
   c = DABS(theta)
   IF (c < 0.00002D0) c = 0.0D0
   IF (c > 89.9998D0 .AND. c < 90.0002D0) c =  90.0D0
   IF (c > 179.998D0 .AND. c < 180.002D0) c = 180.0D0
   IF (c > 269.998D0 .AND. c < 270.002D0) c = 270.0D0
   IF (c > 359.998D0 .AND. c < 360.002D0) c = 360.0D0
   IF (theta < 0.0D0) c = -c
   thetar = c*degrad
   
   c  = DCOS(thetar)
   IF (DABS(c) < epsi) c = 0.0D0
   c2 = c*c
   c4 = c2*c2
   s  = DSIN(thetar)
   IF (DABS(s) < epsi) s = 0.0D0
   s2 = s*s
   s4 = s2*s2
   
   t(1) = c2
   t(2) = s2
   t(3) = c*s
   t(4) = s2
   t(5) = c2
   t(6) =-c*s
   t(7) =-2.0*c*s
   t(8) = 2.0*c*s
   t(9) = c2 - s2
   
!     PROCESSING FOR G3 MATRIX
   
!                       T
!     CALCULATE GDBR = T  X GD X T
   
!     DETERMINE GD MATRIX, WHICH IS EQUAL TO G MATRIX WITH POISSONS
!     RATIO=0.0
!        GD(1) ---- YOUNGS MODULUS IN X-DIRN
!        GD(5) ---- YOUNGS MODULUS IN Y-DIRN
!        GD(9) ---- INPLANE SHEAR MODULUS
   
   DO  ll = 1,9
     gd(ll) = 0.0D0
   END DO
   const = 1.0D0 - (g(2)*g(4))/(g(5)*g(1))
   gd(1) = g(1)*const
   gd(5) = g(5)*const
   gd(9) = g(9)
   
!     MULTIPLY GD X T AND WRITE TO GDT
   
   CALL gmmatd (gd(1),3,3,0, t(1),3,3,0, gdt(1))
   
!               T
!     MULTIPLY T  X GDT AND WRITE TO GDBR
   
   CALL gmmatd (t(1),3,3,1, gdt(1),3,3,0, gdbr(1))
   
!     WRITE GBR TO GDBAR
   
   DO  ll = 1,3
     DO  mm = 1,3
       nn = mm + 3*(ll-1)
       gdbar(ll,mm) = gdbr(nn)
     END DO
   END DO
   
!     *************************************************************
!     *       NOTE TO APPROXIMATE BEAM BEHAVIOUR THE CROSS AND    *
!     *       COUPLING TERMS IN THE GDBAR MATRIX NEED TO BE       *
!     *       DEGRADED I.E SET TO ZERO.                           *
!     *************************************************************
   
   gdbar(1,2) = 0.0D0
   gdbar(2,1) = 0.0D0
   gdbar(1,3) = 0.0D0
   gdbar(3,1) = 0.0D0
   gdbar(2,3) = 0.0D0
   gdbar(3,2) = 0.0D0
   
!     INVERT GDBAR TO DETERMINE EX AND EY
   
   ising = -1
   CALL inverd (3,gdbar,3,dummy,0,determ,ising,indexx)
   
!     THE YOUNGS MODULI EX AND EY IN THE MATERIAL COORD SYSTEM ARE
   
   e(1) = 1.0D0/gdbar(1,1)
   e(2) = 1.0D0/gdbar(2,2)
   
!     PERFORM INTILIZATION
   
   zref = -tlam/2.0D0
   zk1 = zk
   IF (k == 1) zk1 = zref
   IF (itype == 0) zk = zk1 + rz(pidloc+5+4*k)
   IF (itype == 1) zk = zk1 + rz(pidloc+6    )
   IF (itype == 2) zk = zk1 + rz(pidloc+6+2*k)
   
!     BUILD TRANSFORMATION MATRIX U
   
   u(1) = c
   u(2) = s
   u(3) =-s
   u(4) = c
   
!     CALCULATE G3BAR = UT X G3I X U
!     G3I MATRIX  -  LAYER K TRANSFORMED G3, IN MATERIAL COORD-SYS
   
   DO  ll = 1,4
     mm = ll + 9
     g3i(ll) = g(mm)
   END DO
   
!     MULTIPLY G3I X U AND WRITE TO G3IU
   
   CALL gmmatd (g3i(1),2,2,0, u(1),2,2,0, g3iu(1))
   
!     MULTIPLY UT X G3IU AND WRITE TO G3BR
   
   CALL gmmatd (u(1),2,2,1, g3iu(1),2,2,0, g3br(1))
   
!     WRITE G3BR IN TWO DIMENSIONED ARRAY G3BAR
   
   DO  ll = 1,2
     DO  mm = 1,2
       nn = mm + 2*(ll-1)
       g3bar(ll,mm) = g3br(nn)
     END DO
   END DO
   
!     INVERT G3BAR
   
   detrmn = g3bar(1,1)*g3bar(2,2) - g3bar(1,2)*g3bar(2,1)
   IF (detrmn == 0.0D0) GO TO 1230
   
   g3invd(1,1) = g3bar(2,2)/detrmn
   g3invd(1,2) =-g3bar(1,2)/detrmn
   g3invd(2,1) =-g3bar(2,1)/detrmn
   g3invd(2,2) = g3bar(1,1)/detrmn
   
!     G3 MATRIX CALC
   
   zi = (zk + zk1)/2.0D0
   ti =  zk - zk1
   
   DO  ir = 1,2
     ri(ir) = ((fi(ir)/e(ir)) + (zbar(ir)-zk1)*ti - (ti*ti/3.0D0))  &
         * (fi(ir)/e(ir))
     ri(ir) = ri(ir) + zbar(ir)*ti*ti*((zbar(ir)-2.0D0*zk1)/3.0D0  &
         - (ti/4.0D0))
     ri(ir) = ri(ir) + ti*ti*((zk1*zk1)/3.0D0 + (zk1*ti)/4.0D0  &
         + (ti*ti)/20.0D0)
     ri(ir) = ri(ir)*e(ir)*e(ir)*ti
   END DO
   
   DO  ir = 1,2
     DO  ic = 1,2
       gtrflx(ir,ic) = gtrflx(ir,ic) + ri(ir)*g3invd(ir,ic)
     END DO
   END DO
   
   DO  ir = 1,2
     fii(ir) = e(ir)*ti*(zbar(ir)-zi)
     fi(ir)  = fi(ir) + fii(ir)
   END DO
   
!     PROCESS NEXT LAYER
   
 END DO
 
!     FALL HERE IF LAMOPT IS SYMM AND G3 CALCULATION IS REQUIRED
 
 IF (lamopt /= sym) GO TO  810
 DO  kk = 1,nlay
   k = nlay + 1 - kk
   
   IF (itype == 0) matid = z(pidloc+4+4*k)
   IF (itype == 1 .OR. itype == 2) matid = z(pidloc+5)
   IF (k >= 2 .AND. (itype == 0 .AND. mid == matid)) GO TO 710
   IF (k >= 2 .AND. (itype == 1 .OR.  itype == 2)  ) GO TO 710
   
   mid = matid
   CALL mat (elid)
   
!     CALL LPROPD TO GET LAYER PROPERTY MATRICES
   
   CALL lpropd (g)
   
!     BUILD TRANSFORMATION MATRIX T
   
   710 IF (itype == 0) theta = rz(pidloc+6+4*k)
   IF (itype == 1) theta = rz(pidloc+7+  k)
   IF (itype == 2) theta = rz(pidloc+7+2*k)
   c = DABS(theta)
   IF (c < 0.00002D0) c = 0.0D0
   IF (c > 89.9998D0 .AND. c < 90.0002D0) c =  90.0D0
   IF (c > 179.998D0 .AND. c < 180.002D0) c = 180.0D0
   IF (c > 269.998D0 .AND. c < 270.002D0) c = 270.0D0
   IF (c > 359.998D0 .AND. c < 360.002D0) c = 360.0D0
   IF (theta < 0.0D0) c = -c
   thetar = c*degrad
   
   c  = DCOS(thetar)
   IF (DABS(c) < epsi) c = 0.0D0
   c2 = c*c
   c4 = c2*c2
   s  = DSIN(thetar)
   IF (DABS(s) < epsi) s = 0.0D0
   s2 = s*s
   s4 = s2*s2
   
   t(1) = c2
   t(2) = s2
   t(3) = c*s
   t(4) = s2
   t(5) = c2
   t(6) =-c*s
   t(7) =-2.0*c*s
   t(8) = 2.0*c*s
   t(9) = c2 - s2
   
!     PROCESSING FOR G3 MATRIX
   
!                       T
!     CALCULATE GDBR = T  X GD X T
   
!     DETERMINE GD MATRIX, WHICH IS EQUAL TO G MATRIX WITH POISSONS
!     RATIO=0.0
!        GD(1) ---- YOUNGS MODULUS IN X-DIRN
!        GD(5) ---- YOUNGS MODULUS IN Y-DIRN
!        GD(9) ---- INPLANE SHEAR MODULUS
   
   DO  ll = 1,9
     gd(ll) = 0.0D0
   END DO
   const = 1.0D0 - (g(2)*g(4))/(g(5)*g(1))
   gd(1) = g(1)*const
   gd(5) = g(5)*const
   gd(9) = g(9)
   
!     MULTIPLY GD X T AND WRITE TO GDT
   
   CALL gmmatd (gd(1),3,3,0, t(1),3,3,0, gdt(1))
   
!               T
!     MULTIPLY T  X GDT AND WRITE TO GDBR
   
   CALL gmmatd (t(1),3,3,1, gdt(1),3,3,0, gdbr(1))
   
!     WRITE GBR TO GDBAR
   
   DO  ll = 1,3
     DO  mm = 1,3
       nn = mm + 3*(ll-1)
       gdbar(ll,mm) = gdbr(nn)
     END DO
   END DO
   
!     *************************************************************
!     *       NOTE TO APPROXIMATE BEAM BEHAVIOUR THE CROSS AND    *
!     *       COUPLING TERMS IN THE GDBAR MATRIX NEED TO BE       *
!     *       DEGRADED I.E SET TO ZERO.                           *
!     *************************************************************
   
   gdbar(1,2) = 0.0D0
   gdbar(2,1) = 0.0D0
   gdbar(1,3) = 0.0D0
   gdbar(3,1) = 0.0D0
   gdbar(2,3) = 0.0D0
   gdbar(3,2) = 0.0D0
   
!     INVERT GDBAR TO DETERMINE EX AND EY
   
   ising = -1
   CALL inverd (3,gdbar,3,dummy,0,determ,ising,indexx)
   
!     THE YOUNGS MODULI EX AND EY IN THE MATERIAL COORD SYSTEM ARE
   
   e(1) = 1.0D0/gdbar(1,1)
   e(2) = 1.0D0/gdbar(2,2)
   
!     PERFORM INTILIZATION
   
   zref = -tlam/2.0D0
   zk1  = zk
   IF (itype == 0) zk = zk1 + rz(pidloc+5+4*k)
   IF (itype == 1) zk = zk1 + rz(pidloc+6    )
   IF (itype == 2) zk = zk1 + rz(pidloc+6+2*k)
   
!     BUILD TRANSFORMATION MATRIX U
   
   u(1) = c
   u(2) = s
   u(3) =-s
   u(4) = c
   
!     CALCULATE G3BAR = UT X G3I X U
!     G3I MATRIX  -  LAYER K TRANSFORMED G3, IN MATERIAL COORD-SYS
   
   DO  ll = 1,4
     mm = ll + 9
     g3i(ll) = g(mm)
   END DO
   
!     MULTIPLY G3I X U AND WRITE TO G3IU
   
   CALL gmmatd (g3i(1),2,2,0, u(1),2,2,0, g3iu(1))
   
!     MULTIPLY UT X G3IU AND WRITE TO G3BR
   
   CALL gmmatd (u(1),2,2,1, g3iu(1),2,2,0, g3br(1))
   
!     WRITE G3BR IN TWO DIMENSIONED ARRAY G3BAR
   
   DO  ll = 1,2
     DO  mm = 1,2
       nn = mm + 2*(ll-1)
       g3bar(ll,mm) = g3br(nn)
     END DO
   END DO
   
!     INVERT G3BAR
   
   detrmn = g3bar(1,1)*g3bar(2,2) - g3bar(1,2)*g3bar(2,1)
   IF (detrmn == 0.0D0) GO TO 1230
   
   g3invd(1,1) = g3bar(2,2)/detrmn
   g3invd(1,2) =-g3bar(1,2)/detrmn
   g3invd(2,1) =-g3bar(2,1)/detrmn
   g3invd(2,2) = g3bar(1,1)/detrmn
   
!     THE CORRESSPONDING LAYER ON THE OTHER SIDE OF SYMMETRY
   
   zi = (zk + zk1)/2.0D0
   ti =  zk - zk1
   
   DO  ir = 1,2
     ri(ir) = (fi(ir)/e(ir) +(-zk1)*ti-ti*ti/3.0D0 )*fi(ir)/e(ir)  &
         + (zk1*zk1/3.0D0+zk1*ti/4.0D0+ti*ti/20.0D0)*ti*ti
     ri(ir) = ri(ir)*e(ir)*e(ir)*ti
   END DO
   
   DO  ir = 1,2
     DO  ic = 1,2
       gtrflx(ir,ic) = gtrflx(ir,ic) + ri(ir)*g3invd(ir,ic)
     END DO
   END DO
   
   DO  ir = 1,2
     fii(ir) = e(ir)*ti*(zbar(ir)-zi)
     fi(ir)  = fi(ir) + fii(ir)
   END DO
   
!     PROCESS NEXT LAYER
   
 END DO
 
 810 DO  ir = 1,2
   DO  ic = 1,2
     gtrflx(ir,ic) = gtrflx(ir,ic)*tlam/(ei(ir)**2)
   END DO
 END DO
 
!     INVERT GTRFLX
 
 detrmn = gtrflx(1,1)*gtrflx(2,2) - gtrflx(1,2)*gtrflx(2,1)
 IF (detrmn == 0.0D0) GO TO 1230
 
 g3(1,1) = gtrflx(2,2)/detrmn
 g3(1,2) =-gtrflx(1,2)/detrmn
 g3(2,1) =-gtrflx(2,1)/detrmn
 g3(2,2) = gtrflx(1,1)/detrmn
 
!     BECAUSE G3(1,2) IS NOT EQUAL TO G3(2,1) IN GENERAL
!     AN AVERAGE VALUE WILL BE USED FOR THE COUPLING TERMS
 
 g3(1,2) = (g3(1,2) + g3(2,1))/2.0D0
 g3(2,1) = g3(1,2)
 
!    *****************************************************
!    WRITE THE NEWLY GENERATED G1, G2, G3, AND G4 MATRICES
!    TO MPTX IN THE FORM OF MAT2 DATA ENTRIES
!    *****************************************************
 
!      NOTE - THE MID FOR THESE MATRICES ARE AS FOLLOWS-
!         1. MID1  -- PID + 100000000
!         2. MID2  -- PID + 200000000
!         3. MID3  -- PID + 300000000
!         4. MID4  -- PID + 400000000
 
!     INITIALIZE G1, G2, G3, AND G4 MATRICES
 
 830 DO  jj = 1,17
   gmembr(jj) = 0.0D0
   gbendg(jj) = 0.0D0
   gtrshr(jj) = 0.0D0
   gmembd(jj) = 0.0D0
 END DO
 
 imembr(1) = 0
 ibendg(1) = 0
 itrshr(1) = 0
 imembd(1) = 0
 
!     START GENERATING G1 MEMBRANE MATRIX
 
 imembr( 1) = z(pidloc) + 100000000
 gmembr( 2) = g1(1,1)
 gmembr( 3) = g1(1,2)
 gmembr( 4) = g1(1,3)
 gmembr( 5) = g1(2,2)
 gmembr( 6) = g1(2,3)
 gmembr( 7) = g1(3,3)
 gmembr( 8) = rho
 
! ... NEXT 5 TERMS ARE NEW FROM 2/1990 UAI CODE
 
 IF (.NOT.ok uai) GO TO 845
 gmembr( 9) = alfa1
 gmembr(10) = alfa2
 gmembr(11) = alfa12
 gmembr(12) = tref
 gmembr(13) = gsube
 
 845 IF (lamopt == mem .OR. lamopt == symmem) GO TO 850
 
!     START GENERATING G2 BENDING MATRIX
 
 ibendg( 1) = z(pidloc) + 200000000
 gbendg( 2) = g2(1,1)
 gbendg( 3) = g2(1,2)
 gbendg( 4) = g2(1,3)
 gbendg( 5) = g2(2,2)
 gbendg( 6) = g2(2,3)
 gbendg( 7) = g2(3,3)
 
! ... NEXT 3 TERMS ARE NEW FROM 2/1990 UAI CODE
 
 IF (.NOT.ok uai) GO TO 847
!     GBEMDG( 8) = ??
 gbendg( 9) = alfa1
 gbendg(10) = alfa2
 gbendg(11) = alfa12
 
!     START GENERATING G3 TRANSVERSE SHEAR FLEXIBILITY MATRIX
 
 847 itrshr( 1) = z(pidloc) + 300000000
 gtrshr( 2) = g3(1,1)
 gtrshr( 3) = g3(1,2)
 gtrshr( 4) = g3(2,1)
 gtrshr( 5) = g3(2,2)
 
 IF (lamopt == sym) GO TO 850
 
!     START GENERATING G4 MEMBRANE-BENDING COUPLING MATRIX
 
 imembd( 1) = z(pidloc) + 400000000
 gmembd( 2) = g4(1,1)
 gmembd( 3) = g4(1,2)
 gmembd( 4) = g4(1,3)
 gmembd( 5) = g4(2,2)
 gmembd( 6) = g4(2,3)
 gmembd( 7) = g4(3,3)
 
 850 CONTINUE
 
!     ******************************************************
!     GENERATE EQUIVALENT PSHELL BULK DATA ENTIES FOR EVERY
!     PCOMPI BULK DATA ENTRY. THIS IS NECESSARY FOR DEMG TO
!     FUNCTION CORRECTLY WHEN LAMINATED COMPOSITE ELEMENTS
!     ARE PRESENT.
!     ******************************************************
 
 ipshel( 1) = z(pidloc)
 ipshel( 2) = z(pidloc) + 100000000
 rpshel( 3) = tlam
 ipshel( 4) = z(pidloc) + 200000000
 rpshel( 5) = 1.0
 ipshel( 6) = z(pidloc) + 300000000
 rpshel( 7) = 1.0
 rpshel( 8) = rz(pidloc+2)
 rpshel( 9) =-tlam/2.0
 rpshel(10) = tlam/2.0
 ipshel(11) = z(pidloc) + 400000000
 rpshel(12) = 0.0
 ipshel(13) = 0
 ipshel(14) = 0
 rpshel(15) = 0.0
 ipshel(16) = 0
 rpshel(17) = 0.0
 
 zoffs = rz(pidloc+1) + tlam/2.0
 IF (z(pidloc)  ==  BLANK) zoffs = 0.0
 IF (lamopt == mem .OR. lamopt == symmem) zoffs = 0.0
 IF (ABS(zoffs) <= 1.0E-3) zoffs = 0.0
 rpshel(14) = zoffs
 
 IF (lamopt /= mem .AND. lamopt /= symmem) GO TO 860
 ipshel( 4) = 0
 ipshel( 6) = 0
 ipshel(11) = 0
 rpshel(14) = 0.0
 860 IF (lamopt /= sym) GO TO 870
 ipshel(11) = 0
 870 CONTINUE
 
!     UPDATE COUNTER ICOUNT TO INDICATE MAT2 AND PSHELL DATA IS BEING
!     WRITTEN SECOND TIME
 
 icount = icount + 1
 
 IF (icount > 1) GO TO 900
 
 IF (pshlpr /= 1) GO TO 890
 icore = lpcomp + 1 + n2mat
 n = buf5 - icore
 CALL OPEN (*1200,ept,z(buf4),rdrew)
 CALL filpos (ept,pos1)
 CALL READ (*900,*880,ept,z(icore),n,0,eptwds)
 CALL mesage (-8,0,nam)
 880 CALL WRITE (eptx,z(icore),eptwds,0)
 GO TO 900
 890 CALL WRITE (eptx,pshnam,3,0)
 900 CALL WRITE (eptx,ipshel(1),17,0)
 
 IF (icount > 1) GO TO 930
 
 IF (mat2pr /= 1) GO TO 920
 icore = lpcomp + 1 + n2mat
 n = buf5 - icore
 CALL OPEN (*1200,mpt,z(buf2),rdrew)
 CALL filpos (mpt,pos)
 CALL READ (*930,*910,mpt,z(icore),n,0,matwds)
 CALL mesage (-8,0,nam)
 910 CALL WRITE (mptx,z(icore),matwds,0)
 GO TO 930
 920 CALL WRITE (mptx,matnam,3,0)
 930 CALL WRITE (mptx,imembr(1),17,0)
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 940
 CALL WRITE (mptx,ibendg(1),17,0)
 CALL WRITE (mptx,itrshr(1),17,0)
 IF (lamopt == sym) GO TO 940
 CALL WRITE (mptx,imembd(1),17,0)
 940 CONTINUE
 CALL sswtch (40,l40)
 IF (l40 == 0) GO TO 980
 
!     WRITE THE NEWLY GENERATED PROPERTY MATRICES TO THE OUTPUT FILE
 
 CALL page2 (2)
 WRITE (nout,960) imembr(1),(gmembr(ll),ll=2,16)
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 980
 CALL page2 (2)
 WRITE (nout,960) ibendg(1),(gbendg(ll),ll=2,16)
 IF (gtrshr(1) == 0.0) GO TO 950
 CALL page2 (2)
 WRITE (nout,960) itrshr(1),(gtrshr(ll),ll=2,16)
 950 IF (lamopt == sym) GO TO 980
 CALL page2 (2)
 WRITE (nout,960) imembd(1),(gmembd(ll),ll=2,16)
 960 FORMAT(/,' MAT2',7X,i9,7(1X,1P,e11.4),/,9X,8(1X,f11.1))
 
!     UPDATE LOCATION OF NEXT PID
 
 980 pidloc = eoeloc + 1
 istart = pidloc
 
!     WRITE END OF ENTRY (EOE) TO PCOMPS BEFORE PROCESSING
!     NEXT PCOMP ENTRY
 
 CALL WRITE (pcomps,eoe,1,0)
 
!     CHECK IF ALL 'PCOMP' TYPE ENTRIES HAVE BEEN PROCESSED
 
 IF (istart >= ifinis) IF (itype-1) 990,1000,1010
 
!     PROCESS NEXT 'PCOMP' ENTRY
 
 GO TO 250
 
 990 CALL WRITE (pcomps,0,0,1)
 IF (typc1 > 0) GO TO 190
 IF (typc2 > 0) GO TO 200
 GO TO 1020
 
 1000 CALL WRITE (pcomps,0,0,1)
 IF (typc2 > 0) GO TO 200
 GO TO 1020
 
 1010 CALL WRITE (pcomps,0,0,1)
 
!     ALL 'PCOMP' TYPES PROCESSED
!     WRITE EOR ON MPTX AND EPTX
 
 1020 CALL WRITE (mptx,0,0,1)
 CALL WRITE (eptx,0,0,1)
 
!     COPY REMAINDER OF EPT TO EPTX
 
 icore = 1
 n = buf5 - 1
 eptwds = 0
 IF (pshlpr /= 1) CALL OPEN (*1200,ept,z(buf4),rdrew)
 CALL filpos (ept,pos1)
 IF (pshlpr == 1) CALL fwdrec (*1050,ept)
 1030 CALL READ (*1050,*1040,ept,z(icore),n,1,eptwds)
 CALL mesage (-8,0,nam)
 1040 CALL WRITE (eptx,z(icore),eptwds,1)
 eptwds = 0
 GO TO 1030
 
!     READ TRAILER FROM EPT AND WRITE TO EPTX
 
 1050 DO  kk = 1,7
   ieptx(kk) = 0
 END DO
 ieptx( 1) = ept
 
 CALL rdtrl (ieptx)
 ieptx(1) = eptx
 kt721 = andf(pshbit,511)
 k1 = (kt721-1)/16 + 2
 k2 = kt721 - (k1-2)*16 + 16
 ieptx(k1) = orf(ieptx(k1),two(k2))
 CALL wrttrl (ieptx)
 
!     IF EOF ON MPT,THEN ALL MAT2 DATA COPIED TO MPTX
 
 IF (eof == 1) GO TO 1090
 
!     OTHERWISE COPY REMAINDER OF MPT TO MPTX
 
 icore = 1
 n = buf5 - 1
 matwds = 0
 IF (mat2pr /= 1) CALL OPEN (*1200,mpt,z(buf2),rdrew)
 CALL filpos (mpt,pos)
 IF (mat2pr == 1) CALL fwdrec (*1090,mpt)
 1070 CALL READ (*1090,*1080,mpt,z(icore),n,1,matwds)
 CALL mesage (-8,0,nam)
 1080 CALL WRITE (mptx,z(icore),matwds,1)
 matwds = 0
 GO TO 1070
 
!     READ TRAILER FROM MPT AND WRITE TO MPTX
 
 1090 DO  kk = 1,7
   imptx(kk) = 0
 END DO
 imptx( 1) = mpt
 
 CALL rdtrl (imptx)
 imptx(1) = mptx
 kt721 = andf(mt2bit,511)
 k1 = (kt721-1)/16 + 2
 k2 = kt721 - (k1-2)*16 + 16
 imptx(k1) = orf(imptx(k1),two(k2))
 CALL wrttrl (imptx)
 
!     WRITE TO TRAILER OF PCOMPS
 
!     SET TRAILER BIT POSITION TO ZERO IF ENTRY TYPE DOES NOT EXIST
 
 IF (typc  == 0) pcbit(1) = 0
 IF (typc1 == 0) pcbit(2) = 0
 IF (typc2 == 0) pcbit(3) = 0
 
 DO  ll = 1,3
   kt721 = andf(pcbit(ll),511)
   k1 = (kt721-1)/16 + 2
   k2 = kt721 - (k1-2)*16 + 16
   ipcomp(k1) = orf(ipcomp(k1),two(k2))
 END DO
 
!     WHEN ICFIAT IS 11, A 65536 IS LEFT IN IPCOMP(2) ACCIDENTALLY
!     ZERO IT OUT
 
 IF (icfiat == 11) ipcomp(2) = 0
 CALL wrttrl (ipcomp)
 
!     CLOSE ALL FILES
 
 CALL CLOSE (pcomps,1)
 CALL CLOSE (eptx,1)
 CALL CLOSE (mptx,1)
 CALL CLOSE (mpt,1)
 CALL CLOSE (ept,1)
 
 RETURN
 
!     FATAL ERROR MESSAGES
 
 1200 CALL mesage (-1,FILE,nam)
 GO TO 1300
 1210 CALL page2 (2)
 WRITE  (nout,1220)
 1220 FORMAT ('0*** SYSTEM FATAL ERROR.  PCOMP, PCOMP1 OR PCOMP2',  &
     ' DATA NOT FOUND BY SUBROUTINE TA1CPD.')
 nogo = 1
 GO TO 1300
 1230 CALL page2 (4)
 WRITE (nout,1240) matid
 nogo = 1
 1240 FORMAT ('0*** USER FATAL ERROR.  IMPROPER DATA PROVIDED FOR',  &
     ' CALCULATION OF TRANSVERSE SHEAR FLEXIBILITY MATRIX',  &
     /23X,'FOR LAMINA REFERENCING MID ',i8,'.',  &
     /23X,'CHECK DATA ON MAT BULK DATA ENTRY.')
 1300 CONTINUE
 RETURN
END SUBROUTINE ta1cpd
