SUBROUTINE extern (nex,ngrav,gvect,ilist,pg,n1,iharm)
     
!     GENERATES EXTERNAL LOADS
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN OUT)                  :: nex
 INTEGER, INTENT(OUT)                     :: ngrav
 REAL, INTENT(IN OUT)                     :: gvect(1)
 INTEGER, INTENT(IN)                      :: ilist(1)
 INTEGER, INTENT(IN OUT)                  :: pg(1)
 INTEGER, INTENT(IN)                      :: n1
 INTEGER, INTENT(IN OUT)                  :: iharm
 INTEGER :: NAME(2),iz(1),ret
 REAL :: core
 COMMON /tranx / idum(14)
 COMMON /BLANK / nrowsp
 COMMON /zzzzzz/ core(1)
 COMMON /loadx / lcare,slt,bgpdt,OLD,cstm,sil,isil,est,mpt,nn(7),  &
     nobld,idit,icm,ilid
 COMMON /system/ sysbuf
 COMMON /packx / itya,ityb,ii,jj,incur
 COMMON /hmatdd/ iihmat,nnhmat,mptfil,iditfl
 COMMON /pindex/ iest(45)
 COMMON /gpta1 / jdum
 EQUIVALENCE     (core(1),iz(1))
 DATA    casecc, permbd,hcflds,remfls,scr6,hccens, NAME         /  &
     110   , 112   ,304   ,305   ,306 ,307   , 4HEXTE,4HRN  /
 
 iest(1) =-1
 idum(1) = 0
 jopen = 0
 ipre  = 0
 incur = 1
 ii    = 1
 jj    = nrowsp
 ngrav = 0
 OLD   = 0
 icm   = 1
 itya  = 1
 ityb  = 1
 ibuf1 = lcare - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 ibuf4 = ibuf3 - sysbuf
 ibuf5 = ibuf4 - sysbuf
 lcore = ibuf5 - sysbuf
 CALL gopen (slt,core(ibuf1),0)
 CALL gopen (bgpdt,core(ibuf2),0)
 FILE = cstm
 CALL OPEN (*20,cstm,core(ibuf3),0)
 icm  = 0
 CALL skprec (cstm,1)
 20 CALL gopen (sil,core(ibuf4),0)
 FILE = slt
 isil = 0
 IF (lcore < nrowsp) GO TO 1580
 
 iii  = 1
 DO  nloop = 1,n1
   
   ilid = ilist(iii)
   IF (ilid /= 0) GO TO 30
   CALL skprec (slt,1)
   GO TO 1310
   30 DO  i = 1,nrowsp
     core(i) = 0.0
   END DO
   nograv  = 0
   ngrold  = ngrav
   50 CALL READ (*1520,*1300,slt,nobld,1,0,flag)
   CALL fread (slt,ido,1,0)
   IF (nograv ==  1) GO TO 1570
   IF (nobld == -20) GO TO 800
   GO TO (100,100,120,120,140,140,160,200,220,300,  &
       320,340,600,620,630,640,360,700,730,800, 800,800,800,800,400), nobld
   100 DO  j = 1,ido
     CALL DIRECT
   END DO
   GO TO 50
   120 DO  j = 1,ido
     CALL tpont
   END DO
   GO TO 50
   140 DO  j = 1,ido
     CALL fpont
   END DO
   GO TO 50
   160 DO  j = 1,ido
     CALL sload
   END DO
   GO TO 50
   200 IF (nograv == 2) GO TO 1570
   DO  j = 1,ido
     CALL grav (ngrav,gvect(1),nex,ilist(1),nloop)
   END DO
   nograv = 1
   GO TO 50
   220 DO  j = 1,ido
     CALL pload
   END DO
   GO TO 50
   
!     RFORCE CARDS
   
   300 DO  j = 1,ido
     CALL rforce (lcore)
   END DO
   GO TO 50
   
!     PRESAX CARDS
   
   320 DO  j = 1,ido
     CALL presax (iharm)
   END DO
   GO TO 50
   
!     QHBDY CARDS
   
   340 DO  j = 1,ido
     CALL qhbdy
   END DO
   GO TO 50
   
!     PLOAD3 CARDS
   
   360 DO  j = 1,ido
     CALL pload3
   END DO
   GO TO 50
   
!     PLOAD4 CARDS
   
   400 CALL pload4 (ibuf5,ido,jopen)
   GO TO 50
   
!     QVOL CARDS (MODIFIED USER ENTRYS)
   
   600 DO  j = 1,ido
     CALL qvol
   END DO
   GO TO 50
   
!     QBDY1 CARDS (MODIFIED USER ENTRYS)
   
   620 kkkk = 1
   GO TO 650
   
!     QBDY2 CARDS (MODIFIED USER ENTRYS)
   
   630 kkkk = 2
   GO TO 650
   
!     QVECT CARDS (MODIFIED USER ENTRYS)
   
   640 kkkk = 3
   650 DO  j = 1,ido
     CALL qloadl (kkkk)
   END DO
   GO TO 50
   
!     PLOAD1 CARDS
   
   700 IF (ipre == 1) GO TO 710
   ipre  = 1
   lcore = lcore - sysbuf - 1
   mcore = lcore - nrowsp - 1
   IF (lcore < nrowsp) GO TO 1580
   CALL premat (core(nrowsp+1),core(nrowsp+1),core(lcore),mcore,  &
       ncore,mpt,idit)
   710 DO  j = 1,ido
     CALL plbar1 (ido,lcore)
   END DO
   GO TO 50
   
!     PLOADX CARDS
   
   730 DO  j = 1,ido
     CALL ploadx
   END DO
   GO TO 50
   
!     CEMLOOP, SPCFLD, GEMLOOP, MDIPOLE, AND REMFLUX CARDS
   
!     BRING HEAT MATERIALS INTO CORE
   
   800 IF (ipre == 1) GO TO 1230
   ipre = 1
   
!     1ST AND LAST AVAILABLE LOCATIONS IN OPEN CORE
   
   iihmat = nrowsp
   nnhmat = lcore
   mptfil = mpt
   iditfl = idit
   CALL prehma (core)
   
!     NOW NNHMAT CONTAINS LAST LOCATION OF MATERIAL INFO
   
   nextz = nnhmat + 1
   
!     OPEN HCFLDS TO CONTAIN APPLIED MAGNETIC FIELD LOAD
   
   lcore = lcore - sysbuf
   IF (lcore <= nextz) GO TO 1580
   
!     STORE SILS  ON PERMBDY, IF ANY, INTO OPEN CORE
   
   nbdys = 0
   FILE  = permbd
   CALL OPEN (*820,permbd,core(lcore+1),0)
   CALL fwdrec (*1520,permbd)
   CALL READ (*1520,*810,permbd,core(nextz),lcore-nextz+1,0,nbdys)
   GO TO 1580
   810 CALL CLOSE (permbd,1)
   820 CONTINUE
   nextz = nextz + nbdys
   
!     NOW CHECK FOR FORCE REQUESTS ON CASECC(MAGNETIC FIELD REQUESTS)
!     MAKE A UNIQUE LIST OF ELEMENT ID-S CORRESPONDING TO ALL SUBCASES.
!     IF A SUBCASE REQUESTS ALL, NO LIST IS NECESSARY.
   
   all    = 0
   nelout = 0
   ij     = 0
   
!     1ST GET MAXIMUM LENGTH OF CASE CONTROL IN ORDER TO STORE ELEMENT
!     ID-S
   
   ncc = 0
   CALL gopen (casecc,core(lcore+1),0)
   830 CALL READ (*850,*840,casecc,core(nextz),lcore-nextz+1,0,kcc)
   GO TO 1580
   840 ncc = MAX0(ncc,kcc)
   GO TO 830
   850 CALL REWIND (casecc)
   CALL fwdrec (*1520,casecc)
   kset = nextz + ncc
   
   860 CALL READ (*1200,*870,casecc,core(nextz),lcore-nextz+1,0,ncc)
   GO TO 1580
   870 setno = iz(nextz+25)
   IF (setno == 0) GO TO 860
   IF (setno > 0) GO TO 1010
   
!     ALL
   
   1000 all    = 1
   nelout = 0
   GO TO 1200
   
!     CREATE UNIQUE LIST OF ELEMENT ID-S
   
   1010 ilsym  = iz(nextz+165)
   isetno = ilsym  + iz(ilsym+nextz-1) + nextz
   1020 iset   = isetno + 2
   nset   = iz(isetno+1) + iset - 1
   IF (iz(isetno) == setno) GO TO 1030
   isetno = nset + 1
   
!     IF SET CANNOT BE FOUND, SET TO ALL. BUT SHOULD NOT HAPPEN
   
   IF (isetno < ncc+nextz-1) GO TO 1020
   GO TO 1000
   
!     PICK UP ELEMENT ID-S. STORE IN UNIQUE LIST
   
   1030 i = iset
   1040 IF (i    == nset) GO TO 1060
   IF (iz(i+1) > 0) GO TO 1060
   ib = iz(i  )
   n  =-iz(i+1)
   i  = i + 1
   ASSIGN 1050 TO ret
   GO TO 1100
   1050 ib = ib + 1
   IF (ib <= n) GO TO 1100
   GO TO 1070
   1060 ib = iz(i)
   ASSIGN 1070 TO ret
   GO TO 1100
   1070 i = i + 1
   IF (i <= nset) GO TO 1040
   
!     DONE WITH THIS SET. GO BACK FOR ANOTHER
   
   GO TO 860
   
!     SEARCH LIST OF ELEMENT ID-S. ADD ID TO LIST IF NOT A DUPLICATE
   
   1100 IF (ij /= 0) GO TO 1110
   mset = kset
   iz(mset) = ib
   nelout = 1
   ij = mset
   GO TO ret, (1050,1070)
   1110 DO  j = mset,ij
     IF (iz(j) == ib) GO TO ret, (1050,1070)
   END DO
   ij = ij + 1
   IF (ij < lcore) GO TO 1130
   GO TO 1000
   1130 iz(ij) = ib
   nelout = nelout + 1
   GO TO ret, (1050,1070)
   
!     DONE WITH ALL CASES. IF ALL.NE.1, MOVE THE ID-S UP IN CORE
   
   1200 CALL CLOSE (casecc,1)
   IF (all == 1) GO TO 1220
   
   DO  j = 1,nelout
     iz(nextz+j-1) = iz(mset+j-1)
   END DO
   nextz = nextz + nelout
   1220 CONTINUE
   
   CALL gopen (hcflds,core(lcore+1),1)
   i = lcore - sysbuf
   j = i     - sysbuf
   lcore = j - sysbuf
   IF (lcore <= nextz) GO TO 1580
   CALL gopen (remfls,core(i+1),1)
   CALL gopen (hccens,core(j+1),1)
   CALL gopen (scr6,core(lcore+1),1)
   
!     NO DO LOOP ON IDO. IN EANDM WE WILL READ ALL CARDS
   
   1230 CALL eandm (nobld,ido,nextz,lcore,nbdys,all,nelout)
   GO TO 50
   
   
   1300 IF (ngrold /= ngrav) CYCLE
   CALL pack (core,pg,pg(1))
   1310 iii = iii + 1
   
 END DO
 
 CALL CLOSE (bgpdt,1)
 IF (icm == 0) CALL CLOSE (cstm,1)
 CALL CLOSE (slt,1)
 CALL CLOSE (sil,1)
 IF (ipre /= 1) GO TO 1410
 CALL CLOSE (hcflds,1)
 CALL CLOSE (remfls,1)
 CALL CLOSE (hccens,1)
 CALL CLOSE (scr6,1)
 1410 CONTINUE
 RETURN
 
!     FILE ERRORS
 
 1520 ip1 = -2
 GO TO 1600
 1570 ip1 = -7
 GO TO 1600
 1580 ip1 = -8
 1600 CALL mesage (ip1,FILE,NAME(1))
 RETURN
END SUBROUTINE extern
