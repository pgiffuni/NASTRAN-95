SUBROUTINE subph1
     
!     THIS MODULE PERFORMS THE PHASE 1 CONVERSION OF NASTRAN DATA BLOCK
!     TABLES TO THEIR EQUIVALENT SOF ITEMS
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,andf,orf
 LOGICAL :: last
 INTEGER :: buf(10),temp(10),TYPE,sub1(2),icode(32),mcb(7),  &
     ltype1(5),ltype2(5),ltype3(5)
 REAL :: rz(12)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ bsize,out
 COMMON /BLANK / dry,NAME(2),pset,pitm
 COMMON /two   / two(32)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (rz(1),z(1))
DATA    case  , eqex,uset,bgpd,cstm,gpse,ELSE,scrt/  &
     101   , 102 ,103 ,104 ,105 ,106 ,107 ,301 /
 DATA    eqss  / 4HEQSS/,icstm/4HCSTM/,lods /4HLODS/,plts/4HPLTS/,  &
     bgss  / 4HBGSS/
 DATA    iua   / 25    /, sub1/4HSUBP,4HH1   /
 DATA    ltype1/ 4HEXTE,4HRNAL,4H sta,4HTIC ,4HLOAD/
 DATA    ltype2/ 4H    ,4H    ,4HTHER,4HMAL ,4HLOAD/
 DATA    ltype3/ 4H ele,4HMENT,4H def,4HORMA,4HTION/
 DATA    loap  , papp  /4HLOAP,4HPAPP/,  i0 / 0    /
 
 mua = two(iua)
 
!     INITIALLIZE CORE, ETC
 
 IF (dry == 0) RETURN
 nc = korsz(z(1))
 b1 = nc - bsize + 1
 
!     OPEN SCRATCH FILE TO WRITE CONVERTED DATA
 
 b2   =  b1  - bsize
 b3   =  b2  - bsize
 buf1 =  b3  - bsize
 buf2 =  buf1- bsize
 nz   =  buf2- 1
 
!     TEST FOR CORE
 
 IF (nz <= 0) GO TO 4010
 
 CALL sofopn (z(b1),z(b2),z(b3))
 
!     EQSS GENERATION
 
 FILE = uset
 CALL OPEN (*5001,uset,z(buf1),0)
 CALL fwdrec (*5001,uset)
 
!     READ USET INTO CORE
 
 CALL READ (*5001,*20,uset,z(1),nz,0,nu)
 
!     RAN OUT OF CORE
 
 CALL CLOSE (uset,1)
 GO TO 4010
 
 20   CALL CLOSE (uset,1)
 
!     FLAG ELEMENTS IN UA SET  (SET OTHERS TO ZERO)
 
 DO  i = 1,nu
   IF (andf(mua,z(i)) == 0) GO TO 30
   z(i) = 1
   CYCLE
   30   z(i) = 0
 END DO
 
!     READ  SECOND RECORD OF EQEXIN - CONTAINS  G AND SIL PAIRS
 
 FILE = eqex
 CALL OPEN (*5001,eqex,z(buf1),0)
 CALL fwdrec (*5001,eqex)
 CALL fwdrec (*5001,eqex)
 
!     OPEN SCRATCH FILE TO WRITE CONVERTED DATA
 
 CALL OPEN (*5001,scrt,z(buf2),1)
 
!     LOOP ON GRID POINTS
 
 k = 0
 i = 0
 
 50   CALL READ (*5001,*110,eqex,buf,2,0,nwds)
 c = 0
 i = i + 1
 isil = buf(2)/10
 TYPE = buf(2) - 10*isil
 IF (TYPE-2) 60,80,4020
 
!     GRID POINT, DETERMINE UA COMPONENTS, PUT IN BINARY FORM
 
 60   DO  j = 1,6
   iu = isil + j - 1
   IF (z(iu) == 0) CYCLE
   c  = orf(c,lshift(1,j-1))
 END DO
 GO TO 90
 
!     SCALAR POINT
 
 80   IF (z(isil) /= 0) c = 1
 
!     WRITE OUT G AND C
 
 90   IF (c == 0) GO TO 100
 buf(2) = c
 CALL WRITE (scrt,buf,2,0)
 k = k + 1
 100  CONTINUE
 GO TO 50
 
 110  mcb(1) = eqex
 CALL rdtrl (mcb)
 npts = mcb(2)
 CALL REWIND (eqex)
 CALL CLOSE (scrt,1)
 IF (npts*2 > nz) GO TO 4010
 
!     READ FIRST RECORD OF EQEXIN - GET G AND IOLD
!     READ SCRATCH - GET G AND C
!     BUILD TABLE IN CORE
 
 FILE = eqex
 CALL fwdrec (*5001,eqex)
 FILE = scrt
 CALL OPEN (*5001,scrt,z(buf2),0)
 
!     SET CORE TO ZERO
 
 DO  i = 1,npts
   izp = 2*i
   z(izp  ) = 0
   z(izp-1) = 0
 END DO
 nnew = k
 
!     LOOP ON POINTS IN SCRATCH FILE, STORE C IN ITH WORD OF ENTRY
!     POSITION OF ENTRY IS THE INTERNAL SEQUENCE
 
 IF (k <= 0) GO TO 210
 DO  i = 1,k
   FILE = scrt
   CALL READ (*5001,*210,scrt,buf,2,0,nwds)
   FILE = eqex
   180  CALL READ (*5001,*210,eqex,temp,2,0,nwds)
   IF (buf(1)-temp(1) < 0.0) THEN
     GO TO  5001
   ELSE IF (buf(1)-temp(1) == 0.0) THEN
     GO TO   190
   ELSE
     GO TO   180
   END IF
   190  izp = 2*temp(2)
   z(izp) = buf(2)
 END DO
 
!     CORE TABLE IS COMPLETE, FILL IN FIRST ENTRIES
 
 210  CALL CLOSE (scrt,1)
 CALL REWIND (eqex)
 k = 0
 DO  i = 1,npts
   IF (z(2*i) == 0) CYCLE
   k = k + 1
   z(2*i-1) = k
 END DO
 
!     CORE NOW CONTAINS NEW IP VALUES AND C IN OLD IP POSITIONS
 
 FILE = eqss
 
!     CHECK IF SUBSTRUCTURE EXISTS ALREADY
 
 CALL fwdrec (*5001,eqex)
 CALL setlvl (NAME,0,temp,itest,0)
 IF (itest /= 1) WRITE (out,6325) uwm,NAME
 itest = 3
 CALL sfetch (NAME,eqss,2,itest)
 IF (itest == 3) GO TO 340
 WRITE (out,6326) uwm,NAME,eqss
 GO TO 1000
 340  buf(1) = NAME(1)
 buf(2) = NAME(2)
 buf(3) = 1
 buf(4) = nnew
 buf(5) = NAME(1)
 buf(6) = NAME(2)
 
 CALL suwrt (buf,6,2)
 
!     PROCESS EQSS OUTPUT-  G, IP, C - SORTED ON G
 
 DO  i = 1,npts
   
   CALL READ (*5001,*400,eqex,temp,2,0,nwds)
   
   ipt = temp(2)*2 - 1
   IF (z(ipt) == 0) CYCLE
   temp(2) = z(ipt  )
   temp(3) = z(ipt+1)
   CALL suwrt (temp,3,1)
 END DO
 CALL suwrt (temp,0,2)
 
!     BUILD SIL TABLE BY COUNTING C VALUES
 
 nc = 0
 is = 1
 DO  i = 1,npts
   ipt = 2*i - 1
   
   IF (z(ipt) == 0) CYCLE
   is = is + nc
   z(ipt) = is
   
   CALL suwrt (z(ipt),2,1)
   
!     CALCULATE NUMBER OF COMPONENTS FOR NEXT STEP
   
   kcode = z(ipt+1)
   CALL decode (kcode,icode,nc)
 END DO
 CALL suwrt (0,0,2)
 CALL suwrt (temp,0,3)
 1000 CALL CLOSE (eqex,1)
 
!     BGSS GENERATION
 
 FILE = bgpd
 CALL OPEN (*5001,bgpd,z(buf1),0)
 CALL fwdrec (*5001,bgpd)
 itest = 3
 CALL sfetch (NAME,bgss,2,itest)
 IF (itest == 3) GO TO 1100
 WRITE (out,6326) uwm,NAME,bgss
 GO TO 2000
 1100 CONTINUE
 
 buf(1) = NAME(1)
 buf(2) = NAME(2)
 buf(3) = nnew
 CALL suwrt (buf,3,2)
 DO  i = 1,npts
   CALL READ (*5001,*1200,bgpd,buf,4,0,nwds)
   
   IF (z(2*i-1) == 0) CYCLE
   
   CALL suwrt (buf,4,1)
 END DO
 CALL suwrt (0,0,2)
 CALL suwrt (buf,0,3)
 2000 CALL CLOSE (bgpd,1)
 
 
!     CSTM GENERATION
 
 
 CALL OPEN (*2500,cstm,z(buf1),0)
 
!     CSTM EXISTS
 
 CALL fwdrec (*5001,cstm)
 itest = 3
 CALL sfetch (NAME,icstm,2,itest)
 IF (itest == 3) GO TO 2100
 WRITE (out,6326) uwm,NAME,icstm
 GO TO 2400
 
 2100 buf(1) = NAME(1)
 buf(2) = NAME(2)
 CALL suwrt (buf,2,2)
 
!     BLAST COPY
 
 CALL READ (*5001,*2200,cstm,z(1),nz,1,nwds)
 GO TO 4010
 2200 CALL suwrt (z(1),nwds,2)
 CALL suwrt (0,0,3)
 2400 CALL CLOSE (cstm,1)
 
!     LODS GENERATION
 
 2500 nlod = 0
 
 CALL gopen (case,z(buf1),0)
 
 icase = 0
 
 2600 CALL READ (*2800,*2800,case,z(1),9,1,nwds)
 icase = icase + 1
 IF (z(i0+4) == 0) GO TO 2610
 WRITE (out,6327) uim,NAME,icase,ltype1,z(i0+4)
 z(nlod+10) = z(i0+4)
 GO TO 2700
 2610 IF (z(i0+7) == 0) GO TO 2620
 WRITE (out,6327) uim,NAME,icase,ltype2,z(i0+7)
 z(nlod+10) = z(i0+7)
 GO TO 2700
 2620 IF (z(i0+6) == 0) GO TO 2630
 WRITE (out,6327) uim,NAME,icase,ltype3,z(i0+6)
 z(nlod+10) = z(i0+6)
 GO TO 2700
 2630 z(nlod+10) = 0
 2700 nlod = nlod + 1
 GO TO 2600
 2800 itest = 3
 litm  = lods
 IF (pitm == papp) litm = loap
 CALL sfetch (NAME,litm,2,itest)
 IF (itest == 3) GO TO 2810
 WRITE (out,6326) uwm,NAME,litm
 GO TO 2900
 2810 z(   1) = NAME(1)
 z(i0+2) = NAME(2)
 z(i0+3) = nlod
 z(i0+4) = 1
 z(i0+5) = NAME(1)
 z(i0+6) = NAME(2)
 CALL suwrt (z(1),6,2)
 CALL suwrt (nlod,1,1)
 CALL suwrt (z(i0+10),nlod,2)
 CALL suwrt (z(1),0,3)
 2900 CALL CLOSE (case,1)
 
!     PLOT SET DATA (PLTS) GENERATION
 
 IF (pset <= 0) GO TO 4000
 FILE = bgpd
 CALL gopen (bgpd,z(buf1),0)
 
 itest = 3
 CALL sfetch (NAME,plts,2,itest)
 IF (itest == 3) GO TO 3010
 WRITE (out,6326) uwm,NAME,plts
 CALL CLOSE (bgpd,1)
 GO TO 4000
 
 3010 buf(1) = NAME(1)
 buf(2) = NAME(2)
 buf(3) = 1
 buf(4) = NAME(1)
 buf(5) = NAME(2)
 CALL suwrt (buf,5,1)
 DO  i = 1,11
   z(i) = 0
 END DO
 rz( 4) = 1.0
 rz( 8) = 1.0
 rz(12) = 1.0
 CALL suwrt (z,12,2)
 
 CALL READ (*5001,*3020,bgpd,z(1),nz,0,nwds)
 GO TO 4010
 3020 CALL suwrt (z,nwds,2)
 CALL CLOSE (bgpd,1)
 FILE = eqex
 CALL gopen (eqex,z(buf1),0)
 CALL READ (*5001,*3030,eqex,z,nz,1,nwds)
 GO TO 4010
 3030 CALL suwrt (z,nwds,2)
 CALL CLOSE (eqex,1)
 FILE = gpse
 last = .false.
 CALL OPEN (*3500,gpse,z(buf1),0)
 
 CALL fwdrec (*3500,gpse)
 
 CALL READ (*5001,*3050,gpse,z(1),nz,0,nsets)
 GO TO 4010
 
!     FIND PLOT SET ID
 
 3050 IF (nsets == 0) GO TO 3500
 
 DO  i = 1,nsets
   IF (z(i) == pset) GO TO 3070
 END DO
 GO TO 3500
 3070 irec = i - 1
 
 3075 IF (irec == 0) GO TO 3090
 
!     POSITION FILE TO SELECTED SET
 
 DO  i = 1,irec
   CALL fwdrec (*3500,FILE)
 END DO
 3090 CALL READ (*3500,*3100,FILE,z(1),nz,0,nwds)
 GO TO 4010
 3100 CALL suwrt (z(1),nwds,2)
 CALL CLOSE (FILE,1)
 IF (last) GO TO 3300
 last = .true.
FILE = ELSE
CALL OPEN (*3500,ELSE,z(buf1),0)
CALL fwdrec (*3500,ELSE)
 GO TO 3075
 
!     FINISHED
 
 3300 CALL suwrt (z(1),0,3)
 GO TO 4000
 3500 CALL CLOSE (FILE,1)
 WRITE  (out,3510) uwm,pset
 3510 FORMAT (a25,' 6050, REQUESTED PLOT SET NO.',i8, ' HAS NOT BEEN DEFINED')
 
 4000 CALL sofcls
 WRITE  (out,6361) uim,NAME
 6361 FORMAT (a29,' 6361, PHASE 1 SUCCESSFULLY EXECUTED FOR ',  &
     'SUBSTRUCTURE ',2A4)
 RETURN
 
!     INSUFFICIENT CORE
 
 4010 WRITE  (out,4015) ufm,nz
 4015 FORMAT (a23,' 6011, INSUFFICIENT CORE TO LOAD TABLES', /5X,  &
     'IN MODULE SUBPH1, CORE =',1I8)
 dry = -2
 GO TO 6000
 
!     BAD GRID POINT TYPE (IE AXISYMMETRIC OR)
 
 4020 WRITE  (out,4035) ufm,buf(1)
 4035 FORMAT (a23,' 6013 , ILLEGAL TYPE OF POINT DEFINED FOR ',  &
     'SUBSTRUCTURE ANALYSIS.', /5X,'POINT NUMBER =',i9)
 GO TO 6000
 
!     BAD FILE
 
 5001 WRITE  (out,5005) sfm,FILE
 5005 FORMAT (a25,' 6012, FILE =',i4,' IS PURGED OR NULL AND IS ',  &
     'REQUIRED IN PHASE 1 SUBSTRUCTURE ANALYSIS.')
 
 6000 CALL sofcls
 CALL mesage (-61,0,sub1)
 RETURN
 
 
 6325 FORMAT (a25,' 6325, SUBSTRUCTURE PHASE 1, BASIC SUBSTRUCTURE ',  &
     2A4,' ALREADY EXISTS ON SOF.', /32X,  &
     'ITEMS WHICH ALREADY EXIST WILL NOT BE REGENERATED.')
 6326 FORMAT (a25,' 6326, SUBSTRUCTURE ',2A4,', ITEM ',a4,  &
     ' ALREADY EXISTS ON SOF.')
 6327 FORMAT (a29,' 6327, SUBSTRUCTURE ',2A4,' SUBCASE',i9,  &
     ' IS IDENTIFIED BY', /36X,5A4,' SET',i9,' IN LODS ITEM.',  &
     /36X,'REFER TO THIS NUMBER ON LOADC CARDS.')
END SUBROUTINE subph1
