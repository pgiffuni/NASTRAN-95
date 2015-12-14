SUBROUTINE rcovui (ub,lastss,modal)
     
!     THIS ROUTINE CALCULATES THE IMPROVED LOWER LEVEL DISPLACEMENTS
!     ON A REDUCED SUBSTRUCTURE WHICH INCLUDE INERTIA AND DAMPING
!     EFFECTS
 
 
 INTEGER, INTENT(IN)                      :: ub
 INTEGER, INTENT(IN)                      :: lastss(2)
 LOGICAL, INTENT(IN OUT)                  :: modal
 LOGICAL :: reqf
 INTEGER :: sof1       ,sof2       ,sof3      ,  &
     buf1       ,buf2       ,scr2       ,scr3      ,  &
     scr4       ,scr5       ,scr6       ,scr7      ,  &
     scr8       ,scr9       ,mpyz       ,tflag     ,  &
     signab     ,horg       ,bmtx       ,uprt      ,  &
     signc      ,scrm       ,z          ,rc        ,  &
     dry        ,fss        ,rfno       ,rss       ,  &
     ua         ,buf4       ,bgg        ,pid       ,  &
     uao        ,rule       ,typa       ,typb      ,  &
     buf3       ,upart      , gims      ,  &
     dua        ,uad        ,typin      ,typot     , typc       ,NAME(2)
 REAL :: rz(1)
 DOUBLE PRECISION :: dz(1)
 COMMON /BLANK /  dry        ,loop       ,step       ,fss(2)    ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5),  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/  icore      ,lcore      ,buf1       ,buf2      ,  &
     buf3       ,buf4       ,sof1       ,sof2      , sof3
 COMMON /rcovcm/  mrecvr     ,ua         ,pa         ,qa        ,  &
     iopt       ,rss(2)     ,energy     ,uimpro    ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /zzzzzz/  z(1)
 COMMON /mpyadx/  mcba(7)    ,mcbb(7)    ,mcbc(7)    ,mcbd(7)   ,  &
     mpyz       ,tflag      ,signab     ,signc     , mprec      ,scrm
 COMMON /names /  rd         ,rdrew      ,wrt        ,wrtrew    ,  &
     rew        ,norew      ,eofnrw     ,rsp       ,  &
     rdp        ,csp        ,cdp        ,square    ,  &
     rect       ,diag       ,upper      ,lower     , sym
 COMMON /parmeg/  mcb(7)     ,mcb11(7)   ,mcb21(7)   ,mcb12(7)  ,  &
     mcb22(7)   ,mrgz       ,rule
 COMMON /saddx /  nomat      ,lcor       ,mcbaa(7)   ,typa      ,  &
     alpha      ,alp(3)     ,mcbbb(7)   ,typb      ,  &
     beta       ,bet(3)     ,mcbcc(7)   ,typc      ,  &
     gama       ,gam(3)     ,dum(24)    ,mcbxx(7)
 COMMON /packx /  typin      ,typot      ,iro        ,nro       , incrp
 EQUIVALENCE      (dz(1),rz(1),z(1))
 DATA    scr2  ,  scr3,scr4,scr5,scr6,scr7,scr8,scr9 /  &
     302   ,  303 ,304 ,305 ,306 ,307 ,308 ,309  /
 DATA    horg  ,  mmtx,bmtx,uprt /  4HHORG,4HMMTX,4HBMTX,4HUPRT /
 DATA    k4mx  /  4HK4MX/,  k4gg / 110  /
 DATA    gims  ,  nhpdat/ 4HGIMS,4HPDAT /
 DATA    mgg   ,  bgg   / 104,109/
 DATA    NAME  / 4 hrcov, 4HUI   /
 
!     INITILIZE
 
 lcorez = korsz(z) - lreq - icore - 1
 idpcor = icore/2 + 1
 tflag  = 0
 signab = 1
 signc  = 1
 mprec  = 0
 scrm   = 309
 reqf   = .false.
 IF (lastss(1) == fss(1) .AND. lastss(2) == fss(2)) reqf = .true.
 
!     GENERATE THE PARTIAL LOAD VECTOR USING THE NORMAL TRANSFORMATION
 
!     UPARTIAL = HORG*UB
 
 item = horg
 CALL mtrxi (scr2,lastss,horg,0,rc)
 IF (rc /= 1) GO TO 6317
 
 mcba(1) = scr2
 CALL rdtrl (mcba)
 mcbb(1) = ub
 CALL rdtrl (mcbb)
 mcbc(1) = 0
 upart   = scr5
 CALL makmcb (mcbd,upart,mcba(3),rect,mcbb(5))
 mpyz    = lcorez
 CALL sofcls
 CALL mpyad (dz(idpcor),dz(idpcor),dz(idpcor))
 CALL wrttrl (mcbd)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     DETERMINE THE NUMBER OF OMITTED POINTS
 
 nrowo = mcba(3) - mcba(2)
 CALL softrl (lastss,gims,mcba)
 IF (mcba(1) == 1) nrowo = mcba(3)
 
!     GENERATE THE VELOCITIES AND ACCELERATIONS
 
 lcore = buf4 - icore - 1
 CALL rcovva (upart,0,0,0,scr7,scr8,lastss,dz(idpcor),dz(idpcor), dz(idpcor))
 IF (upart <= 0) GO TO 9000
 
!     CALCULATE THE INERTIAL AND DAMPING LOADS
 
!     PID = -M*A - B*V
 
!     CALCULATE THE INERTAIL LOADS
 
 pid = 0
 IF (.NOT.reqf) GO TO 100
 mcba(1) = mgg
 IF (mcba(1) > 0) GO TO 110
 100 CALL mtrxi (scr2,lastss,mmtx,0,rc)
 IF (rc /= 1) GO TO 200
 mcba(1) = scr2
 CALL rdtrl (mcba)
 110 mcbb(1) = scr8
 CALL rdtrl (mcbb)
 mcbc(1) = 0
 CALL makmcb (mcbd,scr6,mcbb(3),rect,mcbb(5))
 signab = -1
 CALL sofcls
 
 CALL mpyad (dz(idpcor),dz(idpcor),dz(idpcor))
 
 DO  i = 1,7
   mcbc(i) = mcbd(i)
 END DO
 pid = scr6
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     CALCULATE THE DAMPING LOADS
 
 200 IF (rfno == 3) GO TO 300
 IF (.NOT.reqf  ) GO TO 201
 mcba(1) = k4gg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 202
 201 CALL mtrxi (scr2,lastss,k4mx,0,rc)
 IF (rc /= 1) GO TO 209
 mcba(1) = scr2
 CALL rdtrl (mcba)
 202 mcbb(1) = scr7
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr8,mcbb(3),rect,mcbb(5))
 signab = -1
 CALL sofcls
 CALL mpyad (dz(idpcor),dz(idpcor),dz(idpcor))
 pid = scr8
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 DO  i = 1,7
   mcbc(i) = mcbd(i)
 END DO
 
 209 IF (.NOT.reqf) GO TO 210
 mcba(1) = bgg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 220
 210 CALL mtrxi (scr2,lastss,bmtx,0,rc)
 IF (rc /= 1) GO TO 300
 mcba(1) = scr2
 CALL rdtrl (mcba)
 220 mcbb(1) = scr7
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr6,mcbb(3),rect,mcbb(5))
 signab = -1
 CALL sofcls
 CALL mpyad(dz(idpcor),dz(idpcor),dz(idpcor))
 pid = scr6
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     PARTITION THE INERTIA AND DAMPING LOADS TO THE OMIT SET
 
!     GET THE PARTITIONING VECTOR FROM THE SOF
 
 300 IF (pid == 0) GO TO 400
 item = uprt
 CALL mtrxi (scr2,lastss,uprt,0,rc)
 IF (rc /= 1) GO TO 6317
 rule = 0
 mrgz = lcorez - 14
 idp  = (icore+14)/2 + 1
 DO  i = 1,7
   mcb(i) = mcbd(i)
 END DO
 pid  = scr4
 CALL makmcb (mcb11,pid,nrowo,rect,mcbd(5))
 mcb11(2) = mcbd(2)
 mcb12(1) = 0
 mcb21(1) = 0
 mcb22(1) = 0
 
!     SET UP A NULL ROW PARTITION VECTOR
 
 z(icore) = scr2
 CALL rdtrl (z(icore))
 CALL makmcb (z(icore+7),0,mcb(2),rect,rsp)
 z(icore+8) = 1
 CALL sofcls
 CALL partn (z(icore+7),z(icore),dz(idp))
 CALL wrttrl (mcb11)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     PERFORM THE FBS TO GET THE LOADS ON THE OMMITTED POINTS.  WE
!     WILL ALSO ADD IN THE EFFECTS OF THE DAMPING AND INERTIAL LOADS
 
 400 CALL rcovuo (pid,uao,lastss)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 IF (iopt < 0) GO TO 9000
 
!     IF RECOVERING A MODAL REDUCED SUBSTRUCTURE, CALCULATE
!     THE MODAL CORRECTION TO THE U PARTIAL
 
 dua = 0
 IF (.NOT.modal) GO TO 900
 
!     IF RF-9, SPLIT THE DISPLACEMENTS FROM THE TOTAL VECTOR
 
 uad = upart
 IF (rfno /= 9) GO TO 500
 uad = scr9
 CALL rcovva (upart,1,0,uad,0,0,lastss,dz(idpcor),dz(idpcor), dz(idpcor))
 
!     PARTITION THE PARTIAL DISPLACEMENTS TO THE OMITTED AND
!     BOUNDARY SIZES
 
 500 item = uprt
 CALL mtrxi (scr2,lastss,uprt,0,rc)
 IF (rc /= 1) GO TO 6317
 rule = 0
 mrgz = lcorez - 14
 idp  = (icore + 14)/2 + 1
 mcb(1) = uad
 CALL rdtrl (mcb)
 CALL makmcb (mcb11,scr3,nrowo,rect,mcb(5))
 CALL makmcb (mcb21,scr4,mcb(3)-nrowo,rect,mcb(5))
 mcb11(2) = mcb(2)
 mcb21(2) = mcb(2)
 mcb12(1) = 0
 mcb22(1) = 0
 
 z(icore) = scr2
 CALL rdtrl (z(icore))
 CALL makmcb (z(icore+7),0,mcb(2),rect,rsp)
 z(icore+8) = 1
 CALL sofcls
 
 CALL bug (nhpdat,500,mcb(1),37)
 CALL partn (z(icore+7),z(icore),dz(idp))
 CALL wrttrl (mcb11)
 CALL wrttrl (mcb21)
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     CALCULATE THE CORRECTION TERMS
 
!     DUO = GI*UB - UO
 
 item = gims
 CALL mtrxi (scr6,lastss,gims,0,rc)
 IF (rc /= 1) GO TO 6317
 mcba(1) = scr6
 CALL rdtrl (mcba)
 DO  i = 1,7
   mcbb(i) = mcb21(i)
   mcbc(i) = mcb11(i)
 END DO
 CALL makmcb (mcbd,scr9,mcba(3),rect,mcbb(5))
 signab = 1
 signc  =-1
 tflag  = 0
 scrm   = 308
 mprec  = 0
 CALL sofcls
 mpyz   = mrgz
 CALL mpyad (dz(idp),dz(idp),dz(idp))
 CALL wrttrl (mcbd)
 
!     MERGE DUO TO -A- SIZE
 
 DO  i = 1,7
   mcb11(i) = mcbd(i)
 END DO
 mcb21(1) = 0
 dua = scr4
 CALL makmcb (mcb,dua,z(icore+2),rect,mcb11(5))
 mcb(2) = mcbd(2)
 IF (rfno == 9) mcb(2) = 3*mcbd(2)
 
!     SET UP A NULL ROW PARTITIONING VECTOR (OR FOR RF-9)
!     SET UP A VECTOR THAT WILL MERGE IN A NULL VELOCITY AND
!     ACCELERATION VECTOR FOR EACH DISPLACEMENT VECTOR
 
 nro = mcb(2)
 CALL makmcb (z(icore+7),scr3,nro,rect,rsp)
 IF (nro+15 > lcorez) GO TO 9008
 DO  i = 1,nro
   rz(icore+14+i) = 0.0
 END DO
 IF (rfno /= 9) GO TO 570
 DO  i = 1,nro,3
   rz(icore+15+i) = 1.0
   rz(icore+16+i) = 1.0
 END DO
 570 CONTINUE
 CALL gopen (scr3,z(buf1),wrtrew)
 typin = 1
 typot = 1
 iro   = 1
 incrp = 1
 CALL pack (z(icore+15),scr3,z(icore+7))
 CALL CLOSE (scr3,rew)
 CALL wrttrl (z(icore+7))
 CALL merge (z(icore+7),z(icore),dz(idp))
 CALL wrttrl (mcb)
 
!     ADD THE PARTIAL DISPLACEMENT VECTOR TO THE DISPLACEMENTS FROM
!     THE OMITS, INERTIAL, DAMPING, AND MODAL CORRECTION EFFECTS
!     TO GET THE FINAL DISPLACEMENT VECTOR FOR THIS SUBSTRUCTURE
 
 900 nomat = 2
 IF (dua /= 0) nomat = 3
 typa  = 1
 alpha = 1.0
 mcbaa(1) = upart
 CALL rdtrl (mcbaa)
 typb = 1
 beta = 1.0
 mcbbb(1) = uao
 CALL rdtrl (mcbbb)
 IF (dua == 0) GO TO 910
 typc = 1
 gama = 1.0
 mcbcc(1) = dua
 CALL rdtrl (mcbcc)
 910 CALL makmcb (mcbxx,ua,mcbaa(3),rect,mcbaa(5))
 mcbxx(2) = mcbaa(2)
 lcor = lcorez
 CALL sofcls
 CALL sadd (dz(idpcor),dz(idpcor))
 CALL wrttrl (mcbxx)
 
!     NORMAL RETURN
 
 signab = 1
 RETURN
 
!     ERROR MESSAGES
 
 6317 IF (rc == 2) rc = 3
 CALL smsg (rc-2,item,lastss)
 9000 iopt = -1
 RETURN
 
 9008 CALL mesage (8,0,NAME)
 GO TO 9000
END SUBROUTINE rcovui
