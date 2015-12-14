SUBROUTINE rcovuo (pid,uao,lastss)
     
!     THIS SUBROUTINE CALCULATES THE FULL SIZE DISPLACEMENT VECTOR ON
!     ANY OMITTED POINTS.  THE OPTIONAL INERTIA AND DAMPING EFFECTS
!     WILL BE INCLUDED IF REQUESTED.
 
!     FILE USAGE IS AS FOLLOWS
 
!     SCR1 AND SCR5 ARE NOT USED
!     SCR4 CONTAINS PID ON INPUT AND IS DESTROYED
!     SCR7 CONTAINS UAO OUTPUT
!     ALL OTHER SCRATCH FILES ARE USED
 
 INTEGER, INTENT(IN)                      :: pid
 INTEGER, INTENT(OUT)                     :: uao
 INTEGER, INTENT(IN OUT)                  :: lastss(2)
 INTEGER :: rule       ,pao        ,buf1       ,              &
            rc         ,typin      ,typot      ,              &
            pove       ,dry        ,step       ,fss        ,  &
            rfno       ,uinms      ,ua         ,              &
!           SOLN       ,SRD        ,SWRT       ,SCHK       ,  &
            iz(1)      ,rd         ,rdrew      ,wrt        ,  &
            wrtrew     ,rew        ,eofnrw     ,rsp        ,  &
            rdp        ,csp        ,cdp        ,square     ,  &
            rect       ,diag       ,upper      ,lower      ,  &
            sym        ,scra       ,scrb       ,sdcmpz     ,  &
            scr4       ,NAME(2)    ,power      ,FILE       ,  &
            mcbpao(7)  ,scr3       ,umcb       ,bmcb       ,  &
            scrc       ,chlsky     ,xmcb       ,fbsz       ,  &
            prec       ,SIGN       ,scr2       ,              &
            uprt       ,scr7       ,scr6       ,scr8       ,  &
            sof1       ,sof2       ,sof3       ,scr9       ,  &
            typa       ,typb
            
 DOUBLE PRECISION :: dz(1)      ,det        ,deti       ,mindia

 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=27) :: swm 
 CHARACTER (LEN=25) :: uwm, sfm
 CHARACTER (LEN=23) :: ufm
 
 COMMON /xmssg /  ufm        ,uwm        ,uim        ,sfm        , swm
 COMMON /BLANK /  dry        ,loop       ,step       ,fss(2)     ,  &
                  rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
                  nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/  icore      ,lcore      ,buf1       ,buf2       ,  &
                  buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/  mrecvr     ,ua         ,pa         ,qa         ,  &
                  iopt       ,rss(2)     ,energy     ,uimpro     ,  &
                  range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/  sysbuf     ,nout
 COMMON /saddx /  nomat      ,lcor       ,mcbaa(7)   ,typa       ,  &
                  alpha      ,alp(3)     ,mcbbb(7)   ,typb       ,  &
                  beta       ,bet(3)     ,dum(36)    ,mcbxx(7)
 COMMON /names /  rd         ,rdrew      ,wrt        ,wrtrew     ,  &
                  rew        ,norew      ,eofnrw     ,rsp        ,  &
                  rdp        ,csp        ,cdp        ,square     ,  &
                  rect       ,diag       ,upper      ,lower      , sym
 COMMON /packx /  typin      ,typot      ,iro        ,nro        , incrp
 COMMON /parmeg/  mcbk(7)    ,mcbk11(7)  ,mcbk21(7)  ,mcbk12(7)  ,  &
                  mcbk22(7)  ,mrgz       ,rule
 COMMON /sfact /  mcba(7)    ,mcbl(7)    ,mcblt(7)   ,scra       ,  &
                  scrb       ,sdcmpz     ,det        ,deti       ,  &
                  power      ,scrc       ,mindia     ,chlsky
 COMMON /fbsx  /  lmcb(7)    ,umcb(7)    ,bmcb(7)    ,xmcb(7)    ,  &
                  fbsz       ,prec       ,SIGN
 COMMON /zzzzzz/  z(1)
 
 EQUIVALENCE      (z(1),iz(1),dz(1))
 
 DATA    NAME  /  4HRCOV,4HUO   /
 DATA    pove  ,  lmtx / 4HPOVE,4HLMTX /
 DATA    uprt  ,  kmtx / 4HUPRT,4HKMTX /
 DATA    scr2  ,  scr3,scr4,scr6,scr7,scr8,scr9 /  &
         302   ,  303 ,304 ,306 ,307 ,308 ,309  /
 
!     SET UP COMMON BLOCKS
 
 lcorez = korsz(z) - lreq - icore - 1
 idpcor = icore/2 + 1
 rule   = 0
 mcbk21(1) = 0
 mcbk12(1) = 0
 mcbk22(1) = 0
 SIGN   = 1
 
!     CALCUATE THE LOADS ON THE OMMITED POINTS
 
 pao = 0
 IF (rfno == 3) GO TO 10
 pao = scr3
 CALL rcovsl (lastss,pove,0,scr6,scr7,scr8,pao,z(icore),z(icore),  &
     sof3-icore-1,.false.,rfno)
 mcbpao(1) = pao
 CALL rdtrl (mcbpao)
 
!     ADD IN OPTIONAL INERTIA AND DAMPING FORCES TO THE LOADS ON THE
!     OMMITED POINTS
 
 10  IF (pid == 0) GO TO 200
 IF (pao == 0) GO TO 120
 nomat = 2
 typa  = 1
 alpha = 1.0
 mcbaa(1) = pid
 CALL rdtrl (mcbaa)
 typb  = 1
 beta  = 1.0
 mcbbb(1) = pao
 CALL rdtrl (mcbbb)
 CALL makmcb (mcbxx,scr6,mcbaa(3),rect,mcbaa(5))
 mcbxx(2) = mcbaa(2)
 lcor  = lcorez
 CALL sofcls
 CALL sadd (dz(idpcor),dz(idpcor))
 CALL wrttrl (mcbxx)
 DO  i = 1,7
   mcbpao(i) = mcbxx(i)
 END DO
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 GO TO 200
 
!     NO STATIC LOADS SO THE ADD IS UNECESSARY
 
 120 mcbpao(1) = pid
 CALL rdtrl (mcbpao)
 
 200 IF (mcbpao(1) <= 0) GO TO 500
 
!     CHECK FOR EXISTENCE OF LMTX ON THE SOF.  IF IT EXISTS
!     SKIP THE PARTN AND DECOMP
 
 CALL softrl (lastss,lmtx,lmcb(1))
 IF (lmcb(1) /= 1) GO TO 395
 
!     BRING IN LMTX FROM SOF AND SET UP FOR FBS DIRECTLY
 
 CALL mtrxi (scr2,lastss,lmtx,0,rc)
 DO  i = 1,7
   bmcb(i) = mcbpao(i)
 END DO
 lmcb(1) = scr2
 CALL sofcls
 GO TO 411
 
!     COMPUTE THE KOO PARTITION OF KMTX FOR LASTSS
 
!     COPY THE PARTITIONING VECTOR TO SCR2
 
 395 CALL mtrxi (scr2,lastss,uprt,0,rc)
 item = uprt
 IF (rc /= 1) GO TO 6317
 
!     COPY KMTX TO SCR5
 
 item = kmtx
 CALL mtrxi (scr8,lastss,kmtx,0,rc)
 IF (rc /= 1) GO TO 6317
 mcbk(1) = scr8
 CALL rdtrl (mcbk)
 
!     PARTITION KMTX INTO KOO.  STORE KOO ON SCR4.
 
 CALL sofcls
 iz(icore) = scr2
 CALL rdtrl (iz(icore))
 CALL makmcb (mcbk11,scr9,mcbpao(3),sym,mcbk(5))
 mcbk11(2) = mcbpao(3)
 mrgz = lcorez - 7
 i    = (icore+7)/2 + 1
 CALL partn (z(icore),z(icore),dz(i))
 CALL wrttrl (mcbk11)
 
!     DECOMPOSE KOO
 
 DO  i = 1,7
   mcba(i) = mcbk11(i)
 END DO
 CALL makmcb (mcbl,scr2,mcba(3),lower,mcba(5))
 mcblt(1) = scr8
 scra = scr3
 IF (scra == mcbpao(1)) scra = scr6
 scrb = scr4
 IF (scrb == mcbpao(1)) scrb = scr6
 scrc = scr7
 sdcmpz = mrgz
 power  = 1
 chlsky = 0
 CALL sdcomp (*6311,dz(idpcor),dz(idpcor),dz(idpcor))
 CALL wrttrl (mcbl)
 
!     FORWARD AND BACKWARD SUBSTITUTION TO SOLVE FOR UAO
 
 DO  i = 1,7
   lmcb(i) = mcbl(i)
   bmcb(i) = mcbpao(i)
 END DO
 411 fbsz   = lcorez
 mattyp = bmcb(5)
 CALL makmcb (xmcb,scr8,bmcb(3),rect,mattyp)
 prec = 2 - (mattyp-2*(mattyp/2))
 CALL fbs (dz(idpcor),dz(idpcor))
 CALL wrttrl (xmcb)
 
!     MERGE UAO INTO THE UA SET
 
!     COPY UPRT BACK TO SCR2
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 item = uprt
 CALL mtrxi (scr2,lastss,uprt,0,rc)
 IF (rc /= 1) GO TO 6317
 CALL sofcls
 iz(icore) = scr2
 CALL rdtrl (iz(icore))
 
!     SETUP MCB-S IN /PARMEG/
 
 DO  i = 1,7
   mcbk11(i) = xmcb(i)
 END DO
 uao = scr7
 CALL makmcb (mcbk,uao,iz(icore+2),rect,mcbk11(5))
 mcbk(2) = xmcb(2)
 IF (rfno == 9) mcbk(2) = 3*xmcb(2)
 
!     SETUP A NULL ROW PARTITIONING VECTOR OR FOR RIGID FORMAT 9 A
!     VECTOR THAT WILL MERGE IN A NULL VELOCITY AND ACCELERATION
!     VECTOR FOR EACH DISPLACEMENT VECTOR
 
 nro = mcbk(2)
 CALL makmcb (z(icore+7),scr6,nro,rect,rsp)
 IF (nro+15 > lcorez) GO TO 9008
 DO  i = 1,nro
   z(icore+14+i) = 0.0
 END DO
 IF (rfno /= 9) GO TO 440
 DO  i = 1,nro,3
   z(icore+15+i) = 1.0
   z(icore+16+i) = 1.0
 END DO
 440 CONTINUE
 CALL gopen (scr6,z(buf1),wrtrew)
 typin = 1
 typot = 1
 iro   = 1
 incrp = 1
 CALL pack (z(icore+15),scr6,iz(icore+7))
 CALL CLOSE (scr6,rew)
 CALL wrttrl (iz(icore+7))
 
 mrgz = lcorez - 14
 i    = (icore+14)/2 + 1
 CALL merge (z(icore+7),z(icore),dz(i))
 CALL wrttrl (mcbk)
 
!     NORMAL RETURN
 
 RETURN
 
!     NO LOADS SO THE DISPLACEMENTS ARE ZERO
 
 500 uao = 0
 CALL sofcls
 RETURN
 
!     ERROR PROCESSING
 
 6311 WRITE  (nout,6312) swm,lastss
 6312 FORMAT (a27,' 6311, SDCOMP DECOMPOSITION FAILED ON KOO MATRIX ',  &
     'FOR SUBSTRUCTURE ',2A4)
 GO TO 9000
 6317 IF (rc == 2) rc = 3
 CALL smsg (rc-2,item,lastss)
 9000 iopt = -1
 RETURN
 
 9008 n    = 8
 iopt = -1
 CALL sofcls
 CALL mesage (n,FILE,NAME)
 CALL CLOSE (pao,rew)
 CALL CLOSE (scr3,rew)
 
 RETURN
END SUBROUTINE rcovuo
