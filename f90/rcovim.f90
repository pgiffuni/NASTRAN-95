SUBROUTINE rcovim (higher)
     
!     THIS SUBROUTINE CALCULATES THE ENERGIES ON THE MODAL COORDINATES
!     IN A SUBSTRUCTURE THAT WAS MODAL REDUCED.  IT WILL ALSO
!     CALCULATE THE TOTAL ENERGY FOR EACH COLUMN.
 
 
 INTEGER, INTENT(IN OUT)                  :: higher(2)
 INTEGER :: fss        ,rss        ,ua         ,rfno      ,  &
     z          ,rc         ,sof1       ,sof2      ,  &
     sof3       ,buf1       ,buf2       ,buf3      ,  &
     tflag      ,signab     ,signc      ,scrm      ,  &
     scr5       ,scr6       ,scr7       ,scr8      ,  &
     scr9       ,NAME(2)    ,buf4       ,FILE      , rsp        , uvec
 REAL :: rz(1)
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /zzzzzz/ z(1)
 COMMON /mpyadx/ mcba(7)    ,mcbb(7)    ,mcbc(7)    ,mcbd(7)    ,  &
     mpyz       ,tflag      ,signab     ,signc      , mprec      ,scrm
 COMMON /unpakx/ itinu      ,iru        ,nru        ,incru
 COMMON /packx / itinp      ,itoutp     ,irp        ,nrp        , incrp
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw     ,rsp        ,  &
     rdp        ,csp        ,cdp        ,square     ,  &
     rect       ,diag       ,upper      ,lower      , sym
 EQUIVALENCE     (z(1),rz(1))
 DATA    uvec  , kmtx,mmtx  / 4HUVEC,4HKMTX,4HMMTX /
 DATA    scr5  , scr6,scr7,scr8,scr9 / 305,306,307,308,309 /
 DATA    NAME  / 4HRCOV,4HIM         /
 
!     INITIALIZE
 
 lcorez = korsz(z)
 mpyz   = lcorez
 tflag  = 0
 signab = 1
 signc  = 1
 mprec  = 0
 
!     GET THE DISPLACEMENT VECTOR FOR THE HIGHER LEVEL REDUCED
!     SUBSTRUCTURE.
 
 item = uvec
 CALL mtrxi (scr5,higher,uvec,0,rc)
 IF (rc /= 1) GO TO 6000
 
!     CALCULATE VELOCITIES IF NOT ALREADY DONE FOR THE OUTPUT PHASE.
 
 intyp = 1
 IF (rfno == 3 .OR. rfno == 8) intyp = 0
 CALL rcovva (scr5,intyp,0,scr8,scr9,0,higher,z(1),z(1),z(1))
 IF (ua <= 0) GO TO 9200
 
!     CALCULATE THE KENETIC ENERTY MULTIPLIER - M * V
 
 item = mmtx
 CALL mtrxi (scr5,higher,mmtx,0,rc)
 IF (rc /= 1) GO TO 6000
 mcba(1) = scr5
 CALL rdtrl (mcba)
 mcbb(1) = scr9
 CALL rdtrl (mcbb)
 ncol    = mcbb(2)
 mcbc(1) = 0
 CALL makmcb (mcbd,scr7,mcbb(3),rect,mcbb(5))
 scrm = scr6
 CALL sofcls
 CALL mpyad (z(1),z(1),z(1))
 CALL wrttrl (mcbd)
 
!     CALCULATE THE KENETIC ENERGIES BY PERFORMING THE SCALAR
!     MULTIPLY IN SINGLE PERCISION.  USE ONLY THE REAL PART IF COMPLEX
!     VECTORS.  APPEND THE TOTAL KINETIC ENERGY TO THE END OF EACH
!     COLUMN.
 
 itinu = rsp
 iru   = 1
 nru   = mcbd(3)
 incru = 1
 itinp = rsp
 itoutp= rsp
 irp   = 1
 nrp   = nru + 1
 incrp = 1
 ivec1 = 1
 ivec2 = ivec1 + nru + 1
 IF (ivec2+nru+1 > sof3) GO TO 9008
 
 FILE = scr9
 CALL gopen (scr7,z(sof1),rdrew)
 CALL gopen (scr9,z(sof2),rdrew)
 CALL gopen (scr6,z(sof3),wrtrew)
 CALL makmcb (mcba,scr6,nrp,rect,rsp)
 
 DO  i = 1,ncol
   isk = 1
   CALL unpack (*130,scr7,rz(ivec1))
   isk = 0
   CALL unpack (*130,scr9,rz(ivec2))
   
   total = 0.0
   DO  j = 1,nru
     k = j - 1
     rz(ivec1+k) = rz(ivec1+k)*rz(ivec2+k)
     total = total + rz(ivec1+k)
   END DO
   rz(ivec1+nru) = total
   GO TO 150
   
   130 DO  j = 1,nrp
     rz(ivec1+j-1) = 0.0
   END DO
   IF (isk /= 0)CALL fwdrec (*9002,scr9)
   
   150 CALL pack (rz(ivec1),scr6,mcba)
   
 END DO
 
 CALL CLOSE (scr7,rew)
 CALL CLOSE (scr9,rew)
 CALL CLOSE (scr6,rew)
 CALL wrttrl (mcba)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     CALCULATE THE POTENTIAL ENERTY MULTPLYIER - K*U
 
 item = kmtx
 CALL mtrxi (scr5,higher,kmtx,0,rc)
 IF (rc /= 1) GO TO 6000
 mcba(1) = scr5
 CALL rdtrl (mcba)
 mcbb(1) = scr8
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr9,mcbb(3),rect,mcbb(5))
 scrm = scr7
 CALL sofcls
 CALL mpyad (z(1),z(1),z(1))
 CALL wrttrl (mcbd)
 
!     CALCULATE THE POTENTIAL ENERGIES BY PERFORMING THE SCALAR
!     MULTIPLY IN SINGLE PERCISION.  USE ONLY THE REAL PART IF COMPLEX
!     VECTORS.  APPEND THE TOTAL POTENTIAL ENERGY TO THE END OF EACH
!     COLUMN.
 
 itinu = rsp
 iru   = 1
 nru   = mcbd(3)
 incru = 1
 itinp = rsp
 itoutp= rsp
 irp   = 1
 nrp   = nru + 1
 incrp = 1
 
 FILE = scr8
 CALL gopen (scr9,z(sof1),rdrew)
 CALL gopen (scr8,z(sof2),rdrew)
 CALL gopen (scr7,z(sof3),wrtrew)
 CALL makmcb (mcba,scr7,nrp,rect,rsp)
 
 DO  i = 1,ncol
   isk = 1
   CALL unpack (*230,scr9,rz(ivec1))
   isk = 0
   CALL unpack (*230,scr8,rz(ivec2))
   total = 0.0
   DO  j = 1,nru
     k = j - 1
     rz(ivec1+k) = rz(ivec1+k)*rz(ivec2+k)
     total = total + rz(ivec1+k)
   END DO
   rz(ivec1+nru) = total
   GO TO 250
   
   230 DO  j = 1,nrp
     rz(ivec1+j-1) = 0.0
   END DO
   IF (isk /= 0)CALL fwdrec (*9002,scr8)
   
   250 CALL pack (rz(ivec1),scr7,mcba)
   
 END DO
 
 CALL CLOSE (scr9,rew)
 CALL CLOSE (scr8,rew)
 CALL CLOSE (scr7,rew)
 CALL wrttrl (mcba)
 
!     NORMAL RETURN
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 RETURN
 
!     ERRORS
 
 6000 CALL smsg (rc-2,item,higher)
 GO TO 9200
 9002 n = 2
 GO TO 9100
 9008 n = 8
 9100 CALL mesage (n,FILE,NAME)
 9200 iopt = -1
 RETURN
END SUBROUTINE rcovim
