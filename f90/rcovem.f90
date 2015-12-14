SUBROUTINE rcovem (noexcl,nrowe)
     
!     THIS SUBROUTINE CALCULATES THE ENERGIES ON THE MODAL COORDINATES
!     THAT WERE EXCLUDED FROM THE MODAL REDUCTION PROCESSING
 
 
 LOGICAL, INTENT(OUT)                     :: noexcl
 INTEGER, INTENT(OUT)                     :: nrowe
 
 INTEGER :: rss        ,rc         ,rule       ,sof1       ,  &
     sof2       ,sof3       ,buf1       ,buf2       ,  &
     z          ,qa         ,pa         ,scr3       ,  &
     scr4       ,scr5       ,scr6       ,scr7       ,  &
     scr8       ,phis       ,soln       ,typa       ,  &
     typb       ,tflag      ,signab     ,signc      ,  &
     buf3       ,rsp        ,csp        ,NAME(2)    , rfno
 REAL :: rz(5)
 COMPLEX :: sc         ,cz(2)      ,sc2        ,dkdc       , dksc
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm        ,uwm
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /condas/ phi        ,twophi
 COMMON /system/ sysbuf     ,nout
 COMMON /zzzzzz/ z(1)
 COMMON /packx / itinp      ,itoutp     ,irp        ,nrp        , incrp
 COMMON /unpakx/ itinu      ,iru        ,nru        ,incru
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw     ,rsp        ,  &
     rdp        ,csp        ,cdp        ,square     , rect
 COMMON /parmeg/ mcbp(7)    ,mcbp11(7)  ,mcbp21(7)  ,mcbp12(7)  ,  &
     mcbp22(7)  ,mrgz       ,rule
 COMMON /saddx / nomat      ,lcorez     ,mcbaa(7)   ,typa       ,  &
     alpha      ,alp(3)     ,mcbbb(7)   ,typb       ,  &
     beta       ,bet(3)     ,dum(36)    ,mcbxx(7)
 COMMON /mpyadx/ mcba(7)    ,mcbb(7)    ,mcbc(7)    ,mcbd(7)    ,  &
     mpyz       ,tflag      ,signab     ,signc      , mprec      ,mscr
 EQUIVALENCE     (z(1),rz(1),cz(1))     ,(cz(2),sc)
 DATA    lams  , phis,soln  /4HLAMS,4HPHIS,4HSOLN /
 DATA    scr3  , scr4,scr5,scr6,scr7,scr8 /303,304,305,306,307,308/
 DATA    NAME  / 4HRCOV,4HEM   /
 
!     INITILIZE
 
 lcorez = korsz(z)
 
!     FROM THE LAST GROUP ON LAMS CREATE A PARTITIONING VECTOR TO
!     DIFFERENTIATE THE INCLUDED AND EXCLUDED MODES
 
 nrowe = 0
 item  = lams
 CALL sfetch (rss,lams,1,rc)
 IF (rc /= 1) GO TO 6000
 n = 2
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 i = 0
 
 10 CALL suread (icode,1,n,rc)
 IF (rc /= 1) GO TO 20
 i = i + 1
 IF (i > buf1) GO TO 9008
 rz(i) = 1.0
 IF (icode == 1) GO TO 10
 rz(i) = 0.0
 nrowe = nrowe + 1
 GO TO 10
 
 20 CONTINUE
 IF (nrowe == 0) GO TO 900
 IF (qa+pa == 0) GO TO 9200
 itinp  = rsp
 itoutp = rsp
 irp   = 1
 nrp   = i
 incrp = 1
 CALL makmcb (mcba,scr8,nrp,rect,rsp)
 CALL gopen  (scr8,z(buf1),wrtrew)
 CALL pack   (rz(1),scr8,mcba)
 CALL CLOSE  (scr8,rew)
 CALL wrttrl (mcba)
 
!     PARTITION THE EIGENVECTOR TO GET THE EXCLUDED MODES OUT
 
 item = phis
 CALL mtrxi (scr7,rss,phis,0,rc)
 IF (rc /= 1) GO TO 6000
 rule = 0
 mcbp(1) = scr7
 CALL rdtrl (mcbp)
 IF (mcbp(2) /= nrp) GO TO 6372
 CALL makmcb (mcbp11,scr6,mcbp(3),rect,mcbp(5))
 mcbp11(2) = nrowe
 mcbp21(1) = 0
 mcbp12(1) = 0
 mcbp22(1) = 0
 
!     SETUP NULL COLUMN PARTITONING VECTOR
 
 CALL makmcb (mcbb,0,mcbp(3),rect,rsp)
 mcbb(2) = 1
 mrgz = lcorez
 CALL sofcls
 
 CALL partn (mcba,mcbb,z(1))
 
 CALL wrttrl (mcbp11)
 
!     IF BOTH LOADS AND SINGLE POINT CONSTRAINT FORCES EXIST, ADD
!     THEM TOGETHER
 
 irh = qa + pa
 IF (qa == 0 .OR. pa == 0) GO TO 100
 nomat = 2
 typa  = 1
 alpha = 1.0
 mcbaa(1) = qa
 CALL rdtrl (mcbaa)
 typb = 1
 beta = 1.0
 mcbbb(1) = pa
 CALL rdtrl (mcbbb)
 CALL makmcb (mcbxx,scr7,mcbaa(3),rect,mcbaa(5))
 mcbxx(2) = mcbaa(2)
 
 CALL sadd (z(1),z(1))
 
 CALL wrttrl (mcbxx)
 irh = scr7
 
!     MULTIPLY   PK = QK(T)*(PA + QA)
 
 100 DO  i = 1,7
   mcba(i) = mcbp11(i)
 END DO
 mcbb(i) = irh
 CALL rdtrl (mcbb)
 mcbc(1) = 0
 mpyz   = lcorez
 tflag  = 1
 signab = 1
 signc  = 1
 mprec  = 0
 mscr   = scr8
 CALL makmcb (mcbd,scr5,nrowe,rect,mcba(5))
 
 CALL mpyad (z(1),z(1),z(1))
 
!     READ MODAL MASS AND TWOPHI*FREQUENCY FOR EACH OF THE EXCLUDED
!     MODES MODES FROM LAMS
!     IF MODE WAS EXCLUDED BECAUSE OF NON-PARTICIPATION, SET ITS
!     FREQUENCY TO ZERO
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 item = lams
 CALL sfetch (rss,lams,1,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n <= 0) GO TO 6200
 imode = 8
 CALL suread (z(imode),-1,n,rc)
 IF (rc /= 2) GO TO 6100
 nmode = imode + n - 1
 IF (nmode > buf3) GO TO 9008
 icode = nmode + 1
 CALL suread (z(icode),-1,n,rc)
 IF (rc /= 2 .AND. rc /= 3) GO TO 6100
 ncode = icode + n - 1
 IF (ncode > buf3) GO TO 9008
 
 i1 = imode - 7
 i2 = imode - 2
 DO  i = icode,ncode
   i1 = i1 + 7
   IF (z(i) == 1) CYCLE
   i2 = i2 + 2
   rz(i2) = rz(i1+3)
   IF (z(i) == 2 .OR. rz(i2) <= 0.001) rz(i2) = 0.0
   rz(i2+1) = rz(i1+5)
 END DO
 nmode = i2 + 1
 
!     POSITION SOLN ITEM TO SOLUTION DATA
 
 item = soln
 CALL sfetch (rss,soln,1,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 
!     SET UP TO LOOP OVER COLUMNS
 
 ncol  = mcbd(2)
 nword = 1
 IF (mcbd(5) >= 3) nword = 2
 ivec1 = (nmode/2)*2 + 3
 icvec1= ivec1/2 + 1
 ivec2 = ivec1 + (nrowe*nword/2) * 2 + 1
 IF (ivec2+nrowe > buf3) GO TO 9008
 
 CALL gopen (scr5,z(buf1),rdrew)
 CALL gopen (scr3,z(buf2),wrtrew)
 CALL gopen (scr4,z(buf3),wrtrew)
 CALL makmcb (mcba,scr3,nrowe,rect,rsp)
 CALL makmcb (mcbb,scr4,nrowe,rect,rsp)
 
 itinu = rsp
 IF (mcbd(5) >= 3) itinu = csp
 iru = 1
 nru = nrowe
 incru = 1
 nrp = nrowe
 
!     LOOP OVER EACH SOLUTION STEP
 
 DO  icol = 1,ncol
   
!     GET FREQUENCY OR POLE FROM SOLN ITEM FOR THIS STEP
   
   IF (rfno > 3) GO TO 310
   CALL suread (rz(1),7,n,rc)
   IF (rc /= 1) GO TO 6100
   IF (mcbd(5) >= 3) GO TO 300
   freq = rz(5)
   s2   = -(twophi*freq)**2
   GO TO 320
   
   300 sc  = cz(2)
   sc2 = sc*sc
   GO TO 320
   
   310 CALL suread (rz(1),1,n,rc)
   IF (rc /= 1) GO TO 6100
   sc  = twophi*rz(1)*(0.0,1.0)
   sc2 = sc*sc
   
!     UNPACK THE NEXT COLUMN
   
   320 CALL unpack (*330,scr5,rz(ivec1))
   GO TO 350
   330 DO  i = 1,nrowe
     j = i - 1
     rz(ivec1+j) = 0.0
     rz(ivec2+j) = 0.0
   END DO
   GO TO 600
   
   350 IF (mcbd(5) >= 3) GO TO 500
   
!     CALCULATE ENERGIES FOR REAL MATRICIES
   
   im = imode - 2
   DO  i = 1,nrowe
     im = im + 2
     j  = i - 1
     IF (rz(im) == 0.0 .OR. (twophi*freq) > rz(im)) GO TO 400
     wk2 = rz(im)**2
     
     dkd =-s2*rz(ivec1+j)/(rz(im+1)*wk2**2*(1.0 + s2/wk2))
     dks = rz(ivec1+j)/(rz(im+1)*wk2)
     
     rz(ivec2+j) = .5*rz(im+1)*wk2*ABS((2.0*dks+dkd)*dkd)
     rz(ivec1+j) = ABS(s2/wk2)*rz(ivec2+j)
     CYCLE
     
     400 rz(ivec1+j) = 0.0
     rz(ivec2+j) = 0.0
     
   END DO
   GO TO 600
   
!     CALCULATE ENERGIES FOR COMPLEX VECTORS
   
   500 im = imode - 2
   DO  i = 1,nrowe
     im = im + 2
     j  = i - 1
     IF (rz(im) == 0.0 .OR. AIMAG(sc) > rz(im)) GO TO 510
     wk2 = rz(im)**2
     
     dkdc =-sc2*cz(icvec1+j)/(rz(im+1)*wk2**2*(1.0+sc2/wk2))
     dksc = cz(icvec1+j)/(rz(im+1)*wk2)
     
     rz(ivec2+j) = .5*rz(im+1)*wk2*cabs((2.0*dksc+dkdc)*dkdc)
     rz(ivec1+j) = cabs(sc**2/wk2)*rz(ivec2+j)
     CYCLE
     
     510 rz(ivec1+j) = 0.0
     rz(ivec2+j) = 0.0
     
   END DO
   
!     PACK OUT THE KENETIC AND POTENTIAL ENERGIES
   
   600 CALL pack (rz(ivec1),scr3,mcba)
   CALL pack (rz(ivec2),scr4,mcbb)
   
 END DO
 
 CALL CLOSE (scr5,rew)
 CALL CLOSE (scr3,rew)
 CALL CLOSE (scr4,rew)
 CALL wrttrl (mcba)
 CALL wrttrl (mcbb)
 
!     NORMAL RETURN
 
 RETURN
 
!     NO EXECLUDED MODES EXIST
 
 900 noexcl = .true.
 RETURN
 
!     ERRORS
 
 6000 CALL smsg (rc-2,item,rss)
 GO TO 9200
 6100 CALL smsg (rc+4,item,rss)
 GO TO 9200
 6200 CALL smsg (7,item,rss)
 GO TO 9200
 6372 WRITE (nout,6373) uwm,rss
 GO TO 9200
 9008 CALL mesage (8,0,NAME)
 9200 WRITE (nout,6371) uwm,rss
 noexcl = .true.
 CALL CLOSE (scr3,rew)
 CALL CLOSE (scr4,rew)
 CALL CLOSE (scr5,rew)
 RETURN
 
!     FORMAT STATEMENTS
 
 6371 FORMAT (a25,' 6371, CALCULATIONS FOR EXCLUDED MODE ENERGIES FOR',  &
     ' SUBSTRUCTURE ',2A4,' ABORTED.')
 6373 FORMAT (a25,' 6372, THE PHIS AND LAMS ITEMS ARE INCONSISTANT FOR',  &
     ' SUBSTRUCTURE ',2A4)
END SUBROUTINE rcovem
