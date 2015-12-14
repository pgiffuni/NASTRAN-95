SUBROUTINE rcove
     
!     THIS SUBROUTINE PRINTS THE ENERGIES ON THE MODAL COORDINATES
!     IN A SUBSTRUCTURE THAT WAS MODAL REDUCED.  IT WILL ALSO PRINT
!     THE ENERGIES ON THOSE MODES EXCLUDED FROM THE REDUCTION
!     PROCESSING.
 
 EXTERNAL        andf
 LOGICAL :: mredu      ,credu      ,noexcl
 INTEGER :: dry        ,fss        ,rss        ,ua        ,  &
     rfno       ,z(3)       ,rc         ,buf(1)    ,  &
     FILE       ,sof1       ,sof2       ,sof3      ,  &
     buf1       ,buf2       ,buf3       ,sysbuf    ,  &
     NAME(2)    ,scr3       ,scr4       ,scr6      ,  &
     scr7       ,eqss       ,rsp        ,grid      ,  &
     soln       ,trigid(2)  ,timode(2)  ,temode(2) ,  &
     TYPE(2)    ,casess     ,casecc(2)  ,mcb(7)    ,  &
     cmask      ,energy     ,buf4       ,BLANK     ,  &
     higher(2)  ,andf       ,names(2)
 REAL :: keng       ,peng
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm        ,uwm
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)    ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5),  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2      ,  &
     buf3       ,buf4       ,sof1       ,sof2      , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa        ,  &
     iopt       ,rss(2)     ,energy     ,uimpro    ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /zzzzzz/ rz(1)
 COMMON /system/ sysbuf     ,nout       ,dum1(6)    ,nlpp      ,  &
     dum2(2)    ,nlines
 COMMON /unpakx/ itinu      ,iru        ,nru        ,incru
 COMMON /output/ ititle(96)
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew    ,  &
     rew        ,norew      ,eofnrw     ,rsp       ,  &
     rdp        ,csp        ,cdp        ,square    , rect       ,diag
 EQUIVALENCE     (buf(1)    ,rz(1))
 EQUIVALENCE     (rz(1)     ,z(1))
 DATA    casecc/ 4HCASE,4HCC   /
 DATA    eqss  , lams, soln    / 4HEQSS,4HLAMS,4HSOLN /
 DATA    casess, scr3, scr4, scr6, scr7           /  &
     101   , 303,  304,  306,  307            /
 DATA    trigid/ 4HINER,4HTIAL /
 DATA    timode/ 4HIN-m,4HODE  /
 DATA    temode/ 4HEX-m,4HODE  /
 DATA    ib    / 1             /
 DATA    mmask / 201326592     /
 DATA    cmask / 67108864      /
 DATA    BLANK / 4H            /
 DATA    NAME  / 4HRCOV,4HE    /
 
!     IF THIS IS A STATICS SOLUTION NO ENERGY CALCULATIONS CAN BE MADE
 
 IF (rfno <= 2) RETURN
 
!     INITIALIZE
 
 sof1 = korsz(z) - sysbuf + 1
 sof2 = sof1 - sysbuf - 1
 sof3 = sof2 - sysbuf
 buf1 = sof3 - sysbuf
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 lcore= buf4 - 1
 IF (lcore <= 0) GO TO 9008
 
!     GET THE NAME OF THE HIGHER LEVEL SUBSTRUCTURE.  IF NONE EXISTS
!     THEN RETURN.
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 names(1) = rss(1)
 names(2) = rss(2)
 CALL fndnxl (rss,higher)
 rc = 4
 IF (higher(1) == BLANK) GO TO 6000
 IF (higher(1) == rss(1) .AND. higher(2) == rss(2)) GO TO 1000
 
!     CHECK IF THE HIGHER LEVEL SUBSTRUCTURE WAS MODAL REDUCED.
!     IF NOT THEN WE HAVE NOTHING TO DO
 
 names(1) = higher(1)
 names(2) = higher(2)
 rc = 4
 CALL fdsub (higher,idit)
 IF (idit < 0) GO TO 6000
 CALL fmdi (idit,imdi)
 mredu = .false.
 credu = .false.
 IF (andf(buf(imdi+ib),mmask) /= 0) mredu = .true.
 IF (andf(buf(imdi+ib),cmask) /= 0) credu = .true.
 IF (.NOT.mredu) GO TO 1000
 
!     READ THE MODAL GROUP OF THE EQSS TO DETERMINE IF THERE ARE ANY
!     RIGID BODY DOF PRESENT.  ALSO GET THE SIL NUMBER OF THE FIRST
!     MODAL CORDINATE.
 
 item = eqss
 CALL sfetch (higher,eqss,1,rc)
 IF (rc /= 1) GO TO 6000
 CALL suread (z(1),3,n,rc)
 IF (rc /= 1) GO TO 6100
 n = z(3)
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 
 nrigid = 0
 10 CALL suread (z(1),3,n,rc)
 IF (rc /= 1) GO TO 6100
 IF (nrigid == 0) ip = z(2)
 IF (z(1) >= 100) GO TO 20
 nrigid = nrigid + 1
 GO TO 10
 
 20 IF (2*ip > sof3) GO TO 9008
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 
 CALL suread (z(1),2*ip,n,rc)
 IF (rc /= 1) GO TO 6100
 i = 2*(ip-1) + 1
 isil = z(i)
 
!     CALCULATE THE ENERGIES ON THE EXCLUDED MODES
 
 noexcl = .true.
 nrowe  = 0
 IF (credu .OR. rfno < 3 .OR. rfno > 8) GO TO 100
 noexcl = .false.
 CALL rcovem (noexcl,nrowe)
 
!     CALCULATE THE ENERGIES ON THE INCLUDED MODE AND THE TOTAL
!     ENERGIES ON EACH VECTOR
 
 100 CALL rcovim (higher)
 IF (iopt < 0) GO TO 9200
 mcb(1) = scr6
 CALL rdtrl(mcb)
 ncol  = mcb(2)
 nrowi = mcb(3)
 nmodei= nrowi - isil + 1
 
!     READ THE MODE DATA FROM LAMS AND SAVE THE MODE NUMBER AND
!     THE FREQUENCY FOR EACH MODE.
 
 names(1) = rss(1)
 names(2) = rss(2)
 item = lams
 CALL sfetch (rss,lams,1,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 imode = 1
 IF (nrigid == 0) GO TO 210
 n = 3*nrigid
 imode = imode + n
 DO  i = 1,n,3
   z(i  ) = 0
   z(i+1) = 0
   z(i+2) = (i-1)/3 + 1
 END DO
 210 nmode = 3*(nmodei-1+nrowe)
 IF (nmode > lcore) GO TO 9008
 
 DO  i = imode,nmode,3
   CALL suread (z(i),7,n,rc)
   IF (rc /= 1) GO TO 6100
   z(i+1) = z(i+4)
   z(i+2) = 0
 END DO
 
!     READ THE LAST GROUP OF LAMS AND GENERATE GRID NUMBERS FOR THE
!     INCLUDED MODES.
 
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 iinc = 100
 DO  i = imode,nmode,3
   CALL suread (icode,1,n,rc)
   IF (rc /= 1) GO TO 6100
   IF (icode > 1) CYCLE
   iinc = iinc + 1
   z(i+2) = iinc
 END DO
 
!     POSITION THE SOLN ITEM TO THE FREQUENCY OR TIME DATA
 
 item = soln
 CALL sfetch (rss,soln,1,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6200
 
!     ALLOCATE INCORE ARRAYS FOR THE ENERGY VECTORS
 
 ivec1 = nmode + 1
 ivec2 = ivec1 + nmodei
 ivec3 = ivec2 + nmodei
 ivec4 = ivec3 + nrowe
 isets = ivec4 + nrowe
 IF (isets > lcore) GO TO 9008
 
!     READ CASESS AND GET THE TITLE AND ANY SET INFORMATION
 
 FILE = casess
 CALL gopen (casess,z(buf1),rdrew)
 450 CALL fread (casess,z(ivec1),2,1)
 IF (z(ivec1) /= casecc(1) .OR. z(ivec1+1) /= casecc(2)) GO TO 450
 
 CALL fread (casess,0,-38,0)
 CALL fread (casess,ititle(1),96,0)
 
 IF (energy <= 0) GO TO 485
 CALL fread (casess,0,-31,0)
 CALL fread (casess,lcc,1,0)
 lskip = 167 - lcc
 CALL fread (casess,0,lskip,0)
 CALL READ (*9002,*480,casess,lseq,1,0,i)
 IF (lseq > 0) CALL fread (casess,0,lseq,0)
 
 460 CALL READ (*9002,*480,casess,iset,1,0,i)
 CALL fread (casess,lset,1,0)
 IF (iset == energy) GO TO 470
 CALL fread (casess,0,-lset,0)
 GO TO 460
 470 IF (isets+lset > lcore) GO TO 9008
 CALL fread (casess,z(isets),lset,0)
 GO TO 485
 
 480 WRITE (nout,63650) uwm,energy
 energy = -1
 
 485 CALL CLOSE (casess,rew)
 
!     LOOP OVER EACH COLUMN AND PRINT THE KINETIC AND POTENTIAL
!     ENERGIES FOR EACH MODAL COORDINATE IF REQUESTED
 
 next = 1
 CALL gopen (scr6,z(buf1),rdrew)
 CALL gopen (scr7,z(buf2),rdrew)
 IF (noexcl) GO TO 490
 CALL gopen (scr3,z(buf3),rdrew)
 CALL gopen (scr4,z(buf4),rdrew)
 
 490 itinu = rsp
 incru = 1
 
 DO  icol = 1,ncol
   
!     SET FLAGS FOR NULL COLUMNS
   
   ikflag = 0
   ipflag = 0
   
!     GET THE FREQUENCY OR TIME FOR THIS VECTOR
   
   IF (rfno > 3) GO TO 500
   
!     NORMAL MODES SOLUTION
   
   CALL suread (z(ivec1),7,n,rc)
   IF (rc /= 1) GO TO 6100
   step = rz(ivec1+4)
   GO TO 505
   
!     DYNAMICS SOLUTION
   
   500 CALL suread (step,1,n,rc)
   IF (rc /= 1) GO TO 6100
   
!     SEE IF THIS COLUMN IS REQUESTED
   
   505 IF (energy <= 0) GO TO 510
   CALL setfnd (*790,z(isets),lset,icol,next)
   
   510 IF (step < range(1) .OR. step > range(2)) GO TO 790
   
!     UNPACK THE KINETIC AND POTENTIAL ENERGIES ON INCLUDED MODES
   
   iru = isil
   nru = nrowi
   CALL unpack (*520,scr6,rz(ivec1))
   GO TO 540
   520 DO  i = 1,nmodei
     rz(ivec1+i-1) = 0.0
   END DO
   rz(ivec1+nmodei-1) = 1.0
   ikflag = 1
   
   540 CALL unpack (*550,scr7,rz(ivec2))
   GO TO 570
   550 DO  i = 1,nmodei
     rz(ivec2+i-1) = 0.0
   END DO
   rz(ivec2+nmodei-1) = 1.0
   ipflag = 1
   
!     UNPACK THE KINETIC AND POTENTIAL ENERGIES ON EXLUDED MODES
   
   570 IF (noexcl) GO TO 580
   iru = 1
   nru = nrowe
   CALL unpack (*580,scr3,rz(ivec3))
   GO TO 600
   580 DO  i = 1,nrowe
     rz(ivec3+i-1) = 0.0
   END DO
   
   600 IF (noexcl) GO TO 610
   CALL unpack (*610,scr4,rz(ivec4))
   GO TO 630
   610 DO  i = 1,nrowe
     rz(ivec4+i-1) = 0.0
   END DO
   
!     INITILIZE FOR THE OUTPUT
   
   630 nlines = nlpp
   
!     GET TOTAL ENERGIES
   
   tkeng = rz(ivec1+nmodei-1)
   tpeng = rz(ivec2+nmodei-1)
   perkt = 1.0
   perpt = 1.0
   
!     LOOP OVER EACH MODAL COORDINATE
   
   iinc = 0
   iexc = 0
   
   DO  i = 1,nmode,3
     
     mode = z(i)
     freq = rz(i+1)
     grid = z(i+2)
     
!     GET ENERGIES FORM THE PROPER VECTOR
     
     IF (noexcl) GO TO 650
     IF (grid == 0) GO TO 660
     650  keng = rz(ivec1+iinc)
     peng = rz(ivec2+iinc)
     iinc = iinc + 1
     TYPE(1) = timode(1)
     TYPE(2) = timode(2)
     
     IF (mode /= 0) GO TO 670
     TYPE(1) = trigid(1)
     TYPE(2) = trigid(2)
     GO TO 670
     
     660 keng = rz(ivec3+iexc)
     peng = rz(ivec4+iexc)
     iexc = iexc + 1
     TYPE(1) = temode(1)
     TYPE(2) = temode(2)
     
!     CALCULATE THE ENERGY PERCENTAGES
     
     670 perk = keng/tkeng
     IF (perk >= 100.0) perk = 99.9999
     perp = peng/tpeng
     IF (perp >= 100.0) perp = 99.9999
     IF (grid /= 0) GO TO 680
     perkt = perkt + perk
     perpt = perpt + perp
     
!     PRINT A LINE OF OUTPUT
     
     680 nlines = nlines + 1
     IF (nlines <= nlpp) GO TO 690
     CALL page1
     WRITE (nout,5000) rss
     IF (rfno == 9) WRITE (nout,5010) step
     IF (rfno /= 9) WRITE (nout,5020) step
     WRITE (nout,5100)
     nlines = 0
     
     690 WRITE (nout,5200) grid,TYPE,mode,freq,keng,perk,peng,perp
     
   END DO
   
!     PRINT THE TOTAL KINETIC AND POTENTIAL ENERGIES FOR THIS COLUMN
   
   IF (perkt >= 100.0) perkt = 99.9999
   IF (perpt >= 100.0) perpt = 99.9999
   IF (ikflag == 0) GO TO 710
   tkeng = 0.0
   perkt = 0.0
   710 IF (ipflag == 0) GO TO 720
   tpeng = 0.0
   perpt = 0.0
   720 WRITE (nout,5300) tkeng,perkt,tpeng,perpt
   CYCLE
   
!     THIS VECTOR IS NOT TO BE PRINTED SO SKIP IT
   
   790 CALL fwdrec (*9002,scr6)
   CALL fwdrec (*9002,scr7)
   IF (noexcl) CYCLE
   CALL fwdrec (*9002,scr3)
   CALL fwdrec (*9002,scr4)
   
 END DO
 
!     CLOSE FILES
 
 CALL CLOSE (scr6,rew)
 CALL CLOSE (scr7,rew)
 IF (noexcl) GO TO 1000
 CALL CLOSE (scr3,rew)
 CALL CLOSE (scr4,rew)
 
!     NORMAL RETURN
 
 1000 CALL sofcls
 RETURN
 
!     ERRORS
 
 6000 CALL smsg (rc-2,item,names)
 GO TO 9200
 6100 CALL smsg (rc+4,item,names)
 GO TO 9200
 6200 CALL smsg (7,item,names)
 GO TO 9200
 9002 n = 2
 GO TO 9100
 9008 n = 8
 9100 CALL mesage (n,FILE,NAME)
 9200 CALL sofcls
 WRITE (nout,63710) uwm,rss
 RETURN
 
!     FORMAT STATEMENTS
 
 5000 FORMAT (//39X,43HMODAL coordinate energies for substructure ,2A4)
 5010 FORMAT (//12X,7HTIME = ,1P,e13.6)
 5020 FORMAT (//12X,12HFREQUENCY = ,1P,e13.6)
 5100 FORMAT (//12X,4HGRID,6X,4HTYPE,6X,4HMODE,7X,9HFREQUENCY,10X,  &
     7HKINETIC,8X,8HKE/total,6X,9HPOTENTIAL,7X,8HPE/total ,/)
 5200 FORMAT (1H ,8X,i8,5X,2A4,2X,i5,5X,1P,e13.6,2(5X,1P,e13.6,5X, 0P,f7.4))
 5300 FORMAT (1H ,55X,2(4X,14H--------------,4X,8H--------), /12X,  &
     28HTOTAL energy for this vector,15X,2(5X,1P,e13.6,5X, 0P,f7.4))
 63710 FORMAT (a25,' 6371, MODAL REDUCTION ENERGY CALCULATIONS FOR ',  &
     'SUBSTRUCTURE ',2A4,' ABORTED.')
 63650 FORMAT (a25,' 6365, REQUESTED OUTPUT SET ID',i6,' IS NOT ',  &
     'DECLARED IN CASE CONTROL. ALL OUTPUT WILL BE PRODUCED')
END SUBROUTINE rcove
