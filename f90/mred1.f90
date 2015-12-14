SUBROUTINE mred1
     
!     THIS SUBROUTINE IS THE MRED1 MODULE WHICH INITIATES THE MODAL
!     SYNTHESIS CALCULATIONS.
 
!     DMAP CALLING SEQUENCE
!     MRED1    CASECC,GEOM4,DYNAMICS/USETX,EEDX,EQST,DMR/*NAMEA*/
!              S,N,DRY/STEP/S,N,NOUS/S,N,SKIPM/S,N,GPARM/TYPE $
 
!     INPUT DATA
!     GINO   - CASECC - CASE CONTROL
!              GEOM4  - BDYC DATA
!                     - BDYS DATA
!                     - BDYS1 DATA
!              DYNAMICS - EIGR DATA
!     SOF    - EQSS   - SUBSTRUCTURE EQUIVALENCE TABLE
!              BGSS   - BASIC GRID POINT IDENTIFICATION TABLE
!              CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRICES DATA
 
!     OUTPUT DATA
!     GINO   - USETX  - S,R,B DEGREES OF FREEDOM
!              EEDX   - EIGR DATA
!              EQST   - TEMPORARY EQSS
!              DMR    - RIGID BODY MATRIX
 
!     PARAMETERS
!     INPUT  - NAMEA  - INPUT SUBSTRUCTURE NAME (BCD)
!              DRY    - OPERATION MODE (INTEGER)
!              STEP   - CONTROL DATA CASECC RECORD (INTEGER)
!              TYPE   - REAL OR COMPLEX (BCD)
!     OUTPUT - DRY    - MODULE OPERATION FLAG (INTEGER)
!              NOUS   - FIXED POINTS FLAG (INTEGER)
!                     = +1 IF FIXED POINTS DEFINED
!                     = -1 IF NO FIXED POINTS DEFINED
!              SKIPM  - MODES FLAG (INTEGER)
!                     =  0 IF MODES NOT PRESENT
!                     = -1 IF MODES PRESENT
!     OTHERS - GBUF   - GINO BUFFERS
!              SBUF   - SOF  BUFFERS
!              KORLEN - CORE LENGTH
!              NEWNAM - NEW SUBSTRUCTURE NAME
!              BNDSET - BOUNDARY SET IDENTIFICATION NUMBER
!              FIXSET - FIXED SET IDENTIFICATION NUMBER
!              IEIG   - EIGENVALUE SET IDENTIFICATION NUMBER
!              IO     - OUTPUT FLAGS
!              RGRID  - FREEBODY MODES FLAGS
!              RNAME  - FREEBODY SUBSTRUCTURE NAME
!              IRSAVE - RSAVE FLAG
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              NCSUBS - NUMBER OF CONTRIBUTING SUBSTRUCTURES
!              NAMEBS - BEGINNING ADDRESS OF CONTRIBUTING SUBSTRUCTURE
!                       NAMES
!              EQSIND - BEGINNING ADDRESS OF EQSS GROUP ADDRESSES
!              NSLBGN - BEGINNING ADDRESS OF SIL DATA
!              NSIL   - NUMBER OF SIL GROUPS
!              BDYCS  - BEGINNING ADDRESS OF BDYC DATA
!              NBDYCC - NUMBER OF BDYC DATA GROUPS
!              USETL  - LENGTH OF USET ARRAY
!              USTLOC - BEGINNING ADDRESS OF USET ARRAY
!              RGRIDX - FREEBODY MODE RELATIVE X COORDINATE
!              RGRIDY - FREEBODY MODE RELATIVE Y COORDINATE
!              RGRIDZ - FREEBODY MODE RELATIVE Z COORDINATE
!              USRMOD - USERMODE  OPTION FLAG
!              BOUNDS - OLDBOUNDS OPTION FLAG
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,andf,orf
 LOGICAL :: usrmod,bounds,ponly,errors
 REAL :: rz(1),range(2),gprm
 DIMENSION       modnam(2),mtrlra(7),mtrlrb(7),mtrlrc(7),mtrlrd(7),  &
     nmonic(16),cctype(2),mtrlre(7),itmnam(2), lstbit(32),errnam(6),letrs(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /BLANK / oldnam(2),dry,step,nous,skipm,TYPE(2),gprm,gbuf1,  &
     gbuf2,sbuf1,sbuf2,sbuf3,korlen,newnam(2),bndset,  &
     fixset,ieig,io,rgrid(2),rname(2),irsave,korbgn,  &
     ncsubs,namebs,eqsind,nslbgn,nsil,bdycs,nbdycc,  &
     usetl,ustloc,rgridx,rgridy,rgridz,usrmod,bounds, ponly
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,iprntr
 EQUIVALENCE     (z(1),rz(1))
 DATA    nmonic/ 4HNAMB,4HBOUN,4HFIXE,4HMETH,4HCMET,4HOUTP,4HRGRI,  &
     4HOLDM,4HOLDB,4HRSAV,4HRNAM,4HRANG,4HNMAX,4HUSER, 4HNAMA,4HGPAR/
 DATA    casecc/ 101          /
 DATA    modnam/ 4HMRED,4H1   /
 DATA    errnam/ 4HLAMS,4HPHIS,4HPHIL,4HGIMS,4HLMTX,4HUPRT/
 DATA    iblank, yes,no,all   /4H    ,4HYES ,4HNO  ,4HALL /
 DATA    letrs / 1HM,1HC      /
 DATA    cctype/ -1,-2        /
 DATA    cred  , nhlods,nhloap,nheqss /4HCRED,4HLODS,4HLOAP,4HEQSS/
 
!     COMPUTE OPEN CORE AND DEFINE GINO, SOF BUFFERS
 
 nozwds = korsz(z(1))
 gbuf1  = nozwds- sysbuf - 2
 gbuf2  = gbuf1 - sysbuf
 sbuf1  = gbuf2 - sysbuf
 sbuf2  = sbuf1 - sysbuf - 1
 sbuf3  = sbuf2 - sysbuf
 korlen = sbuf3 - 1
 korbgn = 1
 IF (korlen <= korbgn) GO TO 430
 
!     INITIALIZE SOF
 
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 
!     INITIALIZE CASE CONTROL PARAMETERS
 
 DO  i = 1, 2
   rgrid(i) = -1
   newnam(i)= iblank
   rname(i) = iblank
 END DO
 bndset = 0
 fixset = 0
 ieig   = 0
 noieig = yes
 io     = 0
 skipm  = 0
 modes  = no
 bounds = .false.
 ponly  = .false.
 ibound = no
 irsave = no
 nous   = 1
 ifree  = no
 nmax   = 2147483647
 imax   = all
 imode  = no
 usrmod = .false.
 iuserm = 1
 module = 1
 gprm   = 0.0
 ibf    = 0
 nrange = 0
 irange = all
 range(1) =-1.0E+35
 range(2) = 1.0E+35
 
!     PROCESS CASE CONTROL
 
 ifile = casecc
 CALL OPEN (*400,casecc,z(gbuf2),0)
 IF (step == 0.0) THEN
   GO TO    40
 END IF
 20 DO  i = 1,step
   CALL fwdrec (*420,casecc)
 END DO
 
!     READ CASECC AND EXTRACT DATA
 
 40 CALL READ (*410,*420,casecc,z(korbgn),2,0,noread)
 IF (z(korbgn) == cred) module = 2
 nowdsc = z(korbgn+1)
 DO  i = 1,nowdsc,3
   CALL READ (*410,*420,casecc,z(korbgn),3,0,noread)
   
!     TEST CASE CONTROL MNEMONICS
   
   DO  j = 1,16
     IF (z(korbgn) == nmonic(j)) GO TO 60
   END DO
   CYCLE
   
!     SELECT DATA TO EXTRACT
   
   60 SELECT CASE ( j )
     CASE (    1)
       GO TO  70
     CASE (    2)
       GO TO  90
     CASE (    3)
       GO TO 100
     CASE (    4)
       GO TO 110
     CASE (    5)
       GO TO 110
     CASE (    6)
       GO TO 120
     CASE (    7)
       GO TO 130
     CASE (    8)
       GO TO 140
     CASE (    9)
       GO TO 150
     CASE (   10)
       GO TO 160
     CASE (   11)
       GO TO 170
     CASE (   12)
       GO TO 102
     CASE (   13)
       GO TO 115
     CASE (   14)
       GO TO 125
     CASE (   15)
       GO TO 132
     CASE (   16)
       GO TO 155
   END SELECT
   
!     EXTRACT NEW SUBSTRUCTURE NAME
   
   70 DO  k = 1,2
     newnam(k) = z(korbgn+k)
   END DO
   CYCLE
   
!     EXTRACT BOUNDARY SET
   
   90 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   bndset = z(korbgn+2)
   ibf = ibf + 2
   CYCLE
   
!     EXTRACT FIXED SET
   
   100 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   fixset = z(korbgn+2)
   ibf = ibf + 1
   CYCLE
   
!     EXTRACT FREQUENCY RANGE
   
   102 IF (z(korbgn+1) /= cctype(2)) GO TO 185
   irange = iblank
   IF (nrange == 1) GO TO 104
   nrange = 1
   range(1) = rz(korbgn+2)
   CYCLE
   104 range(2) = rz(korbgn+2)
   CYCLE
   
!     EXTRACT EIGENVALUE METHOD
   
   110 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   ieig = z(korbgn+2)
   noieig = no
   CYCLE
   
!     EXTRACT MAXIMUM NUMBER OF FREQUENCIES
   
   115 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   IF (z(korbgn+2) == 0) CYCLE
   nmax = z(korbgn+2)
   imax = iblank
   CYCLE
   
!     EXTRACT OUTPUT FLAGS
   
   120 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   io = orf(io,z(korbgn+2))
   CYCLE
   
!     EXTRACT USERMODE FLAG
   
   125 IF (z(korbgn+1) /= cctype(1)) GO TO 185
   imode  = yes
   skipm  = -1
   usrmod = .true.
   IF (z(korbgn+2) == 2) iuserm = 2
   CYCLE
   
!     EXTRACT RIGID BODY GRID POINT ID
   
   130  rgrid(1) = z(korbgn+2)
   IF (z(korbgn+1) /= cctype(1)) rgrid(1) = 0
   ifree = yes
   CYCLE
   
!     EXTRACT OLD SUBSTRUCTURE NAME
   
   132 DO  k = 1,2
     oldnam(k) = z(korbgn+k)
   END DO
   CYCLE
   
!     SET OLDMODES FLAG
   
   140 IF ((z(korbgn+1) == cctype(1)) .OR. (z(korbgn+1) == cctype(2)))  &
       GO TO 185
   IF (z(korbgn+1) /= yes) CYCLE
   skipm = -1
   modes = yes
   CYCLE
   
!     SET OLDBOUND FLAG
   
   150 IF ((z(korbgn+1) == cctype(1)) .OR. (z(korbgn+1) == cctype(2)))  &
       GO TO 185
   IF (z(korbgn+1) /= yes) CYCLE
   bounds = .true.
   ibound = yes
   CYCLE
   
!     EXTRACT GPARAM PARAMETER
   
   155 IF (z(korbgn+1) /= cctype(2)) GO TO 185
   gprm = rz(korbgn+2)
   CYCLE
   
!     SET RSAVE FLAG
   
   160 IF (z(korbgn+1) == no) CYCLE
   irsave = yes
   CYCLE
   
!     EXTRACT RIGID BODY SUBSTRUCTURE NAME
   
   170 DO  k = 1,2
     rname(k) = z(korbgn+k)
   END DO
   IF (rgrid(1) < 0) rgrid(1) = 0
   ifree = yes
   CYCLE
   
!     CASECC COMMAND ERROR
   
   185 WRITE (iprntr,916) uwm,letrs(module),nmonic(j)
 END DO
 CALL CLOSE (casecc,1)
 
!     TEST MODULE OPERATION FLAG
 
 IF (dry < 0.0) THEN
   GO TO   192
 ELSE IF (dry == 0.0) THEN
   GO TO   194
 ELSE
   GO TO   196
 END IF
 192 IF (dry == -2) GO TO 198
 WRITE (iprntr,909) uim
 dry = -2
 GO TO 198
 194 skipm = -1
 itest = 0
 CALL fdsub (newnam,itest)
 IF (itest /= -1) GO TO 510
 WRITE (iprntr,922) ufm,letrs(module),newnam
 GO TO 500
 196 itest = 0
 CALL fdsub (newnam,itest)
 IF (itest == -1) GO TO 198
 IF (bounds .OR. (skipm == -1)) GO TO 198
 CALL sfetch (newnam,nhlods,3,itest)
 IF (itest == 3) GO TO 197
 CALL sfetch (newnam,nhloap,3,itest)
 IF (itest == 3) GO TO 197
 itmnam(1) = newnam(1)
 itmnam(2) = newnam(2)
 GO TO 450
 
!     LOADS ONLY PROCESSING
 
 197 ponly = .true.
 
!     TEST OUTPUT OPTION
 
 198 IF (andf(io,1) == 0) GO TO 200
 CALL page1
 WRITE (iprntr,900) oldnam,newnam
 IF (ibf == 0) WRITE (iprntr,918)
 IF (ibf == 1) WRITE (iprntr,919) fixset
 IF (ibf == 2) WRITE (iprntr,920) bndset
 IF (ibf == 3) WRITE (iprntr,921) bndset,fixset
 IF (rgrid(1) == -1) WRITE (iprntr,906) rname
 IF (rgrid(1) /= -1) WRITE (iprntr,907) rgrid(1),rname
 IF (noieig == no) WRITE (iprntr,908) ibound,modes,ifree,imode, irsave,ieig
 IF (noieig /= no) WRITE (iprntr,908) ibound,modes,ifree,imode, irsave
 IF (imax   == all) WRITE (iprntr,910) imax,gprm
 IF (imax   /= all) WRITE (iprntr,911) nmax,gprm
 IF (irange == all) WRITE (iprntr,912) oldnam,irange
 IF (irange /= all) WRITE (iprntr,913) oldnam,range(1)
 
!     CHECK FOR OLDMODES, OLDBOUND ERRORS
 
 200 errors = .false.
 IF (ponly) GO TO 290
 CALL sfetch (oldnam,errnam(1),3,itest)
 CALL softrl (oldnam,errnam(2),mtrlra)
 CALL softrl (oldnam,errnam(4),mtrlrb)
 CALL softrl (oldnam,errnam(5),mtrlrc)
 CALL softrl (oldnam,errnam(3),mtrlrd)
 CALL softrl (oldnam,errnam(6),mtrlre)
 iflag = 1
 IF (usrmod) GO TO 290
 IF (skipm < 0.0) THEN
   GO TO   210
 ELSE
   GO TO   230
 END IF
 
!     OLDMODES SET - PHIS AND LAMS MUST BE ON SOF
 
 210 IF (itest > 3) GO TO 360
 220 iflag = 2
 IF (mtrlra(1) > 2) GO TO 360
 GO TO 260
 
!     OLDMODES NOT SET - PHIS, PHIL AND LAMS MUST BE DELETED
 
 230 IF (itest < 3) GO TO 370
 240 iflag = 2
 IF (mtrlra(1) < 3) GO TO 370
 250 iflag = 3
 IF (mtrlrd(1) < 3) GO TO 370
 
!     OLDBOUND SET - GIMS AND UPRT MUST BE ON SOF
 
 260 iflag = 4
 IF (.NOT. bounds) GO TO 270
 IF (mtrlrb(1) > 2) GO TO 380
 265 iflag = 6
 IF (mtrlre(1) > 2) GO TO 380
 GO TO 290
 
!     OLDBOUND NOT SET - GIMS AND LMTX MUST BE DELETED
 
 270 IF (mtrlrb(1) < 3) GO TO 390
 280 iflag = 5
 IF (mtrlrc(1) < 3) GO TO 390
 
!     TEST FOR ERRORS
 
 290 IF (errors) GO TO 500
 IF (iuserm == 2) WRITE (iprntr,917) uim
 
!     READ EQSS GROUP 0 DATA AND TEST OPEN CORE LENGTH
 
 itmnam(2) = oldnam(2)
 CALL sfetch (oldnam,nheqss,1,itest)
 IF (itest == 3) GO TO 460
 IF (itest == 4) GO TO 470
 IF (itest == 5) GO TO 480
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 IF (korbgn+nwdsrd >= sbuf3) GO TO 430
 
!     COMPRESS BASIC SUBSTRUCTURE NAMES AND TEST OPEN CORE LENGTH
 
 ncsubs = z(korbgn+2)
 namebs = korbgn
 i = 2*((nwdsrd - 4)/2)
 k = 4
 DO  j = 1,i,2
   z(korbgn+j-1) = z(korbgn+k  )
   z(korbgn+j  ) = z(korbgn+k+1)
   IF (rgrid(1) < 0) GO TO 300
   IF (rname(1) /= iblank) GO TO 298
   rname(1) = z(korbgn+j-1)
   rname(2) = z(korbgn+j  )
   298 CONTINUE
   IF ((z(korbgn+j-1) /= rname(1)) .OR. (z(korbgn+j) /= rname(2))) GO TO 300
   rgrid(2) = (j+1)/2
   300 k = k + 2
 END DO
 eqsind = korbgn + 2*ncsubs
 IF (eqsind >= sbuf3) GO TO 430
 
!     TEST OUTPUT OPTION
 
 IF (andf(io,1) == 0) GO TO 310
 IF (irange  /=  all) GO TO 302
 i = 2*ncsubs
 WRITE (iprntr,901) (z(korbgn+j-1),z(korbgn+j),j=1,i,2)
 GO TO 310
 302 IF (ncsubs >= 5) GO TO 306
 i = 1 + 2*ncsubs
 DO  j = i,10
   z(korbgn+j-1) = iblank
 END DO
 306 k = 10
 WRITE (iprntr,914) (z(korbgn+j-1),z(korbgn+j),j=1,k,2),range(2)
 IF (ncsubs <= 5) GO TO 310
 k = k + 1
 i = 2*ncsubs
 WRITE (iprntr,901) (z(korbgn+j-1),z(korbgn+j),j=k,i,2)
 
!     READ EQSS GROUPS TO END-OF-ITEM
 
 310 korbgn = eqsind + 2*ncsubs
 DO  i = 1, ncsubs
   IF (korbgn >= sbuf3) GO TO 430
   CALL suread (z(korbgn),-1,nwdsrd,itest)
   j = 2*(i - 1)
   z(eqsind+j  ) = korbgn
   z(eqsind+j+1) = nwdsrd
   korbgn = korbgn + nwdsrd
 END DO
 nslbgn = korbgn
 CALL suread (z(korbgn),-2,nwdsrd,itest)
 nsil = nwdsrd/2
 
!     TEST OUTPUT OPTION
 
 IF (andf(rshift(io,3),1) == 0) GO TO 350
 DO  i = 1,ncsubs
   j = 2*(i-1)
   CALL cmiwrt (1,oldnam,z(namebs+j),z(eqsind+j),z(eqsind+j+1),rz,z)
 END DO
 isil = 2*nsil
 CALL cmiwrt (8,oldnam,oldnam,nslbgn,isil,rz,z)
 
!     DETERMINE USET LENGTH
 
 350 korbgn = nslbgn + nwdsrd
 ustloc = korbgn
 icode  = z(korbgn-2)
 CALL decode (icode,lstbit,nwdsd)
 usetl = (z(korbgn-3) + nwdsd) - 1
 
!     PROCESS FIXED SET
 
 CALL mred1a (1)
 CALL mred1b (1)
 
!     PROCESS BOUNDARY SET
 
 CALL mred1a (2)
 CALL mred1b (2)
 
!     CONVERT EQSS DATA TO UB DATA
 
 IF (ponly) GO TO 510
 CALL mred1c
 
!     PROCESS EIGENVALUE DATA
 
 IF (skipm == -1) GO TO 355
 CALL mred1d
 
!     PROCESS FREE BODY MODES
 
 355 CALL mred1e
 GO TO 510
 
!     PHIS, LAMS DO NOT EXIST
 
 360 WRITE (iprntr,902) ufm,errnam(iflag),oldnam
 errors = .true.
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 220
   CASE (    2)
     GO TO 260
 END SELECT
 
!     PHIS, PHIR, LAMS NOT DELETED
 
 370 WRITE (iprntr,903) ufm,errnam(iflag),oldnam
 errors = .true.
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 240
   CASE (    2)
     GO TO 250
   CASE (    3)
     GO TO 260
 END SELECT
 
!     GIMS, UPRT DOES NOT EXIST
 
 380 WRITE (iprntr,904) ufm,errnam(iflag),oldnam
 errors = .true.
 IF (iflag - 5 > 0) THEN
   GO TO   290
 ELSE
   GO TO   265
 END IF
 
!     GIMS, LMTX NOT DELETED
 
 390 WRITE (iprntr,905) ufm,errnam(iflag),oldnam
 errors = .true.
 iflag = iflag - 3
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 280
   CASE (    2)
     GO TO 290
 END SELECT
 
!     PROCESS SYSTEM FATAL ERRORS
 
 400 imsg = -1
 GO TO 440
 410 imsg = -2
 GO TO 440
 420 imsg = -3
 GO TO 440
 430 imsg = -8
 ifile = 0
 440 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 450 imsg = -4
 GO TO 490
 460 imsg = -1
 GO TO 490
 470 imsg = -2
 GO TO 490
 480 imsg = -3
 490 CALL smsg (imsg,nheqss,itmnam)
 RETURN
 
 500 CALL sofcls
 dry = -2
 RETURN
 
!     CLOSE ANY OPEN FILES
 
 510 CALL sofcls
 IF (dry == -2) WRITE (iprntr,915) letrs(module)
 IF (ponly) skipm = -1
 RETURN
 
 900 FORMAT (//38X,46HS u m m a r y    o f    c u r r e n t    p r o,  &
     8H b l e m,//13X,38HNAME of pseudostructure TO be reduced ,  &
     4(2H. ),2A4,6X,40HNAME given TO resultant pseudostructure , 2A4)
 901 FORMAT (16X,2A4,2X,2A4,2X,2A4,2X,2A4,2X,2A4)
 902 FORMAT (a23,' 6617, OLDMODES SET AND REQUESTED SOF ITEM DOES NOT',  &
     ' EXIST.  ITEM ',a4,', SUBSTRUCTURE ',2A4,1H.)
 903 FORMAT (a23,' 6618, OLDMODES NOT SET AND REQUESTED SOF ITEM MUST',  &
     ' BE DELETED.  ITEM ',a4,', SUBSTRUCTURE ',2A4,1H.)
 904 FORMAT (a23,' 6619, OLDBOUND SET AND REQUESTED SOF ITEM DOES NOT',  &
     ' EXIST.  ITEM ',a4,', SUBSTRUCTURE ',2A4,1H.)
 905 FORMAT (a23,' 6620, OLDBOUND NOT SET AND REQUESTED SOF ITEM MUST',  &
     ' BE DELETED.  ITEM ',a4,', SUBSTRUCTURE ',2A4,1H.)
 906 FORMAT (13X,'RIGID BODY GRID POINT IDENTIFICATION NUMBER .',14X,  &
     'RIGID BODY SUBSTRUCTURE NAME ',5(2H. ),2A4)
 907 FORMAT (13X,46HRIGID body grid point identification NUMBER . ,i8,  &
     6X,30HRIGID body substructure NAME  ,5(2H. ),2A4)
 908 FORMAT (13X,18HOLDBOUND flag set ,14(2H. ),a4,10X,12HOLDMODES fla,  &
     6HG set ,11(2H. ),a4,/13X,29HFREE body modes TO be calcula,  &
     5HTED  ,6(2H. ),a4,10X,20HUSER modes flag set ,10(2H. ),a4,  &
     /13X,24HSAVE reduction products ,11(2H. ),a4,10X,7HEIGENVA,  &
     23HLUE extraction method  ,5(2H. ),i8)
 909 FORMAT (a29,' 6630, FOR DRY OPTION IN MODAL REDUCE, INPUT DATA ',  &
     'WILL BE CHECKED', /36X,'BUT NO SOF TABLE ITEMS WILL BE ', 'CREATED.')
 910 FORMAT (13X,42HMAXIMUM NUMBER of frequencies TO be used  ,2(2H. ),  &
     a4,10X,14HGPARAM value  ,13(2H. ),1P,e12.6)
 911 FORMAT (13X,42HMAXIMUM NUMBER of frequencies TO be used  ,2(2H. ),  &
     i8,6X,14HGPARAM value  ,13(2H. ),1P,e12.6)
 912 FORMAT (13X,46HNAMES of component substructures contained in ,2A4,  &
     6X,32HRANGE of frequencies TO be used ,4(2H. ),a4)
 913 FORMAT (13X,46HNAMES of component substructures contained in ,2A4,  &
     6X,32HRANGE of frequencies TO be used ,4(2H. ),1P,e12.6)
 914 FORMAT (16X,5(2A4,2X),47X,1P,e12.6)
 915 FORMAT (10H0  module ,a1,36HREDUCE terminating due TO above erro, 3HRS.)
 916 FORMAT (a25,' 6367, ILLEGAL FORMAT ON THE ',a1,'REDUCE OUTPUT ',  &
     'COMMAND ',a4,'.  COMMAND IGNORED.')
 917 FORMAT (a29,' 6636, NMAX AND RANGE SUB COMMANDS ARE IGNORED ',  &
     'UNDER USERMODES = TYPE 2.')
 918 FORMAT (13X,36HBOUNDARY set identification NUMBER  ,5(2H. ),14X,  &
     32HFIXED set identification NUMBER ,4(2H. ))
 919 FORMAT (13X,36HBOUNDARY set identification NUMBER  ,5(2H. ),14X,  &
     32HFIXED set identification NUMBER ,4(2H. ),i8)
 920 FORMAT (13X,36HBOUNDARY set identification NUMBER  ,5(2H. ),i8,6X,  &
     32HFIXED set identification NUMBER ,4(2H. ))
 921 FORMAT (13X,36HBOUNDARY set identification NUMBER  ,5(2H. ),i8,6X,  &
     32HFIXED set identification NUMBER ,4(2H. ),i8)
 922 FORMAT (a23,' 6220, MODULE ',a1,'REDUCE - RUN EQUALS GO AND ',  &
     'SUBSTRUCTURE ',2A4,' DOES NOT EXIST.')
 
END SUBROUTINE mred1
