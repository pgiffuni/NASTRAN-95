SUBROUTINE mma1 ( zi, zr, zd, zc, zdc )
     
!     MMA1 PERFORMS THE MATRIX OPERATION USING METHODS 10 AND 11
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA1 IS DESIGNED AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  UNPACK AS MANY COLUMNS OF "A" INTO MEMORY AS POSSIBLE
!           LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  CALL UNPACK TO READ MATRICES "A" AND "C".
!       5.  FOR METHOD 10, CALL UNPACK TO READ COLUMNS OF MATRIX "B".
!       6.  FOR METHOD 11, CALL MMARC1,2,3,4 TO READ COLUMNS OF MATRIX "B"
!           INTO MEMORY IN COMPACT FORM.
 
 
 
 INTEGER, INTENT(OUT)                     :: zi(2)
 REAL, INTENT(IN OUT)                     :: zr(2)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(2)
 COMPLEX, INTENT(IN OUT)                  :: zc(2)
 REAL, INTENT(IN OUT)                     :: zdc
 INTEGER :: module(3),sysbuf,scrtch
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: ofile      ,filea    ,fileb ,filec , filed
 
 
 
 DOUBLE COMPLEX    zdc(2)
 COMMON / names  / rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 COMMON / TYPE   / iprc(2)    ,nwords(4),irc(4)
 COMMON / mpyadx / filea(7)   ,fileb(7) ,filec(7)  &
     ,                 filed(7)   ,nz       ,t     ,signab,signc ,prec1  &
     ,                 scrtch     ,time
 INCLUDE           'MMACOM.COM'
 COMMON / system / ksystm(152)
 COMMON / unpakx / typeu      ,iurow1   ,iurown, incru
 COMMON / packx  / typei      ,typep    ,iprow1, iprown , incrp
 EQUIVALENCE       (ksystm( 1),sysbuf)  , (ksystm( 2),nout  )
 EQUIVALENCE       (ksystm(58),ksys58)
 EQUIVALENCE       (filea(2)  ,nac   )  , (filea(3)  ,nar   )  &
     ,                 (filea(4)  ,naform)  , (filea(5)  ,natype)  &
     ,                 (filea(6)  ,nanzwd)  , (filea(7)  ,nadens)
 EQUIVALENCE       (fileb(2)  ,nbc   )  , (fileb(3)  ,nbr   )  &
     ,                 (fileb(4)  ,nbform)  , (fileb(5)  ,nbtype)  &
     ,                 (fileb(6)  ,nbnzwd)  , (fileb(7)  ,nbdens)
 EQUIVALENCE       (filec(2)  ,ncc   )  , (filec(3)  ,ncr   )  &
     ,                 (filec(4)  ,ncform)  , (filec(5)  ,nctype)  &
     ,                 (filec(6)  ,ncnzwd)  , (filec(7)  ,ncdens)
 EQUIVALENCE       (filed(2)  ,ndc   )  , (filed(3)  ,ndr   )  &
     ,                 (filed(4)  ,ndform)  , (filed(5)  ,ndtype)  &
     ,                 (filed(6)  ,ndnzwd)  , (filed(7)  ,nddens)
 
 DATA    module / 4HMMA1  , 4H     ,4H    /
 DATA    kzero  / 1H0   /
 DATA    kone   / 1H1   /
DATA    jbegn  / 4HBEGN/ , jend  / 3HEND /

module( 3 ) = jbegn
IF ( nbstor == 1 ) module( 2 ) = kzero
IF ( nbstor == 2 ) module( 2 ) = kone
CALL conmsg ( module, 3, 0 )
incru  = 1
typei  = ndtype
typep  = ndtype
SIGN   = signab
nwdd   = nwords( ndtype )
nwdb   = nwords( nbtype )
irfile = fileb( 1 )

!   OPEN CORE ALLOCATION AS FOLLOWS:
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "B" MATRIX
!     Z( IDX      ) = ARRAY FOR ONE COLUMN OF "D" MATRIX
!     Z( IAX      ) = ARRAY FOR MULTIPLE COLUMNS OF "A" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE

IF ( nbstor == 1 .OR. ksys58 == 10 ) idx = 1 + nwdd*nbr
IF ( nbstor /= 2 .AND. ksys58 /= 11 ) GO TO 90

! REDEFINE IDX AND INSURE A QUAD WORD BOUNDARY FOR COMPLEX DOUBLE

idx = 1 + nwdd*nbr + nbr
itest = MOD( idx, 4 )
IF ( itest == 1 ) GO TO 90
IF ( itest == 0 ) idx = idx + 1
IF ( itest == 2 ) idx = idx + 3
IF ( itest == 3 ) idx = idx + 2
90    CONTINUE
idx2   = ( ( idx+1 ) / 2 ) - 1
idx4   = ( idx+1 ) / 4
iax    = idx  + nwdd*ndr
ibuf1  = nz    - sysbuf
ibuf2  = ibuf1 - sysbuf
ibuf3  = 0
IF ( filec( 1 ) == 0 .OR. signc == 0 ) GO TO 100
ibuf3  = ibuf2 - sysbuf
ibuf4  = ibuf3 - sysbuf
GO TO 200
100   CONTINUE
ibuf4  = ibuf2 - sysbuf
200   CONTINUE
lasmem = ibuf4 - 1
iprow1 = 1
iprown = ndr
incrp  = 1
CALL gopen  ( filea, zr( ibuf1 ), rdrew )
CALL gopen  ( fileb, zr( ibuf2 ), rdrew )

!   DETERMINE HOW MANY COLUMNS OF A CAN BE READ INTO MEMORY IN ONE PASS

220   CONTINUE
iavail = lasmem - iax + 1

!   NCOLPP  -  NUMBER OF COLUMNS OF "A" THAT CAN BE READ IN ONE PASS
!   NPASS   -  NUMBER OF PASSES NEEDED TO READ ENTIRE "A" MATRIX

ncolpp = iavail / ( 2+nwdd*nar )
IF ( ncolpp > nac ) ncolpp = nac
IF ( ncolpp <= 0 ) CALL mesage ( -8, iavail+nwdd*nar, module )
npass  = ( (nac-1) / ncolpp ) + 1
IF ( npass == 1 .OR. ibuf3 /= 0 ) GO TO 250

! MUST ALLOCATE TWO BUFFERS FOR MULTIPLE PASSES

ibuf3  = ibuf2 - sysbuf
ibuf4  = ibuf3 - sysbuf
lasmem = ibuf4 - 1
GO TO 220
250   CONTINUE
DO  m = 1, npass
  ipass = m
  ibrow = ( m-1 ) * ncolpp
  IF ( m == npass ) GO TO 400
  
! MULTIPLE PASSES REQUIRED, DETERMINE PROPER FILE FOR OUTPUT SO THAT
! REQUESTED OUTPUT FILE IS USED ON THE LAST PASS
  
  itest = npass - m
  itest = MOD( itest, 2 )
  IF ( itest /= 0 ) GO TO 350
  ifile = scrtch
  ofile = filed( 1 )
  GO TO 380
  350   ifile = filed( 1 )
  ofile = scrtch
  380   CONTINUE
  IF ( m == 1 ) GO TO 300
  CALL REWIND( fileb )
  CALL skprec( fileb, 1 )
  CALL gopen ( ifile, zr( ibuf3 ), rdrew )
  CALL gopen ( ofile, zr( ibuf4 ), wrtrew)
  GO TO 490
! FIRST PASS, OPEN "C" FILE IF IT EXISTS
  300   CONTINUE
  CALL gopen ( ofile, zr( ibuf4 ), wrtrew)
  310   ifile = filec( 1 )
  IF ( signc == 0 ) ifile = 0
  IF ( ifile == 0 ) GO TO 490
  CALL gopen  ( ifile, zr( ibuf3 ), rdrew )
  GO TO 490
! LAST PASS, CREATE OUTPUT FILE
  400   CONTINUE
  ncolpp = nac - ncolpp*(npass-1)
  ofile = filed( 1 )
  ifile = scrtch
  CALL REWIND( fileb )
  CALL skprec( fileb, 1 )
  CALL gopen ( filed, zr( ibuf4 ), wrtrew)
  filed( 2 ) = 0
  filed( 6 ) = 0
  filed( 7 ) = 0
  IF ( m == 1 ) GO TO 310
  CALL gopen  ( ifile, zr( ibuf3 ), rdrew )
  490   CONTINUE
  indx   = iax
  typeu  = ndtype
  DO  i = 1, ncolpp
    iurow1 = -1
    CALL unpack ( *500, filea, zr( indx+2 ) )
    zi( indx   ) = iurow1
    zi( indx+1 ) = iurown
    indx    = indx + 2 + nwdd*nar
    CYCLE
    500   CONTINUE
! NULL COLUMN READ
    zi( indx   ) = 0
    zi( indx+1 ) = 0
    indx    = indx + 2 + nwdd*nar
  END DO
  IF ( ksys58 == 10 ) GO TO 950
  IF ( ksys58 == 11 ) GO TO 1000
  IF ( nbstor == 2  ) GO TO 1000
! PROCESS ALL OF THE COLUMNS OF "B", ADD "C" DATA ON FIRST PASS
  950   CONTINUE
  IF ( ndtype == 1 ) CALL mma101( zi, zr )
  IF ( ndtype == 2 ) CALL mma102( zi, zd )
  IF ( ndtype == 3 ) CALL mma103( zi, zc )
  IF ( ndtype == 4 ) CALL mma104( zi, zd, zdc )
  GO TO 60000
  1000  CONTINUE
  IF ( ndtype == 1 ) CALL mma111( zi, zr )
  IF ( ndtype == 2 ) CALL mma112( zi, zd )
  IF ( ndtype == 3 ) CALL mma113( zi, zc )
  IF ( ndtype == 4 ) CALL mma114( zi, zd, zdc )
  60000 CONTINUE
  CALL CLOSE ( ifile, clsrew )
  CALL CLOSE ( ofile, clsrew )
END DO
CALL CLOSE ( filea, clsrew )
CALL CLOSE ( fileb, clsrew )
module( 3 ) = jend
CALL conmsg ( module, 3, 0 )
RETURN
END SUBROUTINE mma1

