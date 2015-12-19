SUBROUTINE mma2 ( zi, zr, zd, zc, zdc )
     
!     MMA2 PERFORMS THE MATRIX OPERATION USING METHODS 20 AND 21
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA2 IS DESIGNED AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  UNPACK AS MANY COLUMNS OF "B" INTO MEMORY AS POSSIBLE
!           LEAVING SPACE FOR A COLUMN OF "D" FOR EVERY COLUMN "B" READ.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  CALL UNPACK TO READ MATRICES "B" AND "C".
!       5.  FOR METHOD 20, CALL UNPACK TO READ COLUMNS OF MATRIX "A".
!       6.  FOR METHOD 21, CALL MMARC1,2,3,4 TO READ COLUMNS OF MATRIX "A"
!           INTO MEMORY IN COMPACT FORM.
 
 
 INTEGER, INTENT(OUT)                     :: zi(2)
 REAL, INTENT(OUT)                        :: zr(2)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(2)
 COMPLEX, INTENT(IN OUT)                  :: zc(2)
 INTEGER :: module(3),sysbuf,scrtch
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: filea      ,fileb ,filec , filed
 
 
 
 DOUBLE COMPLEX    zdc(2)
 INCLUDE           'MMACOM.COM'
 COMMON / names  / rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 COMMON / TYPE   / iprc(2)    ,nwords(4),irc(4)
 COMMON / mpyadx / filea(7)   ,fileb(7) ,filec(7)  &
     ,                 filed(7)   ,nz       ,t     ,signab,signc ,prec1  &
     ,                 scrtch     ,time
 COMMON / system / ksystm(152)
 COMMON / unpakx / typeu      ,iurow1   ,iurown, incru
 COMMON / packx  / typei      ,typep    ,iprow1, iprown , incrp
 EQUIVALENCE       (ksystm( 1),sysbuf)  , (ksystm( 2),nout  )  &
     ,                 (ksystm(58),ksys58)
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
 
 DATA    module / 4HMMA2  , 4H     ,4H    /
 DATA    kzero  / 1H0    /
 DATA    kone   / 1H1    /
DATA    jbegn  / 4HBEGN/ , jend  / 3HEND /

IF ( nastor == 1 .OR. ksys58 == 20 ) module( 2 ) = kzero
IF ( nastor == 2 .OR. ksys58 == 21 ) module( 2 ) = kone
module( 3 ) = jbegn
CALL conmsg ( module, 3, 0 )
incru  = 1
typei  = ndtype
typep  = ndtype
nwdd   = nwords( ndtype )
irfile = filea( 1 )

!   OPEN CORE ALLOCATION AS FOLLOWS:
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "A" MATRIX
!     Z( IDX      ) = ARRAY FOR MULTIPLE COLUMNS OF "D" MATRIX
!     Z( IBX      ) = ARRAY FOR MULTIPLE COLUMNS OF "B" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE

idx    = 1 + nwdd*nar
IF ( nastor /= 2 .AND. ksys58 /= 21 ) GO TO 90

! REDEFINE IDX AND INSURE A QUAD WORD BOUNDARY FOR COMPLEX DOUBLE

idx = 1 + nwdd*nar + nar
itest  = MOD ( idx, 4 )
IF ( itest == 1 ) GO TO 90
IF ( itest == 0 ) idx = idx + 1
IF ( itest == 2 ) idx = idx + 3
IF ( itest == 3 ) idx = idx + 2
90    CONTINUE
idx2   = ( ( idx+1 ) / 2 ) - 1
idx4   = ( idx+1 ) / 4
ibuf1  = nz    - sysbuf
ibuf2  = ibuf1 - sysbuf
IF ( filec( 1 ) == 0 ) GO TO 100
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
SIGN   = 1.0
CALL gopen  ( filea, zr( ibuf1 ), rdrew )
CALL gopen  ( fileb, zr( ibuf2 ), rdrew )
IF ( filec( 1 ) /= 0 ) CALL gopen  ( filec, zr( ibuf3 ), rdrew )
CALL gopen  ( filed, zr( ibuf4 ), wrtrew )
filed( 2 ) = 0
filed( 6 ) = 0
filed( 7 ) = 0

!   DETERMINE HOW MANY COLUMNS OF "B" CAN BE READ INTO MEMORY AND HOW
!   MANY COLUMNS OF "D" CAN BE HELD IN MEMORY FOR ONE PASS

iavail = lasmem - idx + 1

!   NCOLPP  -  NUMBER OF COLUMNS OF "B" THAT CAN BE READ IN ONE PASS
!   NPASS   -  NUMBER OF PASSES NEEDED TO READ ENTIRE "B" MATRIX

nwddndr = nwdd * ndr
nwddnbr = nwdd * nbr
ncolpp  = iavail / ( 2+ nwddnbr + nwddndr )
IF ( ncolpp <= 0 ) CALL mesage ( -8, iavail+nwddnbr+nwddndr, module)
IF ( ncolpp > nbc ) ncolpp = nbc
npass   = ( (nbc-1) / ncolpp ) + 1
ibx     = idx + ncolpp*nwddndr
DO  m = 1, npass
  ipass   = m
  IF ( m == npass ) ncolpp = nbc - ( ncolpp*(npass-1) )
  CALL REWIND ( filea )
  CALL skprec ( filea, 1 )
  indxb  = ibx
  indxd  = idx
  typeu  = ndtype * signab
  DO  i = 1, ncolpp
    iurow1 = -1
    CALL unpack ( *500, fileb, zr( indxb+2 ) )
    zi( indxb   ) = iurow1
    zi( indxb+1 ) = iurown
    GO TO 550
    500   CONTINUE
! NULL COLUMN READ ON "B"
    zi( indxb   ) = 0
    zi( indxb+1 ) = 0
    550   CONTINUE
    indxb   = indxb + nwddnbr + 2
  END DO
  IF ( filec( 1 ) == 0 .OR. signc == 0 ) GO TO 800
  typeu   = ndtype * signc
  iurow1 = 1
  iurown = ncr
  DO  i = 1, ncolpp
    CALL unpack ( *650, filec, zr( indxd ) )
    GO TO 680
    
! NULL COLUMN READ ON "C"
    
    650   CONTINUE
    LEN     = indxd + nwddndr - 1
    DO  k = indxd, LEN
      zr( k ) = 0.0
    END DO
    680   CONTINUE
    indxd   = indxd + nwddndr
  END DO
  GO TO 900
  
! "C" MATRIX IS NULL OR "SIGNC" IS ZERO
  
  800   CONTINUE
  LEN = idx + ncolpp*nwddndr - 1
  DO  k = idx, LEN
    zr( k ) = 0.
  END DO
  900   CONTINUE
  
! PROCESS ALL OF THE COLUMNS OF "A"
  
  IF ( ksys58 == 21 ) GO TO 1000
  IF ( ksys58 == 20 ) GO TO 950
  IF ( nastor == 2  ) GO TO 1000
  950   CONTINUE
  IF ( ndtype == 1 ) CALL mma201 ( zi, zr )
  IF ( ndtype == 2 ) CALL mma202 ( zi, zd )
  IF ( ndtype == 3 ) CALL mma203 ( zi, zc )
  IF ( ndtype == 4 ) CALL mma204 ( zi, zd, zdc )
  CYCLE
  1000  CONTINUE
  IF ( ndtype == 1 ) CALL mma211 ( zi, zr )
  IF ( ndtype == 2 ) CALL mma212 ( zi, zd )
  IF ( ndtype == 3 ) CALL mma213 ( zi, zc )
  IF ( ndtype == 4 ) CALL mma214 ( zi, zd, zdc )
END DO
CALL CLOSE ( filea, clsrew )
CALL CLOSE ( fileb, clsrew )
CALL CLOSE ( filec, clsrew )
CALL CLOSE ( filed, clsrew )
module( 3 ) = jend
CALL conmsg ( module, 3, 0 )
RETURN
END SUBROUTINE mma2

