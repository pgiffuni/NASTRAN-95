SUBROUTINE mma4 ( zi, zr, zd, zc, zdc )
     
!     MMA4 PERFORMS THE MATRIX OPERATION USING METHODS 40 AND 41
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA4 IS DESIGNED AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  PACK (IN COMPACT FORM) AS MANY COLUMNS OF THE "B" MATRIX INTO
!           MEMORY AS POSSIBLE LEAVING SPACE FOR A FULL COLUMN OF THE
!           "D" MATRIX FOR EACH COLUMN OF THE "B" MATRIX READ.
!           SEE SUBROUTINES MMARM1,2,3,4 FOR FORMAT OF COMPACT FORM.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  CALL UNPACK TO READ MATRIX "C".
!       5.  FOR METHOD 40, CALL UNPACK TO READ COLUMNS OF MATRIX "A".
!       6.  FOR METHOD 41, CALL MMARC1,2,3,4 TO READ COLUMNS OF "A" INTO
!           MEMORY IN COMPACT FORM.
 
 
 INTEGER, INTENT(IN OUT)                  :: zi(2)
 REAL, INTENT(OUT)                        :: zr(2)
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
 
 DATA    module / 4HMMA4  , 4H     ,4H    /
 DATA    kzero  / 1H0   /
 DATA    kone   / 1H1   /
DATA    jbegn  / 4HBEGN/ , jend  / 3HEND /

module( 3 ) = jbegn
IF ( nastor == 1 ) module( 2 ) = kzero
IF ( nastor == 2 ) module( 2 ) = kone
CALL conmsg ( module, 3, 0 )
incru  = 1
typei  = ndtype
typep  = ndtype
nwdd   = nwords( ndtype )
nwdb   = nwords( nbtype )

!   OPEN CORE ALLOCATION AS FOLLOWS:
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "A" MATRIX
!     Z( IBX      ) = ARRAY FOR MULTIPLE COLUMNS OF "B" MATRIX
!     Z( IDX      ) = ARRAY FOR MULTIPLE COLUMNS OF "D" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE

ibx = 1 + nwdd*nar
IF ( nastor == 2 .OR.  ksys58 == 41 ) ibx = 1 + nwdd*nar + nar
IF ( MOD ( ibx , 2 ) == 0 ) ibx = ibx + 1
ibx2   = ( ( ibx+1 ) / 2 )
ibuf1  = nz    - sysbuf
ibuf2  = ibuf1 - sysbuf
ibuf3  = ibuf2 - sysbuf
ibuf4  = ibuf3 - sysbuf
lasmems= ibuf4 - 1

! INSURE THAT LASMEM IS ON A QUAD WORD BOUNDARY TO ALLOW FOR COMPLEX
! DOUBLE DATA TO BE REFERENCED FOR "D" MATRIX COLUMNS

itest = MOD( lasmems, 4 )
IF ( itest == 1 ) GO TO 90
IF ( itest == 0 ) lasmems = lasmems - 3
IF ( itest == 2 ) lasmems = lasmems - 5
IF ( itest == 3 ) lasmems = lasmems - 6
90    CONTINUE
iprow1 = 1
iprown = ndr
incrp  = 1
CALL gopen ( filea, zr( ibuf1 ), rdrew )
CALL gopen ( fileb, zr( ibuf2 ), rdrew )
IF ( filec( 1 ) /= 0 .AND. signc /= 0 )  &
    CALL gopen ( filec, zr( ibuf3 ), rdrew )
ipass  = 0
ircoln = 0
nwddndr= nwdd*ndr
CALL gopen ( filed, zr( ibuf4 ), wrtrew)
filed( 2 ) = 0
filed( 6 ) = 0
filed( 7 ) = 0
100   ipass  = ipass + 1
ircol1 = ircoln + 1
ircoln = nbc
irfile = fileb( 1 )
SIGN   = signab
lasmem = lasmems - ibx
IF ( ipass  /= 1 ) CALL dsspos ( irfile, irpos( 1 ), irpos( 2 ),irpos( 3 ) )
IF ( ndtype == 1 ) CALL mmarm1 ( zi( ibx ), zr( ibx  ), nwddndr)
IF ( ndtype == 2 ) CALL mmarm2 ( zi( ibx ), zd( ibx2 ), nwddndr)
IF ( ndtype == 3 ) CALL mmarm3 ( zi( ibx ), zc( ibx2 ), nwddndr)
IF ( ndtype == 4 ) CALL mmarm4 ( zi( ibx ), zd( ibx2 ), nwddndr)
ncolpp = ircoln - ircol1 + 1
ibrow  = ircol1 - 1
idx    = lasmem + ibx
idx2   = ( ( idx+1 ) / 2 ) - 1
idx4   = ( idx+1 ) / 4
IF ( ipass == 1 ) GO TO 300
CALL REWIND( filea )
CALL skprec( filea, 1 )
300   CONTINUE

! NOW READ INTO MEMORY THE "C" FILE.
! READ AS MANY COLUMNS OF THIS AS WERE READ OF THE "B" MATRIX

IF ( filec( 1 ) == 0 .OR. signc == 0 ) GO TO 800
indxd  = idx
typeu  = ndtype * signc
iurow1 = 1
iurown = ncr
DO  i = 1, ncolpp
  CALL unpack( * 650, filec, zr( indxd ) )
  GO TO 680
  
! NULL COLUMN READ ON "C"
  
  650   CONTINUE
  LEN = indxd + nwddndr - 1
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
  zr( k ) = 0.0
END DO
900   CONTINUE

! PROCESS ALL OF THE COLUMNS OF "A" MATRIX

SIGN = 1
IF ( ksys58 == 40 ) GO TO 950
IF ( ksys58 == 41 ) GO TO 1000
IF ( nastor == 2  ) GO TO 1000
! PROCESS ALL OF THE COLUMNS OF "B", ADD "C" DATA ON FIRST PASS
950   CONTINUE
IF ( ndtype == 1 ) CALL mma401( zi, zr )
IF ( ndtype == 2 ) CALL mma402( zi, zd )
IF ( ndtype == 3 ) CALL mma403( zi, zc )
IF ( ndtype == 4 ) CALL mma404( zi, zd, zdc )
GO TO 60000
1000  CONTINUE
IF ( ndtype == 1 ) CALL mma411( zi, zr )
IF ( ndtype == 2 ) CALL mma412( zi, zd )
IF ( ndtype == 3 ) CALL mma413( zi, zc )
IF ( ndtype == 4 ) CALL mma414( zi, zd, zdc )
60000 CONTINUE
IF ( ircoln < nbc ) GO TO 100

! ALL COLUMNS OF A HAVE BEEN PROCESSED, MULTIPLICATION COMPLETE

CALL CLOSE ( filea, clsrew )
CALL CLOSE ( fileb, clsrew )
CALL CLOSE ( filec, clsrew )
CALL CLOSE ( filed, clsrew )
module( 3 ) = jend
CALL conmsg ( module, 3, 0 )
RETURN
END SUBROUTINE mma4
