SUBROUTINE mma3 ( zi, zr, zd, zc, zdc )
     
!     MMA3 PERFORMS THE MATRIX OPERATION USING METHODS 30, 31 AND 32
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA3 IS DESIGNED AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  PACK (IN COMPACT FORM) AS MANY COLUMNS OF "A" INTO MEMORY
!           AS POSSIBLE LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!           SEE SUBROUTINES MMARM1,2,3,4 FOR FORMAT OF COMPACT FORM.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  FOR METHODS 30 AND 31, CALL UNPACK TO READ MATRIX "C".
!       5.  FOR METHOD 30, CALL UNPACK TO READ COLUMNS OF MATRIX "B".
!       6.  FOR METHOD 31, CALL MMARC1,2,3,4 TO READ COLUMNS OF "B" INTO
!           MEMORY IN COMPACT FORM.
!       7.  FOR METHOD 32, CALL MMARC1,2,3,4 TO READ COLUMNS OF "B" AND
!           "C" INTO MEMORY IN COMPACT FORM.
!       8.  FOR METHODS 30 AND 31, CALL PACK TO WRITE "D" MATRIX.
!       9.  FOR METHOD 32, CALL BLDPK TO WRITE "D" MATRIX.
 
 
 INTEGER, INTENT(IN OUT)                  :: zi(2)
 REAL, INTENT(IN OUT)                     :: zr(2)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(2)
 COMPLEX, INTENT(IN OUT)                  :: zc(2)
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
 
 DATA    module / 4HMMA3  , 4H     ,4H    /
 DATA    kzero  / 1H0   /
 DATA    kone   / 1H1   /
 DATA    ktwo   / 1H2   /
DATA    jbegn  / 4HBEGN/ , jend  / 3HEND /

module( 3 ) = jbegn
IF ( method == 30 ) module( 2 ) = kzero
IF ( method == 31 ) module( 2 ) = kone
IF ( method == 32 ) module( 2 ) = ktwo
CALL conmsg ( module, 3, 0 )
incru  = 1
typei  = ndtype
typep  = ndtype
nwdd   = nwords( ndtype )
nwdb   = nwords( nbtype )

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

idx = 1 + nwdd*nbr
IF ( method /= 31 .AND. method /= 32) GO TO 90

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
IF ( method /= 32 ) GO TO 96

! FOR METHOD 32, INSURE IAX IS ON QUAD WORD BOUNDARY FOR COMPLEX DOUBLE

iax    = idx  + nwdd*ndr + ndr
itest = MOD( iax, 4 )
IF ( itest == 1 ) GO TO 96
IF ( itest == 0 ) iax = iax + 1
IF ( itest == 2 ) iax = iax + 3
IF ( itest == 3 ) iax = iax + 2
96    CONTINUE
iax2   = ( ( iax+1 ) / 2 )
ibuf1  = nz    - sysbuf
ibuf2  = ibuf1 - sysbuf
ibuf3  = ibuf2 - sysbuf
ibuf4  = ibuf3 - sysbuf
lasmem = ibuf4 - 1
lasmem = lasmem - iax
iprow1 = 1
iprown = ndr
incrp  = 1
CALL gopen  ( filea, zr( ibuf1 ), rdrew )
CALL gopen  ( fileb, zr( ibuf2 ), rdrew )
ipass  = 0
ircoln = 0
100   ipass  = ipass + 1
ircol1 = ircoln + 1
ircoln = nac
irfile = filea( 1 )
SIGN   = signab
IF ( ipass  /= 1 ) CALL dsspos ( irfile, irpos( 1 ), irpos( 2 ),irpos( 3 ) )
IF ( ndtype == 1 ) CALL mmarm1 ( zi( iax ), zr( iax  ), 0 )
IF ( ndtype == 2 ) CALL mmarm2 ( zi( iax ), zd( iax2 ), 0 )
IF ( ndtype == 3 ) CALL mmarm3 ( zi( iax ), zc( iax2 ), 0 )
IF ( ndtype == 4 ) CALL mmarm4 ( zi( iax ), zd( iax2 ), 0 )
ncolpp = ircoln - ircol1 + 1
ibrow  = ircol1 - 1
IF ( ircoln == nac ) GO TO 400
itest = MOD( ipass, 2 )
IF ( itest == 0 ) GO TO 350
ifile = scrtch
ofile = filed( 1 )
GO TO 380
350   ifile = filed( 1 )
ofile = scrtch
380   CONTINUE
IF ( ipass == 1 ) GO TO 300
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
IF ( ifile == 0 ) ifile = scrtch
IF ( ofile == filed( 1 ) .AND. ipass /= 1 ) CALL filswi( ifile, ofile )
ofile = filed( 1 )
ifile = scrtch
CALL REWIND( fileb )
CALL skprec( fileb, 1 )
CALL gopen ( filed, zr( ibuf4 ), wrtrew)
filed( 2 ) = 0
filed( 6 ) = 0
filed( 7 ) = 0
IF ( ipass == 1 ) GO TO 310
CALL gopen  ( ifile, zr( ibuf3 ), rdrew )
490   CONTINUE
SIGN = 1
IF ( method == 30 ) GO TO 950
IF ( method == 31 ) GO TO 1000
IF ( method == 32 ) GO TO 2000
! PROCESS ALL OF THE COLUMNS OF "B", ADD "C" DATA ON FIRST PASS
950   CONTINUE
IF ( ndtype == 1 ) CALL mma301( zi, zr )
IF ( ndtype == 2 ) CALL mma302( zi, zd )
IF ( ndtype == 3 ) CALL mma303( zi, zc )
IF ( ndtype == 4 ) CALL mma304( zi, zd, zdc )
GO TO 60000
1000  CONTINUE
IF ( ndtype == 1 ) CALL mma311( zi, zr )
IF ( ndtype == 2 ) CALL mma312( zi, zd )
IF ( ndtype == 3 ) CALL mma313( zi, zc )
IF ( ndtype == 4 ) CALL mma314( zi, zd, zdc )
GO TO 60000
2000  CONTINUE
IF ( ndtype == 1 ) CALL mma321( zi, zr )
IF ( ndtype == 2 ) CALL mma322( zi, zd )
IF ( ndtype == 3 ) CALL mma323( zi, zc )
IF ( ndtype == 4 ) CALL mma324( zi, zd )
60000 CONTINUE
CALL CLOSE ( ifile, clsrew )
CALL CLOSE ( ofile, clsrew )
IF ( ircoln < nac ) GO TO 100

! ALL COLUMNS OF A HAVE BEEN PROCESSED, MULTIPLICATION COMPLETE

CALL CLOSE ( filea, clsrew )
CALL CLOSE ( fileb, clsrew )
module( 3 ) = jend
CALL conmsg ( module, 3, 0 )
RETURN
END SUBROUTINE mma3
