SUBROUTINE mma ( zi, zr, zd )
     
!     MMA PERFORMS THE MATRIX OPERATION
!       (+/-)A    * B (+/-)C = D   OR
!       (+/-)A(T) * B (+/-)C = D
 
!     USING METHODS 10, 11, 20, 21, 30, 31, 32, 40, 41
 
 
!  IN REGARDS TO THE METHODS BELOW, WHEN MULTIPLE COLUMNS OF A MATRIX
!  ARE STORED AND READ BY GETSTR, THEN THE MATRIX IS STORED IN MEMORY IN
!  COMPACT FORM.  SEE SUBROUTINES 'MMARM1,2,3,4' FOR A DESCRIPTION OF
!  THIS COMPACT FORM.  WHEN ONLY A SINGLE COLUMN OF A MATRIX IS STORED
!  AND IT IS BEING READ BY GETSTR, IT IS STORED IN COMPACT FORM IN MEMORY.
!  SEE SUBROUTINES 'MMARC1,2,3,4' FOR A DESCRIPTION OF THIS FORM.
 
!   ------------------------------------------------------------------------
!   METHOD     METHOD OF READING MATRIX    MULTIPLE COLUMNS OF MATRIX STORED
!                 A        B       C           A         B        D
!   ------------------------------------------------------------------------
!     10        UNPACK  UNPACK   UNPACK       YES        NO       NO
!     11        UNPACK  GETSTR   UNPACK       YES        NO       NO
!     20        UNPACK  UNPACK   UNPACK       NO         YES      YES
!     21        GETSTR  UNPACK   UNPACK       NO         YES      YES
!     30        GETSTR  UNPACK   UNPACK       YES        NO       NO
!     31        GETSTR  GETSTR   UNPACK       YES        NO       NO
!     32        GETSTR  GETSTR   GETSTR       YES        NO       NO
!     40        UNPACK  GETSTR   UNPACK       NO         YES      YES
!     41        GETSTR  GETSTR   UNPACK       NO         YES      YES
!   ------------------------------------------------------------------------
 
!   TO DETERMINE WHICH METHOD TO USE, THE FOLLOWING RATIONAL IS USED.
 
!   1.  DETERMINE THE METHOD FOR READING MATRICES "A" AND "B".  THIS IS
!       DETERMINED BY EXAMINING THE FOLLOWING PERCENTAGE:
 
!            (MEMORY TO CONTAIN ENTIRE MATRIX)
!            ------------------------------------  = PERCENTAGE
!            (MEMORY TO CONTAIN COMPACTED MATRIX)
 
!       IF THE PERCENTAGE IS .GE. THE VARIABLE "TESTPCT", THEN UNPACK IS
!       USED.  OTHERWISE, GETSTR IS USED.
 
!            NiSTOR (i=A or B) = 1, CALL UNPACK TO READ MATRIX
!                              = 2, CALL GETSTR TO READ MATRIX
 
!    2. THE RESULTS OF THE FIRST TEST WILL NARROW THE OPTIONS TO TWO
!       DIFFERENT METHODS AS FOLLOWS:
 
!                                  CANDIDATE METHOD
!                10     11     20     21     30     31     32     40     41
!       NASTOR =  1      1      1      2      2      2      2      1      2
!       NBSTOR =  1      2      1      1      1      2      2      2      2
 
!          FOR NASTOR = 1 AND NBSTOR = 1, METHODS 10 AND 20 ARE CONSIDERED
!          FOR NASTOR = 1 AND NBSTOR = 2, METHODS 11 AND 40 ARE CONSIDERED
!          FOR NASTOR = 2 AND NBSTOR = 1, METHODS 21 AND 30 ARE CONSIDERED
!          FOR NASOTR = 2 AND NBSTOR = 2, METHODS 31,32 AND 41 ARE CONSIDERED
!            (NOTE, METHOD 32 IS ONLY AVAILABLE WITH "A" TRANSPOSED)
 
!    3. LASTLY, DETERMINE THE ESTIMATED NUMBER OF PASSES FOR EACH OF THE
!       TWO CANDIDATE METHODS.  THE METHOD WITH THE FEWER NUMBER OF PASSES
!       IS CHOSEN.
 
!       MPASSii (ii=10,11,20,21,30,31,32,40,41) = ESTIMATED NUMBER OF PASSES
!                                                 FOR METHOD ii.
 
!       NiTOTAL (i=A,B,C) = MEMORY WORDS TO CONTAIN ENTIRE FULL MATRIX
!       NiPACK  (i=A,B,C) = MEMORY WORDS TO CONTAIN ENTIRE MATRIX IN COMPACT
!                           FORM.
!       NWDD              = NUMBER OF WORDS FOR EACH ELEMENT OF THE "D" MATRIX
 
 
!     THE FOLLOWING SUBROUTINES ARE CALLED FOR THE DIFFERENT METHODS AND
!     MATRIX "D" TYPES (RS,RD,CS,CD).
 
!          METHODS  MAIN      OTHER SUBROUTINES DEPENDING ON TYPE
!                 SUBROUTINE    RS     RD     CS     CD
!            10     MMA1      MMA101 MMA102 MMA103 MMA104
!            11     MMA1      MMA111 MMA112 MMA113 MMA114
!            20     MMA2      MMA201 MMA202 MMA203 MMA204
!            21     MMA2      MMA211 MMA212 MMA213 MMA214
!            30     MMA3      MMA301 MMA302 MMA303 MMA304
!            31     MMA3      MMA311 MMA312 MMA313 MMA314
!            32     MMA3      MMA321 MMA322 MMA323 MMA324 (TRANSPOSE ONLY)
!            40     MMA4      MMA401 MMA402 MMA403 MMA404
!            41     MMA4      MMA411 MMA412 MMA413 MMA414
! ---------------------------------------------------------------------------
 
 INTEGER, INTENT(IN OUT)                  :: zi(2)
 REAL, INTENT(IN OUT)                     :: zr(2)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(2)
 INTEGER :: namea(2) ,nameb(2)  ,namec(2) , NAMED(2)
 INTEGER :: prntyp(4) ,module(3),prec1
 INTEGER :: blk1(15) ,blk2(15)  ,eol      ,eor
 INTEGER :: signab   ,signc     ,t        ,scrtch
 INTEGER :: filea    ,fileb     ,filec    ,filed
 INTEGER :: sysbuf   ,typei     ,typep    ,typeu
 INTEGER :: isave(9)
 
 DOUBLE PRECISION :: ad(2)   , dd(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 CHARACTER (LEN=6) :: upmeth(2)
 CHARACTER (LEN=2) :: ct
 INCLUDE           'MMACOM.COM'
 COMMON / names  / rd       ,rdrew     ,wrt     ,wrtrew,clsrew,cls
 COMMON / xmssg  / ufm      ,uwm       ,uim
 COMMON / logout / lout
 COMMON / mpyadx / filea(7) ,fileb(7)  ,filec(7)  &
     ,                 filed(7) ,nz        ,t       ,signab,signc,prec1  &
     ,                 scrtch   ,time
 COMMON / system / ksystm(152)
 COMMON / packx  / typei    ,typep     ,irow1p  ,irownp, incrp
 COMMON / unpakx / typeu    ,irowu     ,irownu  ,incru
 COMMON / zblpkx / d(4)     ,irowbk
 COMMON / zntpkx / a(4)     ,irowin    ,eol     ,eor
 EQUIVALENCE      (ad(1)     ,a(1)  ) , (dd(1)     ,d(1) )
 EQUIVALENCE      (ksystm( 1),sysbuf) , (ksystm( 2),nout )  &
     ,                (ksystm(58),ksys58) , (ksystm(40),nbpw )  &
     ,                (ksystm(55),iprec )
 EQUIVALENCE      (filea(2)  ,nac   ) , (filea(3)  ,nar   )  &
     ,                (filea(4)  ,naform) , (filea(5)  ,natype)  &
     ,                (filea(6)  ,nanzwd) , (filea(7)  ,nadens)
 EQUIVALENCE      (fileb(2)  ,nbc   ) , (fileb(3)  ,nbr   )  &
     ,                (fileb(4)  ,nbform) , (fileb(5)  ,nbtype)  &
     ,                (fileb(6)  ,nbnzwd) , (fileb(7)  ,nbdens)
 EQUIVALENCE      (filec(2)  ,ncc   ) , (filec(3)  ,ncr   )  &
     ,                (filec(4)  ,ncform) , (filec(5)  ,nctype)  &
     ,                (filec(6)  ,ncnzwd) , (filec(7)  ,ncdens)
 EQUIVALENCE      (filed(2)  ,ndc   ) , (filed(3)  ,ndr   )  &
     ,                (filed(4)  ,ndform) , (filed(5)  ,ndtype)  &
     ,                (filed(6)  ,ndnzwd) , (filed(7)  ,nddens)
 
 DATA    module / 4HMMA , 2*4H    /
 DATA    jbegn  /  4HBEGN/, jend  / 3HEND/
 DATA    upmeth / 'UNPACK', 'STRING' /
 DATA prntyp / 2HRS, 2HRD, 2HCS, 2HCD /
 DATA testpct / .8 /
 
 isave( 1 ) = typei
 isave( 2 ) = typep
 isave( 3 ) = irow1p
 isave( 4 ) = irownp
 isave( 5 ) = incrp
 isave( 6 ) = typeu
 isave( 7 ) = irowu
 isave( 8 ) = irownu
 isave( 9 ) = incru
 CALL sswtch ( 19, l19 )
 module( 3 ) = jbegn
 CALL conmsg ( module, 3, 0 )
 ndr = nar
 ndc = nbc
 IF ( t   /= 0   ) ndr = nac
 IF ( ndform /= 0 ) GO TO 50
 ndform = 2
 IF ( ndr == ndc ) ndform = 1
 50    CONTINUE
 IF ( filea( 6 ) == 0 .OR. fileb( 6 ) == 0 ) GO TO 5000
 IF ( signab     == 0 ) GO TO 5000
 IF ( t        /= 0   ) GO TO 100
 IF ( nac      /= nbr ) GO TO 7001
 IF ( filec(1) == 0   ) GO TO 200
 IF ( nar      /= ncr ) GO TO 7001
 IF ( nbc      /= ncc ) GO TO 7001
 GO TO 200
 100   CONTINUE
 IF ( nar      /= nbr ) GO TO 7001
 IF ( filec(1) == 0   ) GO TO 200
 IF ( nac      /= ncr ) GO TO 7001
 IF ( nbc      /= ncc ) GO TO 7001
 200   CONTINUE
 nwdc = 0
 CALL dssize ( filea, ncols, naterms, nastrgs, nwda )
 CALL dssize ( fileb, ncols, nbterms, nbstrgs, nwdb )
 IF ( filec( 1 ) /= 0 ) CALL dssize ( filec, ncols, ncterms, ncstrgs, nwdc )
 nwdd    = MAX0 ( nwda, nwdb, nwdc )
 ndtype  = 2
 IF ( nwdd == 4 ) ndtype = 4
 IF ( nwdd == 1 ) ndtype = 1
 IF ( ndtype == 1 .OR. ndtype == 4 ) GO TO 250
 itest1  = MIN0 ( natype, nbtype, nctype )
 itest2  = MAX0 ( natype, nbtype, nctype )
 ndtype  = 3
 IF ( itest2 == 3 .AND.  &
     ( natype == 2 .OR. nbtype == 2 .OR. nctype == 2 ) ) ndtype = 4
 IF ( itest2 <= 2 ) ndtype = 2
 250   CONTINUE
 natotal = nac * nar * nwdd
 nbtotal = nbc * nbr * nwdd
 IF ( filec(1) /= 0 ) nctotal = ncc * ncr * nwdd
 ndtotal = ndc * ndr * nwdd
 napack  = 2*nac + 2*nastrgs + naterms*nwdd
 nbpack  = 2*nbc + 2*nbstrgs + nbterms*nwdd
 IF ( filec(1) /= 0 ) ncpack  = 2*ndc + 2*ncstrgs + ncterms*nwdd
 denstya = ( nadens*1.) / 10000.
 denstyb = ( nbdens*1.) / 10000.
 IF ( filec(1) /= 0 ) denstyc = ( ncdens*1.) / 10000.
 nastor  = 2
 nbstor  = 2
 ncstor  = 2
 x = natotal
 y = napack
 percnta = y / x
 x = nbtotal
 y = nbpack
 percntb = y / x
 IF ( filec( 1 ) == 0 ) GO TO 300
 x = nctotal
 y = ncpack
 percntc = y / x
 300   CONTINUE
 IF ( percnta >= testpct ) nastor = 1
 IF ( percntb >= testpct ) nbstor = 1
 IF ( filec(1) /= 0 .AND. percntc >= testpct ) ncstor = 1
 memavl  = nz - 4*sysbuf
 mpass10 = (natotal / ( memavl - (nbr + ndr)*nwdd       ) ) + 1
 mpass11 = (natotal / ( memavl - ndr*nwdd - (nbpack/nbc)) ) + 1
 mpass20 = ((nbtotal + ndtotal) / (memavl - nar*nwdd    ) ) + 1
 mpass21 = ((nbtotal + ndtotal) / (memavl - (napack/nac)) ) + 1
 mpass30 = (napack  / ( memavl - (nbr + ndr)*nwdd       ) ) + 1
 mpass31 = (napack  / ( memavl - ndr*nwdd - (nbpack/nbc)) ) + 1
 mpass32 = (napack  / ( memavl - (ncpack/ndc) - (nbpack/nbc)) ) + 1
 mpass40 = ((nbpack + ndtotal) / (memavl - nar*nwdd     ) ) + 1
 mpass41 = ((nbpack + ndtotal) / (memavl - (napack/nac) ) ) + 1
 IF ( nastor == 1 .AND. nbstor == 1 ) GO TO 1000
 IF ( nastor == 2 .AND. nbstor == 1 ) GO TO 1100
 IF ( nastor == 1 .AND. nbstor == 2 ) GO TO 1200
 IF ( nastor == 2 .AND. nbstor == 2 ) GO TO 1300
 1000  CONTINUE
!---------USE UNPACK FOR MATRICES "A" AND "B"  (CHOOSE METHOD 10 OR 20)
 method = 10
 IF ( mpass10 == 1       ) GO TO 2000
 IF ( mpass10 <= mpass20 ) GO TO 2000
 method = 20
 GO TO 2000
 1100  CONTINUE
!---------USE GETSTR FOR MATRIX "A"; UNPACK FOR MATRIX "B"
!         (CHOOSE METHOD 21 OR 30)
 method = 21
 IF ( mpass21 == 1       ) GO TO 2000
 IF ( mpass21 <= mpass30 ) GO TO 2000
 method = 30
 GO TO 2000
 1200  CONTINUE
!---------USE UNPACK FOR MATRIX "A"; GETSTR FOR MATRIX "B"
!         (CHOOSE METHOD 11 OR 40)
 method = 11
 IF ( mpass11 == 1       ) GO TO 2000
 IF ( mpass11 <= mpass40 ) GO TO 2000
 method = 40
 GO TO 2000
 1300  CONTINUE
!---------USE GETSTR FOR MATRICES "A" AND "B" (CHOOSE METHOD 31, 32 OR 41)
 method = 31
 IF ( mpass31 == 1       ) GO TO 1310
 IF ( mpass31 <= mpass41 ) GO TO 1310
 method = 41
 GO TO 2000
 1310  CONTINUE
 IF ( ncstor == 2 .AND. t /= 0 ) method = 32
 2000  CONTINUE
 IF(l19 == 0) GO TO 3000
 CALL fname ( filea, namea )
 CALL fname ( fileb, nameb )
 CALL fname ( filec, namec )
 CALL fname ( filed, NAMED )
 WRITE( lout,2001, IOSTAT=ierr )  &
     namea, nar, nac, naterms, denstya, prntyp( natype )  &
     ,        nameb, nbr, nbc, nbterms, denstyb, prntyp( nbtype )
 2001  FORMAT(  &
     '  /-----------------------------------------------------------/' ,/  &
     ,'  /     MATRIX      ROWS   COLS     TERMS   DENS    TYPE      /' ,/  &
     ,'  /-----------------------------------------------------------/' ,/  &
     ,'     A- ',2A4,i8,i7,i10,f7.4, 5X, a2 ,/  &
     ,'     B- ',2A4,i8,i7,i10,f7.4, 5X, a2 )
 IF (filec(1) == 0) GO TO 2010
 WRITE( lout,2002, IOSTAT=ierr )  &
     namec, ncr, ncc, ncterms, denstyc, prntyp( nctype )
 2002  FORMAT( '     C- ',2A4,i8,i7,i10, f7.4, 5X, a2 )
 2010  WRITE( lout, 2003 ) NAMED, ndr, ndc, prntyp(ndtype)
 2003  FORMAT('     D- ',2A4, i8, i7, 10X, 7X,   5X, a2 )
 WRITE( lout, 2004 ) signab, signc, nz, ksys58
 2004  FORMAT('     SIGNAB =',i2,'  SIGNC =',i2,'  MEMORY =',i10  &
     ,'  SYSTEM(58)=',i3 )
 WRITE( lout, 2005 ) upmeth( nastor ), natotal, napack  &
     ,                   upmeth( nbstor ), nbtotal, nbpack
 IF ( filec( 1 ) /= 0 ) WRITE( lout, 20051) upmeth( ncstor ), nctotal, ncpack
 WRITE( lout, 20052) t, method, prntyp( ndtype )
 2005  FORMAT(  &
     '  /-----------------------------------------------------------/' ,/  &
     ,'  /    READ METHOD   MEMORY (FULL MATRIX)    MEMORY (STRINGS) /' ,/  &
     ,'  /-----------------------------------------------------------/' ,/  &
     ,'     A-  ',a6,i21,i21 ,/  &
     ,'     B-  ',a6,i21,i21 )
 20051 FORMAT( '     C-  ',a6,i21,i21 )
 20052 FORMAT( '     T =',i2,'    SUGGESTED METHOD =',i2  &
     ,'    "D" MATRIX TYPE:',1X,a2)
 WRITE( lout, 2006 ) mpass10,mpass11,mpass20,mpass21,mpass30  &
     ,                           mpass31,mpass32,mpass40,mpass41
 2006  FORMAT(  &
     '  /-----------------------------------------------------------/' ,/  &
     '  /       ESTIMATED NUMBER OF PASSES REQUIRED PER METHOD      /' ,/  &
     ,'  /         10   11   20   21   30   31   32   40   41        /' ,/  &
     ,'  /-----------------------------------------------------------/' ,/  &
     ,'         ',9I5 ,/  &
     ,'  /-----------------------------------------------------------/' )
 3000  CONTINUE
 IF ( filed( 1 ) < 0 ) GO TO 7777
 IF (  ksys58 /= 0 .AND. (ksys58 >= 10 .AND. ksys58 <= 11) .OR.  &
     (ksys58 >= 20 .AND. ksys58 <= 21) .OR.  &
     (ksys58 >= 30 .AND. ksys58 <= 31) .OR.  &
     (ksys58 >= 40 .AND. ksys58 <= 41) )  method = ksys58
 IF ( ksys58 == 32 .AND. t /= 0 ) method = ksys58
 IF ( method == 10 ) nbstor = 1
 IF ( method == 11 ) nbstor = 2
 IF ( method == 20 ) nastor = 1
 IF ( method == 21 ) nastor = 2
 IF ( method == 30 ) nbstor = 1
 IF ( method == 31 ) nbstor = 2
 IF ( method == 32 ) nbstor = 2
 IF ( method == 40 ) nastor = 1
 IF ( method == 41 ) nastor = 2
 IF ( method == 10 .OR. method == 11 ) CALL mma1 ( zi, zr, zd, zr, zd )
 IF ( method == 20 .OR. method == 21 ) CALL mma2 ( zi, zr, zd, zr, zd )
 IF ( method >= 30 .AND. method <= 32 ) CALL mma3 ( zi, zr, zd, zr, zd )
 IF ( method == 40 .OR. method == 41 ) CALL mma4 ( zi, zr, zd, zr, zd )
 ct = 'NT'
 IF ( t /= 0 ) ct = 't '
 WRITE ( lout, 2007 ) method, ct, ipass
 2007  FORMAT('   METHOD USED = ',i2,a2,'  ACTUAL NUMBER OF PASSES =',i4)
 GO TO 7777
 
! "A" AND "B" MATRICES ARE NULL, MOVE "C" TO "D" IF "C" EXISTS
 
 5000  CONTINUE
 IF ( filed( 1 ) < 0 ) GO TO 7777
 ndtype = nctype
 WRITE ( lout, 9002 )
 9002  FORMAT('       MMA - NULL MATRIX PRODUCT')
 ibuf1 = nz    - sysbuf
 ibuf2 = ibuf1 - sysbuf
 IF ( filec( 1 ) == 0 ) GO TO 5900
 IF ( signc      == 0 ) GO TO 5900
 IF ( signc      < 0 ) GO TO 5500
 
! USE CPYSTR TO COPY "C" TO "D"
 
 blk1( 1 ) = filec( 1 )
 blk2( 1 ) = filed( 1 )
 CALL gopen ( filec, zr( ibuf1 ), rdrew )
 CALL gopen ( filed, zr( ibuf2 ), wrtrew)
 DO  i = 1, ncc
   CALL cpystr ( blk1, blk2, 0, 0 )
 END DO
 CALL CLOSE ( filed, clsrew )
 CALL CLOSE ( filec, clsrew )
 filed( 2 ) = filec( 2 )
 filed( 3 ) = filec( 3 )
 filed( 4 ) = filec( 4 )
 filed( 5 ) = filec( 5 )
 filed( 6 ) = filec( 6 )
 filed( 7 ) = filec( 7 )
 GO TO 7777
 
! USE INTPK/BLDPK TO COPY C TO D BECAUSE SIGNS CONFLICT
 
 5500  CONTINUE
 filed( 2 ) = 0
 filed( 6 ) = 0
 filed( 7 ) = 0
 CALL gopen ( filec, zr( ibuf1 ), rdrew )
 CALL gopen ( filed, zr( ibuf2 ), wrtrew)
 DO  i = 1, ncc
   CALL bldpk ( ndtype, ndtype, filed, blk1, 1 )
   CALL intpk ( *5550 , filec , 0, ndtype*signc, 0 )
   5510  CALL zntpki
   CALL bldpki ( a, irowin, filed, blk1 )
   IF ( eol == 0 ) GO TO 5510
   5550  CALL bldpkn ( filed, blk1, filed )
 END DO
 filed( 3 ) = filec( 3 )
 filed( 4 ) = filec( 4 )
 filed( 5 ) = filec( 5 )
 CALL CLOSE ( filec, clsrew )
 CALL CLOSE ( filed, clsrew )
 GO TO 7777
 
! CREATE NULL MATRIX BECAUSE "C" MATRIX IS NULL
 
 5900  CONTINUE
 ndr = 0
 ndc = 0
 CALL gopen ( filed, zr( ibuf1 ) , wrtrew )
 ndc = nbc
 ndr = nar
 IF ( nar == nbc ) ndr = nac
 dd( 1 ) = 0.0D0
 incrp   = 1
 irow1p  = 1
 irownp  = 1
 typei   = prec1
 IF ( typei == 0 ) typei = 1
 typep   = typei
 numc    = ndc
 filed( 2 ) = 0
 filed( 3 ) = ndr
 filed( 5 ) = iprec
 filed( 6 ) = 0
 filed( 7 ) = 0
 DO  i = 1, numc
   CALL pack ( dd, filed, filed )
 END DO
 CALL CLOSE ( filed, clsrew )
 GO TO 7777
! MATRICES ARE INCOMPATIBLE FOR MULTIPLICATION
 7001  CONTINUE
 WRITE ( nout, 9001 ) ufm
 9001  FORMAT( a23, ' MATRICES FOR MULTIPLICATION HAVE INCOMPATIBLE SIZES',/)
 CALL fname ( filea, namea )
 CALL fname ( fileb, nameb )
 CALL fname ( filec, namec )
 CALL fname ( filed, NAMED )
 WRITE( nout,2001, IOSTAT=ierr )  &
     namea, nar, nac, naterms, denstya, prntyp( natype )  &
     ,        nameb, nbr, nbc, nbterms, denstyb, prntyp( nbtype )
 IF ( filec(1) == 0) GO TO 7002
 WRITE( nout,2002, IOSTAT=ierr )  &
     namec, ncr, ncc, ncterms, denstyc, prntyp( nctype )
 7002  CALL mesage ( -61, 0, 0 )
 7777  CONTINUE
 module( 3 ) = jend
 CALL conmsg ( module, 3, 0 )
 typei  = isave( 1 )
 typep  = isave( 2 )
 irow1p = isave( 3 )
 irownp = isave( 4 )
 incrp  = isave( 5 )
 typeu  = isave( 6 )
 irowu  = isave( 7 )
 irownu = isave( 8 )
 incru  = isave( 9 )
 RETURN
END SUBROUTINE mma
