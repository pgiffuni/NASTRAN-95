SUBROUTINE exfowr  &
        ( iun, iprec, FORM, indata, nwds )
!********************************************************************
!    EXPECTED TYPES OF FORMAT CODES ARE AS FOLLOWS
!        NH------       NENN.N       NDNN.N         NX
!        NFNN.N         NINN         NGNN.N         NAN
!        NPENN.N        NPFNN.N      NPN(----)
!        SPECIAL CHARACTERS:  /(),
!     ICHAR = CURRENT CHARACTER NUMBER BEING PROCESSED IN "FORM"
!     ICOL  = CURRENT CHARACTER COLUMN POSITION WITHIN THE LINE
!     NCNT  = NUMBER OF VALUES OF IDATA AND DATA THAT HAVE BEEN PROCESSE
!********************************************************************
 CHARACTER (LEN=1) :: FORM(1000)
 CHARACTER (LEN=1) :: slash , BLANK
 CHARACTER (LEN=1) :: lparen, rparen, period, comma, NUMBER(10)
 CHARACTER (LEN=1) :: h, e, d, x, f, i, g, a, p
 CHARACTER (LEN=2) :: pfact
 CHARACTER (LEN=4) :: cdata(200)
 CHARACTER (LEN=132) :: line
 CHARACTER (LEN=132) :: tform
 INTEGER*4       idata(200)
 REAL*4          DATA(200)
 REAL*8          ddata(100)
 INTEGER*4       indata(200)
 COMMON /system/ isysbf, iwr
 EQUIVALENCE     ( idata, DATA, ddata, cdata )
 DATA            h/'H'/, e/'E'/, d/'D'/, x/'X'/, f/'F'/
 DATA            i/'I'/, g/'G'/, a/'A'/, p/'P'/
 DATA            lparen /'('/, rparen/')'/, period/'.'/
 DATA            comma  /','/, slash /'/'/, BLANK /' '/
 DATA            NUMBER /'0','1','2','3','4','5','6','7','8','9'/
 IF ( nwds <= 200 ) GO TO 2
 PRINT *,' word limit exceeded in exfowr, limit=200'
 RETURN
 2     DO  kb = 1, nwds
   idata( kb ) = indata( kb )
 END DO
 iloop = 0
 ICHAR = 1
 ncnt  = 1
 icol  = 1
 line  = BLANK
 pfact = BLANK
 icycle= 0
 5     IF ( FORM(ICHAR) == lparen ) GO TO 75
 ICHAR = ICHAR + 1
 IF ( ICHAR <= 1000 ) GO TO 5
 GO TO 7702
 70    IF ( ICHAR > 1000 ) GO TO 7702
 IF ( ncnt  > nwds ) GO TO 1200
 IF ( FORM(ICHAR) == BLANK ) GO TO 75
 IF ( FORM(ICHAR) == slash ) GO TO 100
 IF ( FORM(ICHAR) >= NUMBER(1) .AND. FORM(ICHAR) <= NUMBER(10) ) GO TO 200
 IF ( FORM(ICHAR) == a ) GO TO 300
 IF ( FORM(ICHAR) == i ) GO TO 400
 IF ( FORM(ICHAR) == h ) GO TO 500
 IF ( FORM(ICHAR) == x ) GO TO 600
 IF ( FORM(ICHAR) == p ) GO TO 700
 IF ( FORM(ICHAR) == f ) GO TO 800
 IF ( FORM(ICHAR) == g ) GO TO 800
 IF ( FORM(ICHAR) == d ) GO TO 800
 IF ( FORM(ICHAR) == e ) GO TO 800
 IF ( FORM(ICHAR) == lparen ) GO TO 1000
 IF ( FORM(ICHAR) == rparen ) GO TO 1100
 IF ( FORM(ICHAR) /= comma  ) GO TO 7702
 IF ( icycle == 0 ) pfact = BLANK
 75    ICHAR = ICHAR + 1
 GO TO 70
! PROCESS SLASH
 100   CONTINUE
 IF ( line /= BLANK ) WRITE ( iwr,900 ) line
 900   FORMAT(a132)
 IF ( line == BLANK ) WRITE ( iwr,901 )
 901   FORMAT(/)
 line   = BLANK
 IF ( icycle == 0 ) pfact = BLANK
 icol = 1
 GO TO 75
! GET MULTIPLIER FOR FIELD CONVERSION
 200   CALL fornum ( FORM, ICHAR, imult )
 GO TO 70
! PROCESS ALPHA FIELD--FORMAT(NNANNN) (NN=IMULT,NNN=IFIELD)
 300   ICHAR = ICHAR + 1
 IF ( ncnt > nwds ) GO TO 1200
 CALL fornum ( FORM, ICHAR, ifield )
 IF ( imult == 0 ) imult = 1
 WRITE ( tform, 902 ) imult, ifield
 902   FORMAT('(',i2,'A',i2,')')
 i1 = icol
 length = imult*ifield
 nend   = ncnt + imult - 1
 IF ( nend > nwds ) nend = nwds
 last   = icol + length - 1
 WRITE( line(icol:last), tform ) (cdata(kk),kk=ncnt,nend)
 icol   = icol + length
 ncnt   = ncnt + imult
 imult = 1
 GO TO 70
! PROCESS INTEGER FIELD -- FORMAT(NNINNN) (NN=IMULT,NNN=IFIELD)
 400   ICHAR = ICHAR + 1
 IF ( ncnt > nwds ) GO TO 1200
 CALL fornum ( FORM, ICHAR, ifield )
 IF ( imult == 0 ) imult = 1
 WRITE ( tform, 903 ) imult, ifield
 903   FORMAT('(',i2,'I',i2,')')
 i1 = icol
 length = imult*ifield
 nend   = ncnt + imult - 1
 last   = icol + length - 1
 IF ( nend > nwds ) nend = nwds
 WRITE( line(icol:last), tform ) (idata(kk),kk=ncnt,nend)
 icol   = icol + length
 ncnt   = ncnt + imult
 imult  = 1
 GO TO 70
! PROCESS HOLERITH FIELD -- FORMAT(NNH----) (NN=IMULT)
 500   last   = icol  + imult - 1
 ICHAR  = ICHAR + 1
 lchar  = ICHAR + imult - 1
 WRITE ( line(icol:last), 904 ) (FORM(kk),kk=ICHAR,lchar)
 904   FORMAT(133A1)
 icol   = icol  + imult
 ICHAR  = lchar
 imult  = 1
 GO TO 75
! PROCESS X FIELD -- FORMAT(NNX) (NN=IMULT)
 600   WRITE ( tform, 905 ) imult
 905   FORMAT('(',i2,'X',')')
 last   = icol + imult - 1
 WRITE( line(icol:last), tform )
 icol   = icol + imult
 imult  = 1
 GO TO 75
! PROCESS P FACTOR FOR FLOATING FORMAT
 700   WRITE ( pfact,904 ) FORM(ICHAR-1), FORM(ICHAR)
 IF ( ncnt > nwds ) GO TO 1200
 GO TO 75
! PROCESS FLOATING FIELD -- FORMAT(NPNNXNNN.NNNN)  WHERE
!          (NP = PFACT, NN=IMULT, NNN=IFIELD, NNNN=IDEC)
 800   itype = ICHAR
 IF ( ncnt > nwds ) GO TO 1200
 ICHAR = ICHAR + 1
 CALL fornum ( FORM, ICHAR, ifield )
 810   IF ( FORM( ICHAR ) == period ) GO TO 820
 ICHAR = ICHAR + 1
 GO TO 810
 820   ICHAR = ICHAR + 1
 CALL fornum ( FORM, ICHAR, idec )
 IF ( imult == 0 ) imult = 1
 WRITE ( tform, 906 ) pfact, imult, FORM(itype),ifield, idec
 906   FORMAT('(',a2,i2,a1,i2,'.',i2,')')
 i1 = icol
 length = imult*ifield
 nend   = ncnt + imult - 1
 last   = icol + length - 1
 IF ( nend > nwds ) nend = nwds
 IF ( iprec == 2 ) WRITE( line(icol:last), tform ) (ddata(kk),kk=ncnt,nend)
 IF ( iprec /= 2 ) WRITE( line(icol:last), tform ) (DATA(kk),kk=ncnt,nend)
 icol   = icol + length
 ncnt   = ncnt + imult
 imult  = 1
 GO TO 70
! PROCESS LEFT PAREN (NOT THE FIRST LEFT PAREN BUT ONE FOR A GROUP)
! IMULT HAS THE MULTIPLIER TO BE APPLIED TO THE GROUP
 1000  icycle = imult-1
 icsave = ICHAR+1
 iloop  = 1
 imult  = 1
 GO TO 75
! PROCESS RIGHT PAREN ( CHECK IF IT IS THE LAST OF THE FORMAT)
! IF IT IS PART OF A GROUP, THEN ICYCLE WILL BE NON-ZERO
 1100  IF ( icycle > 0 ) GO TO 1110
 IF ( iloop  /= 0 ) GO TO 1120
 IF ( ncnt > nwds ) GO TO 1200
! NO GROUP, THEREFORE MUST RE CYCLE THROUGH FORMAT
! UNTIL LIST IS SATISFIED
 WRITE ( iun,900 ) line
 ICHAR  = 2
 line   = BLANK
 pfact  = BLANK
 icol   = 1
 GO TO 70
! GROUP BEING PROCESSED, DECREMENT COUNT AND RESET ICHAR TO BEGINNING
! OF THE GROUP
 1110  icycle = icycle - 1
 ICHAR  = icsave
 GO TO 70
! FINISHED WITH LOOP, CONTINUE WITH FORMAT
 1120  iloop  = 0
 icycle = 0
 GO TO 75
 1200  WRITE ( iun,900 ) line
 IF ( ncnt > nwds ) GO TO 7000
 line = BLANK
 GO TO 70
 7000  CONTINUE
 RETURN
 7702  WRITE( iwr, 9901 ) ICHAR, FORM
 9901  FORMAT(///' SUBROUTINE EXFOWR UNABLE TO DECIPHER THE FOLLOWING'  &
     ,' FORMAT AT CHARACTER ',i4,/,' FORMAT GIVEN WAS THE FOLLOWING:'  &
     ,/,(1X,131A1))
END SUBROUTINE exfowr
