SUBROUTINE strqd2( npts, ti )
     
!     ****PHASE II OF STRESS DATA RECOVERY*********
 
!     NPTS = 3 IMPLIES STRIA1 OR STRIA2  (PHASE II)
!     NPTS = 4 IMPLIES SQUAD1 OR SQUAD2  (PHASE II)
 
 
 INTEGER, INTENT(IN)                      :: npts
 REAL, INTENT(IN)                         :: ti(6)
 LOGICAL :: flag
 LOGICAL :: strain
 
 REAL :: sdelta(3),sstrss(3),frlast(2)
 
 INTEGER :: ist(10)
 INTEGER :: tloads   ,eject    ,ifrvec(12),ished(7),tr       ,qu  &
     ,       ontw(2)  ,istyp(2)
 
 DIMENSION nsil(4), str(18), nph1ou(2), si(36)
 
 COMMON /BLANK / idummy(10), strain
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35), ivec, ivecn, ldtemp, deform,dum8(8), tloads
 COMMON /sdr2x7/ ph1out(200),forvec(6)
 COMMON /sdr2x8/ temp,delta,npoint,i,j,npt1,vec(5),tem  &
     ,            z1ovri, z2ovri,stress(3),chpout(30),cvec(5)  &
     ,            cfrvec(12)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 COMMON /sdr2de/ skp2de(8),ieltyp
 
 EQUIVALENCE (nsil(1),ph1out(2))  &
     ,(nph1ou(1),ph1out(1)) ,(si(1),ph1out(9))  &
     ,(str(1),ph1out(101)) ,(ldtemp,ftemp)  &
     ,(f1,n1) ,  (cfrvec(1),ifrvec(1)) ,  (ished(6),frlast(1))  &
     ,  (ished(1),lsub) , (ished(2),lld)
 EQUIVALENCE (str(1), ist(1))
 
 DATA tr,qu,ontw / 4HTRIA , 4HQUAD , 1H1 , 1H2  / , lld  / -100 /
 DATA lsub,frlast / -100 , -1.0E30 , -1.0E30   /
 DATA iblank /4H    /
! **********************************************************************
! **********************************************************************
 
!     PHASE I OUTPUT FROM THE PLATE IS THE FOLLWOING
 
!     PH1OUT(1)                        ELEMENT ID
!     PH1OUT(2 THRU 5)                 3 SILS AND DUMMY OR 4 SILS
!     PH1OUT(6)                        I
!     PH1OUT(7 THRU 8)                 Z1 AND Z2
!     PH1OUT(9 THRU 30*NPTS+8)         3 OR 4 S SUB I  5X6 ARRAYS
!     PH1OUT(30*NPTS+9 THRU 30*NPTS+11)  S SUB T MATRIX
 
! **********************************************************************
 
!     PHASE I OUTPUT FROM THE MEMBRANE IS THE FOLLOWING
!     NOTE..BEGIN = 30*NPTS+11
 
!     PH1OUT(BEGIN + 1)                ELEMENT ID
!     PH1OUT(BEGIN + 2 THRU BEGIN + 5) 3 SILS AND DUMMY OR 4 SILS
!     PH1OUT(BEGIN + 6)                T SUB 0
!     PH1OUT(BEGIN + 7 THRU BEGIN + 9) S SUB T  3X1 ARRAY
!     PH1OUT(BEGIN + 10 THRU BEGIN + 9*NPTS+9) 3 OR 4 S SUB I 3X3 ARRAYS
 
! **********************************************************************
! **********************************************************************
 
!     THE ABOVE ELEMENTS ARE COMPOSED OF PLATES AND MEMBRANES...
!     SOME MAY ONLY CONTAIN PLATES WHILE OTHERS MAY ONLY CONTAIN
!     MEMBRANES.
 
!     A CHECK FOR A ZERO FIRST SIL IN THE PHASE I OUTPUT, WHICH
!     INDICATES WHETHER ONE OR THE OTHER HAS BEEN OMITTED, IS MADE BELOW
 
 
 
!     FIRST GET FORCE VECTOR FOR THE PLATE CONSIDERATION
 
!     M ,  M ,  M  ,  V ,  V
!      X    Y    XY    X    Y
 
!                                NPTS
!     THE  5X1 FORCE VECTOR = SUMMATION  (S )(U )
!                                I=1       I   I
 
! *********************************************************************
 
!  . ZERO FORVEC AND PRECISION CHECK STORAGE...
 
 DO  i = 1,6
   forvec(i) = 0.0E0
   cfrvec(i) = 0.0E0
   cfrvec(i+6) = 0.0E0
 END DO
 forvec(1) = ph1out(1)
 
!     ZERO OUT LOCAL STRESSES
 
 sigx1 = 0.0E0
 sigy1 = 0.0E0
 sigxy1 = 0.0E0
 sigx2 = 0.0E0
 sigy2 = 0.0E0
 sigxy2 = 0.0E0
 
 IF( nsil(1) == 0 ) GO TO 30
 
!     FORM SUMMATION
 
 DO  i=1,npts
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nsil(i) - 1
   
   CALL smmats (si(30*i-29),5,6,0, z(npoint),6,1,0, vec(1),cvec(1) )
   DO  j=2,6
     cfrvec(j) = cfrvec(j) + cvec(j-1)
     forvec(j) = forvec(j) + vec(j-1)
   END DO
   
 END DO
 IF (strain) GO TO 220
 IF( tloads == 0 ) GO TO 25
 jst = 30*npts+8
 flag = .false.
 f1 = ti(6)
 IF( n1 == 1 ) GO TO 22
 forvec(2) = forvec(2) - ti(2)
 forvec(3) = forvec(3) - ti(3)
 forvec(4) = forvec(4) - ti(4)
 IF( ti(5) == 0.0 .AND. ti(6) == 0.0 ) flag = .true.
 GO TO 25
 22 forvec(2) = forvec(2) + ti(2)*ph1out(jst+1)
 forvec(3) = forvec(3) + ti(2)*ph1out(jst+2)
 forvec(4) = forvec(4) + ti(2)*ph1out(jst+3)
 IF( ti(3) == 0.0 .AND. ti(4) == 0.0 ) flag = .true.
 25 CONTINUE
 
!     FORCE VECTOR IS NOW COMPLETE
 
 z1 = ph1out(7)
 z2 = ph1out(8)
 
 z1ovri = - ph1out(7) / ph1out(6)
 z2ovri = - ph1out(8) / ph1out(6)
 z1i = ABS (z1ovri)
 z2i = ABS (z2ovri)
 
 k1 = 0
 ASSIGN 26 TO iretrn
 GO TO 170
 
 26 sigx1 = forvec(2) * z1ovri - sdelta(1)
 sigy1 = forvec(3) * z1ovri - sdelta(2)
 sigxy1 = forvec(4) * z1ovri - sdelta(3)
 cfrvec(7) = cfrvec(2) * z1i
 cfrvec(8) = cfrvec(3) * z1i
 cfrvec(9) = cfrvec(4) * z1i
 
 k1 = 1
 ASSIGN 27 TO iretrn
 GO TO 170
 
 27 sigx2 = forvec(2) * z2ovri - sdelta(1)
 sigy2 = forvec(3) * z2ovri - sdelta(2)
 sigxy2 = forvec(4) * z2ovri - sdelta(3)
 cfrvec(10) = cfrvec(2) * z2i
 cfrvec(11) = cfrvec(3) * z2i
 cfrvec(12) = cfrvec(4) * z2i
!     *******************************
 
 GO TO 40
 30 z1 = 0.0E0
 z2 = 0.0E0
 
!     FIND SIG X, SIG Y, SIG XY, FOR MEMBRANE CONSIDERATION
 40 IF( nph1ou(30*npts+13) == 0 ) GO TO 90
 
!     ZERO STRESS VECTOR STORAGE
 
 stress(1) = 0.0E0
 stress(2) = 0.0E0
 stress(3) = 0.0E0
 sstrss(1) = 0.0E0
 sstrss(2) = 0.0E0
 sstrss(3) = 0.0E0
 
!                            I=NPTS
!        STRESS VECTOR = (  SUMMATION(S )(U )  ) - (S )(LDTEMP - T )
!                            I=1       I   I         T            0
 
 DO  i=1,npts
   
!     POINTER TO I-TH SIL IN PH1OUT
   npoint = 30*npts + 12 + i
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   npoint = ivec + nph1ou (npoint) - 1
   
!     POINTER TO S SUB I 3X3
   npt1 = 30*npts + 12 + 9*i
   
   CALL smmats (ph1out(npt1),3,3,0, z(npoint),3,1,0, vec(1),cvec(1))
   DO  j=1,3
     sstrss(j) = sstrss(j) + cvec(j)
     stress(j) = stress(j) + vec(j)
   END DO
   
 END DO
 
 IF (strain) GO TO 230
 IF(ldtemp == (-1) ) GO TO 80
 
!     POINTER TO T SUB 0 = 30*NPTS + 17
 
 tem = ftemp - ph1out(30*npts + 17)
 DO  i=1,3
   npoint = 30*npts + 17 + i
   stress(i) = stress(i) -ph1out(npoint) *tem
 END DO
 
!     ADD MEMBRANE STRESSES TO PLATE STRESSES
 
 80 sigx1 = sigx1 + stress(1)
 sigy1 = sigy1 + stress(2)
 sigxy1 = sigxy1 + stress(3)
 sigx2 = sigx2 + stress(1)
 sigy2 = sigy2 + stress(2)
 sigxy2 = sigxy2 + stress(3)
 cfrvec( 7) = cfrvec( 7) + sstrss(1)
 cfrvec( 8) = cfrvec( 8) + sstrss(2)
 cfrvec( 9) = cfrvec( 9) + sstrss(3)
 cfrvec(10) = cfrvec(10) + sstrss(1)
 cfrvec(11) = cfrvec(11) + sstrss(2)
 cfrvec(12) = cfrvec(12) + sstrss(3)
 
!     STRESS OUTPUT VECTOR IS THE FOLLOWING
 
!      1) ELEMENT ID
!      2) Z1 = FIBER DISTANCE 1
!      3) SIG X  1
!      4) SIG Y  1
!      5) SIG XY 1
!      6) ANGLE OF ZERO SHEAR AT Z1
!      7) SIG P1 AT Z1
!      8) SIG P2 AT Z1
!      9) TAU MAX = MAXIMUM SHEAR STRESS AT Z1
 
!     10) ELEMENT ID
!     11) Z2 = FIBER DISTANCE 2
!     12) SIG X  2
!     13) SIG Y  2
!     14) SIG XY 2
!     15) ANGLE OF ZERO SHEAR AT Z2
!     16) SIG P1 AT Z2
!     17) SIG P2 AT Z2
!     S7) SIG P2 AT Z2
!     18) TAU MAX = MAXIMUM SHEAR STRESS AT Z2
 
 
 90 IF( nph1ou(2) == 0 .AND. nph1ou(30*npts+13) == 0 ) GO TO 120
 
!     COMPUTE PRINCIPAL STRESSES
 
 str( 1) = ph1out(1)
 str( 2) = z1
 str( 3) = sigx1
 str( 4) = sigy1
 str( 5) = sigxy1
 str(10) = ph1out(1)
 str(11) = z2
 str(12) = sigx2
 str(13) = sigy2
 str(14) = sigxy2
 
 DO  i=3,12,9
   temp = str(i) - str(i+1)
   str(i+6) = SQRT( (temp/2.0E0)**2 + str(i+2)**2 )
   delta = (  str(i)  +  str(i+1)  )  /  2.0E0
   str(i+4) = delta + str(i+6)
   str(i+5) = delta - str(i+6)
   delta = 2.0E0 * str(i+2)
   IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15)GO TO 100
   str(i+3) = ATAN2( delta,temp ) * 28.6478898E0
   CYCLE
   100 str(i+3) = 0.0E0
 END DO
 
 GO TO 140
 120 DO  i=2,18
   str(i) = 0.0E0
 END DO
 140 str(1) = ph1out(1)
 str(10) = ph1out(1)
 
 
!     ADDITION TO ELIMINATE 2ND ELEMENT ID IN OUTPUT
 
 DO  i=10,17
   str(i) = str(i+1)
 END DO
 IF (.NOT.strain) GO TO 152
 ist( 2) = iblank
 str( 5) = 2.0*str(5)
 str( 9) = 2.0*str(9)
 ist(10) = iblank
 str(13) = 2.0*str(13)
 str(17) = 2.0*str(17)
 152 CONTINUE
 
!  . STRESS CHECK...
 
 IF (nchk <= 0 ) GO TO 999
 cfrvec(1) = ph1out(1)
 k = 0
!  . FORCES...
 CALL sdrchk (forvec(2),cfrvec(2),5,k)
!  . STRESSES...
 CALL sdrchk (str(3),cfrvec( 7),3,k)
 CALL sdrchk (str(11),cfrvec(10),3,k)
 
 IF (k == 0) GO TO 999
 
!  . LIMITS EXCEEDED...
 j = 0
 istyp(1) = tr
 istyp(2) = ontw(1)
 IF (ieltyp > 17) istyp(1) = qu
 IF (IABS (ieltyp-17) < 2) istyp(2) = ontw(2)
 
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 156
 
 lsub = isub
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 1
 CALL page1
 154 CALL sd2rhd (ished,j)
 line = line + 1
 WRITE(nout,155)
 155 FORMAT (7X,51HTYPE     eid     mx     my    mxy     vx     vy    ,  &
     38HSX1    sy1   sxy1    sx2    sy2   sxy2)
 GO TO 157
 
 156 IF (eject(2) /= 0) GO TO 154
 
 157 WRITE(nout,158) istyp,ifrvec(1),(cfrvec(ii),ii=2,12)
 158 FORMAT (1H0,5X,a4,a2,i7,11F7.1)
 
 GO TO 999
 
!     INTERNAL SUBROUTINE
 
 170 IF( tloads == 0 .OR. flag ) GO TO 200
 jst = 30*npts + 8
 IF( n1 == 1 ) GO TO 190
 ff = ti(k1+5) - ti(1)
 sdelta(1) = (ph1out(jst+1)*ff + ti(2)*ph1out(k1+7)) / ph1out(6)
 sdelta(2) = (ph1out(jst+2)*ff + ti(3)*ph1out(k1+7)) / ph1out(6)
 sdelta(3) = (ph1out(jst+3)*ff + ti(4)*ph1out(k1+7)) / ph1out(6)
 GO TO 210
 190 ff = (ti(k1+3) - ph1out(k1+7)*ti(2) - ti(1))/ph1out(6)
 sdelta(1) = ph1out(jst+1)*ff
 sdelta(2) = ph1out(jst+2)*ff
 sdelta(3) = ph1out(jst+3)*ff
 GO TO 210
 200 sdelta(1) = 0.0
 sdelta(2) = 0.0
 sdelta(3) = 0.0
 210 GO TO iretrn,(26,27)
 
!     SPECIAL CALCULATIONS FOR STRAINS.
 
 220 sigx2 = forvec(2)
 sigy2 = forvec(3)
 sigxy2= forvec(4)
 GO TO 40
 
 230 sigx1 = stress(1)
 sigy1 = stress(2)
 sigxy1= stress(3)
 GO TO 90
 999 RETURN
END SUBROUTINE strqd2
