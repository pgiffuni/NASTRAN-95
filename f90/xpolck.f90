SUBROUTINE xpolck (dbn1,dbn2,fn,l)
     
 
 INTEGER, INTENT(IN OUT)                  :: dbn1
 INTEGER, INTENT(IN)                      :: dbn2
 INTEGER, INTENT(OUT)                     :: fn
 INTEGER, INTENT(OUT)                     :: l
 INTEGER          :: outtap
 EXTERNAL        lshift   ,  andf    ,   orf
 DIMENSION       npolck(2),  ddbn( 1),   dfnu( 1),  fcum( 1),  &
     fcus(  1),  fdbn( 1),   fequ( 1),  FILE( 1),  &
     fknd(  1),  fmat( 1),   fntu( 1),  fpun( 1),  &
     fon (  1),  ford( 1),   minp( 1),  mlsn( 1),  &
     mout(  1),  mscr( 1),   sal ( 1),  sdbn( 1), sntu(  1),  sord(  1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ksystm(65)
 COMMON /xfiat / fiat(7)
 COMMON /xfist / fist
 COMMON /xdpl  / dpd(6)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(5)
 EQUIVALENCE    (ksystm(2),outtap  )
 EQUIVALENCE               (dpd  (1),dnaf    ),(dpd  (2),dmxlg   ),  &
     (dpd  (3),dculg   ),(dpd  (4),ddbn (1)),(dpd  (6),dfnu (1)),  &
     (fiat (1),funlg   ),(fiat (2),fmxlg   ),(fiat (3),fculg   ),  &
     (fiat (4),fequ (1)),(fiat (4),FILE (1)),(fiat (4),ford (1)),  &
     (fiat (5),fdbn (1)),(fiat (7),fmat (1)),(md   (1),mlgn    ),  &
     (md   (2),mlsn (1)),(md   (3),minp (1)),(md   (4),mout (1)),  &
     (md   (5),mscr (1)),(sos  (1),slgn    ),(sos  (2),sdbn (1)),  &
     (sos  (4),sal  (1)),(sos  (4),sntu (1)),(sos  (4),sord (1)),  &
     (xf1at(1),fntu (1)),(xf1at(1),fon  (1)),(xf1at(2),fpun (1)),  &
     (xf1at(3),fcum (1)),(xf1at(4),fcus (1)),(xf1at(5),fknd (1))
 EQUIVALENCE               (comm (1),almsk   ),(comm (2),apndmk  ),  &
     (comm (3),cursno  ),(comm (4),entn1   ),(comm (5),entn2   ),  &
     (comm (6),entn3   ),(comm (7),entn4   ),(comm (8),flag    ),  &
     (comm (9),fnx     ),(comm(10),lmsk    ),(comm(11),lxmsk   ),  &
     (comm(12),macsft  ),(comm(13),rmsk    ),(comm(14),rxmsk   ),  &
     (comm(15),s       ),(comm(16),scornt  ),(comm(17),tapmsk  ),  &
     (comm(18),thcrmk  ),(comm(19),zap     )
 DATA   pool    /4HPOOL       /
 DATA   npolck  /4HXPOL,4HCK  /
 
 
!     XPOLCK CHECKS THE DATA POOL DICT FOR A DATA BLOCK NAME
 
 lmt1 = dculg*entn4
 DO  i = 1,lmt1,entn4
   IF (dbn1 /= ddbn(i) .OR. dbn2 /= ddbn(i+1)) CYCLE
   fn = andf(rmsk,dfnu(i))
   l  = i
   RETURN
   
 END DO
 fn = 0
 RETURN
 
 
 ENTRY xfilps (NEW)
!     ==================
 
!     XFILPS POSITIONS THE POOL TAPE FORWARD OR BACKWARD
 
!     NEW IS DESIRED POSITION
!     FNX IS CURRENT POSITION
 
 fdif = NEW - fnx
 IF (fdif < 0.0) THEN
   GO TO    10
 ELSE IF (fdif == 0.0) THEN
   GO TO    30
 ELSE
   GO TO    20
 END IF
 10 fdif = fdif - 1
 20 CALL skpfil (pool,fdif)
 IF (fdif < 0 .AND. NEW /= 1) CALL skpfil (pool,+1)
 30 RETURN
 
 
 ENTRY xpleqk (nx,ny)
!     ====================
 
!     XPLEQK MOVES SECONDARY EQUIVALENCED DATA BLOCK NAMES FROM THE
!     POOL DICT. TO FIAT.   NTU-LTU DATA ARE ALSO STORED IN FIAT FOR THE
!     EQUIV. D.B.   NTU-LTU DATA IS EXTRACTED FROM SOS IF FOUND, IF NOT,
!     IT IS COPIED FROM THE CALLING PRIMARY D.B.
 
!     NX IS THE POOL DICT. INDEX
!     NY IS THE FIAT INDEX FOR PRIMARY D.B.
 
 fequ(ny) = orf(s,fequ(ny))
 kfil = andf(rmsk,dfnu(nx))
 lmt1 = dculg*entn4
 lmt2 = slgn *entn2
 lmt3 = fculg*entn1
 nfculg = lmt3 + 1
 
!     SEARCH FOR EQUIV FILES IN DICT
 
 DO  i = 1,lmt1,entn4
   IF (ddbn(i) == 0 .AND. ddbn(i+1) == 0) CYCLE
   IF (kfil /= andf(rmsk,dfnu(i))) CYCLE
   IF (i == nx) CYCLE
   
!     SEE IF NAME IS IN FIAT
   
   DO  j = 1,lmt3,entn1
     IF (ddbn(i) == fdbn(j) .AND. ddbn(i+1) == fdbn(j+1)) GO TO 115
   END DO
   fdbn(nfculg  ) = ddbn(i  )
   fdbn(nfculg+1) = ddbn(i+1)
   FILE(nfculg) = FILE(ny)
   fntu(nfculg) = fntu(ny)
   ford(nfculg) = orf(lshift(1000,16),andf(rmsk,FILE(nfculg)))
   fequ(nfculg) = orf(s,fequ(nfculg))
   DO  j = 1,lmt2,entn2
     IF (ddbn(i) == sdbn(j) .AND. ddbn(i+1) == sdbn(j+1)) GO TO 120
   END DO
   GO TO 140
   
!     FILE ALREADY ALLOCATED  BE SURE EQUIVED
   
   115 FILE(j) = orf(andf(rmsk,FILE(ny)),andf(lmsk,FILE(j)))
   fequ(j) = orf(s,fequ(j))
   CYCLE
   120 ford(nfculg) = orf(andf(lmsk,sord(j)),andf(rmsk,FILE(nfculg)))
   fequ(nfculg) = orf(s,fequ(nfculg))
   fntu(nfculg) = sntu(j)
   140 nfculg = nfculg+ entn1
   fculg  = fculg + 1
   
!     FLAG INDICATES D.B. S HAVE BEEN ADDED TO FIAT
   
   flag = -1
   IF (fculg > fmxlg) GO TO 900
 END DO
 RETURN
 
 900 WRITE  (outtap,901) sfm
 901 FORMAT (a25,' 1051, FIAT OVERFLOW')
 CALL mesage (-37,0,npolck)
 RETURN
END SUBROUTINE xpolck
