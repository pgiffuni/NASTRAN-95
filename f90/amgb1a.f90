SUBROUTINE amgb1a (INPUT,matout,ajj,ajjt,tsonx,tamach,tredf)
     
    !     COMPUTE AJJ MATRIX FOR COMPRESSOR BLADES
 
 
    INTEGER, INTENT(IN OUT)                  :: INPUT
    INTEGER, INTENT(IN OUT)                  :: matout
    COMPLEX, INTENT(IN OUT)                  :: ajj(nstns,1)
    COMPLEX, INTENT(IN OUT)                  :: ajjt(nstns)
    INTEGER, INTENT(OUT)                     :: tsonx(1)
    REAL, INTENT(OUT)                        :: tamach(1)
    REAL, INTENT(OUT)                        :: tredf(1)
    LOGICAL :: tsonic,debug
    INTEGER :: sln,NAME(2)
    REAL :: minmac,maxmac,mach
 
 
    COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq
    COMMON /condas/ pi,twopi,radeg,degra,s4pisq
    COMMON /packx / iti,ito,ii,nn,incr
    COMMON /bamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
        refmac,refden,refvel,refflo,sln,nstnsx,stager,  &
        chord,radius,bspace,mach,den,vel,flowa,amach, redf,blspc,amachr,tsonic
    COMMON /amgbug/ debug
    DATA    NAME  / 4HAMGB,4H1A  /
 
    !     LOOP ON STREAMLINES, COMPUTE AJJ FOR EACH STREAMLINE AND THEN
    !     PACK AJJ INTO AJJL MATRIX AT CORRECT POSITION
 
    ii = 0
    nn = 0
    nstns3 = 3*nstns
    DO  line = 1,nlines
   
        !     READ STREAMLINE DATA (SKIP COORDINATE DATA)
   
        CALL READ (*999,*999,INPUT,sln,10,0,nwar)
        CALL READ (*999,*999,INPUT,0,-nstns3,0,nwar)
   
        !     COMPUTE PARAMETERS
   
        amach = mach*COS(degra*(flowa-stager))
        redf  = rfreq*(chord/refcrd)*(refvel/vel)*(mach/amach)
        blspc = bspace/chord
        IF (debug) CALL bug1 ('BAMG1L    ',5,iref,26)
   
        !     COMPUTE POINTER FOR LOCATION INTO AJJ MATRIX
   
        iajjc = 1
        IF (tsonic) iajjc = nstns*(line-1) + 1
   
        !     BRANCH TO SUBSONIC, SUPERSONIC OR TRANSONIC CODE
   
        tamach(line) = amach
        tredf(line)  = redf
        IF (amach <= maxmac) GO TO 10
        IF (amach >= minmac) GO TO 20
   
        !     TRANSONIC STREAMLINE. STORE DATA FOR TRANSONIC INTERPOLATION
   
        tsonx(line) = iajjc
        CYCLE
   
        !     SUBSONIC STREAMLINE
   
10      CALL amgb1b (ajj(1,iajjc))
        GO TO 30
   
        !     SUPERSONIC STREAMLINE
   
20      CALL amgb1c (ajj(1,iajjc))
30  CONTINUE
   
    !     IF THERE ARE NO TRANSONIC STREAMLINES OUTPUT THIS AJJ SUBMATRIX
   
    IF (tsonic) GO TO 60
    ii = nn + 1
    nn = nn + nstns
   
    !     OUTPUT AJJ MATRIX
   
    DO  i = 1,nstns
        IF (debug) CALL bug1 ('SS-AJJL   ',40,ajj(1,i),nstns*2)
        CALL pack (ajj(1,i),matout,mcb)
    END DO
    CYCLE
60  tsonx(line) = 0
END DO
 
!     PERFORM TRANSONIC INTERPOLATION, IF NECESSARY
 
IF (.NOT.tsonic) GO TO 300
IF (debug) CALL bug1 ('TSONX     ',102,tsonx,nlines)
IF (debug) CALL bug1 ('TAMACH    ',103,tamach,nlines)
IF (debug) CALL bug1 ('TREDF     ',104,tredf,nlines)
CALL amgb1d (ajj,tsonx,tamach,tredf)
 
!     OUTPUT AJJ FOR EACH STREAMLINE
 
DO  nline = 1,nlines
    ii = nn + 1
    nn = nn + nstns
    DO  i = ii,nn
        IF (debug) CALL bug1 ('STS-AJJL  ',110,ajj(1,i),nstns*2)
        CALL pack (ajj(1,i),matout,mcb)
    END DO
END DO
300 RETURN
 
!     ERROR MESSAGES
 
!     INPUT NOT POSITIONED PROPERLY OR INCORRECTLY WRITTEN
 
999 CALL mesage (-7,0,NAME)

    RETURN
END SUBROUTINE amgb1a
