SUBROUTINE amgt1a (INPUT,matout,ajj,tsonx,tamach,tredf,nstns2)
     
    !     COMPUTE AJJ MATRIX FOR SWEPT TURBOPROP BLADES.
 
 
    INTEGER, INTENT(IN OUT)                  :: INPUT
    INTEGER, INTENT(IN OUT)                  :: matout
    COMPLEX, INTENT(IN OUT)                  :: ajj(nstns2,1)
    INTEGER, INTENT(OUT)                     :: tsonx(1)
    REAL, INTENT(OUT)                        :: tamach(1)
    REAL, INTENT(OUT)                        :: tredf(1)
    INTEGER, INTENT(IN)                      :: nstns2
    LOGICAL :: tsonic,debug
    INTEGER :: sln,NAME(2)
    REAL :: minmac,maxmac,mach
 
 
    COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq
    COMMON /condas/ pi,twopi,radeg,degra,s4pisq
    COMMON /packx / iti,ito,ii,nn,incr
    COMMON /tamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
        refmac,refden,refvel,refswp,sln,nstnsx,stager,  &
        chord,dcbdzb,bspace,mach,den,vel,sweep,amach, redf,blspc,amachr,tsonic
    COMMON /amgbug/ debug
    DATA    NAME  / 4HAMGT,4H1A  /
 
    !     LOOP ON STREAMLINES, COMPUTE AJJ FOR EACH STREAMLINE AND THEN
    !     PACK AJJ INTO AJJL MATRIX AT CORRECT POSITION
 
    ii = 0
    nn = 0
    nstns3 = 3*nstns
    DO  line = 1,nlines
   
        !     READ STREAMLINE DATA (SKIP COORDINATE DATA)
   
        CALL READ (*400,*400,INPUT,sln,10,0,nwar)
   
        !     COMPUTE PARAMETERS
   
        amach = mach
        redf  = rfreq*(chord/refcrd)*(refvel/vel)
        blspc = bspace/chord
   
        !     COMPUTE C3 AND C4 FOR THIS STREAMLINE.
   
        !     INPUT IS POSITIONED AT THE FIRST 10 WORDS OF THE NEXT
        !     STREAMLINE WHEN IT RETURNS FROM AMGT1T
   
        CALL amgt1t (nlines,line,INPUT,nstns,c3,c4)
   
        IF (debug) CALL bug1 ('TAMG1L    ',5,iref,26)
   
        !     COMPUTE POINTER FOR LOCATION INTO AJJ MATRIX
   
        iajjc = 1
        IF (tsonic) iajjc = nstns2*(line-1) + 1
   
        !     BRANCH TO SUBSONIC, SUPERSONIC OR TRANSONIC CODE
   
        tamach(line) = amach
        tredf(line)  = redf
        IF (amach <= maxmac) GO TO 10
        IF (amach >= minmac) GO TO 20
   
        !     TRANSONIC STREAMLINE. STORE DATA FOR TRANSONIC INTERPOLATION
   
        tsonx(line) = iajjc
        CYCLE
   
        !     SUBSONIC STREAMLINE
   
10      CALL amgt1b (ajj(1,iajjc),nstns2,c3,c4)
        GO TO 30
   
        !     SUPERSONIC STREAMLINE
   
20      CALL amgt1c (ajj(1,iajjc),nstns2,c3,c4)
30  CONTINUE
   
    !     IF THERE ARE NO TRANSONIC STREAMLINES OUTPUT THIS AJJ SUBMATRIX
   
    IF (tsonic) GO TO 60
    ii = nn + 1
    nn = nn + nstns2
   
    !     OUTPUT AJJ MATRIX
   
    DO  i = 1,nstns2
        IF (debug) CALL bug1 ('SS-AJJL   ',40,ajj(1,i),nstns2*2)
        CALL pack (ajj(1,i),matout,mcb)
    END DO
    CYCLE
60  tsonx(line) = 0
END DO
 
!     PERFORM TRANSONIC INTERPOLATION, IF NECESSARY
 
IF (.NOT.tsonic) GO TO 300
IF (debug) CALL bug1 ('TSONX     ', 80,tsonx,nlines)
IF (debug) CALL bug1 ('TAMACH    ', 90,tamach,nlines)
IF (debug) CALL bug1 ('TREDF     ',100,tredf,nlines)
CALL amgt1d (ajj,tsonx,tamach,tredf,nstns2)
 
!     OUTPUT AJJ FOR EACH STREAMLINE
 
DO  nline = 1,nlines
    ii = nn + 1
    nn = nn + nstns2
    DO  i = ii,nn
        IF (debug) CALL bug1 ('STS-AJJL  ',110,ajj(1,i),nstns2*2)
        CALL pack (ajj(1,i),matout,mcb)
    END DO
END DO
300 RETURN
 
!     ERROR MESSAGES
 
!     INPUT NOT POSITIONED PROPERLY OR INCORRECTLY WRITTEN
 
400 CALL mesage (-7,0,NAME)

    RETURN
END SUBROUTINE amgt1a
