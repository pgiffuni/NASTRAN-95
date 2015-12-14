SUBROUTINE combo (cdata,nx,extra,nnam,NAME,nn,var,ier)
     
    !     THIS ROUTINE  PROCESSES THE  COMBINE INPUT.
    !        THE  INPUT/ OUTPUTS  ARE
 
    !                  CDATA  -  XRCARD  IMAGE OF  COMBINE CARD  (IN)
    !                  NX     -  NUMBER OF EXTRAS                (IN)
    !                  EXTRA  -  3 BY NX ARRAY OF EXTRAS         (IN)
    !                  NNAM   -  NUMBER OF CURRENT SUBS NAMES (IN/OUT)
    !                  NAMES  -  ARRAY OF  CURRENT SUBS NAMES (IN/OUT)
    !                  NN     -  NUMBER OF SUBS TO BE COMBINED  (OUT)
    !                  VAR    -  3 BY NVAR ARRAY OF VARIABLES   (OUT)
    !                            ARRANGED AS- KEY WORD + 2 DATA WORDS
 
 
    INTEGER, INTENT(IN OUT)                  :: cdata(5)
    INTEGER, INTENT(IN)                      :: nx
    INTEGER, INTENT(IN OUT)                  :: extra(3,1)
    INTEGER, INTENT(IN OUT)                  :: nnam
    INTEGER, INTENT(IN OUT)                  :: NAME(2,1)
    INTEGER, INTENT(OUT)                     :: nn
    INTEGER, INTENT(OUT)                     :: var( 3,2)
    INTEGER, INTENT(OUT)                     :: ier
    EXTERNAL rshift   ,complf
    INTEGER :: rshift   ,complf        ,eqsn
 
    DIMENSION   inum(7)   ,numbs(7) ,mopt(3)   ,msort(3)  ,nai(7)
 
    DATA   inum   / 4HN1  , 4HN2  ,4HN3  ,4HN4  ,4HN5  ,4HN6  ,4HN7  /
    DATA   lprn,    nopt,   nsort, mopt                ,msort        /  &
        4H(    , 4HOPTS, 4HSORT,4HAUTO,4HMAN ,4HREST,4HX   ,4HY   , 4HZ    /
    DATA   manu   / 4HMANU/
    DATA   nno    / 4HNAME/ ,nnc /4HNAMC/, nams /4HNAMS/
    DATA   nai    / 4HNA1 ,4HNA2 ,4HNA3  ,4HNA4 ,4HNA5 ,4HNA6 ,4HNA7 /
    DATA   ncno   / 4HNCNO/
 
    DATA   eqsn   / 4H=   /
 
 
    lword  =  rshift( complf(0),1)
    ier    = 0
    !     COMBINE OPERATION
    !         PROCESS PRIMARY CARD -COMBINE( OPTS,SORT) = NAME1,NAME2, ETC
    !     SET DEFAULTS
    DO  i =1,150
        var(i,1) = 0
    END DO
    jnam = 6
    var(1,1) = nopt
    var(2,1) = mopt(1)
    var(1,2) = nsort
    var(2,2) = msort(1)
    IF( cdata(5) /= lprn) GO TO 1220
    k = 6
 
    !     PROCESS  AUTO/MAN  OR XYZ
 
    1211 DO  i =1,3
        IF (  cdata(k) /= mopt(i)) GO TO 1212
        var(2,1) = mopt(i)
        GO TO 1216
1212    IF ( cdata(k) /= manu) GO TO  1213
        var(2,1) = mopt(2)
        GO TO 1216
1213    IF (  cdata(k) /= msort(i)) CYCLE
        var(2,2) = msort(i)
        GO TO 1216
    END DO
    !     NOT VALID    ASSUME EQ SIGN OR NAME
 
    GO TO 1222
1216 k  = k+2
    GO TO 1211
 
    !     NO  OPTION
1220 k = 4
 
    !     CHECK FOR  EQ SIGN
1222 IF ( cdata( k+1) == eqsn)  k =k+2
 
    !     PROCESS NAMES
    nn = 0
    DO  i = 1,7
        kn = k + 2*i -2
        IF ( cdata( kn) == lword) GO TO 1236
   
        var(1,i+2) = nams
        var(2,i+2) = cdata(kn)
        var(3,i+2) = cdata(kn+1)
   
        !     FIND STRUCTURE NUMBER
        IF ( nnam == 0 ) GO TO 1231
        DO  j =1, nnam
            IF ( cdata(kn) /= NAME(1,j) .OR. cdata(kn+1) /= NAME(2,j)) CYCLE
            numbs(i) = j
            GO TO 1232
        END DO
   
        !     NEW NAME
   
1231    nnam  = nnam +1
        numbs(i) = nnam
        NAME(1,nnam) = cdata(kn)
        NAME(2,nnam) = cdata(kn+1)
1232    nn= nn+1
    END DO
 
 
    !     MOVE  EXTRAS INTO PLACE  CHANGE NAME TO NAMC
1236 ic = 0
    DO   j = 1,nx
        ix  = j +3*nn  +2
        IF ( extra(1,j) /= nno ) GO TO 1238
        extra(1,j) = nnc
        ic = ix
        1238 DO  k = 1,3
            var( k,ix) = extra(k,j)
        END DO
    END DO
 
    !     SET  STRUCTURE NUMBER KEYS
 
    IF( nn == 0) GO TO 1248
 
    DO   i = 1, nn
   
        ix =  i + nn  +2
        var(1,ix) = inum(i)
        var(2,ix) = -1
        var(3,ix) = numbs(i)
        iy = ix+nn
        var(1,iy) = nai(i)
        var(2,iy) = var(2,i+2)
        var(3,iy) = var(3,i+2)
    END DO
    GO TO 1250
1248 ier = 1
 
    !     CHECK  FOR NAMC AS A PREVIOUS NAME  OR MISSING
1250 IF ( ic == 0) GO TO 1265
    DO  j =1,nnam
        IF (var(2,ic) /= NAME(1,j).OR.var(3,ic) /= NAME(2,j)) CYCLE
        GO TO 1265
    END DO
 
    !     OK -NEW NAME , ADD TO LIST
 
    nnam = nnam+1
    NAME(1,nnam) = var(2,ic)
    NAME(2,nnam) = var(3,ic)
    ix = nx+3*nn+3
    var(1,ix) = ncno
    var(2,ix) = -1
    var(3,ix) = nnam
    RETURN
1265 ier = ier +2

    RETURN
END SUBROUTINE combo
