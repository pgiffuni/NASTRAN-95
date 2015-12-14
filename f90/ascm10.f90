SUBROUTINE ascm10 (NAME,iphase,isol,nogo)
     
    !     SUBSTRUCTURE UTILITY COMMANDS DMAP DATA
 
 
    INTEGER, INTENT(IN OUT)                  :: NAME
    INTEGER, INTENT(IN OUT)                  :: iphase
    INTEGER, INTENT(IN OUT)                  :: isol
    INTEGER, INTENT(OUT)                     :: nogo
    INTEGER :: comnd(6,6),xtra(1),subnam(2),rdmap(18,2),ptbs(7,10)
    COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(109)
    DATA    comnd/ 4HDEST    ,  2    ,  0    ,  0    ,  3    ,  0  ,  &
        4HEDIT    ,  2    ,  0    ,  0    ,  3    ,  0  ,  &
        4HEQUI    ,  2    ,  1    ,  0    ,  5    ,  0  ,  &
        4HSOFP    ,  2    ,  0    ,  0    , 10    ,  0  ,  &
        4HDELE    ,  2    ,  0    ,  0    , 10    ,  0  ,  &
        4HRENA    ,  2    ,  0    ,  0    ,  5    ,  0  /
    DATA    slash/ 1H/       /
    DATA    rdmap/  &
        4HSOFU,4HT   ,4H  //,4HDRY/,4H*nam,4HE   ,4H *!*,4HOPER,4H*/OP,  &
    4HT!*N,4HAME0,4H002*,4H!*PR,4HEF*/,4H*ITM,4H1*!*,4HITM2,4H*/  ,  &
    4H    ,4H    ,4H  *i,4HTM3*,4H!*IT,4HM4*/,4H*ITM,4H5* $,4H    ,  &
    4H    ,8*4H          /
    DATA    xtra / 4HPREF /
    DATA    ptbs / 1  , 16  , 18  ,  8  ,4HNAME  ,         0  ,  0  ,  &
        1  , 27  , 29  ,  4  ,4HOPER  ,         0  ,  0  ,  &
        1  , 34  , 35  ,  3  ,4HOPTI  ,         0  ,  0  ,  &
        1  , 38  , 40  ,  8  ,4HNEW   ,         0  ,  0  ,  &
        1  , 49  , 51  ,  4  ,4HPREF  ,         0  ,  0  ,  &
        1  , 56  , 58  ,  4  ,4HITM1  ,         0  ,  0  ,  &
        1  , 63  , 65  ,  4  ,4HITM2  ,         0  ,  0  ,  &
        2  , 11  , 12  ,  4  ,4HITM3  ,         0  ,  0  ,  &
        2  , 17  , 19  ,  4  ,4HITM4  ,         0  ,  0  ,  &
        2  , 24  , 26  ,  4  ,4HITM5  ,         0  ,  0  /
    DATA   subnam / 4HASCM,2H10    /
 
    !     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
    !     (SEE ASCM01 FOR EXPLANATION))
 
    rdmap(7, 1) = khrfn1(rdmap(7, 1),3,slash,1)
    rdmap(10,1) = khrfn1(rdmap(10,1),2,slash,1)
    rdmap(13,1) = khrfn1(rdmap(13,1),1,slash,1)
    rdmap(16,1) = khrfn1(rdmap(16,1),3,slash,1)
    rdmap(5, 2) = khrfn1(rdmap(5, 2),1,slash,1)
 
    !     VALIDATE COMMAND AND SET POINTERS
 
    DO  i = 1,6
        IF (NAME == comnd(1,i)) GO TO 20
    END DO
    GO TO 70
20  icomnd = i
    irdm   = 1
    nrdm   = comnd(2,icomnd)
    ixtra  = irdm  + 18*nrdm
    nxtra  = comnd(3,icomnd)
    ioct   = ixtra + nxtra
    noct   = comnd(4,icomnd)
    iptbs  = ioct  + 3*noct
    nptbs  = comnd(5,icomnd)
    iph    = iptbs + 7*nptbs
    nph    = comnd(6,icomnd)
 
    !     MOVE RDMAP DATA
 
    k = 0
    IF (nrdm == 0) GO TO 35
    DO  j = 1,nrdm
        DO  i = 1,18
            k = k + 1
            idat(k) = rdmap(i,j)
        END DO
    END DO
35 CONTINUE
 
   !     MOVE XTRA DATA
 
   IF (nxtra == 0) GO TO 45
   DO  i = 1,nxtra
       k = k + 1
       idat(k) = xtra(i)
   END DO
45 CONTINUE
 
   !     MOVE PTBS DATA
 
   IF (nptbs == 0) GO TO 65
   DO  j = 1,nptbs
       DO  i = 1,7
           k = k + 1
           idat(k) = ptbs(i,j)
       END DO
   END DO
65 CONTINUE
 
   RETURN
 
   !     INPUT ERROR
 
70 CALL mesage (7,0,subnam)
   nogo = 1

   RETURN
 END SUBROUTINE ascm10
