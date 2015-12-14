SUBROUTINE ascm11 (NAME,iphase,isol,nogo)
     
!     EXIO COMMANDS DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,7),subnam(2),rdmap(18,2),xtra(4), ptbs(7,12),isave(21)
 COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(126)
 DATA    comnd/ 4HSOFI    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HSOFO    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HREST    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HDUMP    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HCHEC    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HCOMP    ,  2    ,  4    ,  0    , 12    ,  0  ,  &
     4HAPPE    ,  2    ,  4    ,  0    , 12    ,  0  /
 DATA    slash/ 1H/       /
 DATA    isave/ 1, 7,1,  1,11,3,  1,13,2,  1,15,1,  2, 6,1,  2,11,3,  2,14,2/
 DATA    rdmap/  &
     4HEXIO,4H    ,4H  //,4HS,n,,4HDRY/,4HMACH,4H!*DE,4HVI*/,4H*UNI,  &
     4HTNAM,4HE*!*,4HFORM,4H*!*M,4HODE*,4H!*PO,4HSI*/,4H*ITE,4HM*/ ,  &
     4H    ,4H    ,4H  *n,4HAME0,4H001*,4H!*NA,4HME00,4H02*/,4H*NAM,  &
     4HE000,4H3*!*,4HNAME,4H0004,4H*!*N,4HAME0,4H005*,4H $  ,4H    /
 DATA    xtra / 4HMACH,4HPOSI,4HITEM,4HNAME /
 DATA    ptbs / 1  , 21  , 21  ,  4  ,   101  ,         0  ,  0  ,  &
     1  , 27  , 27  ,  4  ,   102  ,         0  ,  0  ,  &
     1  , 34  , 34  ,  8  ,   103  ,         0  ,  0  ,  &
     1  , 45  , 45  ,  4  ,   104  ,         0  ,  0  ,  &
     1  , 52  , 52  ,  4  ,   105  ,         0  ,  0  ,  &
     1  , 59  , 59  ,  4  ,   106  ,         0  ,  0  ,  &
     1  , 66  , 66  ,  4  ,   107  ,         0  ,  0  ,  &
     2  , 12  , 12  ,  8  ,   108  ,         0  ,  0  ,  &
     2  , 23  , 23  ,  8  ,   109  ,         0  ,  0  ,  &
     2  , 34  , 34  ,  8  ,   110  ,         0  ,  0  ,  &
     2  , 45  , 45  ,  8  ,   111  ,         0  ,  0  ,  &
     2  , 56  , 56  ,  8  ,   112  ,         0  ,  0  /
 DATA subnam  / 4HASCM,2H11  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))
 
 DO  l = 1,21,3
   i = isave(l+1)
   j = isave(l  )
   k = isave(l+2)
   rdmap(i,j) = khrfn1(rdmap(i,j),k,slash,1)
 END DO
 
!     VALIDATE COMMAND AND SET POINTERS
 
 DO  i = 1,7
   IF (NAME == comnd(1,i)) GO TO 20
 END DO
 GO TO 70
 20 icomnd = i
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
 
END SUBROUTINE ascm11
