SUBROUTINE ascm12 (NAME,iphase,isol,nogo)
     
    !     PLOT COMMAND DMAP DATA
 
 
    INTEGER, INTENT(IN OUT)                  :: NAME
    INTEGER, INTENT(IN OUT)                  :: iphase
    INTEGER, INTENT(IN OUT)                  :: isol
    INTEGER, INTENT(OUT)                     :: nogo
    INTEGER :: comnd(6,1),subnam(2),rdmap(18,6),ptbs(7,15)
    COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(213)
    DATA    comnd/ 4HPLOT,  6,  0,  0, 15,  0   /
    DATA    rdmap/  &
        4HPLTM,4HRG  ,4H  ca,4HSECC,4H,pcd,4HB/pl,4HTSTP,4H,gps,4HTP,e,  &
        4HLSTP,4H,bgs,4HTP,c,4HASST,4HP,EQ,4HSTP/,4H*nam,4HE   ,4H */ ,  &
        4H    ,4H    ,4H  s,,4HN,ng,4HP/s,,4HN,ls,4HIL/s,4H,n,n,4HPSET,  &
        4H $  , 8*4H    ,  &
        4HSETV,4HAL  ,4H  //,4HS,n,,4HPLTF,4HLG/1,4H/s,n,4H,pfi,4HL/0 ,  &
        4H$   , 8*4H    ,  &
        4HPLOT,4H    ,4H  pl,4HTSTP,4H,gps,4HTP,e,4HLSTP,4H,cas,4HSTP,,  &
        4HBGST,4HP,EQ,4HSTP,,4H,,,,,4H,,/p,4HMSTP,4H/ngp,4H/lsi,4HL/  ,  &
        4H    ,4H    ,4H  s,,4HN,np,4HSET/,4HS,n,,4HPLTF,4HLG/s,4H,n,p,  &
        4HFIL ,4H$   , 7*4H    , 4HPRTM,4HSG  ,4H  pm,4HSTP/,4H/ $ ,13*4H        /
    DATA    ptbs / 1  , 26  , 26  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 32  , 32  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 38  , 38  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 44  , 44  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 51  , 51  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 57  , 57  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        1  , 62  , 62  ,  8  ,4HNAME  ,         0  ,  0  ,  &
        4  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 20  , 20  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 26  , 26  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 33  , 33  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 39  , 39  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 45  , 45  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        4  , 58  , 58  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
        6  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  /
    DATA  subnam  / 4HASCM,2H12  /
 
    !     VALIDATE COMMAND AND SET POINTERS
 
    IF (NAME /= comnd(1,1)) GO TO 70
    icomnd = 1
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
END SUBROUTINE ascm12
