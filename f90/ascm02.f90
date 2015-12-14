SUBROUTINE ascm02 (NAME,iphase,isol,nogo)
     
    !     RUN COMMAND DATA
 
 
    INTEGER, INTENT(IN OUT)                  :: NAME
    INTEGER, INTENT(IN OUT)                  :: iphase
    INTEGER, INTENT(IN OUT)                  :: isol
    INTEGER, INTENT(OUT)                     :: nogo
    INTEGER :: comnd(6,2),subnam(2),rdmap(18,6),ptbs(7,1)
    COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(117)
    DATA    comnd/ 4HRUN     ,  1    ,  0    ,  0    ,  1    ,  0  ,  &
        4HENDD    ,  6    ,  0    ,  0    ,  0    ,  0  /
    DATA    rdmap/  &
        4HPARA,4HM   ,4H  //,4H*add,4H*/dr,4HY/-1,4H /0 ,4H$   ,4H    ,  &
        4H    , 8*4H    , 4HLABE,4HL   ,4H  lb,4HSEND,4H $  ,13*4H    ,  &
        4HPARA,4HM   ,4H  //,4H*add,4H*/dr,4HY/dr,4HY/1 ,4H$   ,4H    ,  &
        4H    , 8*4H    , 4HCOND,4H    ,4H  fi,4HNIS,,4HDRY ,4H$   ,12*4H    ,  &
        4HREPT,4H    ,4H  lb,4HSBEG,4H,1 $,13*4H    ,  &
        4HJUMP,4H    ,4H  fi,4HNIS ,4H$   ,13*4H           /
    DATA    ptbs / 1,  22, 23, 3,  4HRUN ,  0,  0       /
 
    DATA   subnam/ 4HASCM,2H02  /
 
    !     VALIDATE COMMAND AND SET POINTERS
 
    DO  i = 1, 2
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
    IF (nrdm == 0) GO TO 40
    DO  j = 1,nrdm
        DO  i = 1,18
            k = k + 1
            idat(k) = rdmap(i,j)
        END DO
    END DO
40 CONTINUE
 
   !     MOVE PTBS DATA
 
   IF (nptbs == 0) GO TO 60
   DO  j = 1,nptbs
       DO  i = 1,7
           k = k + 1
           idat(k) = ptbs(i,j)
       END DO
   END DO
60 CONTINUE
 
   RETURN
 
   !     INPUT ERROR
 
70 CALL mesage (7,0,subnam)
   nogo = 1

   RETURN
END SUBROUTINE ascm02
