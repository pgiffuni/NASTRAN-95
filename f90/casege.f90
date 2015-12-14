SUBROUTINE casege
     
    ! GENERATES IDENTICAL SUBCASES LMODES*NDIR TIMES FOR DDAM
 
    !     CASEGEN  CASECC/CASEDD/C,Y,LMODES/V,N,NDIR/V,N,NMODES $
    !    EQUIV CASEDD,CASECC  $
 
    INTEGER :: buf1,buf2,sysbuf,casecc,casedd
    DIMENSION iz(1),nam(2),mcb(7)
    COMMON/system/sysbuf
    COMMON/BLANK/lmodes,ndir,nmodes
    COMMON/zzzzzz/z(1)
    EQUIVALENCE (z(1),iz(1))
    DATA casecc,casedd/101,201/
    DATA nam/4HCASE,4HGE  /
 
    lcore=korsz(z)
    buf1=lcore-sysbuf+1
    buf2=buf1-sysbuf
    lcore=buf2-1
    IF(lcore <= 0)GO TO 1008
 
    CALL gopen(casecc,z(buf1),0)
    CALL gopen(casedd,z(buf2),1)
    CALL READ (*1002,*10,casecc,z,lcore,0,iwords)
    GO TO 1008
10  IF(lmodes > nmodes)lmodes=nmodes
    itot=lmodes*ndir
    DO  i=1,itot
        iz(1)=i
        CALL WRITE(casedd,z,iwords,1)
    END DO
    CALL CLOSE(casecc,1)
    CALL CLOSE(casedd,1)
    mcb(1)=casecc
    CALL rdtrl(mcb)
    mcb(1)=casedd
    mcb(2)=itot
    CALL wrttrl(mcb)
    RETURN
 
1002 CALL mesage(-2,casecc,nam)
1008 CALL mesage(-8,0,nam)

    RETURN
END SUBROUTINE casege
