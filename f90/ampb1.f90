SUBROUTINE ampb1(ipvct,noh,noe)
     
    !     THIS ROUTINE BUILDS A PARTITIONING VECTOR WHICH WILL APPEND NOE
    !       TERM(OR COLUMNS)
 
 
    INTEGER, INTENT(IN OUT)                  :: ipvct
    INTEGER, INTENT(IN)                      :: noh
    INTEGER, INTENT(IN)                      :: noe
    INTEGER :: sysbuf,mcb(7)
 
    COMMON  /zblpkx/a(4),ii
    COMMON /system/sysbuf
    COMMON /zzzzzz/ z(1)
 
    !-----------------------------------------------------------------------
 
    ibuf1=korsz(z)-sysbuf+1
    CALL gopen(ipvct,z(ibuf1),1)
    CALL makmcb(mcb,ipvct,noh+noe,2,1)
    CALL bldpk(1,1,ipvct,0,0)
    ii=noh
    DO  i=1,noe
        a(1)=1.0
        ii=ii+1
        CALL zblpki
    END DO
    CALL bldpkn(ipvct,0,mcb)
    CALL CLOSE(ipvct,1)
    CALL wrttrl(mcb)

    RETURN
END SUBROUTINE ampb1
