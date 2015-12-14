SUBROUTINE ampc2 (inp,outp,scrf)
     
    !     THE PURPOSE OF THIS ROUTINE IS TO COPY SCR5 ONTO THE BOTTOM OF
    !     OUTPUT
 
 
    INTEGER, INTENT(IN)                      :: inp
    INTEGER, INTENT(IN)                      :: outp
    INTEGER, INTENT(IN OUT)                  :: scrf
    INTEGER :: sysbuf,mcbi(7),mcbo(7)
    COMMON /packx / it1,it2,ii,jj,incr
    COMMON /unpakx/ it3,ii1,jj1,incr1
    COMMON /system/ sysbuf
    COMMON /zzzzzz/ iz(1)
    COMMON /TYPE  / isk(2),iword(4)
 
 
    mcbi(1) = inp
    CALL rdtrl (mcbi)
    mcbo(1) = outp
    CALL rdtrl (mcbo)
 
    !     IS THIS THE FIRST ENTRY
 
    IF (mcbo(2) /= 0) GO TO 10
 
    !     SWITCH SCRATCH FILES
 
    CALL filswi (inp,outp)
    RETURN
 
    !     MUST DO COPY
 
10  CALL filswi (outp,scrf)
    ibuf1 = korsz(iz) - sysbuf + 1
    ibuf2 = ibuf1 - sysbuf
    ibuf3 = ibuf2 - sysbuf
    CALL gopen (inp,iz(ibuf1),0)
    CALL gopen (scrf,iz(ibuf2),0)
    CALL gopen (outp,iz(ibuf2),1)
    ncol  = mcbi(2)
    nrowo = mcbi(3) + mcbo(3)
    it1   = mcbi(5)
    it2   = it1
    it3   = it1
    incr  = 1
    incr1 = 1
    nterm = nrowo*iword(it1)
    ii    = 1
    jj    = nrowo
    nrowis= mcbo(3)*iword(it1) + 1
    ii1   = 1
    nri   = mcbi(3)
    nro   = mcbo(3)
    mcbo(2) = 0
    mcbo(6) = 0
    mcbo(7) = 0
    mcbo(3) = nrowo

    DO  i = 1,ncol
        DO  j = 1,nterm
            iz(j) = 0
        END DO
        jj1 = nro
        CALL unpack (*40,scrf,iz)
40      CONTINUE
        jj1 = nri
        CALL unpack (*50,inp,iz(nrowis))
50      CALL pack (iz,outp,mcbo)
    END DO

    CALL CLOSE (scrf,1)
    CALL CLOSE (inp,1)
    CALL CLOSE (outp,1)
    CALL wrttrl (mcbo)

    RETURN
END SUBROUTINE ampc2
