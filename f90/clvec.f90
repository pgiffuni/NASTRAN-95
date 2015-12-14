SUBROUTINE clvec (lamd,nvect,phidl,ih,ibuf,ibuf1)
    !*****
    !     CLVEC CACLULATES THE LEFT EIGENVECTORS FOR THE DETERMINANT AND
    !     UPPER HESSENBERG APPROACHES TO THE COMPLEX EIGENVALUE PROBLEM
    !*****
 
    INTEGER, INTENT(IN OUT)                  :: lamd
    INTEGER, INTENT(IN)                      :: nvect
    INTEGER, INTENT(IN OUT)                  :: phidl
    INTEGER, INTENT(IN OUT)                  :: ih(7)
    INTEGER, INTENT(IN OUT)                  :: ibuf
    INTEGER, INTENT(IN)                      :: ibuf1
    DOUBLE PRECISION :: di1,dnrow,dz(1),lambda,mindia
    INTEGER :: clsrew,flag, rdrew,switch,sysbuf
    INTEGER :: filek,filem,fileb,scr
    DIMENSION NAME(2),buf(6)
    COMMON / cdcmpx / dumdcp(30),mindia
    COMMON / zzzzzz / z(1)
    COMMON / cinvpx / filek(7),filem(7),fileb(7),dum(15),scr(11)
    COMMON / cinvxx / lambda(2),switch
    COMMON / names  / rd,rdrew,wrt,wrtrew,clsrew,norew
    COMMON / packx  / it1,it2,ii,jj,inc
    COMMON / system / sysbuf
    EQUIVALENCE (nrow,filek(3))
    EQUIVALENCE (dz(1),z(1))
    DATA NAME   / 4HCLVE,4HC    /
    !*****
    !     INITIALIZATION
    !*****
    ibuf2 = ibuf1 - sysbuf
    IF (fileb(1) < 0) fileb(1) = 0
    IF (fileb(6) == 0) fileb(1) = 0
    DO  i=1,11
        scr(i) = 300 + i
    END DO
    switch = -204
    fnrow = FLOAT(nrow)
    dnrow = fnrow
    !*****
    !     OPEN SORTED EIGENVALUE FILE
    !*****
    CALL gopen (lamd,z(ibuf),rdrew)
    CALL skprec (lamd,1)
    !*****
    !     LOOP TO CALCULATE LEFT EIGENVECTORS
    !*****
    DO  i=1,nvect
        ! READ EIGENVALUE
        CALL READ(*9002,*9003,lamd,buf,6,0,flag)
        lambda(1) = buf(3)
        lambda(2) = buf(4)
        ! CREATE DYNAMIC MATRIX
100     CALL cinvp1
        ! DECOMPOSE DYNAMIC MATRIX
        CALL cinvp2(*900)
        ! BUILD LOAD FOR FBS
        fi1 = FLOAT(i-1)
        di1 = fi1
        j2 = 2*nrow
        DO  j=1,j2,2
            f = FLOAT((j+1)/2)
            dz(j) = mindia/(1.0D0 + (1.0D0 - f/dnrow)*di1)
            dz(j+1) = 0.0D0
        END DO
        ! PERFORM FORWARD-BACKWARD SUBSTITUTION - U(T)*L(T)*PHI
        CALL cdifbs (dz(1),z(ibuf2))
        ! NORMALIZE LEFT EIGENVECTOR
        CALL cnorm1 (dz(1),nrow)
        ! PACK LEFT EIGENVECTOR ONTO PHIDL
        it1 = 4
        it2 = 3
        ii = 1
        jj = nrow
        inc = 1
        CALL pack (dz(1),phidl,ih)
        CYCLE
        ! ADJUST CURRENT EIGENVALUE
900     lambda(1) = 1.01D0*lambda(1)
        lambda(2) = 1.01D0*lambda(2)
        GO TO 100
    ! END OF LOOP
    END DO
    CALL CLOSE (lamd,clsrew)
    RETURN
    !*****
    !     ERRORS
    !*****
9002 n = -2
    GO TO 9999
9003 n = -3
9999 CALL mesage (n,lamd,NAME)

    RETURN
END SUBROUTINE clvec
