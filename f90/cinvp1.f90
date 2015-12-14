SUBROUTINE cinvp1

    !*******
    !     CINVP1 INITIALIZES AND CALLS SUBROUTINE ADD FOR CINVPR
    !*******

    DOUBLE PRECISION :: alpha(2) ,beta(2)  ,lambda
    INTEGER :: scr1     ,scr2     ,scr11     ,sqr
    INTEGER :: filek    ,filem    ,fileb     ,cdp
    INTEGER :: ifila(7) ,ifilb(7) ,ifilc(7)  ,csp
    INTEGER :: sysbuf   ,switch
 
    COMMON   /cinvpx/  filek(7)  ,filem(7) ,fileb(7) ,dum(15)  ,  &
        scr1      ,scr2     ,scr(8)   ,scr11
    COMMON   /cinvxx/  lambda(2) ,switch
    COMMON   /saddx /  nomat     ,nz       ,mcbs(67)
    COMMON   /zzzzzz/  z(1)
    COMMON   /names /  dumm(9)   ,csp      ,cdp      ,sqr
    COMMON   /system/  sysbuf
 
    EQUIVALENCE        ( mcbs( 1), ifila(1) ) ,( mcbs( 8), itypal   )  &
        ,( mcbs(61), ifilc(1) ) ,( mcbs(13), ifilb(1) )  &
        ,( mcbs(20), itypbt   ) ,( mcbs(21), beta(1)  )  &
        ,( mcbs( 9), alpha(1) )

    !*******
    !     FORM -(B+LAMBDA*M) ON SCR2
    !*******

    nomat = 2
    DO  i = 1,7
        ifila(i) = filem(i)
        ifilb(i) = fileb(i)
    END DO
    alpha(1) = -lambda(1)
    alpha(2) = -lambda(2)
    beta(1)  = -1.d0
    beta(2)  = 0.d0
    itypal   = cdp
    itypbt   = cdp
    nz = korsz(z)
    IF (switch == -204) nz = nz - 2*sysbuf
    ifilc(1) = scr2
    IF (switch /= 0) ifilc(1) = scr11
    ifilc(2) = filek(2)
    ifilc(3) = filek(3)
    ifilc(4) = 1
    ifilc(5) = cdp
    CALL sadd (z,z)
    !*******
    !     FORM (LAMBDA**2*M+LAMBDA*B+K) ON SCR1
    !*******
    DO  i = 1,7
        ifila(i) = filek(i)
    END DO
    ifilb(1) = ifilc(1)
    ifilb(2) = filek(2)
    ifilb(3) = filek(3)
    ifilb(4) = sqr
    alpha(2) = 0.d0
    beta(1)  = -lambda(1)
    beta(2)  = -lambda(2)
    ifilb(5) = cdp
    alpha(1) = 1.d0
    ifilc(1) = scr1
    CALL sadd (z,z)

    RETURN
END SUBROUTINE cinvp1
