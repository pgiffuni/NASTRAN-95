SUBROUTINE mtimsu (y,x,buf)
     
    !     MTIMSU  FORMS THE  PRODUCT  X = M*Y
 
    DOUBLE PRECISION, INTENT(IN)             :: y(1)
    DOUBLE PRECISION, INTENT(OUT)            :: x(1)
    REAL, INTENT(IN OUT)                     :: buf(1)
    INTEGER :: diag      ,eol       ,filem    ,filek
    DOUBLE PRECISION :: da
 
    COMMON   /invpwx/  filek(7)  ,filem(7)
    COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
        rew       ,norew    ,eofnrw   ,rsp      ,  &
        rdp       ,csp      ,cdp      ,sqr      ,  &
        rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
    COMMON   /invpxx/  dumm(13)  ,nzero
    COMMON   /zntpkx/  a(4)      ,ii       ,eol
    !     COMMON   /DESCRP/  LENGTH    ,MAJOR(1)
    EQUIVALENCE        (a(1),da)
 
 
    ncol = filek(2)
    DO  i = 1,ncol
        x(i) = 0.d0
    END DO
 
    !     MASS MATRIX IS NOT DIAGONAL
 
    nzero = 0
    DO  i = 1,ncol
        IF (y(i) == 0.d0) GO TO 30
        CALL intpk (*40,filem(1),0,rdp,0)
        nzero = nzero + 1
20      CALL zntpki
        x(ii) = da*y(i) + x(ii)
        IF (eol == 0.0) THEN
            GO TO    20
        ELSE
            GO TO    40
        END IF
30      CALL skprec (filem,1)
        nzero = nzero + 1
40  CONTINUE
    END DO
    GO TO 90
90  CALL REWIND (filem(1))
    CALL skprec (filem,1)
    nzero = ncol - nzero

    RETURN
END SUBROUTINE mtimsu
