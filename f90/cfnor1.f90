SUBROUTINE cfnor1 (right,left,size2,option,ri)
     
    !     CFNOR1 IS A SINGLE-PRECISION ROUTINE (CREATED FOR USE BY
    !     THE COMPLEX FEER METHOD) WHICH NORMALIZES A COMPLEX PAIR
    !     OF VECTORS TO MAGNITUDE UNITY
 
    !     DEFINITION OF INPUT PARAMETERS
 
    !     RIGHT    = ORIGINAL RIGHT-HANDED COMPLEX SINGLE PRECISION VECTOR
    !     LEFT     = ORIGINAL LEFT -HANDED COMPLEX SINGLE PRECISION VECTOR
    !     SIZE2    = LENGTH OF EITHER VECTOR IN SINGLE PRECISION WORDS
    !                (I.E., TWICE THE LENGTH OF THE COMPLEX VECTORS)
    !     OPTION   = 0  NORMALIZE THE INPUT VECTORS, AND OUTPUT THE
    !                   SQUARE ROOT OF THE INNER PRODUCT IN RI(2)
    !              = 1  ONLY OUTPUT INNER-PRODUCT, IN RI(2)
    !              = 2  ONLY OUTPUT SQUARE ROOT OF INNER-PRODUCT, IN RI(2)
 
    !     DEFINITION OF OUTPUT PARAMETERS
 
    !     RIGHT    = NORMALIZED RIGHT-HANDED VECTOR
    !     LEFT     = NORMALIZED LEFT -HANDED VECTOR
    !     RI       = INNER-PRODUCT, OR SQUARE ROOT OF INNER-PRODUCT (SEE
    !                OPTION)
 
 
    REAL, INTENT(IN OUT)                     :: right(1)
    REAL, INTENT(IN OUT)                     :: left(1)
    INTEGER, INTENT(IN)                      :: size2
    INTEGER, INTENT(IN OUT)                  :: option
    REAL, INTENT(OUT)                        :: ri(2)
    LOGICAL :: skip     ,qpr
 
 
    DIMENSION  rj(2)
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm
    COMMON /feerxc/ dumxc(21),qpr
    COMMON /system/ ksystm(65)
    EQUIVALENCE     (ksystm(2),nout)
 
    skip = .false.
 
    !     COMPUTE INNER PRODUCT (LEFT*RIGHT)
 
5   ri(1) = 0.
    ri(2) = 0.
    DO  i = 1,size2,2
        j = i + 1
        ri(1) = ri(1) + left(i)*right(i) - left(j)*right(j)
        ri(2) = ri(2) + left(j)*right(i) + left(i)*right(j)
    END DO
    IF (option == 1) GO TO 50
    IF (skip) GO TO 200
 
    !     COMPUTE MAGNITUDE OF SQUARE ROOT
 
    rsqrt = SQRT(SQRT(ri(1)**2+ri(2)**2))
    IF (rsqrt > 0.) GO TO 30
    WRITE  (nout,20) uwm
20  FORMAT (a25,' 3162', //5X,'ATTEMPT TO NORMALIZE NULL VECTOR. ',  &
        'NO ACTION TAKEN.'//)
    GO TO 50
 
    !     COMPUTE MODULUS OF SQUARE ROOT
 
30  theta2 = .5*ATAN2(ri(2),ri(1))
 
    !     COMPUTE REAL AND IMAGINARY PARTS OF SQUARE ROOT OF INNER PRODUCT
 
    ri(1) = rsqrt*COS(theta2)
    ri(2) = rsqrt*SIN(theta2)
    IF (option == 2) GO TO 50
    rj(1) = ri(1)
    rj(2) = ri(2)
 
    !     INVERT THE ABOVE COMPLEX NUMBER (THETA2 IS DUMMY)
 
    theta2= 1./(ri(1)**2+ri(2)**2)
    ri(1) = ri(1)*theta2
    ri(2) =-ri(2)*theta2
 
    !     NORMALIZE THE INPUT VECTORS
 
    DO  i = 1,size2,2
        j = i + 1
        theta2   = right(i)
        right(i) = ri(1)*right(i) - ri(2)*right(j)
        right(j) = ri(2)*theta2   + ri(1)*right(j)
        theta2   = left(i)
        left(i)  = ri(1)*left(i) - ri(2)*left(j)
        left(j)  = ri(2)*theta2  + ri(1)*left(j)
    END DO
 
    !     ----------- SPECIAL PRINT ----------------------------------------
    IF (.NOT.qpr) GO TO 45
    skip = .true.
    GO TO 5
200 theta2 = SQRT(ri(1)**2 + ri(2)**2)
    WRITE  (nout,300) theta2,ri
300 FORMAT (3H --,32(4H----), /,7H cfnor1,6X,  &
        16HOUTPUT magnitude,e16.8,8X,2E16.8, /,3H --,32(4H----))
    WRITE  (nout,400) (right(i),i=1,size2)
400 FORMAT ((1H ,4E25.16))
    WRITE  (nout,500)
500 FORMAT (3H --,32(4H----))
    WRITE  (nout,400) (left(i),i=1,size2)
    WRITE  (nout,500)
    !     ------------------------------------------------------------------
 
45  ri(1) = rj(1)
    ri(2) = rj(2)

50  RETURN
END SUBROUTINE cfnor1
