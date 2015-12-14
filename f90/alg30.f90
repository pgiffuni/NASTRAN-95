SUBROUTINE alg30 (a,b)
 
    REAL, INTENT(IN OUT)                     :: a(9,9)
    REAL, INTENT(IN OUT)                     :: b(9)
    DIMENSION  INDEX(9)
 
    n=9
    DO  j=1,n
        INDEX(j)=0
    END DO
110 amax=-1.0
    DO  j=1,n
        IF(INDEX(j) /= 0)CYCLE
        DO  l=1,n
            IF(INDEX(l) /= 0)CYCLE
            pv=ABS(a(j,l))
            IF(pv <= amax)CYCLE
            ir=j
            ic=l
            amax=pv
        END DO
    END DO
    IF(amax <= 0.0)RETURN
    INDEX(ic)=ir
    IF(ic == ir)GO TO 150
    DO  l=1,n
        pv=a(ir,l)
        a(ir,l)=a(ic,l)
        a(ic,l)=pv
        IF(l > 1)CYCLE
        pv=b(ir)
        b(ir)=b(ic)
        b(ic)=pv
    END DO
150 pv=1.0/a(ic,ic)
    a(ic,ic)=1.0
    DO  l=1,n
        a(ic,l)=a(ic,l)*pv
        IF(l > 1)CYCLE
        b(ic)=b(ic)*pv
    END DO
    DO  l1=1,n
        IF(l1 == ic)CYCLE
        pv=a(l1,ic)
        a(l1,ic)=0.0
        DO  l=1,n
            a(l1,l)=a(l1,l)-a(ic,l)*pv
        END DO
        b(l1)=b(l1)-b(ic)*pv
    END DO
    GO TO 110

END SUBROUTINE alg30
