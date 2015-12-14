SUBROUTINE alg29 (y,x,fxy,n)
 
    REAL, INTENT(IN)                         :: y(3)
    REAL, INTENT(IN)                         :: x(3)
    REAL, INTENT(OUT)                        :: fxy(3)
    INTEGER, INTENT(IN)                      :: n
 
 
    x1=(x(3)+x(2))*(y(3)-y(2))/(x(3)-x(2))
    fxy(2)=x1/(x(3)-x(1))
    n2=n-2
    DO  j=3,n2
        x2=(x(j+1)+x(j))*(y(j+1)-y(j))/(x(j+1)-x(j))
        fxy(j)=(x2-x1)/(x(j+1)-x(j-1))
        x1=x2
    END DO
    fxy(n-1)=-x1/(x(n)-x(n-2))
    fxy(1)=fxy(2)-(fxy(3)-fxy(2))/(x(3)-x(2))*(x(2)-x(1))
    fxy(n)=fxy(n-1)+(fxy(n-1)-fxy(n-2))/(x(n-1)-x(n-2))*(x(n)-x(n-1))

    RETURN
END SUBROUTINE alg29
