SUBROUTINE a82int (*,a,n,b,INT)
 
    REAL, INTENT(IN OUT)                     :: a(2)
    INTEGER, INTENT(IN OUT)                  :: n
    REAL, INTENT(IN OUT)                     :: b
    INTEGER, INTENT(OUT)                     :: INT
    CHARACTER (LEN=8) :: c
 
    COMMON /xreadx/ nout
 
    !     THESE ROUTINES PERFORM IN THE OPPOSITE DIRECTION AS THOSE OF THE
    !     INT2A8 GROUP OF ROUTINES
    !     THIS ROUTINE IS MACHINE INDEPENDENT
 
    !     ENTRY POINTS   A8 2 INT  (BCD-INTEGER VERSION)
    !                    K8 2 INT  (CHARACTER-INTEGER VERSION)
    !                    A8 2 FP   (BCD-REAL VERSION)
    !                    K8 2 FP   (CHARACTER-REAL VERSION)
 
    nt = +1
    GO TO 20
 
    ENTRY k82int (*,c,n,b,int)
    !     ****************************
 
    nt = +1
    GO TO 30
 
    ENTRY a82fp (*,a,n,b,int)
    !     ***************************
 
    nt = -1
 
20  IF (n > 8) GO TO 50
    INT = nt
    CALL na12if (*80,a,n,b,int)
    RETURN
 
    ENTRY k82fp (*,c,n,b,int)
    !     ***************************
 
    nt = -1
 
30  IF (n > 8) GO TO 50
    INT = nt
    CALL nk12if (*80,c,n,b,int)
    RETURN
 
50  WRITE  (nout,60) n,nt
60  FORMAT ('  N.GT.8/A82INT',i5,7X,'NT=',i2)

80  RETURN 1
END SUBROUTINE a82int
