SUBROUTINE a8 2 INT (*,a,n,b,INT)
     
 
 , INTENT(IN OUT)                         :: *
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
 
 ENTRY k8 2 INT (*,c,n,b,INT)
!     ****************************
 
 nt = +1
 GO TO 30
 
 ENTRY a8 2 fp (*,a,n,b,INT)
!     ***************************
 
 nt = -1
 
 20   IF (n > 8) GO TO 50
 INT = nt
 CALL na12if (*80,a,n,b,INT)
 RETURN
 
 ENTRY k8 2 fp (*,c,n,b,INT)
!     ***************************
 
 nt = -1
 
 30   IF (n > 8) GO TO 50
 INT = nt
 CALL nk12if (*80,c,n,b,INT)
 RETURN
 
 50   WRITE  (nout,60) n,nt
 60   FORMAT ('  N.GT.8/A82INT',i5,7X,'NT=',i2)
 80   RETURN 1
END SUBROUTINE a8 2 INT
