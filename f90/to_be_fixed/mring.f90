SUBROUTINE mring (points)
     
!     HEAT CONDUCTIVITY SMA2 ROUITNE FOR TRIANGULAR (POINTS=3) AND
!     TRAPEZOIDAL (POINTS=4) RING ELEMENTS.
!     THIS ROUTINE IS SEPARATE FROM MTRAPR AND MTRIRG SO AS TO BE
!     IN OVERLAY WITH MTRMEM AND MQDMEM.
 
 
 INTEGER, INTENT(IN)                      :: points
 LOGICAL :: nogo
 INTEGER :: outpt    ,sysbuf   ,tint     ,map(15)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON /system/ sysbuf   ,outpt    ,nogo
 COMMON /sma2et/ ecpt(100)
 EQUIVALENCE     (t,tint)
 DATA    pi23  / 2.0943951024  /
 DATA    map   / 1,2,3, 1,2,3, 2,3,4, 3,4,1, 4,1,2 /
 
!     ECPT LISTS
 
!     ECPT     TRIRG -------- TRMEM          TRAPRG ------- QDMEM
!     ===========================================================
!      1       EL-ID          EL-ID          EL-ID          EL-ID
!      2       SIL-1          SIL-1          SIL-1          SIL-1
!      3       SIL-2          SIL-2          SIL-2          SIL-2
!      4       SIL-3          SIL-3          SIL-3          SIL-3
!      5       THETA          THETA          SIL-4          SIL-4
!      6       MATID          MATID          THETA          THETA
!      7       CSID-1         T              MATID          MATID
!      8       X1             NS-MASS        CSID-1         T
!      9       Y1             CSID-1         X1             NS-MASS
!     10       Z1             X1             Y1             CSID-1
!     11       CSID-2         Y1             Z1             X1
!     12       X2             Z1             CSID-2         Y1
!     13       Y2             CSID-2         X2             Z1
!     14       Z2             X2             Y2             CSID-2
!     15       CSID-3         Y2             Z2             X2
!     16       X3             Z2             CSID-3         Y2
!     17       Y3             CSID-3         X3             Z2
!     18       Z3             X3             Y3             CSID-3
!     19       AVG-TEMP       Y3             Z3             X3
!     20                      Z3             CSID-4         Y3
!     21                      AVG-TEMP       X4             Z3
!     22                                     Y4             CSID-4
!     23                                     Z4             X4
!     24                                     AVG-TEMP       Y4
!     25                                                    Z4
!     26                                                    AVG-TEMP
 
!     GEOMETRY CHECKS X  MUST BE .GT.0, AND Y = 0  FOR I = 1,2,..,POINTS
!                      I                     I
 
 i1 = points + 4
 i2 = i1 + 4*points - 1
 DO  i = i1,i2,4
   IF (ecpt(i+1) < 0.0) THEN
     GO TO   140
   END IF
   10 IF (ecpt(i+2) == 0.0) THEN
     GO TO    20
   ELSE
     GO TO   140
   END IF
 END DO
 
!     POINT ORDERING CHECK.
 
 IF (points == 4) GO TO 30
 i1 = 1
 i2 = 3
 GO TO 40
 30 i1 = 4
 i2 = 15
 40 jpoint = points + 1
 DO  i = i1,i2,3
   ir = map(i  )*4 + jpoint
   is = map(i+1)*4 + jpoint
   it = map(i+2)*4 + jpoint
   temp = (ecpt(is) - ecpt(ir))*(ecpt(it+2) - ecpt(is+2))  -  &
       (ecpt(it) - ecpt(is))*(ecpt(is+2) - ecpt(ir+2))
   IF (temp > 0.0) THEN
     GO TO    50
   ELSE
     GO TO   140
   END IF
 END DO
 
!     TRAPEZOID TEST.
 
 IF (points /= 4) GO TO 100
 IF (ecpt(11) - ecpt(15) == 0.0) THEN
   GO TO    60
 ELSE
   GO TO    70
 END IF
 60 IF (ecpt(19) - ecpt(23) == 0.0) THEN
   GO TO    90
 END IF
 70 CALL page2 (-2)
 WRITE  (outpt,80) swm,ecpt(1)
 80 FORMAT (a27,' 3091, A TRAPRG ELEMENT =',i14,' DOES NOT HAVE ',  &
     'SIDE 1-2 PARALLEL TO SIDE 3-4.')
 
!     THICKNESS OF TRMEM OR QDMEM TO BE CALLED BELOW.
!     QDMEM WILL SUBDIVIDE THICKNESS FOR SUB-TRIANGLES AND THUS
!     T IS SET = INTEGER 1 AS A FLAG TO QDMEM ROUTINE WHICH WILL
!     COMPUTE T FOR EACH.
 
!     TEMP. PATH FOR APPROX. THICKNESS
 
 90 t = pi23*(ecpt(9) + ecpt(13) + ecpt(17) + ecpt(21))*3.0/4.0
 GO TO 110
 100 t = pi23*(ecpt(8) + ecpt(12) + ecpt(16))
 
!  CONVERT ECPT TO THAT OF A TRMEM OR QDMEM.
 
 110 j = 5*points + 6
 k = 4*points + 1
 DO  i = 1,k
   ecpt(j) = ecpt(j-2)
   j = j - 1
 END DO
 ecpt(points+4) = t
 ecpt(points+5) = 0.0
 IF (points == 4) GO TO 130
 
!     MTRMEM CALL
 
 CALL masstq (4)
 RETURN
 
!     MQDMEM CALL
 
 130 CALL masstq (1)
 RETURN
 
!     BAD GEOMETRY FATAL ERROR.
 
 140 WRITE  (outpt,150) ufm,ecpt(1)
 150 FORMAT (a23,' 3092, TRIRG OR TRAPRG ELEMENT =',i14,' POSSESSES ',  &
     'ILLEGAL GEOMETRY.')
 nogo = .true.
 RETURN
END SUBROUTINE mring
