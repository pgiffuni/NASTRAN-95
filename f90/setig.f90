SUBROUTINE setig (kg1,kg2,ig,norig)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS ROUTINE SETS IG(KG1,-)=KG2 AND IG(KG2,-)=KG1 IF THIS
!     CONNECTION HAS NOT ALREADY BEEN SET.
!     NEDGE = NUMBER OF UNIQUE EDGES.
 
 
 INTEGER, INTENT(IN)                      :: kg1
 INTEGER, INTENT(IN)                      :: kg2
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN OUT)                  :: norig(1)
 INTEGER :: bunpk
 DIMENSION  sub(2)
 COMMON /bands /  nn,       mm,       dum2(2),  maxgrd,   maxdeg,  &
     dum3(3),  nedge
 COMMON /system/  ibuf,     nout
 DATA             sub /     4HSETI,   4HG       /
 
 IF (kg1 == 0 .OR. kg2 == 0 .OR. kg1 == kg2) GO TO 80
 l=kg1
 k=kg2
 DO  loop=1,2
   IF (loop == 1) GO TO 20
   l=kg2
   k=kg1
   20 m=0
   30 m=m+1
   IF (m > maxdeg) GO TO 60
   is=bunpk(ig,l,m)
   IF (is == 0) GO TO 40
   IF (is /= k) GO TO 30
   GO TO 80
   40 CALL bpack (ig,l,m,k)
   mm=MAX0(mm,m)
   IF (loop == 1) nedge = nedge + 1
 END DO
 GO TO 80
 
 60 WRITE (nout,70) norig(l),maxdeg
 70 FORMAT (34H0***  fatal error - - - grid point,i10,  &
     48H  has degree exceeding the nodal degree limit of,i8)
 CALL mesage (-8,0,sub)
 80 RETURN
END SUBROUTINE setig
