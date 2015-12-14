SUBROUTINE redu (cdata,nx,ix,nas,ias,nvar,var,ipre,ier)
     
 
 INTEGER, INTENT(IN)                      :: cdata(5)
 INTEGER, INTENT(IN)                      :: nx
 INTEGER, INTENT(IN)                      :: ix(3,1)
 INTEGER, INTENT(IN OUT)                  :: nas
 INTEGER, INTENT(IN OUT)                  :: ias(2,1)
 INTEGER, INTENT(OUT)                     :: nvar
 INTEGER, INTENT(OUT)                     :: var(3,6)
 INTEGER, INTENT(IN OUT)                  :: ipre
 INTEGER, INTENT(OUT)                     :: ier
 INTEGER :: BLANK,eqs
 DIMENSION  keys(6)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,iout
 DATA    keys  / 4HNAMA,4HNAMB,4HNONA,4HNONB,4HPREC,4HBOUN /
 
 DATA    NAME  / 4HNAME/, eqs /4H=   /,BLANK/4H            /
 
!     INITIALLIZE
 
 DO  i = 1,6
   var(1,i) = keys(i)
 END DO
 
 nvar = 18
 DO  i = 1,2
   var(2,i) = BLANK
   var(3,i) = BLANK
 END DO
 DO  i = 3,6
   var(2,i) = -1
   var(3,i) = 0
 END DO
 
!     DECODE COMMAND
 
 i2 = 4
 IF (cdata(5) ==  eqs) i2 = 6
 IF (cdata(1)*2 < i2) GO TO 100
 
 var(2,1) = cdata(i2  )
 var(3,1) = cdata(i2+1)
 
 nvx = 6
 
!     FIND NAME
 
 DO  i = 1,nx
   IF (ix(1,i) /= NAME) GO TO 35
   var(2,2) = ix(2,i)
   var(3,2) = ix(3,i)
   CYCLE
   35  IF (ix(1,i) /= keys(6)) GO TO 37
   var(2,6) = ix(2,i)
   var(3,6) = ix(3,i)
   CYCLE
   37 nvx = nvx + 1
   DO  j = 1,3
     var(j,nvx) = ix(j,i)
   END DO
 END DO
 IF (var(2,2) == BLANK) GO TO 100
 IF (var(3,6) <=     0) GO TO 120
 IF (ipre <= 0 .OR. ipre > 2) ipre = 1
 
 var(3,5) = ipre
 
!     FIND STRUCTURE NUMBERS, B MAY NOT PRE-EXIST
 
 IF (nas == 0)  GO TO  80
 DO  i = 1,nas
   IF (var(2,1) /= ias(1,i) .OR. var(3,1) /= ias(2,i)) GO TO 55
   var(3,3) = i
   CYCLE
   55 IF (var(2,2) == ias(1,i) .AND. var(3,2) == ias(2,i)) GO TO 100
 END DO
 80 nas = nas + 1
 var(3,4) = nas
 ias(1,nas) = var(2,2)
 ias(2,nas) = var(3,2)
 IF (var(3,3) /= 0) GO TO 90
 nas = nas + 1
 var(3,3) = nas
 ias(1,nas) = var(2,1)
 ias(2,nas) = var(3,1)
 90 ier  = 0
 nvar = nvx*3
 RETURN
 
 100 WRITE  (iout,101) ufm
 101 FORMAT (a23,' 6614, ILLEGAL OR NON-EXISTANT STRUCTURE NAME USED ',  &
     'ABOVE')
 GO TO 130
 120 WRITE  (iout,121) ufm
 121 FORMAT (a23,' 6615, ILLEGAL BOUNDARY SET IDENTIFICATION NUMBER')
 130 ier = 1
 RETURN
END SUBROUTINE redu
