SUBROUTINE sofi
     
!     MODULE USED TO COPY SELECTED ITEMS FROM SELECTED SUBSTRUCTURES
!     ONTO NASTRAN MATRIX FILES.   THE CALLING SEQUENCE TO THE MODULE
!     IS
!     SOFI   /A,B,C,D,E/V,N,DRY/C,N,NAME/C,N,IA/C,N,IB/C,N,IC/C,N,ID/
!                       C,N,IE $
 
 INTEGER :: dry,FILE,sysbuf,xxxx
 DIMENSION       FILE(5),modnam(2),mcb(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / dry,NAME(2), items(2,5)
 COMMON /zzzzzz/ iz(1)
 COMMON /system/ sysbuf,nout
 DATA    FILE  / 201,202,203,204,205 /
 DATA    iblnk , xxxx/4H    ,4HXXXX  /
 DATA    modnam/ 4HSOFI,4H           /
 
 DO  i = 1,5
   IF (items(1,i) == xxxx .OR. items(1,i) == 0) items(1,i) = iblnk
 END DO
 
 nz = korsz(iz)
 IF (3*sysbuf > nz) CALL mesage (-8,0,modnam(1))
 ib1 = nz  - sysbuf + 1
 ib2 = ib1 - sysbuf - 1
 ib3 = ib2 - sysbuf
 CALL sofopn (iz(ib1),iz(ib2),iz(ib3))
 IF (dry >= 0) GO TO 60
 
!     CHECK THE EXISTENCE OF THE SOF FILE.
 
 DO  i = 1,5
   IF (items(1,i) == iblnk) CYCLE
   mcb(1) = FILE(i)
   CALL rdtrl (mcb)
   IF (mcb(1) < 0) CYCLE
   CALL softrl (NAME(1),items(1,i),mcb)
   itest = mcb(1)
   SELECT CASE ( itest )
     CASE (    1)
       GO TO 50
     CASE (    2)
       GO TO 50
     CASE (    3)
       GO TO 20
     CASE (    4)
       GO TO 30
     CASE (    5)
       GO TO 40
   END SELECT
   20 WRITE (nout,1020) uwm,items(1,i),NAME(1),NAME(2)
   GO TO 45
   30 WRITE (nout,1030) uwm,NAME(1),NAME(2)
   dry = -2
   GO TO 130
   40 WRITE (nout,1040) uwm,items(1,i)
   45 dry = -2
 END DO
 GO TO 130
 
!     COPY SOF DATA INTO NASTRAN DATA BLOCKS
 
 60 DO  i = 1,5
   IF (items(1,i) == iblnk) CYCLE
   mcb(1) = FILE(i)
   CALL rdtrl (mcb)
   IF (mcb(1) < 0) CYCLE
   CALL mtrxi (FILE(i),NAME(1),items(1,i),0,itest)
   SELECT CASE ( itest )
     CASE (    1)
       GO TO 120
     CASE (    2)
       GO TO 70
     CASE (    3)
       GO TO 80
     CASE (    4)
       GO TO 90
     CASE (    5)
       GO TO 100
     CASE (    6)
       GO TO 120
   END SELECT
   70 WRITE (nout,1050) uwm,items(1,i),NAME(1),NAME(2)
   GO TO 110
   80 WRITE (nout,1020) uwm,items(1,i),NAME(1),NAME(2)
   CYCLE
   90 WRITE (nout,1030) uwm,NAME(1),NAME(2)
   dry = -2
   EXIT
   100 WRITE (nout,1040) uwm,items(1,i)
   110 dry = -2
 END DO
 130 CALL sofcls
 RETURN
 
!     ERROR MESSAGES.
 
 1020 FORMAT (a25,' 6216, MODULE SOFI - ITEM ',a4,' OF SUBSTRUCTURE ',  &
     2A4,' DOES NOT EXIST.')
 1030 FORMAT (a25,' 6212, MODULE SOFI - THE SUBSTRUCTURE ',2A4,  &
     ' DOES NOT EXIST.')
 1040 FORMAT (a25,' 6213, MODULE SOFI - ',a4,' IS AN ILLEGAL ITEM NAME')
 1050 FORMAT (a25,' 6215, MODULE SOFI - ITEM ',a4,' OF SUBSTRUCTURE ',  &
     2A4,' PSEUDO-EXISTS ONLY.')
END SUBROUTINE sofi
