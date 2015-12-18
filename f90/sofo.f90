SUBROUTINE sofo
     
!     MODULE USED TO TRANSFER NASTRAN DATA BLOCKS TO THE SOF FILE FOR
!     PURPOSES OF SAVING THE DATA FOR SUBSEQUENT RUNS OR SUBSEQUENT
!     EXECUTION STEPS.  THE CALLING SEQUENCE TO THE MODULE IS
 
!     SOFO     A,B,C,D,E//V,N,DRY/C,N,NAME/C,N,IA/C,N,IB/C,N,IC/
!                         C,N,ID/C,N,IE $
 
 INTEGER :: sysbuf,dry,FILE,xxxx
 DIMENSION       FILE(5),modnam(2),mcb(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / dry,NAME(2),items(2,5)
 COMMON /zzzzzz/ iz(1)
 COMMON /system/ sysbuf,nout
 DATA    FILE  / 101,102,103,104,105/
 DATA    modnam/ 4HSOFO,4H    /
 DATA    iblnk , xxxx / 4H    ,4HXXXX/
 
 DO  i = 1,5
   IF (items(1,i) == xxxx .OR. items(1,i) == 0) items(1,i) = iblnk
 END DO
 
 IF (dry < 0) RETURN
 nz = korsz(iz)
 IF (3*sysbuf > nz) CALL mesage (-8,0,modnam(1))
 ib1 = nz  - sysbuf + 1
 ib2 = ib1 - sysbuf - 1
 ib3 = ib2 - sysbuf
 CALL sofopn (iz(ib1),iz(ib2),iz(ib3))
 
!     COPY MATRICES FROM NASTRAN DATA BLOCKS TO SOF.
 
 DO  i = 1,5
   IF (items(1,i) == iblnk) CYCLE
   mcb(1) = FILE(i)
   CALL rdtrl (mcb)
   IF (mcb(1) < 0) CYCLE
   CALL mtrxo (FILE(i),NAME(1),items(1,i),0,itest)
   SELECT CASE ( itest )
     CASE (    1)
       GO TO 10
     CASE (    2)
       GO TO 50
     CASE (    3)
       GO TO 50
     CASE (    4)
       GO TO 20
     CASE (    5)
       GO TO 30
     CASE (    6)
       GO TO 50
   END SELECT
   10 WRITE (nout,1010) uwm,items(1,i),NAME(1),NAME(2)
   dry = -2
   CYCLE
   20 WRITE (nout,1020) uwm,NAME(1),NAME(2)
   dry = -2
   EXIT
   30 WRITE (nout,1030) uwm,items(1,i)
   dry = -2
   50 CONTINUE
 END DO
 60 CALL sofcls
 RETURN
 
!     ERROR MESSAGES.
 
 1010 FORMAT (a25,' 6211, MODULE SOFO - ITEM ',a4,' OF SUBSTRUCTURE ',  &
     2A4,' HAS ALREADY BEEN WRITTEN.')
 1020 FORMAT (a25,' 6212, MODULE SOFO - THE SUBSTRUCTURE ',2A4,  &
     ' DOES NOT EXIST.')
 1030 FORMAT (a25,' 6213, MODULE SOFO - ',a4,' IS AN ILLEGAL ITEM NAME')
END SUBROUTINE sofo
