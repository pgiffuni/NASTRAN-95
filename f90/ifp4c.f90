SUBROUTINE ifp4c (FILE,scrt,buf1,buf2,eof)
     
!     THIS ROUTINE, CALLED BY IFP4, OPENS THE 2 FILES AND COPIES THE
!     HEADER RECORD FROM -FILE- TO -SCRT-.
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: scrt
 INTEGER, INTENT(IN OUT)                  :: buf1(10)
 INTEGER, INTENT(IN OUT)                  :: buf2(10)
 LOGICAL, INTENT(OUT)                     :: eof
 
 INTEGER :: work(10),flag,NAME(2), name2(2),eor,trail(7)
 COMMON /names/ rd,rdrew,wrt,wrtrew,clsrew,cls
 DATA    NAME / 4HIFP4,4HC   /, eor,noeor/1,0/
 
 trail(1) = FILE
 DO  i = 2,7
   trail(i) = 0
 END DO
 CALL rdtrl (trail)
 DO  i = 2,7
   IF (trail(i) == 0.0) THEN
     GO TO    70
   ELSE
     GO TO    60
   END IF
   70 CONTINUE
 END DO
 GO TO 1000
 60 CALL OPEN (*1002,FILE,buf1,rdrew)
 eof = .false.
 CALL OPEN (*2000,scrt,buf2,wrtrew)
 80 CALL READ (*1001,*100,FILE,work,10,noeor,flag)
 CALL WRITE (scrt,work,10,noeor)
 GO TO 80
 100 CALL WRITE (scrt,work,flag,eor)
 RETURN
 
!     FILE IS NULL
 
 1000 eof = .true.
 CALL OPEN (*2000,scrt,buf2,wrtrew)
 CALL fname (FILE,name2)
 CALL WRITE (scrt,name2,2,eor)
 RETURN
 
 2000 CALL mesage (-1,scrt,NAME)
 1001 CALL mesage (-2,FILE,NAME)
 1002 CALL mesage (-1,FILE,NAME)
 RETURN
END SUBROUTINE ifp4c
