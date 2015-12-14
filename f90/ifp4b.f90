SUBROUTINE ifp4b (FILE,scrt,any,SPACE,lspace,recid,eof)
     
!     THIS ROUTINE, CALLED BY IFP4, COPIES DATA FROM -FILE- TO -SCRT-
!     UP TO THE -RECID- SPECIFIED IF IT EXISTS AND COPIES THE -RECID-
!     IN ANY EVENT.  -ANY- IS SET TRUE IF THE -RECID- WAS FOUND.
!     -EOF- IS SET TRUE AS SOON AS AN END OF FILE IS HIT ON -FILE-.
!     -SPACE- IS A WORKING AREA OF LENGTH -LSPACE-.  IF RECID(1) = -1,
!     THE BALANCE OF -FILE- IS COPIED TO -SCRT- AND THEN -FILE- IS
!     REWOUND AND -SCRT- IS COPIED BACK ONTO -FILE-.  BOTH FILES ARE
!     THEN CLOSED.
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: scrt
 LOGICAL, INTENT(OUT)                     :: any
 INTEGER, INTENT(IN OUT)                  :: SPACE(5)
 INTEGER, INTENT(IN)                      :: lspace
 INTEGER, INTENT(IN)                      :: recid(2)
 LOGICAL, INTENT(OUT)                     :: eof
 LOGICAL :: bit,nogo
 INTEGER :: output,NAME(2),sysbuf, flag,buf1,buf2,REC(3),ilimit(3),eor
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,output,nogo
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 DATA    NAME  , noeor,eor  /4HIFP4,4HB   ,0,1/, ilimit/ 3*2147483647 /
 
 IF (recid(1) == -1) GO TO 3000
 any = .false.
 
!     CHECK TRAILER BIT TO SEE IF RECORD EXISTS
 
 CALL ifp4f (recid(2),FILE,bit)
 IF (.NOT.bit) GO TO 2000
 IF (eof) GO TO 5001
 
!     READ A 3-WORD RECORD ID
 
 any = .true.
 ifile = FILE
 650 CALL READ  (*5001,*6002,FILE,REC(1),3,noeor,flag)
 CALL WRITE (scrt,REC,3,noeor)
 IF (REC(1) == recid(1)) RETURN
 
!     NOT THE CORRECT RECORD, THUS COPY BALANCE OF RECORD OVER.
 
 750 CALL READ  (*6001,*800,FILE,SPACE,lspace,noeor,flag)
 CALL WRITE (scrt,SPACE,lspace,noeor)
 GO TO 750
 800 CALL WRITE (scrt,SPACE,flag,eor)
 GO TO 650
 
!     RECORD DOES NOT CURRENTLY EXIST, THUS START ONE
 
 2000 CALL WRITE (scrt,recid,2,noeor)
 CALL WRITE (scrt,0,1,noeor)
 
!     PUT BIT IN TRAILER
 
 CALL ifp4g (recid(2),FILE)
 RETURN
 
!     WRAP UP FILES
 
 3000 IF (eof) GO TO 3400
 3100 CALL READ  (*3500,*3200,FILE,SPACE,lspace,noeor,flag)
 CALL WRITE (scrt,SPACE,lspace,noeor)
 GO TO 3100
 3200 CALL WRITE (scrt,SPACE,flag,eor)
 GO TO 3100
 3400 CALL WRITE (scrt,ilimit,3,eor)
 
!     FILE IS ALL COPIED TO SCRT.  REWIND AND RETURN.
 
 3500 eof = .true.
 CALL CLOSE (scrt,clsrew)
 CALL CLOSE (FILE,clsrew)
 
!     COPY DATA FROM SCRT TO FILE.
 
 buf1 = 1
 buf2 = sysbuf + 2
 i    = 2*sysbuf + 4
 j    = lspace - i
 IF (i > lspace) CALL mesage (-8,0,NAME)
 ifile = FILE
 CALL OPEN (*6003,FILE,SPACE(buf1),wrtrew)
 ifile = scrt
 CALL OPEN  (*6003,scrt,SPACE(buf2),rdrew)
 3800 CALL READ  (*4000,*3900,scrt,SPACE(i),j,noeor,flag)
 CALL WRITE (FILE,SPACE(i),j,noeor)
 GO TO 3800
 3900 CALL WRITE (FILE,SPACE(i),flag,eor)
 GO TO 3800
 4000 CALL CLOSE (scrt,clsrew)
 CALL CLOSE (FILE,clsrew)
 RETURN
 
!     ERROR CONDITIONS
 
 5001 nogo = .true.
 WRITE  (output,5002) ufm,recid(1),recid(2),FILE
 5002 FORMAT (a23,' 4056, RECORD ID =',2I10,' IS OUT OF SYNC ON DATA ',  &
     'BLOCK NUMBER',i10, /5X,'AN IFP4 SYSTEM ERROR.')
 eof = .true.
 RETURN
 
 6001 CALL mesage (-2,ifile,NAME)
 6002 CALL mesage (-3,ifile,NAME)
 6003 CALL mesage (-1,ifile,NAME)
 RETURN
END SUBROUTINE ifp4b
