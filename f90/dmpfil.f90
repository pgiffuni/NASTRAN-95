SUBROUTINE dmpfil (ifile,z,lz)
     
!     DUMPS A FILE ON DIAG 20 SETTING.
 
 
 INTEGER, INTENT(IN OUT)                  :: ifile
 INTEGER, INTENT(OUT)                     :: z(2)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: sysbuf,outpt,buf, FILE,NAME(2)
 COMMON /system/ sysbuf,outpt
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /unpakx/ iout,irow,nrow,incr
 DATA    NAME  / 4HDMPF,2HIL  /
 
 1 FORMAT (1H0,i5,10(1X,i10,1X))
 2 FORMAT (1H ,5X,10(1X,1P,e11.4))
 3 FORMAT (1H ,5X,10(4X,a4,4X))
 
 CALL sswtch (20,l)
 IF (l == 0) RETURN
 
 FILE = IABS(ifile)
 buf  = lz - sysbuf + 1
 IF (buf < 6) GO TO 91
 lcore = (buf-1)/5
 lcore = lcore*5
 CALL OPEN (*90,FILE,z(buf),rdrew)
 WRITE  (outpt,10) FILE
 10 FORMAT (14H1DUMP of FILE ,i3)
 IF (ifile <= 0) GO TO 100
 
 irec = 0
 20 WRITE  (outpt,30) irec
 30 FORMAT (8H0RECORD ,i6,6X,100(1H-))
 40 CALL READ (*70,*60,FILE,z,lcore,0,iwords)
 
 i1 = -9
 50 i1 = i1 + 10
 i2 = MIN0(i1+9,lcore)
 WRITE (outpt,1) i1,(z(i),i=i1,i2)
 WRITE (outpt,2)    (z(i),i=i1,i2)
 WRITE (outpt,3)    (z(i),i=i1,i2)
 IF (lcore-i2 > 0) THEN
   GO TO    50
 ELSE
   GO TO    40
 END IF
 
 60 i1 = -9
 65 i1 = i1 + 10
 i2 = MIN0(i1+9,iwords)
 WRITE (outpt,1) i1,(z(i),i=i1,i2)
 WRITE (outpt,2)    (z(i),i=i1,i2)
 WRITE (outpt,3)    (z(i),i=i1,i2)
 IF (iwords > i2) GO TO 65
 irec = irec + 1
 GO TO 20
 
 70 z(1) = FILE
 CALL CLOSE (FILE,clsrew)
 CALL rdtrl (z)
 WRITE  (outpt,80) (z(i),i=1,7)
 80 FORMAT (4H0EOF ,//,8H0TRAILER  ,/,7(1X,i12 /))
 90 RETURN
 
 91 CALL mesage (8,0,NAME)
 GO TO 90
 
 100 CALL READ (*70,*101,FILE,z,2,1,iwords)
 101 WRITE  (outpt,102) z(1),z(2)
 102 FORMAT (14H0HEADER record  ,/1H0,2A4)
 z(1) = FILE
 CALL rdtrl (z)
 ncols = z(2)
 IF (ncols > 300) ncols = 100
 iout = 1
 incr = 1
 IF (ncols > 0) THEN
   GO TO   110
 ELSE
   GO TO    70
 END IF
 110 DO  j = 1,ncols
   WRITE  (outpt,115) j
   115 FORMAT (7H0COLUMN  ,i5)
   irow = 0
   nrow = 0
   CALL unpack (*140,FILE,z)
   WRITE  (outpt,118) irow,nrow
   118 FORMAT (1H+,20X,3HROW  ,i4,11H   thru row   ,i5)
   IF (nrow > 300) nrow = 100
   nels = nrow - irow + 1
   IF (nels <= 0) CYCLE
   WRITE  (outpt,119) (z(k),k=1,nels)
   119 FORMAT (1P,10E13.4)
   CYCLE
   140 WRITE  (outpt,141)
   141 FORMAT (13H null column  )
 END DO
 GO TO 70
 
END SUBROUTINE dmpfil
