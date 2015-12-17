SUBROUTINE nasopn ( *, lu, dsn )
     
 INTEGER, INTENT(OUT)                     :: lu
 CHARACTER (LEN=80), INTENT(OUT)          :: dsn
 CHARACTER (LEN=80) :: ifile
 INCLUDE 'NASNAMES.COM'
 LOGICAL :: iexist
 
 klen = INDEX( rfdir, ' ' )
 ifile = rfdir(1:klen-1) // '/NASINFO'
 dsn = ifile
 INQUIRE ( FILE=ifile, EXIST=iexist )
 IF ( .NOT. iexist ) GO TO 100
 OPEN ( UNIT=lu, FILE=ifile, STATUS='OLD', ERR=100 )
 RETURN
 100   RETURN 1
END SUBROUTINE nasopn
