SUBROUTINE rfopen (member,lu)
     
!     THIS .MIS ROUTINE OPENS THE RIGID FORMAT FILE, AS AN ORDINARY
!     FORTRAN FILE. USE REGULAR FORTRAN READ TO READ THE FILE
 
!     ENTRY POINT RFCLSE TO CLOSE IT
 
!     IF RIGID FORMAT FILE OPENS OK, LU IS THE FORTRAN UNIT NUMBER
!     OTHERWISE, LU = 0
 
!     THIS ROUTINE REPLACES ALL THE MACHINE DEPENDENT DSXOPN, DSXCLS,
!     DSXREA, AND DSXFRE ROUTINES. PLUS DSXRDS, DSXIO, AND DSXSIO IN
!     IBM VERSION, AND DSXRET AND DSXZER IN CDC
 
!     NOTE - FORTRAN UNIT 'IN' IS USED TO READ THE RIGID FORMAT FILE.
!            UNIT 'IN' IS SYNCHRONOUS WITH ANY READFILE OR NESTED
!            READFILE OPERATION.
 
!     WRITTEN BY G.CHAN/UNISYS.   10/1990
 
 
 INTEGER, INTENT(IN OUT)                  :: member(2)
 INTEGER, INTENT(OUT)                     :: lu
 INTEGER :: facsf
 CHARACTER (LEN=1) :: bk,mb1(8)
 CHARACTER (LEN=6) :: mb6
 CHARACTER (LEN=5) :: mb5
 CHARACTER (LEN=8) :: mb8,free8,add(3)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBI
 CHARACTER (LEN=44) :: rfdir, dsn
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /xxread/ in
 COMMON /system/ ibuf,nout,nogo
 EQUIVALENCE     (mb1(1),mb5,mb6,mb8)
 DATA    bk,     add(1),   add(3), free8     /  &
     ' '   , '@ADD,E ',' .  ', '@FREE   '/
 
 CALL a42k8 (member(1),member(2),mb8)
 IF (mach == 3) GO TO 30
 in = in + 1
 IF (in < 60) in = 60
 j  = 5
 IF (mb1(6) /= bk) j = 6
 
!           DUMMY  IBM  UNVC  CDC  VAX  ULTRIX  SUN   AIX   HP
!             S/G  MAC  CRAY CNVX  NEC  FUJTSU   DG  AMDL PRIME
!             486 DUMMY ALFA RESV
!            ---- ----  ---- ---- ----  ------ ----  ---- -----
 GO TO (  60,  50,   30,  50,  50,     50,  50,   50,  50,  &
     50,  70,   70,  70,  70,     70,  70,   70,  70, 50,  60,   50,  70), mach
 
!     UNIVAC ONLY -
!     ADD FILE TO INPUT STREAM
 
 30   add(2) = mb8
 j = facsf(add)
 lu = 5
 GO TO 130
 50    CONTINUE
 rfdir = ' '
 CALL getenv ( 'RFDIR', rfdir )
 DO  i = 44, 1, -1
   IF ( rfdir( i:i ) == ' ' ) CYCLE
   lenr = i
   GO TO 56
 END DO
 lenr = 44
 56    dsn = ' '
 dsn = rfdir(1:lenr) // '/' // mb6
!WKBR IF (J .EQ. 6) OPEN (UNIT=IN,FILE=MB6,ACCESS='SEQUENTIAL',ERR=100,
 OPEN (UNIT=in,FILE=dsn,ACCESS='SEQUENTIAL',ERR=100,  &
     FORM='FORMATTED',STATUS='OLD')
 GO TO 80
 
!     OTHERS -
 
 60   GO TO 100
 
 70   OPEN (UNIT=in,FILE=mb8,ACCESS='SEQUENTIAL',ERR=100,STATUS='OLD',  &
     FORM='FORMATTED')
 
!     VERIFY FILE EXISTANCE
 
 80   READ (in,90,ERR=100,END=100) j
 90   FORMAT (a1)
 REWIND in
 lu = in
 GO TO 130
 
!WKBR100  WRITE  (NOUT,110) SFM,MB8
 100  WRITE  (nout,110) sfm,dsn
!WKBR 110  FORMAT (A25,', RFOPEN CAN NOT OPEN ',A8)
 110  FORMAT (a25,', RFOPEN CAN NOT OPEN ',a44)
 
 IF (mach > 7 .AND. mach /= 21) WRITE (nout,120) mach
 120  FORMAT (5X,'MACHINE',i4,' IS NOT AVAILABLE/RFOPEN')
 lu   = 0
 nogo = 1
 
 130  RETURN
 
 
 ENTRY rfclse (lu)
!     =================
 
 IF (mach == 3) GO TO 150
 IF (lu  < 60) WRITE (nout,140) sfm,lu
 140  FORMAT (a25,'. RFCLSE/RFOPEN ERROR.  LU =',i4)
 CLOSE (UNIT=lu)
 in = in - 1
 IF (in < 60) in = 0
 GO TO 160
 
 150  add(1) = free8
 j = facsf(add)
 160  lu = 0
 RETURN
END SUBROUTINE rfopen
