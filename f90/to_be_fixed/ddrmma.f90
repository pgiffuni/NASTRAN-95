SUBROUTINE ddrmma( setup )
!*****
!  UNPACKS DATA FROM A TRANSIENT OR FREQUENCY RESPONSE SOLUTION
!  COLUMN AS REQUIRED TO FORM ONE OFP OUTPUT LINE ENTRY.
 
!  BEFORE CALLING FOR ENTRY CONSTRUCTION ONE SETUP CALL IS REQUIRED
!  FOR EACH COLUMN. (SETUP = .TRUE.)
!*****
 
 LOGICAL, INTENT(IN)                      :: setup
 REAL :: lambda   ,rbufa(75),rbufb(75)
 
 INTEGER :: buf(150), bufa(75), bufb(75), elwork(300), phase, complx
 INTEGER :: scrt,buff,FILE,outfil,setid,dhsize,entrys,filnam,passes
 INTEGER :: sets,device,FORM,uvsol
 INTEGER :: typout
 INTEGER :: savdat,savpos,bufsav
 
 LOGICAL :: trnsnt   ,sort2    ,col1     ,frstid
 LOGICAL :: lminor
 
 COMMON/stdata/    lminor    ,nstxtr   ,npos     ,savdat(75)  &
     ,savpos(25)         ,bufsav(10)
 COMMON/ddrmc1/    idrec(146),buff(6)  ,passes   ,outfil   ,jfile  &
     ,mcb(7)   ,entrys   ,sets(5,3),infile   ,lambda  &
     ,FILE     ,sort2    ,col1     ,frstid   ,ncore  &
     ,nsols    ,dhsize   ,filnam(2),rbuf(150),idout  &
     ,icc      ,ncc      ,ilist    ,nlist    ,nwds  &
     ,setid    ,trnsnt   ,i1       ,i2       ,phase  &
     ,itype1   ,itype2   ,nptsf    ,lsf      ,nwdsf  &
     ,scrt(7)  ,ierror   ,itemp    ,device   ,FORM  &
     ,istlst   ,lstlst   ,uvsol    ,nlambs   ,nwords ,omega    ,ipass
 COMMON/clstrs/    complx(1)
 COMMON/zntpkx/ a(4), irow, ieol, IEOR
 
 EQUIVALENCE(buf(1),rbuf(1),bufa(1),rbufa(1)) ,(rbufb(1),bufb(1),buf(76))
!*****
!  PERFORM SOLUTION COLUMN SETUP WHEN SETUP = .TRUE.
!*****
 IF( .NOT. setup ) GO TO 10
 typout = 3
 IF( trnsnt ) typout = 1
 icomp = 1
 CALL intpk(*5,scrt(6),0,typout,0)
 CALL zntpki
 RETURN
 5 irow = 0
 RETURN
!*****
!  FILL BUFFER WITH REAL AND OR COMPLEX VALUES.
!*****
 10 k = i1 - 1
 DO  i = 1,k
   bufb(i) = bufa(i)
 END DO
 DO  i = i1,i2
   IF( icomp == irow ) GO TO 30
   rbufa(i) = 0.0
   rbufb(i) = 0.0
   GO TO 60
   
!     NON-ZERO COMPONENT AVAILABLE.
   
   30 IF( .NOT. trnsnt ) SELECT CASE ( ipass )
     CASE (    1)
       GO TO 31
     CASE (    2)
       GO TO 32
     CASE (    3)
       GO TO 33
   END SELECT
   
!     TRANSIENT RESPONSE
   
   rbufa(i) = a(1)
   IF( ieol  > 0) THEN
     GO TO    50
   ELSE
     GO TO    40
   END IF
   
!     FREQUENCY RESPONSE FOR DISPLACEMENTS OR SPCFS PASS
   
   31 rbufa(i) = a(1)
   rbufb(i) = a(2)
   IF( ieol  > 0) THEN
     GO TO    50
   ELSE
     GO TO    40
   END IF
   
!     FREQUENCY RESPONSE VELOCITYS PASS
   
   32 rbufa(i) = -omega * a(2)
   rbufb(i) =  omega * a(1)
   IF( ieol  > 0) THEN
     GO TO    50
   ELSE
     GO TO    40
   END IF
   
!     FREQUENCY RESPONSE ACCELERATIONS PASS
   
   33 rbufa(i) = omega * a(1)
   rbufb(i) = omega * a(2)
   IF( ieol  > 0) THEN
     GO TO    50
   END IF
   40 CALL zntpki
   GO TO 60
   
   50 irow = 0
   
   60 icomp = icomp + 1
   
 END DO
!*****
!  IF TRANSIENT (REAL) THEN RETURN. FOR FREQUENCY (COMPLEX) COMBINE DATA
!  FOR OUTPUT AND CONVERT TO MAGNITUDE PHASE IF NECESSARY.
 
!  BUFA CONTAINS THE REAL PART
!  BUFB CONTAINS THE IMAGINARY PART
!*****
 IF (trnsnt) GO TO 81
 IF (itype1 == 4)  GO TO 90
 IF (itype1 == 5)  GO TO 81
 
!     POINT DATA
 
 DO  k = 1,6
   IF( FORM == 3 ) CALL magpha( bufa(k+2), bufb(k+2) )
   bufa(k+8) = bufb(k+2)
 END DO
 nwdsf = 14
 RETURN
 
!     ELEMENT STRESS OR FORCE DATA
 
 81 IF (lminor)  GO TO 90
 DO  k=1,nstxtr
   j=savpos(npos+k-1)
   buf(j) = bufsav(k)
 END DO
 90 IF (trnsnt) RETURN
 iout = 0
 i = nptsf
 100 npt = complx(i)
 IF( npt  < 0) THEN
   GO TO   110
 ELSE IF ( npt  == 0) THEN
   GO TO   140
 ELSE
   GO TO   130
 END IF
 110 npt = -npt
 IF( FORM /= 3 ) GO TO 130
 
!     COMPUTE MAGNITUDE PHASE
 
 CALL magpha( bufa(npt), bufb(npt) )
 120 iout = iout + 1
 elwork(iout) = bufa(npt)
 i = i + 1
 GO TO 100
 130 IF( npt <= lsf ) GO TO 120
 npt = npt - lsf
 iout = iout + 1
 elwork(iout) = bufb(npt)
 i = i + 1
 GO TO 100
 
!     MOVE OUTPUT DATA
 
 140 DO  i = 1,iout
   buf(i) = elwork(i)
 END DO
 nwdsf = iout
 RETURN
END SUBROUTINE ddrmma
