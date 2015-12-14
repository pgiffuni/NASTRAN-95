SUBROUTINE read3 (novect,ncol,sr1fil,sr2fil,filc,kdblm)
     
!     READ3 PACKS THE EIGENVECTORS AND EIGENVALUES AND PUTS THEM OUT IN
!     ASCENDING ORDER
 
!     LAST REVISED  1/92, BY G.CHAN/UNISYS
!     ZERO OUT RIGID BODY FREQUENCIES IF METHOD IS 'FEER' (NOT 'FEER-X'
!     NOR 'FEER-Q')
 
 
 INTEGER, INTENT(IN)                      :: novect
 INTEGER, INTENT(IN)                      :: ncol
 INTEGER, INTENT(IN)                      :: sr1fil
 INTEGER, INTENT(IN)                      :: sr2fil
 INTEGER, INTENT(IN)                      :: filc
 INTEGER, INTENT(IN)                      :: kdblm
 INTEGER :: sysbuf    ,iz(1)    ,rsp      ,rdp      , filelm    ,filevc   ,  &
      option   ,optn2    ,feer     , dashz     ,sturm
 INTEGER :: rdrew
 DOUBLE PRECISION :: dxx(2)
 DIMENSION          nam(2)    ,filevc(7),filelm(7)
 COMMON   /zzzzzz/  z(1)
 COMMON   /sturmx/  sturm     ,shftpt
 COMMON   /reigkr/  option    ,optn2
 COMMON   /system/  sysbuf    ,nout     ,systm(52),iprec
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp
 COMMON   /packx /  itypa     ,itypb    ,ipak     ,jpak     , incr
 COMMON   /unpakx/  itypu     ,iunp     ,junp     ,incru
 EQUIVALENCE        (iz(1),z(1))
 DATA      feer  ,  dashz /4HFEER,   4H-x   /
 DATA      nam   /  4HREAD,4H3   /,  i2 / 2 /
 
!     FILELM (=KDBLM=LAMA=201) WILL HOLD THE EIGENVALUES  UPON RETURN
!     FILEVC (=FILC =PHIA=202) WILL HOLD THE EIGENVECTORS UPON RETURN
 
 filelm(1) = kdblm
 filevc(1) = filc
 itypa = rsp
 itypb = rsp
 incr  = 1
 ipak  = 1
 jpak  = ncol
 ncol2 = iprec*ncol
 itypu = rsp
 incru = 1
 nocl  = 2*ncol + 2
 nz    = korsz(z)
 ibuf1 = nz    - sysbuf
 ibuf2 = ibuf1 - sysbuf
 
!     READ IN ALL EIGENVALUES
 
 ifile = sr1fil
 CALL gopen (sr1fil,z(ibuf1),rdrew)
 i = 1
 10 CALL fread (sr1fil,dxx,iprec,1)
 z(i+1) = dxx(1)
 i = i + 1
 IF (i <= novect) GO TO 10
 CALL CLOSE (sr1fil,rew)
 
!     SET UP AN INDEX VECTOR AND SORT THE EIGENVALUES
 
 j = ncol + 2
 k = j + ncol - 1
 ii = 1
 DO  i = j,k
   iz(i) = ii
   ii = ii + 1
 END DO
 z(1) = z(i2)
 j = 2
 k = j + novect - 1
 DO  i = j,k
   IF (z(i) < z(1)) z(1) = z(i)
 END DO
 DO  i = 1,novect
   k = i
   30 IF (z(k+1) >= z(k)) CYCLE
   zz     = z(k  )
   z(k  ) = z(k+1)
   z(k+1) = zz
   j  = k + ncol
   ii = iz(j)
   iz(j  ) = iz(j+1)
   iz(j+1) = ii
   k = k - 1
   GO TO 30
 END DO
 
!     ZERO OUT RIGID BODY EIGENVALUES IF THEY ARE PRESENT AND METHOD IS
!     'FEER-Z'
!     I.E. ZERO FREQUENCIES BELOW PTSHFT AND KEEP, AS CHECKED BY STURM
!     SEQUENCE
 
 IF (sturm < 0) GO TO 45
 DO  i = 2,novect
   ik = i + sturm
   IF (z(ik) >= shftpt .OR. ik > novect) GO TO 45
   IF (z(i) < 0. .AND. option == feer .AND. optn2 == dashz) z(i)= 0.
 END DO
 
!     READ THE EIGENVECTORS AND PACK THEM IN ASCENDING ORDER
 
 45 CALL gopen (filevc,z(ibuf1),1)
 ifile = sr2fil
 CALL gopen (sr2fil,z(ibuf2),rdrew)
 ipos = 1
 CALL makmcb (filevc(1),filc,ncol,2,rsp)
 
 DO  i = 1,novect
   k  = i + ncol + 1
   no = iz(k)
   IF (no-ipos < 0) THEN
     GO TO    50
   ELSE IF (no-ipos == 0) THEN
     GO TO    80
   ELSE
     GO TO    70
   END IF
   50 CALL REWIND (sr2fil)
   ipos = no
   IF (no <= 0) GO TO 120
   60 CALL skprec (sr2fil,no)
   GO TO 80
   70 no   = no - ipos
   ipos = ipos + no
   GO TO 60
   80 iunp = 0
   CALL unpack (*90,sr2fil,z(nocl))
   ipos = ipos + 1
   ipak = iunp
   jpak = junp
   GO TO 100
   90 ipak = 1
   jpak = 1
   z(nocl) = 0.0
   100 CALL pack (z(nocl),filevc,filevc)
 END DO
 
 CALL CLOSE  (filevc(1),rew)
 CALL CLOSE  (sr2fil,rew)
 CALL wrttrl (filevc)
 
!     OUTPUT THE EIGENVALUES, 1ST DATA RECORD
 
 CALL gopen (filelm,z(ibuf1),1)
 CALL WRITE (filelm,z(i2),novect,1)
 
!     SAVE ORDER FOUND IN 2ND DATA RECORD
 
 CALL WRITE (filelm,iz(ncol+2),novect,1)
 CALL CLOSE (filelm(1),rew)
 filelm(2) = novect
 CALL wrttrl (filelm)
 RETURN
 
 120 CALL mesage (-7,FILE,nam)
 RETURN
 
END SUBROUTINE read3
