SUBROUTINE emg1b (buf,sil,ii,FILE,dampc)
     
!     THIS ROUTINE REPLACES SMA1B AND GROUPS TOGETHER THE
!     SUB-PARTITIONS OF A PIVOT-PARTITION.
 
!     THE SUB-PARTIONS ARE ARRANGED IN CORE BY ASCENDING SILS OF THE
!     ELEMENT INVOLVED.
 
 
 DOUBLE PRECISION, INTENT(IN)             :: buf(1)
 INTEGER, INTENT(IN OUT)                  :: sil
 INTEGER, INTENT(IN OUT)                  :: ii
 INTEGER, INTENT(IN OUT)                  :: FILE
 DOUBLE PRECISION, INTENT(IN)             :: dampc
 LOGICAL :: anycon, error, DOUBLE, last, heat
 INTEGER :: z, zbase, posvec, precis, rowsiz,  &
     eltype, elid, dict, outpt, estid, filtyp, sils, subr(2), flags
 REAL :: rz(1)
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm, uwm, uim, sfm
 COMMON /emgprm/  icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat, icmbar, lcstm, lmat, lhmat
 COMMON /emg1bx/  nsils, posvec(10), ibloc, nbloc, irows, dict(15),  &
     filtyp, sils(10), last
 COMMON /emgdic/  eltype, ldict, nlocs, elid, estid
 COMMON /sma1io/  smaio(36)
 COMMON /zzzzzz/  z(1)
 COMMON /system/  ksystm(65)
 COMMON /iemg1b/  icall, ilast
 
 EQUIVALENCE      (ksystm(2), outpt)
 EQUIVALENCE      (z(1), dz(1), rz(1)),   (c, ic)
 EQUIVALENCE      (smaio(13), if4gg)
 
 DATA    subr  /  4HEMG1,4HB   /
 
 IF (error) RETURN
 IF (sil == -1111111) GO TO 155
 DOUBLE = .false.
 IF (precis == 2) DOUBLE = .true.
 icall = icall + 1
 
!     IF -FILE- EQUALS IF4GG FOR THE OLD ELEMENTS, THEN THE ELEMENT
!     DAMPING CONSTANT SENT IS PLACED IN THE DICTIONARY AND A SIMPLE
!     RETURN IS MADE.
 
 IF (heat) GO TO  10
 IF (FILE /= if4gg) IF (ii) 20,20,10
 c = dampc
 dict(5) = ic
 icall = icall - 1
 RETURN
 
 10 irows = 1
 dict(4) = 1
 GO TO 30
 20 irows = 6
 dict(4) = 63
 30 IF (icall > 1) GO TO 70
 rowsiz = nsils*irows
 dict(3) = rowsiz
 ibloc = jcore + MOD(jcore+1,2)
 IF (dict(2) == 2) GO TO 40
 nbloc = ibloc + rowsiz*irows*precis - 1
 GO TO 50
 40 nbloc = ibloc + irows*precis - 1
 50 IF (nbloc > ncore) CALL mesage (-8,nbloc-ncore,subr)
 IF (DOUBLE) GO TO 60
 DO    i = ibloc, nbloc
   rz(i) = 0.0E0
 END DO
 GO TO 70
 
 60 idbloc = ibloc/2 + 1
 ndbloc = nbloc/2
 DO  i = idbloc,ndbloc
   dz(i) = 0.0D0
 END DO
 
!     INSERT SUB-PARTITION OF PARTITION IN POSITION OF SIL ORDER.
 
!     BUF IS ASSUMED DOUBLE PRECISION.
 
 70 DO  i = 1,nsils
   IF (sil == sils(i)) GO TO 100
 END DO
 WRITE  (outpt,90) sfm,elid
 90 FORMAT (a25,' 3116, ELEMENT ID',i10,' SENDS BAD SIL TO ROUTINE ',  &
     'EMG1B.')
 CALL mesage (-37,0,subr)
 
 100 IF (dict(2) == 2) GO TO 130
 zbase = irows*(i-1)
 kmat = 1
 IF (DOUBLE) GO TO 125
 
!     SINGLE PRECISION ADDITION OF DATA
 
 j1 = ibloc + zbase
 j2 = j1 + irows - 1
 DO    i = 1,irows
   DO    j = j1,j2
     rz(j) = rz(j) + SNGL(buf(kmat))
     kmat = kmat + 1
   END DO
   j1 = j1 + rowsiz
   j2 = j2 + rowsiz
 END DO
 GO TO 150
 
!     DOUBLE PRECISION ADDITION OF MATRIX DATA.
 
 125 j1 = idbloc + zbase
 j2 = j1 + irows - 1
 DO    i = 1,irows
   DO    j = j1,j2
     dz(j) = dz(j) + buf(kmat)
     kmat = kmat + 1
   END DO
   j1 = j1 + rowsiz
   j2 = j2 + rowsiz
 END DO
 GO TO 150
 
!     SIMPLE DIAGONAL MATRIX INSERTION
 
 130 kmat = 1
 IF (DOUBLE) GO TO 145
 j1 = ibloc
 DO    i = 1,irows
   rz(j1) = rz(j1) + SNGL(buf(kmat))
   j1 = j1 + 1
   kmat = kmat + 14
 END DO
 GO TO 150
 
 145 j1 = idbloc
 DO    i = 1,irows
   dz(j1) = dz(j1) + buf(kmat)
   j1 = j1 + 1
   kmat = kmat + 14
 END DO
 
 150 RETURN
 
!     OUTPUT PIVOT-ROWS-PARTITION
 
 155 IF (icall <= 0) GO TO 161
 IF (.NOT. last) GO TO 160
 ilast = 1
 160 CALL emgout (z(ibloc),z(ibloc),(nbloc-ibloc+1)/precis,ilast,  &
     dict,filtyp,precis)
 161 ilast = 0
 icall = 0
 RETURN
END SUBROUTINE emg1b
