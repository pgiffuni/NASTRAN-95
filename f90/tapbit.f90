LOGICAL FUNCTION tapbit (FILE)
     
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 EXTERNAL  andf
 INTEGER :: fist,xfiat,fiat, andf,nam(2)
 COMMON   /xxfiat/ xfiat(1) /xpfist/ npfist  &
     /xfist / nfist,lfist,fist(1) /xfiat / mfiat,nfiat,lfiat,fiat(1)
 COMMON   /system/ ib(45)
 COMMON   /two   / itwo(32)
 DATA      nam   / 4HTAPB,4HIT   /
 
 tapbit = .true.
 DO  j = 1,npfist
   IF (fist(2*j-1) == FILE) GO TO 20
 END DO
 npf1 = npfist + 1
 DO  j = npf1,lfist
   IF (fist(2*j-1) == FILE) GO TO 30
 END DO
 CALL mesage (-21,FILE,nam)
 
 20 j = -fist(2*j)
 IF (andf(itwo(32-j),ib(45)) /= 0) RETURN
 IF (andf(xfiat(j+1),32768)  == 0) tapbit = .false.
 RETURN
 
 30 j = fist(2*j)
 IF (andf(fiat(j+1),32768) == 0) tapbit = .false.
 RETURN
END FUNCTION tapbit
