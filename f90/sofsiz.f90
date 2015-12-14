INTEGER FUNCTION sofsiz (dum)
!*****
!     RETURNS THE REMAINING NUMBER OF AVAILABLE WORDS ON THE SOF.
!*****
 
 REAL, INTENT(IN OUT)                     :: dum
 INTEGER :: blksiz,avblks
 DIMENSION    nmsbr(2)
 COMMON /sys/ blksiz,dirsiz,supsiz,avblks
 DATA  nmsbr/ 4HSOFS,4HIZ  /
!*****
 CALL chkopn (nmsbr(1))
 sofsiz = blksiz*avblks
 RETURN
END FUNCTION sofsiz
