SUBROUTINE permut(ia,id,n,isw)
     
 
 INTEGER, INTENT(IN)                      :: ia(1)
 INTEGER, INTENT(OUT)                     :: id(10)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: isw
 DIMENSION  ib(32),ic(32)
 DO  i=1,n
   ic(i)=ia(i)
   ib(i)=i
 END DO
 n1=n-1
 DO  i=1,n1
   i1=i+1
   DO  j=i1,n
     IF(ic(j)-ic(i) < 0) THEN
       GO TO    40
     ELSE
       GO TO    30
     END IF
     40 is1=ib(j)
     ib(j)=ib(i)
     ib(i)=is1
     is1 = ic(j)
     ic(j)=ic(i)
     ic(i) = is1
     30 CONTINUE
   END DO
   50 CONTINUE
 END DO
 DO  i = 1,n
   IF(ic(i)-isw < 0) THEN
     GO TO    50
   ELSE
     GO TO    60
   END IF
 END DO
 k=1
 GO TO 71
 60 DO  j=i,n
   k=j-i+1
   id(k)=ib(j)
 END DO
 IF(k == n) GO TO 90
 k=k+1
 71 DO  j=k,n
   l=j-k+1
   id(j)=ib(l)
 END DO
 90 RETURN
END SUBROUTINE permut
