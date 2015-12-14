SUBROUTINE scalex(ilval,code,l)
     
 
 INTEGER, INTENT(IN)                      :: ilval
 INTEGER, INTENT(IN)                      :: code
 INTEGER, INTENT(OUT)                     :: l(1)
 INTEGER :: expnd(6)
 DO  i=1,6
   l(i)=0
 END DO
 IF(code > 0.0) THEN
   GO TO   103
 END IF
 102 l(1) = ilval
 GO TO 110
 103 id=code
 DO  i=1,6
   inv=7-i
   expnd(inv)=MOD(id,10)
   id=id/10
 END DO
 j=0
 loop107:  DO  i=1,6
   IF(expnd(i) == 0) CYCLE loop107
   IF(i < 2) GO TO 106
   ii=i-1
   DO  k=1,ii
     IF(expnd(k) == expnd(i)) CYCLE loop107
   END DO
   106 j=j+1
   l(j)=expnd(i)
 END DO loop107
 i=0
 108 i=i+1
 l(i)=ilval+l(i)-1
 IF(i-j < 0) THEN
   GO TO   108
 END IF
 110 RETURN
END SUBROUTINE scalex
