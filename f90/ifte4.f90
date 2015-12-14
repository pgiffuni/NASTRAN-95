SUBROUTINE ifte4(tha,rp,cp)
     
 
 REAL, INTENT(IN)                         :: tha
 REAL, INTENT(OUT)                        :: rp
 REAL, INTENT(OUT)                        :: cp
 DATA thao,epsi /.1,1.e-9 /
 IF(ABS(tha) < thao) GO TO 100
 d = tha**4/24.
 rp = ((.5*(tha*tha))-1.+ COS(tha))/d
 cp = ((tha**3/6.)-tha+SIN(tha))/d
 RETURN
 
!  EVALUATE SERIES
 
 100 CONTINUE
 rn = 1.0
 d = 1.0
 SIGN = -1.
 rps = 1.
 tsq = tha*tha
 t1 = 5.
 t2 = 6.
 it = 1
 101 CONTINUE
 DO  i=1,50
   rn = rn*tsq
   d = d*t1*t2
   trm = rn/d*SIGN
   rps = rps + trm
   IF(ABS(trm) < epsi) EXIT
   SIGN = -SIGN
   t1 = t1+2.
   t2 = t2+2.
 END DO
 120 CONTINUE
 IF(it == 2) GO TO 125
 rp = rps
 rn = tha
 d = 5.0
 SIGN = -1.
 rps = tha/5.
 t1 = 6.
 t2 = 7.
 it = 2
 GO TO 101
 125 CONTINUE
 cp = rps
 RETURN
END SUBROUTINE ifte4
