SUBROUTINE sdcom1 (p,ac,wa,wb)
     
 
 REAL, INTENT(IN OUT)                     :: p(1)
 INTEGER, INTENT(IN OUT)                  :: ac(1)
 REAL, INTENT(IN)                         :: wa(1)
 REAL, INTENT(OUT)                        :: wb(1)
 INTEGER :: row,c,start
 
 COMMON /sdcomx/ row,c,spflg,start,frstpc,lastpl,lasti
 
 j  = 1
 l  = 1
 k1 = lastpl + 1
 iend   = MIN0(lastpl,lasti)
 istart = MAX0(k1,start)
 IF (c == lastpl) GO TO 200
 IF (start > lastpl) GO TO 100
 DO  i = start,iend
   pi   =-p(i)/p(1)
   ijmk = j - i
   ilmk = l - i
   DO  k = i,lastpl
     wb(k+ijmk) = pi*p(k) + wa(k+ilmk)
   END DO
   l = ilmk + k1
   DO  k = k1,c
     IF (ac(k) > 0) GO TO 12
     wb(k+ijmk) = pi*p(k)
     CYCLE
     12 wb(k+ijmk) = pi*p(k) + wa(l)
     l = l + 1
   END DO
   j = ijmk + c + 1
   p(i) = pi
 END DO
 IF (lastpl >= lasti) RETURN
 100 DO  i = istart,lasti
   pi   = -p(i)/p(1)
   ijmk = j - i
   IF (ac(i) < 0) GO TO 120
   DO  k = i,c
     IF (ac(k) > 0) GO TO 112
     wb(k+ijmk) = pi*p(k)
     CYCLE
     112 wb(k+ijmk) = pi*p(k) + wa(l)
     l = l + 1
   END DO
   GO TO 140
   120 DO  k = i,c
     wb(k+ijmk) = pi*p(k)
   END DO
   140 j = ijmk + c + 1
   p(i) = pi
 END DO
 RETURN
 
 200 IF (start > lastpl) GO TO 300
 DO  i = start,iend
   pi   = -p(i)/p(1)
   ijmk = j - i
   ilmk = l - i
   DO  k = i,lastpl
     wb(k+ijmk) = pi*p(k) + wa(k+ilmk)
   END DO
   j = ijmk + k1
   l = ilmk + k1
   p(i) = pi
 END DO
 IF (lastpl >= lasti) RETURN
 300 DO  i = istart,lasti
   pi   = -p(i)/p(1)
   ijmk = j - i
   IF (ac(i) < 0) GO TO 320
   DO  k = i,c
     IF (ac(k) > 0) GO TO 312
     wb(k+ijmk) = pi*p(k)
     CYCLE
     312 wb(k+ijmk) = pi*p(k) + wa(l)
     l = l + 1
   END DO
   GO TO 340
   320 DO  k = i,c
     wb(k+ijmk) = pi*p(k)
   END DO
   340 j = ijmk + c + 1
   p(i) = pi
 END DO
 RETURN
END SUBROUTINE sdcom1
