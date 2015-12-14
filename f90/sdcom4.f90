SUBROUTINE sdcom4( p, ac, wa, wb )
!******
 
! SDCOM4 COMPUTES THE CONTRIBUTIONS OF THE PIVOT ROW FOR SDCOMP IN CDP
 
!******
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: p(2)
 INTEGER, INTENT(IN)                      :: ac(1)
 DOUBLE PRECISION, INTENT(IN)             :: wa(1)
 DOUBLE PRECISION, INTENT(OUT)            :: wb(1)
 INTEGER :: row, c, start, sc
 
 DOUBLE PRECISION :: pir, pii
 DOUBLE PRECISION :: p1,  p2
 
 COMMON/ sdcomx / row, c, spflg, start, frstpc, lastpl, lasti, sc
 
 j = 1
 l = 1
 
! FOR THE OUTER LOOP I RUNS FROM START TO LASTI.
! BEGIN BY FORMING -P(I)/P(1). THEN DECIDE WHICH INNER LOOP TO EXECUTE
 
 p2 = p(1)**2 + p(2)**2
 p1 = p(1) / p2
 p2 = p(2) / p2
 DO  i=start,lasti
   pir = -p(2*i-1) * p1 - p(2*i) * p2
   pii =  p(2*i-1) * p2 - p(2*i) * p1
   IF( i <= lastpl ) GO TO 30
   IF( ac(i) < 0 ) GO TO 20
   k1 = i
   
! LOOP 1 -- L IS INCREMENTED WHENEVER AC(K) .GT. 0
   
   10 DO  k=k1,c
     IF( ac(k) > 0 ) GO TO 12
     wb(j  ) = pir*p(2*k-1) - pii*p(2*k  )
     wb(j+1) = pir*p(2*k  ) + pii*p(2*k-1)
     GO TO 14
     12 wb(j  ) = pir*p(2*k-1) - pii*p(2*k  ) + wa(l  )
     wb(j+1) = pir*p(2*k  ) + pii*p(2*k-1) + wa(l+1)
     l = l + 2
     14 j = j + 2
   END DO
   GO TO 40
   
! LOOP 2 -- L IS NEVER INCREMENTED
   
   20 DO  k=i,c
     wb(j  ) = pir*p(2*k-1) - pii*p(2*k  )
     wb(j+1) = pir*p(2*k  ) + pii*p(2*k-1)
     j = j + 2
   END DO
   GO TO 40
   
! LOOP 3 -- K RUNS FROM I TO LASTPL AND L IS INCREMENTED EVERY TIME
!           THEN, IF LASTPL .LT. C, LOOP 1 IS EXECUTED TO FINISH IT UP
   
   30 DO  k=i,lastpl
     wb(j  ) = pir*p(2*k-1) - pii*p(2*k  ) + wa(l  )
     wb(j+1) = pir*p(2*k  ) + pii*p(2*k-1) + wa(l+1)
     l = l + 2
     j = j + 2
   END DO
   IF( lastpl == c ) GO TO 40
   k1 = lastpl + 1
   GO TO 10
   
! END OUTER LOOP BY STORING -P(I)/P(1) AT P(1).
   
   40 p(2*i-1 ) = pir
   p(2*i   ) = pii
 END DO
 RETURN
END SUBROUTINE sdcom4
