SUBROUTINE intert (nl,nl1,nl2,nm,ajj,ta)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: nl
 INTEGER, INTENT(IN OUT)                  :: nl1
 INTEGER, INTENT(IN OUT)                  :: nl2
 INTEGER, INTENT(IN)                      :: nm
 REAL, INTENT(OUT)                        :: ajj(1)
 REAL, INTENT(IN)                         :: ta(1)
 
 
 t1 = ta(nl1)
 t2 = ta(nl2)
 t  = ta(nl)
 n1 = nm *(nl1-1)
 n2 = nm *(nl2-1)
 n  = nm*(nl -1)
 fract = (t-t1) / (t2 -t1)
 DO  i=1,nm
   ajj(i+n) = ajj(i+n1) + fract *(ajj(i+n2) - ajj(i+n1))
 END DO
 RETURN
END SUBROUTINE intert
