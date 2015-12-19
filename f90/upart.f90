SUBROUTINE upart (uset,scr1,major,sub0,sub1)
     
!     UPART ALONG WITH MPART WILL PERFORM A SYMMETRIC PARTITION OF A
!     MATRIX
 
 
 INTEGER, INTENT(IN)                      :: uset
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN OUT)                  :: major
 INTEGER, INTENT(IN OUT)                  :: sub0
 INTEGER, INTENT(IN OUT)                  :: sub1
 INTEGER :: rule,pvect,uset1
 
 COMMON /parmeg/ ia(7),ia11(7),ia12(7),ia21(7),ia22(7),lcore,rule
 COMMON /patx  / lc,n1,n2,n3,uset1,pvect(7)
 COMMON /zzzzzz/ core(1)
 
 
 uset1 = uset
 
!     TRANSFER OF PVECT TRAILER AS LOADED BY CALCV IS NOW BY /PATX/
 
 rule  = 0
 lc    = korsz(core)
 lcore = lc
 CALL calcv (scr1,major,sub0,sub1,core)
 n4 = n2 + n3
 ia11(2) = n1
 ia11(3) = n1
 ia21(2) = n4
 ia21(3) = n1
 ia21(4) = 2
 ia12(2) = n1
 ia12(3) = n4
 ia12(4) = 2
 ia22(2) = n4
 ia22(3) = n4
 10 RETURN
 
 
 ENTRY mpart (ia1,ia111,ia121,ia211,ia221)
!     =========================================
 
 ia(1) = ia1
 CALL rdtrl (ia)
 IF (ia(1) < 0) THEN
   GO TO    10
 END IF
 20 ia11(1) = ia111
 ia12(1) = ia121
 ia21(1) = ia211
 ia22(1) = ia221
 ia11(4) = ia(4)
 ia11(5) = ia(5)
 ia21(5) = ia(5)
 ia12(5) = ia(5)
 ia22(4) = ia(4)
 ia22(5) = ia(5)
 CALL partn (pvect,pvect,core)
 DO  i = 1,4
   j = (i-1)*7 + 1
   IF (ia11(j) == 0) THEN
     GO TO    40
   END IF
   30 CALL wrttrl (ia11(j))
   40 CONTINUE
 END DO
 GO TO 10
END SUBROUTINE upart
