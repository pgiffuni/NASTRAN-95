SUBROUTINE calcv (pvact,set1,sub1,sub2,core)
 
 INTEGER, INTENT(IN OUT)                  :: pvact
 INTEGER, INTENT(IN OUT)                  :: set1
 INTEGER, INTENT(IN OUT)                  :: sub1
 INTEGER, INTENT(IN OUT)                  :: sub2
 INTEGER, INTENT(IN OUT)                  :: core(1)
 EXTERNAL        andf
 INTEGER :: a1,pvect, uset,sysbuf, andf
 COMMON /system/ sysbuf
 COMMON /two   / two1(32)
 COMMON /patx  / lc,n,no,n3,uset,pvect(7)
 COMMON /zblpkx/ b1(4),n1
 
 n  = 0
 n3 = 0
 no = 0
 n1 = 0
 CALL makmcb (pvect,pvact,0,2,1)
 lcore = lc - sysbuf
 CALL gopen (uset,core(lcore+1),0)
 lcore = lcore - sysbuf
 CALL gopen (pvact,core(lcore+1),1)
 CALL bldpk (1,1,pvact,0,0)
 20 CALL READ (*90,*90,uset,a1,1,0,flag)
 IF (andf(two1(set1),a1) == 0) GO TO 20
 n1 = n1 + 1
 IF (andf(two1(sub1),a1) == 0) GO TO 50
 n  = n  + 1
 GO TO 20
 50 IF (andf(two1(sub2),a1) == 0) GO TO 60
 no = no + 1
 b1(1) = 1.0
 GO TO 70
 60 b1(1) = 2.0
 n3 = n3 + 1
 70 CONTINUE
 CALL zblpki
 GO TO 20
 90 CONTINUE
 CALL bldpkn (pvact,0,pvect)
 pvect(3) = n1
 CALL wrttrl (pvect)
 CALL CLOSE  (uset,1)
 CALL CLOSE  (pvact,1)
  
 RETURN
END SUBROUTINE calcv
