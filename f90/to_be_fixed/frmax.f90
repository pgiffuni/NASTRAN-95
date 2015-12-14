SUBROUTINE frmax(ifk,ifm,n,ipr,rsn,rsm)
     
 INTEGER, INTENT(IN OUT)                  :: ifk
 INTEGER, INTENT(IN OUT)                  :: ifm
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: ipr
 DOUBLE PRECISION, INTENT(OUT)            :: rsn
 DOUBLE PRECISION, INTENT(OUT)            :: rsm
 DIMENSION          zk(1)  ,zm(1)
 DOUBLE PRECISION :: ratio  ,ratinv  , dzk(1) ,dzm(1)
 EQUIVALENCE       (dzk(1) ,zk(1) ),(dzm(1) ,zm(1) )
 COMMON  /unpakx/   iprc   ,ip     ,np     ,incr
 
 iprc = ipr
 incr = 1
 rsn = 0.d0
 rsm = 0.d0
 DO  i = 1,n
   ip = i
   np = i
   CALL unpack(*30,ifk,dzk(1))
   CALL unpack(*30,ifm,dzm(1))
   IF(ipr == 2) GO TO 10
   IF (zk(1) == 0.OR.zm(1) == 0) CYCLE
   ratio = zk(1) / zm(1)
   GO TO 20
   10 IF (dzk(1) == 0.0D0.OR.dzm(1) == 0.0D0) CYCLE
   ratio = dzk(1)/dzm(1)
   20 ratinv = 1.d0 /ratio
   IF(ratio > rsm) rsm = ratio
   IF(ratinv > rsn) rsn = ratinv
 END DO
 rsn = 1.d0 / rsn
 RETURN
END SUBROUTINE frmax
