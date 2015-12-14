SUBROUTINE fndpnt (iary,id)
     
 
 INTEGER, INTENT(OUT)                     :: iary(4)
 INTEGER, INTENT(IN)                      :: id
 INTEGER :: NAME(2),OLD,bgpdt,sil,edt
 DIMENSION  isave(4),arry(3),iry(3),iedt(2),icore(1), ifed(2)
 COMMON /system/ ibuf,nout
 COMMON /fpt   / dum(3),nrow1,lcore
 COMMON /zzzzzz/ core(1)
 COMMON /loadx / i1(2),bgpdt,OLD,cstm,sil,isil,i2,mpt,gptt,edt, impt,igptt,ied
 EQUIVALENCE     (iry(1),arry(1)), (core(1),icore(1))
 DATA    NAME  / 4HFNDP,4HNT  /
 DATA    iedt  / 4HEDT ,4HFEDT/, ifed/4HFEDT,4HST  /
 
!     FIND POINT ON BGPDT
 
 IF (id < 0) GO TO 90
 IF (id < 268435455 .AND. OLD >= 0) GO TO 10
!               268435455 = 2**28 - 1
 WRITE  (nout,5) id,OLD
 5 FORMAT (//,' BAD DATA PASSED TO FNDPNT, ID,OLD =',2I14)
 CALL mesage (-37,0,NAME)
 10 ns = 4*(id-OLD)
 IF (ns-4 < 0) THEN
   GO TO    70
 ELSE IF (ns-4 == 0) THEN
   GO TO    30
 END IF
 20 CALL READ (*90,*90,bgpdt,isave(1),-ns+4,0,flag)
 30 CALL READ (*90,*90,bgpdt,isave(1),    4,0,flag)
 OLD = id
 40 DO  i = 1,4
   iary(i) = isave(i)
 END DO
 60 RETURN
 
 70 IF (ns == 0) THEN
   GO TO    40
 END IF
 80 CALL bckrec (bgpdt)
 OLD = 0
 GO TO 10
 
 90 ipm = bgpdt
 100 CALL mesage (-2,ipm,NAME)
 110 ipm = sil
 GO TO 100
 120 ipm = edt
 GO TO 100
 
 
 ENTRY fndsil (ip)
!     =================
 
!     FIND SIL VALUE
 
 130 ns = ip - isil
 IF (ns-1 < 0) THEN
   GO TO   140
 ELSE IF (ns-1 == 0) THEN
   GO TO   170
 ELSE
   GO TO   160
 END IF
 140 IF (ns == 0) THEN
   GO TO   180
 END IF
 150 CALL bckrec (sil)
 isil = 0
 GO TO 130
 160 CALL READ (*110,*110,sil,i,-ns+1,0,flag)
 170 CALL READ (*110,*110,sil,IF,   1,0,flag)
 isil = ip
 180 ip   = IF
 GO TO 60
 
 
 ENTRY fedtst (idef)
!     ===================
 
!     FIND ENFORCED DISPLACEMENT
 
!     PUT DEFORM EID S AND VALUES INTO CORE FOR THIS SET
 
 icp = nrow1 + 1
 k   = 0
 CALL READ (*120,*120,edt,arry(1),-3,0,flag)
 200 CALL READ (*120,*210,edt,arry(1), 3,0,flag)
 IF (idef /= iry(1) .AND. k == 0) GO TO 200
 IF (idef /= iry(1)) GO TO 210
 k = k + 2
 core(icp+k   ) = arry(3)
 icore(icp+k-1) = iry(2)
 IF (lcore-nrow1+k <= 0) CALL mesage (-8,ipm,ifed)
 GO TO 200
 210 IF (k == 0) CALL mesage (-32,idef,iedt)
 CALL bckrec (edt)
 GO TO 60
 
 
 ENTRY fedt (ied1,delta,idef)
!     ============================
 
!     FIND VALUE FOR EID IF IT EXISTS
 
 DO  i = 1,k,2
   IF (ied1 /= icore(icp+i)) CYCLE
   icore(icp+i) = -icore(icp+i)
   delta = core(icp+i+1)
   GO TO 60
 END DO
 delta = 0.0
 GO TO 60
 
 
 ENTRY fedted (idef)
!     ===================
 
!     CHECK TO SEE IF ALL ELEMENTS IN THE SET WERE USED
 
 ifound = 0
 DO  i = 1,k,2
   IF (icore(icp+i) < 0) CYCLE
   iedt(1) = icore(icp+i)
   iedt(2) = idef
   CALL mesage (30,139,iedt)
   ifound = 1
 END DO
 IF (ifound == 1) CALL mesage (-61,0,0)
 GO TO 60
END SUBROUTINE fndpnt
