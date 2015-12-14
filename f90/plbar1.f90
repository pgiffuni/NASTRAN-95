SUBROUTINE plbar1 (ido,lcore)
     
!     THIS SUBROUTINE SETS UP THE DATA NEEDED TO CALL PLOAD1
!     TO GET THE APPLIED CONCENTRATED, UNIFORMLY OR LINEARLY DISTRIBUTED
!     LOADS, ON A BAR ELEMENT FROM A PLOAD1 CARD
!     AND INSERTS THE VECTOR INO PV
 
 
 INTEGER, INTENT(IN OUT)                  :: ido
 INTEGER, INTENT(IN OUT)                  :: lcore
 INTEGER :: bar,nam(2),oldid,est,iz(1),islt(7)
 DIMENSION       ta(9),tb(9),pa(6),pb(6),pg(42)
 COMMON /system/ ibuff,nout
 COMMON /ssga1x/ pv(1)
 COMMON /loadx / lc,slt,d1(5),est,d2(11),ilid
 COMMON /matin / matid,inflag,temp
 EQUIVALENCE     (pg(1),iz(1))
 DATA    nam   / 4HPLBA,4HR1  /, n,oldid,islt / 9*0 /
 DATA    iect  , iept,ibg,nwds,bar / 1,16,34,42,34  /
 
!     INITIALIZE AND OPEN EST
 
 IF (n /= 0) GO TO 30
 CALL gopen (est,iz(lcore),0)
 10 CALL READ (*110,*20,est,i,1,0,flag)
 20 IF (i == bar) GO TO 30
 CALL fwdrec (*110,est)
 GO TO 10
 
!     READ SLT THEN FIND BAR ELEMENT
 
 30 CALL READ (*100,*100,slt,islt,7,0,flag)
 IF (islt(1) == oldid) GO TO 60
 40 CALL READ (*110,*110,est,iz(iect),nwds,0,flag)
 oldid = iz(iect)
 IF (iz(iect)-islt(1) < 0) THEN
   GO TO    40
 ELSE IF (iz(iect)-islt(1) == 0) THEN
   GO TO    50
 ELSE
   GO TO   110
 END IF
 
!     CONVERT COORD. SYSTEMS
 
 50 IF (iz(iect+6) /= 0) CALL glbbas (pg(iect+ 3),pg(iect+ 3),  &
     pg(ibg + 1),iz(iect+ 6))
 IF (iz(ibg   ) /= 0) CALL glbbas (pg(iect+ 9),pg(iect+ 9),  &
     pg(ibg + 1),iz(ibg    ))
 IF (iz(ibg +4) /= 0) CALL glbbas (pg(iect+12),pg(iect+12),  &
     pg(ibg + 5),iz(ibg + 4))
 CALL gbtran (iz(ibg  ),iz(ibg+1),ta)
 CALL gbtran (iz(ibg+4),iz(ibg+5),tb)
 
!     DATA READY
 
 inflag = 1
 temp   = pg(ibg+8)
 matid  = iz(iept)
 CALL mat (oldid)
 60 CALL pload1 (1,islt,pg(iect+3),pg(iect+9),pg(iect+12),pg(ibg+1),  &
     pg(ibg+5),pa,pb,ta,tb,islt,iz(iect))
 
!     INSERT INTO PV
 
 ipg = iz(iect+1) - 1
 DO  i = 1,6
   pv(ipg+i) = pv(ipg+i) + pa(i)
 END DO
 ipg = iz(iect+2) - 1
 DO  i = 1,6
   pv(ipg+i) = pv(ipg+i) + pb(i)
 END DO
 n = n + 1
 IF (n /= ido) GO TO 150
 n = 0
 oldid = 0
 CALL CLOSE (est,1)
 GO TO 150
 
!     ERROR
 
 100 CALL mesage (-1,slt,nam)
 110 WRITE  (nout,120) islt(1),ilid
 120 FORMAT ('0*** USER FATAL MESSAGE 2286, CBAR ELEMENT',i9,  &
     ' REFERENCED ON PLOAD1',i9,' NOT FOUND')
 CALL mesage (-61,0,nam)
 
 150 RETURN
END SUBROUTINE plbar1
