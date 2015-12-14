SUBROUTINE gravl2(nvect,fild,pg)
     
 
 INTEGER, INTENT(IN)                      :: nvect
 INTEGER, INTENT(IN OUT)                  :: fild
 INTEGER, INTENT(IN)                      :: pg(7)
 INTEGER :: sysbuf, sil
 INTEGER :: NAME(2)
 
 COMMON /BLANK/nrowsp
 COMMON /system/ sysbuf
 COMMON  /zntpkx/ a(4),ll,ieol
 COMMON /zblpkx/ b(4),ii
 COMMON /zzzzzz/ core(1)
 COMMON  /loadx/ n(2),bgpdt,OLD,cstm,sil,istl,nn(8),mass
 
 DATA NAME/4HGRAV,4HL2  /
 
! ----------------------------------------------------------------------
 
 lcore=korsz(core)
 nz = lcore
 lcore=lcore-sysbuf
 CALL OPEN(*170,pg(1),core(lcore+1),0)
 CALL skpfil (pg,1)
 CALL skpfil (pg,-1)
 CALL CLOSE (pg,2)
 CALL OPEN(*170,pg(1),core(lcore+1),3)
 lcore = lcore-sysbuf
 CALL gopen(fild,core(lcore+1),0)
 lcore=lcore-sysbuf
 CALL gopen( sil,core(lcore+1),0)
 ibuf=lcore
 isil=0
 DO  iloop=1,nvect
   50 CALL READ(*210,*130,sil,isil1,1,0,flag)
   IF(isil1 < 0) THEN
     GO TO    50
   END IF
   60 ASSIGN 100 TO iout
   CALL bldpk(1,1,pg(1),0,0)
   CALL intpk(*150,fild,0,1,0)
   70 CALL READ(*210,*130,sil,isil2,1,0,flag)
   IF(isil2 < 0) THEN
     GO TO    70
   END IF
   80 IF (isil2-isil1-1 == 0) THEN
     GO TO    90
   ELSE
     GO TO   140
   END IF
   90 GO TO iout,(100,150)
   100 IF (ieol /= 0) GO TO 150
   CALL zntpki
   IF (ll-isil1 < 0) THEN
     GO TO   120
   ELSE IF (ll-isil1 == 0) THEN
     GO TO    90
   ELSE
     GO TO    70
   END IF
   120 b(1)=a(1)
   ii=ll
   CALL zblpki
   GO TO 90
   130 ASSIGN 150 TO iout
   IF(nrowsp-isil1 == 0) THEN
     GO TO   150
   END IF
   140 isil1 = 999999
   GO TO 90
   150 CALL REWIND(sil)
   CALL bldpkn(pg(1),0,pg)
   CALL skprec(sil,1)
   isil=0
 END DO
 CALL CLOSE (sil,1)
 CALL CLOSE (fild,1)
 CALL wrttrl (pg)
 CALL CLOSE (pg,1)
 RETURN
 
 170 ipm=pg(1)
 CALL mesage (-1,ipm,NAME)
 
 210 CALL mesage (-3,sil,NAME)
 RETURN
 
END SUBROUTINE gravl2
