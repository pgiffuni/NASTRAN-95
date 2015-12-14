SUBROUTINE gravl1(nvect,gvect,sr1,iharm)
     
 
 INTEGER, INTENT(IN)                      :: nvect
 REAL, INTENT(IN)                         :: gvect(1)
 INTEGER, INTENT(IN OUT)                  :: sr1
 INTEGER, INTENT(IN OUT)                  :: iharm
 INTEGER :: gravt(7),OLD,sysbuf, bgpdt,sil,cstm
 INTEGER :: NAME(2)
 
 DIMENSION  igpco(4), vect(3)
 
 COMMON /BLANK/nrowsp
 COMMON /system/sysbuf
 COMMON /zblpkx/ b(4),ii
 COMMON /zzzzzz/ core(1)
 COMMON  /loadx/ n(2),bgpdt,OLD,cstm,sil,istl,nn(8),mass
 
 DATA NAME/4HGRAV,4HL1  /
 
! ----------------------------------------------------------------------
 
 IF (iharm == 0) GO TO 5
 CALL gravl3(nvect,gvect,sr1,iharm)
 RETURN
 5 CONTINUE
 lcore=korsz(core)
 icm = 1
 nz = lcore
 lcore=lcore-sysbuf
 CALL gopen(sr1,core(lcore+1),1)
 lcore =lcore - sysbuf
 CALL gopen(bgpdt,core(lcore+1),0)
 OLD =0
 lcore =lcore -sysbuf
 CALL OPEN(*10,cstm,core(lcore+1),0)
 icm = 0
 CALL skprec(cstm,1)
 lcore =lcore-sysbuf
 10 CALL gopen(sil,core(lcore+1),0)
 isil=0
 CALL makmcb(gravt,sr1,nrowsp,2,1)
 DO  iloop=1,nvect
   20 CALL READ(*200,*120,sil,isil1,1,0,flag)
   IF(isil1 < 0) THEN
     GO TO    20
   END IF
   30 il=(iloop-1)*3
   ASSIGN 60 TO iout
   ipont=1
   CALL bldpk(1,1,gravt(1),0,0)
   40 CALL READ(*200,*120,sil,isil2,1,0,flag)
   IF(isil2 < 0) THEN
     GO TO    40
   END IF
   50 IF(isil2 -isil1-1 == 0) THEN
     GO TO    60
   ELSE
     GO TO    70
   END IF
   60 isil1 = isil2
   ipont = ipont+1
   GO TO 40
   70 CALL fndpnt (igpco(1),ipont)
   DO  i=1,3
     in= i+il
     vect(i) = gvect(in)
   END DO
   IF (igpco(1) /= 0) CALL basglb (vect(1),vect(1),igpco(2),igpco(1))
   DO  i=1,3
     b(1)=vect(i)
     ii = isil1-1+i
     CALL zblpki
   END DO
   GO TO iout,(60,130)
   
!     END SIL
   
   120 ASSIGN 130 TO iout
   IF(nrowsp-isil1 == 0) THEN
     GO TO   130
   ELSE
     GO TO    70
   END IF
   130 CALL REWIND(bgpdt)
   CALL REWIND(sil)
   CALL bldpkn(gravt(1),0,gravt)
   CALL skprec(sil,1)
   isil=0
   CALL skprec(bgpdt,1)
   OLD=0
 END DO
 CALL CLOSE(bgpdt,1)
 IF(icm == 0) CALL CLOSE(cstm,1)
 CALL CLOSE (sil,1)
 CALL CLOSE (gravt(1),1)
 CALL wrttrl (gravt)
 RETURN
 
 200 CALL mesage (-3,ipm,NAME)
 RETURN
 
END SUBROUTINE gravl1
