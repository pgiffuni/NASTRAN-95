SUBROUTINE sdr1b (ipvect,im1,im2,iout,major,sub1,sub2,iuset,  &
        iopt,iout1)
     
 
 INTEGER, INTENT(IN)                      :: ipvect
 INTEGER, INTENT(IN)                      :: im1
 INTEGER, INTENT(IN)                      :: im2
 INTEGER, INTENT(IN)                      :: iout
 INTEGER, INTENT(IN OUT)                  :: major
 REAL, INTENT(IN OUT)                     :: sub1
 REAL, INTENT(IN OUT)                     :: sub2
 INTEGER, INTENT(IN)                      :: iuset
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN)                      :: iout1
 INTEGER :: NAME(2),core(7),sysbuf,ipv1(7)
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ kore(1)
 COMMON /patx  / nz,nsub1,nsub2,nsub3,iuset1
 COMMON /parmeg/ ia(7),ia11(7),ia12(7),ib11(7),ib12(7),nz1,irule
 COMMON /unpakx/ itu1,iiu1,jju1,incr1
 COMMON /packx / itp1,itp2,iip1,jjp1,incr
 EQUIVALENCE     (core(1),kore(1))
 DATA    NAME  / 4HSDR1,4HB   /
 
 
 nz  = korsz(core)
 nz1 = nz
 iuset1 = iuset
 DO  i = 2,7
   ia11(i) = 0
   ia12(i) = 0
   ia(i)   = 0
 END DO
 ia11(1) = im1
 IF (im1 == 0) GO TO 20
 CALL rdtrl (ia11)
 20 ia12(1) = im2
 CALL rdtrl (ia12)
 IF (ia11(1) < 0 .AND. ia12(1) < 0) RETURN
 CALL calcv (ipvect,major,sub1,sub2,core)
 IF (iopt /= 0) GO TO 60
 IF (ia12(1) <= 0) ia12(1) = 0
 30 ib11(1) = 0
 ib12(1) = 0
 ia(3) = nsub1 + nsub2 + nsub3
 ia(2) = MAX0(ia11(2),ia12(2))
 ia(4) = 2
 IF (im2 == 0) ia12(5) = ia11(5)
 iprec = MIN0(1-MOD(ia11(5),2),1-MOD(ia12(5),2))
 itype = 1
 IF (ia11(5) > 2 .OR. ia12(5) > 2) itype = 3
 ia(5) = iprec + itype
 40 irule = 0
 ia(1) = iout
 ipv1(1) = ipvect
 CALL rdtrl (ipv1)
 core(1) = 0
 core(2) = 1
 core(3) = ia(2)
 core(4) = 2
 core(5) = 1
 core(6) = 0
 core(7) = 0
 CALL merge (core,ipv1,core)
 CALL wrttrl (ia)
 RETURN
 
!     EXPAND YS
 
 60 nz = nz - sysbuf
 CALL OPEN (*130,im2,core(nz+1),0)
 nz = nz - sysbuf
 CALL OPEN (*150,iout1,core(nz+1),1)
 CALL fname (im2,core)
 CALL WRITE (iout1,core,2,1)
 ia(1) = im2
 CALL rdtrl (ia)
 noys  = ia(2)
 ia(2) = 0
 ia(1) = iout1
 ia(6) = 0
 ia(7) = 0
 CALL fwdrec (*130,im2)
 nload = ia11(2)
 itu1  = ia(5)
 incr  = 1
 itp1  = itu1
 itp2  = itp1
 incr1 = 1
 DO  i = 1,nload
   IF (i > noys) GO TO 81
   iiu1 = 0
   CALL unpack (*80,im2,core)
   iip1 = iiu1
   jjp1 = jju1
   81 CALL pack( core,iout1,ia)
   CYCLE
   80 core(1) = 0
   core(2) = 0
   core(3) = 0
   core(4) = 0
   iip1 = 1
   jjp1 = 1
   GO TO 81
 END DO
 CALL CLOSE (iout1,1)
 CALL CLOSE (im2,1)
 CALL wrttrl (ia)
 ia12(1) = iout1
 CALL rdtrl (ia12)
 GO TO 30
 
 
 ENTRY sdr1c (ipvect,im1,iout)
!     =============================
 
!     EXPAND ROWS OF IM1 TO D SET SIZE
 
 DO  i = 1,7
   ia12(i) = 0
   ib11(i) = 0
   ib12(i) = 0
 END DO
 ia11(1) = im1
 CALL rdtrl (ia11)
 ia(1) = im1
 CALL rdtrl (ia)
 ia(3) = nsub1 + nsub2 + nsub3
 GO TO 40
 
!     ERROR MESAGES
 
 130 ip1 = -1
 ip2 = im2
 140 CALL mesage (ip1,ip2,NAME)
 150 ip1 = -1
 ip2 = iout1
 GO TO 140
END SUBROUTINE sdr1b
