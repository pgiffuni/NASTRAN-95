SUBROUTINE gkad1a (usetd,GO,god,scr1,ue,ua,ud)
     
!     GKAD1A WILL EXPAND GO BY NULL MATRIX TO MAKE GOD, AND
!     AA-S TO D-S ADDING ZEROS FOR E-S
 
 
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN)                      :: GO
 INTEGER, INTENT(IN)                      :: god
 INTEGER, INTENT(IN)                      :: scr1
 REAL, INTENT(IN OUT)                     :: ue
 REAL, INTENT(IN OUT)                     :: ua
 REAL, INTENT(IN OUT)                     :: ud
 INTEGER :: uset1, ipv1(7), core,baa,b1dd
 COMMON /patx  / lc,n,no,n4,uset1
 COMMON /zzzzzz/ core(1)
 COMMON /parmeg/ ia(7),ia11(7),ia12(7),ib11(7),ib12(7),nz,irule
 COMMON /system/ idum(54),iprec
 
 
 ient = 0
 
!     COMPUTE CORE FOR CALCV AND MERGE
 
 20 lc = korsz(core)
 
!     BUILD PART VECTOR
 
 uset1 = usetd
 CALL calcv (scr1,ud,ua,ue,core(1))
 
!     SET UP FOR MERGE
 
 nz = lc
 irule = 0
 DO  i = 1,7
   ia11 (i) = 0
   ia   (i) = 0
   ia12 (i) = 0
   ib11 (i) = 0
   ib12 (i) = 0
 END DO
 ipv1(1) = scr1
 CALL rdtrl (ipv1)
 IF (ient /= 0) GO TO 30
 
!     SET UP FOR 2 WAY MERGE
 
 ia11(1) = GO
 CALL rdtrl (ia11)
 ia(1) = god
 ia(2) = n+no+n4
 ia(3) = ia11(3)
 ia(4) = ia11(4)
 ia(5) = ia11(5)
!     BUILD NULL COLUMN IN CORE
 i = 0
 core(  1) = 0
 core(i+2) = 1
 core(i+3) = ia(3)
 core(i+4) = 2
 core(i+5) = 1
 core(i+6) = 0
 core(i+7) = 0
 CALL merge (ipv1(1),core(1),core(1))
 CALL wrttrl (ia)
 40 RETURN
 
 
 ENTRY gkad1b (usetd,kaa,maa,baa,k4aa,k1dd,m1dd,b1dd,k41dd,ua,ue, ud,scr1)
!     ================================================================
 
 ient = 1
 GO TO 20
 
!     VECTOR MADE, SET UP MCB-S
 
 30 ia(2) = n+no+n4
 ia(3) = ia(2)
 ia(4) = 6
 ia(5) = iprec
 ia11(1) = kaa
 ia(1) = k1dd
 iout  = 1
 CALL rdtrl (ia11)
 IF (ia11(1) > 0) GO TO 35
 k1dd = 0
 GO TO 31
 35 CALL merge (ipv1(1),ipv1(1),core(1))
 CALL wrttrl (ia)
 31 SELECT CASE ( iout )
   CASE (    1)
     GO TO 32
   CASE (    2)
     GO TO 33
   CASE (    3)
     GO TO 34
   CASE (    4)
     GO TO 40
 END SELECT
 32 iout  = 2
 ia(1) = b1dd
 ia11(1) = baa
 CALL rdtrl (ia11)
 IF (ia11(1) > 0) GO TO 35
 b1dd = 0
 GO TO 31
 33 iout  = 3
 ia(1) = m1dd
 ia11(1) = maa
 CALL rdtrl (ia11)
 IF (ia11(1) > 0) GO TO 35
 m1dd = 0
 GO TO 31
 34 iout  = 4
 ia(1) = k41dd
 ia11(1) = k4aa
 CALL rdtrl (ia11)
 IF (ia11(1) > 0) GO TO 35
 k41dd = 0
 GO TO 31
END SUBROUTINE gkad1a
