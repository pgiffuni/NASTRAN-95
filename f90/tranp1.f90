SUBROUTINE tranp1 (in,iout,nscrth,is1,is2,is3,is4,is5,is6,is7,is8)
     
!     DRIVER OF THE OUT-OF-CORE MATRIX TRANSPOSE ROUTINE TRNSP
!     (DTRANP IS THE TRNSP MODULE DRIVER)
 
 
 
 INTEGER, INTENT(IN)                      :: in
 INTEGER, INTENT(IN)                      :: iout
 INTEGER, INTENT(IN)                      :: nscrth
 INTEGER, INTENT(IN)                      :: is1
 INTEGER, INTENT(IN)                      :: is2
 INTEGER, INTENT(IN)                      :: is3
 INTEGER, INTENT(IN)                      :: is4
 INTEGER, INTENT(IN)                      :: is5
 INTEGER, INTENT(IN)                      :: is6
 INTEGER, INTENT(IN)                      :: is7
 INTEGER, INTENT(IN)                      :: is8
 INTEGER :: scr,nam(2)
 COMMON /zzzzzz/ core(1)
 COMMON /trnspx/ ia(7),iat(7),lcore,nscrh,scr(8)
 DATA    nam   / 4HTRNS,4HP1  /
 
 IF (nscrth > 8) CALL mesage (-37,0,nam)
 ia(1)  = in
 CALL rdtrl (ia)
 iat(1) = iout
 iat(2) = ia(3)
 iat(3) = ia(2)
 iat(5) = ia(5)
 iat(4) = ia(4)
 
!     REVERSE THE FORM OF THE LOWER AND UPPER TRIANGULAR MATRIX
 
 IF (ia(4) == 4) iat(4) = 5
 IF (ia(4) == 5) iat(4) = 4
 lcore  = korsz(core)
 nscrh  = nscrth
 scr(1) = is1
 scr(2) = is2
 scr(3) = is3
 scr(4) = is4
 scr(5) = is5
 scr(6) = is6
 scr(7) = is7
 scr(8) = is8
 CALL trnsp  (core)
 CALL wrttrl (iat)
 RETURN
END SUBROUTINE tranp1
