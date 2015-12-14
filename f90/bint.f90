FUNCTION bint(i,j,a,b,iv,iw,r,z)
     
 
    INTEGER, INTENT(IN OUT)                  :: i
    INTEGER, INTENT(IN OUT)                  :: j
    REAL, INTENT(IN)                         :: a
    REAL, INTENT(IN)                         :: b
    INTEGER, INTENT(IN)                      :: iv
    INTEGER, INTENT(IN)                      :: iw
    REAL, INTENT(IN)                         :: r(1)
    REAL, INTENT(IN OUT)                     :: z(1)
 
    bint = 0.0
    iw1 = iw + 1
    c1p = b
    c2p = a
    c1 = c1p
    c2 = c2p
    aw = 0.0
    IF( r(i) /= 0.0E0 .AND. r(j) /= 0.0E0 ) aw = ALOG(r(j)/r(i))
    DO  it = 1,iw1
        ic = iw - it + 1
        IF (ic == 0) c1 = 1.0
        IF (it == 1) c2 = 1.0
   
        ! THE FOLLOWING CODE REPLACES REAL FUNCTION COEF
   
        IF(it == 1) GO TO 20
        in = 1
        id = 1
        DO  k=2,it
            in = in*(iw-k+2)
            id = id*(k-1)
        END DO
        coef = in/id
        GO TO 30
20      coef = 1.0
30  CONTINUE
   
   
    ! THE FOLLOWING CODE REPLACES REAL FUNCTION AJ
   
    is1 = ic+iv+1
    IF(is1 == 0) GO TO 60
    sp1 = is1
    aj = (r(j)**is1-r(i)**is1) / sp1
    GO TO 70
60  aj = aw
70 CONTINUE
   
   bint = bint + c1 ** ic * aj * c2 ** (it - 1) * coef
   c1 = c1p
   c2 = c2p
   END DO
   aw = iw
   bint = bint / aw

   RETURN
END FUNCTION bint
