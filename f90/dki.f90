DOUBLE PRECISION FUNCTION dki(i,j,k,l,m,n,ip,iq,r,z)
     
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN)                      :: j
 INTEGER, INTENT(IN OUT)                  :: k
 INTEGER, INTENT(IN OUT)                  :: l
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN)                      :: ip
 INTEGER, INTENT(IN)                      :: iq
 DOUBLE PRECISION, INTENT(IN)             :: r(1)
 DOUBLE PRECISION, INTENT(IN)             :: z(1)
 DOUBLE PRECISION :: ai, rd, abs1, amkl, akkl, ammn, akmn, arr
 DOUBLE PRECISION :: dkint, dk89, dk100, dk211
 DOUBLE PRECISION :: xx, amm
 
 
 IF (r(i) == r(j)) GO TO 20
 rd = r(j)
 IF (r(j) == 0.0D0) rd = r(i)
 abs1 =  DABS( (r(i) - r(j)) / rd )
 IF (abs1 <= 0.1D-3) GO TO 20
 amkl = (r(l)*z(k)-r(k)*z(l)) / (r(l)-r(k))
 akkl = (z(l)-z(k)) / (r(l)-r(k))
 ammn = (r(n)*z(m)-r(m)*z(n)) / (r(n)-r(m))
 akmn = (z(n)-z(m)) / (r(n)-r(m))
 IF (akmn /= akkl .OR. ammn /= amkl) GO TO 30
 20 ai = 0.0D0
 GO TO  510
 30 CONTINUE
 iss = IABS(ip)
 irr = IABS(iq)
 IF (iq + 1 < 0) THEN
   GO TO   100
 ELSE IF (iq + 1 == 0) THEN
   GO TO   300
 END IF
 50 CONTINUE
 mm = ip
 nn = iq + 1
 ai =dkint(i,j,ammn,akmn,mm,nn,r,z) -dkint(i,j,amkl,akkl,mm,nn,r,z)
 GO TO  510
 100 CONTINUE
 IF (ip < 0) GO TO 200
 mm = ip
 nn = irr - 1
 ai =dk89(i,amkl,akkl,mm,nn,r)  -  dk89(i,ammn,akmn,mm,nn,r)  &
     -dk89(j,amkl,akkl,mm,nn,r)  +  dk89(j,ammn,akmn,mm,nn,r)
 arr = irr
 ai = (1.0D0 / (1.0D0 - arr)) * ai
 GO TO  510
 200 CONTINUE
 mm = iss
 nn = irr - 1
 ai =dk100(i,amkl,akkl,mm,nn,r) -dk100(i,ammn,akmn,mm,nn,r)  &
     -dk100(j,amkl,akkl,mm,nn,r) +dk100(j,ammn,akmn,mm,nn,r)
 arr = irr
 ai = (1.0D0 / (1.0D0 - arr)) * ai
 GO TO  510
 300 CONTINUE
 IF (ip + 1 < 0) THEN
   GO TO   400
 ELSE IF (ip + 1 == 0) THEN
   GO TO   500
 END IF
 301 CONTINUE
 mm = ip + 1
 amm=mm
 xx=r(i)**mm/amm
 ai=   ( +xx*DLOG(DABS(amkl+akkl*r(i)))-akkl/amm*dk89(i,amkl,akkl,mm,1,r)  &
     -xx*DLOG(DABS(ammn+akmn*r(i)))+akmn/amm*dk89(i,ammn,akmn,mm,1,r) )
 xx=r(j)**mm/amm
 ai=   ( -xx*DLOG(DABS(amkl+akkl*r(j)))+akkl/amm*dk89(j,amkl,akkl,mm,1,r)  &
     +xx*DLOG(DABS(ammn+akmn*r(j)))-akmn/amm*dk89(j,ammn,akmn,mm,1,r) ) + ai
 GO TO  510
 400 CONTINUE
 mm = iss - 1
 amm=mm
 xx=amm*r(i)**mm
 ai=   ( -DLOG(DABS(amkl+akkl*r(i)))/xx+akkl/amm*dk100(i,amkl,akkl,mm,1,r)  &
     +DLOG(DABS(ammn+akmn*r(i)))/xx-akmn/amm*dk100(i,ammn,akmn,mm,1,r) )
 xx=amm*r(j)**m
 ai=   ( +DLOG(DABS(amkl+akkl*r(j)))/xx-akkl/amm*dk100(j,amkl,akkl,mm,1,r)  &
     -DLOG(DABS(ammn+akmn*r(j)))/xx+akmn/amm*dk100(j,ammn,akmn,mm,1,r) ) + ai
 GO TO  510
 500 CONTINUE
 ai = dk211(i,amkl,akkl,r) - dk211(i,ammn,akmn,r)  &
     - dk211(j,amkl,akkl,r) + dk211(j,ammn,akmn,r)
 510 CONTINUE
 dki = ai
 RETURN
END FUNCTION dki
