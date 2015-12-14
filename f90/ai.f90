FUNCTION ai (i,j,k,l,m,n,ip,iq,r,z)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN)                      :: j
 INTEGER, INTENT(IN OUT)                  :: k
 INTEGER, INTENT(IN OUT)                  :: l
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN)                      :: ip
 INTEGER, INTENT(IN)                      :: iq
 REAL, INTENT(IN)                         :: r(1)
 REAL, INTENT(IN)                         :: z(1)
 
 
 IF (r(i) == r(j)) GO TO 20
 rd   = r(j)
 IF (r(j) == 0.0) rd = r(i)
 abs1 = ABS((r(i)-r(j))/rd)
 IF (abs1 <= .0001) GO TO 20
 amkl = (r(l)*z(k)-r(k)*z(l))/(r(l)-r(k))
 akkl = (z(l)-z(k))/(r(l)-r(k))
 ammn = (r(n)*z(m)-r(m)*z(n))/(r(n)-r(m))
 akmn = (z(n)-z(m))/(r(n)-r(m))
 IF (akmn /= akkl .OR. ammn /= amkl) GO TO 30
 20 ai  = 0.0
 GO TO 510
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
 ai = bint(i,j,ammn,akmn,mm,nn,r,z) - bint(i,j,amkl,akkl,mm,nn,r,z)
 GO TO 510
 100 CONTINUE
 IF (ip < 0) GO TO 200
 mm = ip
 nn = irr - 1
 ai = f89(i,amkl,akkl,mm,nn,r) - f89(i,ammn,akmn,mm,nn,r)  &
     - f89(j,amkl,akkl,mm,nn,r) + f89(j,ammn,akmn,mm,nn,r)
 arr= irr
 ai = (1.0/(1.0 - arr))*ai
 GO TO 510
 200 CONTINUE
 mm = iss
 nn = irr - 1
 ai = ff100(i,amkl,akkl,mm,nn,r) -ff100(i,ammn,akmn,mm,nn,r)  &
     - ff100(j,amkl,akkl,mm,nn,r) +ff100(j,ammn,akmn,mm,nn,r)
 arr= irr
 ai = (1.0/(1.0-arr))*ai
 GO TO 510
 300 CONTINUE
 IF (ip + 1 < 0) THEN
   GO TO   400
 ELSE IF (ip + 1 == 0) THEN
   GO TO   500
 END IF
 301 CONTINUE
 mm = ip + 1
 amm= mm
 xx = r(i)**mm/amm
 ai = ( +xx*ALOG(ABS(amkl+akkl*r(i)))-akkl/amm*f89(i,amkl,akkl,mm,1,r)  &
     -xx*ALOG(ABS(ammn+akmn*r(i)))+akmn/amm*f89(i,ammn,akmn,mm,1,r) )
 xx = r(j)**mm/amm
 ai = ( -xx*ALOG(ABS(amkl+akkl*r(j)))+akkl/amm*f89(j,amkl,akkl,mm,1,r)  &
     +xx*ALOG(ABS(ammn+akmn*r(j)))-akmn/amm*f89(j,ammn,akmn,mm,1,r) ) + ai
 GO TO 510
 400 CONTINUE
 mm = iss - 1
 amm= mm
 xx = amm*r(i)**mm
 ai = ( -ALOG(ABS(amkl+akkl*r(i)))/xx+akkl/amm*ff100(i,amkl,akkl,mm,1,r)  &
     +ALOG(ABS(ammn+akmn*r(i)))/xx-akmn/amm*ff100(i,ammn,akmn,mm,1,r) )
 xx = amm*r(j)**m
 ai = ( +ALOG(ABS(amkl+akkl*r(j)))/xx-akkl/amm*ff100(j,amkl,akkl,mm,1,r)  &
     -ALOG(ABS(ammn+akmn*r(j)))/xx+akmn/amm*ff100(j,ammn,akmn,mm,1,r) ) + ai
 GO TO 510
 500 CONTINUE
 ai = f6211(i,amkl,akkl,r) - f6211(i,ammn,akmn,r)  &
     - f6211(j,amkl,akkl,r) + f6211(j,ammn,akmn,r)
 510 CONTINUE
 RETURN
END FUNCTION ai
