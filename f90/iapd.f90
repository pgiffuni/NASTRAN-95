FUNCTION iapd(i,j,nc,ncrd)
 IF(j /= 1) GO TO 10
 iapd=ncrd+1
 IF(i == 1) RETURN
 iapd=iapd+1
 IF(i == 2) RETURN
 iapd=3+3*(i-2)+ncrd
 RETURN
 10 IF(j /= 2) GO TO 20
 iapd=3+ncrd
 IF(i == 1) RETURN
 iapd=4+ncrd
 IF(i == 2) RETURN
 iapd=4+3*(i-2)+ncrd
 RETURN
 20 iapd=j+nc*(2*j-3)+ncrd
 IF(i == 1) RETURN
 iapd=iapd+1
 IF(i == 2) RETURN
 iapd=iapd+2*(i-2)
 RETURN
END FUNCTION iapd
