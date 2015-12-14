SUBROUTINE amgrod(d,beta)
     
 REAL, INTENT(IN OUT)                     :: d(1)
 REAL, INTENT(IN OUT)                     :: beta
 INTEGER :: NAME(2),sysbuf
 INTEGER :: scr1,scr2,ecore
!  D IS REALLY A 2-D ARRAY D(2,NTZS)
 
 COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
     iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
     icg,ixij,ix,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria  &
     ,inasb,ifla1,ifla2,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5
 COMMON /zzzzzz / z(1)
 COMMON /system/ sysbuf
 DATA NAME /4HAMGR,4HOD  /
 
 CALL sswtch(30,iprnt)
 nfzb = 1
 nlzb = nbz
 nfyb = nb+1-nby
 nlyb = nb
 ibuf1 = ecore - sysbuf
 
!     CALCULATE DZ ON SCR1
 
 IF(ntzs == 0) GO TO 100
 IF(next+2*ntzs > ibuf1) CALL mesage(-8,0,NAME)
 CALL gopen(scr1,z(ibuf1),1)
 idzdy = 0
 CALL dzymat(d,nfzb,nlzb,ntzs,idzdy,scr1,z(ix),beta,iprnt,z(inb),  &
     z(inc),z(iys),z(izs),z(isg),z(icg),z(iyb),z(izb),z(inbea1))
 CALL CLOSE(scr1,1)
 100 IF(ntys == 0) GO TO 200
 IF(next+2*ntys > ibuf1) CALL mesage(-8,0,NAME)
 CALL gopen(scr2,z(ibuf1),1)
 idzdy = 1
 CALL dzymat(d,nfyb,nlyb,ntys,idzdy,scr2,z(ix),beta,iprnt ,z(inb),  &
     z(inc),z(iys),z(izs),z(isg),z(icg),z(iyb),z(izb),z(inbea1))
 CALL CLOSE(scr2,1)
 200 RETURN
END SUBROUTINE amgrod
