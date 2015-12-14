SUBROUTINE cmtrce (iertab,iwds,itomny)
     
!     THIS ROUTINE TRACES BACK IMPROPER CONNECTIONS FINDING
!     GIRD POINT IDS FOR INTERNAL POINT  NUMBERS
 
 
 INTEGER, INTENT(IN)                      :: iertab(1)
 INTEGER, INTENT(IN)                      :: iwds
 INTEGER, INTENT(IN OUT)                  :: itomny
 INTEGER :: combo, z,of,iout(6),nam(2)
 
 COMMON /cmb003/ combo(7,5)
 COMMON /zzzzzz/ z(1)
 COMMON /system/ junk,of,ijunk(6),nlpp,ij2(2),line
 DATA    nheqss/ 4HEQSS /
 
 CALL sort(0,0,4,2,iertab(1),iwds)
 ib = 1
 CALL page1
 WRITE(of,2000)
 WRITE(of,2100)
 nline = nline + 5
 
 50 ips = iertab(ib+1)
 nam(1) = combo(ips,1)
 nam(2) = combo(ips,2)
 CALL sfetch(nam,nheqss,1,itest)
 CALL suread(z(1),-1,nout,itest)
 ipt = nout
 
!     READ EQSS FOR EACH COMPONENT
 
 ncomp = 3
 ncomp = z(ncomp)
 ist = ipt + ncomp + 2
 z(ipt+1) = ist
 DO  i=1,ncomp
   CALL suread(z(ist),-1,nout,itest)
   z(ipt+1+i) = nout + ist
   CALL sort(0,0,3,2,z(ist),nout)
   ist = ist + nout
 END DO
 DO  i=ib,iwds,4
   IF( iertab(i+1) /= ips ) GO TO 1000
   loop220:  DO  j=1,2
     ip = iertab(i+1+j)
     DO  jj=1,ncomp
       ii = z(ipt+jj)
       nwds = z(ipt+jj+1) - z(ipt+jj)
       CALL bisloc(*210,ip,z(ii+1),3,nwds/3,iloc)
       iout(3*j) = z(ii+iloc-1)
       iout(3*j-2) = z(2*jj+3)
       iout(3*j-1) = z(2*jj+4)
       CYCLE loop220
     END DO
   END DO loop220
   line = line + 1
   IF( line <= nlpp ) GO TO 230
   CALL page1
   WRITE(of,2100)
   line = line + 2
   230 CONTINUE
   WRITE(of,2200) iertab(i),iout
 END DO
 GO TO 1100
 
!     GET NEXT PSEUDOSTRUCUTRE
 
 1000 ib = i
 GO TO 50
 1100 IF( itomny == 0 ) RETURN
 WRITE(of,2300)
 RETURN
 2000 FORMAT(/1X,  &
     61HTHE following connections have been found TO be inconsistant.,  &
     /1X,57HATTEMPTS have been made TO connect internal points within,  &
     /1X,57HTHE same pseudostructure due TO split degrees of freedom.,  &
     /1X,79HTHESE errors must be resolved by the user via reles DATA o  &
     r manual connections. /)
 2100 FORMAT(5X,3HDOF,5X,12HSUBSTRUCTURE,5X,8H grid id,5X,  &
     12HSUBSTRUCTURE,5X,8H grid id   /)
 2200 FORMAT(6X,i1,10X,2A4,5X,i8,9X,2A4,5X,i8)
 2300 FORMAT(/5X,93HTHE NUMBER of fatal messages exceeded the available  &
     storage. some messages have been deleted. )
END SUBROUTINE cmtrce
