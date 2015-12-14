SUBROUTINE magbdy
     
!     THIS ROUTINE PICKS UP THE GRIDS ON THE AIR/IRON INTERFACES
!     FROM A PERMBDY CARD,CONVERTS EXTERNAL TO INTERNAL SILS, AND
!     STORES RESULTS ON PERMBD WHICH IS READ IN SSG1. SSG1 WILL NEED TO
!     COMPUTE MAGNETIC LOADS ONLY AT THESE POINTS.
 
!     MAGBDY   GEOM1,HEQEXIN/PERMBD/V,N,IPG $
 
 INTEGER :: buf1,FILE,geom1,eqexin,permbd,sysbuf,permby(2)
 DIMENSION       iz(1),nam(2),mcb(7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,iout
 COMMON /BLANK / ipg
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    nam   / 4HMAGB,4HDY  /
 DATA    geom1 , eqexin,permbd/101,102,201/
 DATA    permby/ 4201,42 /
 
 lcore = korsz(z)
 buf1  = lcore - sysbuf
 lcore = buf1 - 1
 IF (lcore <= 0) GO TO 1008
 
!     SEE IF A PERMBDY CARD EXISTS
 
 ipg  = -1
 FILE = geom1
 CALL preloc (*1001,z(buf1),geom1)
 CALL locate (*10,z(buf1),permby,idx)
 ipg  = 1
 GO TO 20
 
!     NO PERMBDY CARD - RETURN
 
 10 CALL CLOSE (geom1,1)
 RETURN
 
!     READ PERMBDY INTO CORE
 
 20 CALL READ (*1002,*30,geom1,z,lcore,0,npts)
 GO TO 1008
 30 CALL CLOSE (geom1,1)
 
!     READ IN 1ST RECORD OF EQEXIN
 
 lcore = lcore - npts
 ieqex = npts
 CALL gopen (eqexin,z(buf1),0)
 FILE  = eqexin
 CALL READ (*1002,*40,eqexin,z(ieqex+1),lcore,0,neq)
 GO TO 1008
 40 CALL CLOSE (eqexin,1)
 ngrids = neq/2
 lcore  = lcore - neq
 
!     GET THE INTERNAL NUMBER (=SIL NUMBER FOR HEAT TRAMSFER)FOR EACH
!     POINT ON PERMBDY AND STORE IT BACK ONTO EXTERNAL NUMBER,SINCE THE
!     EXTERNAL IS NO LONGER NEEDED
 
 DO  i = 1,npts
   CALL bisloc (*60,iz(i),iz(ieqex+1),2,ngrids,jloc)
   iz(i) = iz(ieqex+jloc+1)
 END DO
 GO TO 70
 
 60 WRITE  (iout,65) ufm,iz(i)
 65 FORMAT (a23,', GRID',i9,' ON PERMBDY CARD DOES NOT EXIST')
 CALL mesage(-61,0,0)
 
!     WRITE THESE INTERNAL ID-S ONTO PERMBD
 
 70 CALL gopen (permbd,z(buf1),1)
 CALL WRITE (permbd,iz(1),npts,1)
 CALL CLOSE (permbd,1)
 mcb(1) = permbd
 mcb(2) = npts
 DO  i = 3,7
   mcb(i) = 0
 END DO
 CALL wrttrl(mcb)
 
 RETURN
 
 1001 n =-1
 GO TO 1010
 1002 n =-2
 GO TO 1010
 1008 n =-8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE magbdy
