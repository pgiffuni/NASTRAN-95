SUBROUTINE modacc
     
!     THIS IS THE MODULE MODACC
 
!     DMAP CALL
 
!     MODACC  CASECC,TOL,UDV1T,PPT,PDT,PST/TOL1,UDV3T,PP3,PDT3,PST3/
!             C,N,TRAN $
 
!     THE PURPOSE OF THIS MODULE IS TO REDUCE THE COLUMN LENCTHS OF
!     UDV1T,PPT,PDT,PST  TO THE  LENGTH SPECIFIED BY OFREQ IN CASECC.
!     THE CURRENT LIST OF TIMES IS ON  TOL
 
 INTEGER :: casecc, tol,udv1t,ppt,pdt,pst,tol1,udv3t,pp3,pdt3,pst3
 COMMON /BLANK / iop(2)
 COMMON /modac3/ nfo,nfn,nz,id
 COMMON /zzzzzz/ iz(1)
 DATA    casecc, tol,udv1t,ppt,pdt,pst,tol1,udv3t,pp3,pdt3,pst3 /  &
     101   , 102,103  ,104,105,106,201 ,202  ,203,204 ,205  /
 DATA    itran / 4HTRAN/  ,iceign /4HCEIG /
 DATA    ireig / 4HREIG/
 DATA    istat / 4HSTAT/
 
 id = 1
 IF (iop(1) == iceign) id = 2
 IF (iop(1) ==  itran) id = 3
 IF (iop(1) ==  ireig) id = 4
 IF (iop(1) ==  istat) id = 5
 
!     FOR EIGENVALUES STOP LIST AT NUMBER OF VECTORS
 
 nfo = 0
 iz(1) = udv1t
 CALL rdtrl(iz)
 j   = 2
 nfo = 2 * iz(j)
 nz  = korsz(iz(1))
 
!     BUILD LIST OF NEW TIMES, KEEP/REMOVE LIST
 
 CALL modac1 (casecc,tol,tol1,pp3,ppt)
 
!     COPY DISPLACEMENTS
 
 id1 = 1
 IF (id == 3) id1 = 3
 CALL modac2 (id1,udv1t,udv3t)
 IF (id == 2 .OR. id == 4) RETURN
 
!     COPY P LOAD S  (+ HEAD STUFF FOR NOW)
 
 CALL modac2 (-1,ppt,pp3)
 
!     COPY D LOADS
 
 CALL modac2 (1,pdt,pdt3)
 
!     COPY S LOADS
 
 CALL modac2 (1,pst,pst3)
 RETURN
END SUBROUTINE modacc
