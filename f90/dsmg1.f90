SUBROUTINE dsmg1
     
!     THIS ROUTINE IS THE DRIVER FOR THE DIFFERENTIAL STIFFNESS MATRIX
!     GENERATOR MODULE OF THE NASTRAN SYSTEM.  SUBROUTINE DS1 APPENDS
!     TEMPERATURE, ELEMENT DEFORMATION AND DISPLACEMENT INFORMATION TO
!     THE ECPT DATA BLOCK AND A SCRATCH FILE, ECPTDS, OF THIS MERGED
!     INFORMATION IS CREATED.  SUBROUTINE DS1A IS STRUCTURED IDENTICALLY
!     TO SMA1A. IT READS THE ECPTDS FILE AND CREATES A SECOND ORDER
!     APPROXIMATION TO THE KGG, WHICH IS CALLED KDGG.
 
!     DMAP CALL -
 
!     DSMG1    CASECC,GPTT,SIL,EDT,UGV,CSTM,MPT,ECPT,GPCT,DIT/KDGG/
 
 CALL ds1 (iarg)
 IF (iarg > 0) GO TO 10
 
!     ECPTDS IS EMPTY. WRITE MESSAGE AND CALL EXIT.
 
 CALL mesage (30,81,0)
 CALL mesage (-61,0,0)
 10 CALL ds1a
 RETURN
END SUBROUTINE dsmg1
