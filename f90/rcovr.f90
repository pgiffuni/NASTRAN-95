SUBROUTINE rcovr
     
!     MAIN DRIVER FOR PHASE 2 SUBSTRUCTURING RECOVER OPERATION
 
!     THIS MODULE WILL CALCULATE THE DISPLACEMENT AND REACTION MATRICES
!     FOR ANY OF THE SUBSTRUCTURES COMPOSING THE FINAL SOLUTION STRUC-
!     TURE.  OUTPUT DATA MAY BE PLACED ON OFP PRINT FILES OR SAVED ON
!     THE SOF FOR SUBSEQUENT PROCESSING.
 
!     DMAP CALLING SEQUENCES
 
!     RIGID FORMATS 1 AND 2  (STATIC ANALYSIS)
 
!     RCOVR   CASESS,GEOM4,KGG,MGG,PG,UGV,,,,,/OUGV1,OPG1,OQG1,U1,
!             U2,U3,U4,U5/DRY/ILOOP/STEP/FSS/RFNO/0/LUI/U1NM/U2NM/
!             U3NM/U4NM/U5NM/S,N,NOSORT2/V,Y,UTHRESH/V,Y,PTHRESH/
!             V,Y,QTHRESH $
 
!     RIGID FORMAT 3  (MODAL ANALYSIS)
 
!     RCOVR   CASESS,LAMA,KGG,MGG,,PHIG,,,,,/OPHIG,,OQG1,U1,U2,U3,
!             U4,U5/DRY/ILOOP/STEP/FSS/RFNO/NEIGV/LUI/U1NM/U2NM/
!             U3NM/U4NM/U5NM/S,N,NOSORT2/V,Y,UTHRESH/V,Y,PTHRESH/
!             V,Y,QTHRESH $
 
!     RIGID FORMAT 8  (FREQUENCY ANALYSIS)
 
!     RCOVR   CASESS,GEOM4,KGG,MGG,PPF,UPVC,DIT,DLT,BGG,K4GG,PPF/
!             OUGV1,OPG1,OQG1,U1,U2,U3,U4,U5/DRY/ILOOP/STEP/FSS/
!             RFNO/0/LUI/U1NM/U2NM/U3NM/U4UN/U5NM/S,N,NOSORT2/
!             V,Y,UTHRESH/V,Y,PTHRESH/V,Y,QTHRESH $
 
!     RIGID FORMAT 9  (TRANSIENT ANALYSIS)
 
!     RCOVR   CASESS,GEOM4,KGG,MGG,PPT,UPV,DIT,DLT,BGG,K4GG,TOL/
!             OUGV1,OPG1,OQG1,U1,U2,U3,U4,U5/DRY/ILOOP/STEP/FSS/
!             RFNO/0/LUI/U1NM/U2NM/U3NM/U4UN/U5NM/S,N,NOSORT2/
!             V,Y,UTHRESH/V,Y,PTHRESH/V,Y,QTHRESH $
 
!     MRECOVER  (ANY RIGID FORMAT)
 
!     RCOVR   ,,,,,,,,,,/OPHIG,,OQG1,U1,U2,U3,U4,U5/DRY/ILOOP/
!             STEP/FSS/3/NEIGV/LUI/U1NM/U2NM/U3NM/U4NM/U5NM/
!             S,N,NOSORT2/V,Y,UTHRESH/V,Y,PTHRESH/V,Y,QTHRESH $
 
!     MAJOR SUBROUTINES FOR RCOVR ARE -
 
!     RCOVA - COMPUTES THE SOLN ITEM FOR THE FINAL SOLUTION STRUCTURE
!     RCOVB - PERFORMS BACK-SUBSTITUTION TO RECOVER DISPLACEMENTS OF
!             LOWER LEVEL SUBSTRUCTURES FROM THOSE OF THE FINAL SOLUTION
!             STRUCTURE
!     RCOVC - COMPUTES REACTION MATRICES AND WRITES OUTPUT DATA BLOCKS
!             FOR THE OFP
!     RCOVO - PROCESS CASESS FOR THE RCOVER COMMAND AND ANY OUTPUT
!             REQUESTS SPECIFIED
!     RCOVE - COMPUTES MODAL ENERGIES AND ERRORS FOR A MODAL REDUCED
!             SUBSTRUCTURE
 
!     JUNE 1977
 
 INTEGER :: energy
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 
 nosort = -1
 CALL rcovo
 CALL rcova
 IF (iopt < 0) GO TO 10
 CALL rcovb
 IF (iopt <= 0) GO TO 10
 CALL rcovc
 IF (energy /= 0) CALL rcove
 10 RETURN
END SUBROUTINE rcovr
