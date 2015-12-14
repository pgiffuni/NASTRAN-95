SUBROUTINE insert(ncol,nrow,ndof,ngrid,jcore,z,dz,temp,dtemp,ipr)
     
! INSERT INSERTS MATRIX PARTITONS INTO OPEN CORE FOR IS2D8
 
 
 INTEGER, INTENT(IN OUT)                  :: ncol
 INTEGER, INTENT(IN OUT)                  :: nrow
 INTEGER, INTENT(IN)                      :: ndof
 INTEGER, INTENT(IN)                      :: ngrid
 INTEGER, INTENT(IN)                      :: jcore
 REAL, INTENT(OUT)                        :: z(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dz(1)
 REAL, INTENT(IN)                         :: temp(9)
 DOUBLE PRECISION, INTENT(IN)             :: dtemp(9)
 INTEGER, INTENT(IN OUT)                  :: ipr
 
 
 
 is1=ngrid*ndof**2
 
! COMPUTE STARTING POINTS INTO OPEN CORE FOR THIS PARTITION AND ITS TRAN
 
 iz1=is1*(nrow-1)+ndof*(ncol-1)+jcore-1
 iz2=is1*(ncol-1)+ndof*(nrow-1)+jcore-1
 
! IZ1 GETS TEMP,  IZ2 GETS THE TRANSPOSE
 
 i1=iz1
 i2=iz2
 
 IF (ipr == 2) GO TO 20
 
 IF(ndof == 1)GO TO 10
 
! 3 X 3 PARTITION
! I1 GETS TEMP. I2 GETS THE TRANSPOSE
! IF I1=I2, THEN HALF OF THE ENTRIES WILL BE DUPLICATED
! THAT-S OK SINCE THERE ARE NO ADDITIONS
 
 z(i1+1)=temp(1)
 z(i2+1)=temp(1)
 z(i1+2)=temp(2)
 z(i2+25)=temp(2)
 z(i1+3)=temp(3)
 z(i2+49)=temp(3)
 z(i1+25)=temp(4)
 z(i2+2)=temp(4)
 z(i1+26)=temp(5)
 z(i2+26)=temp(5)
 z(i1+27)=temp(6)
 z(i2+50)=temp(6)
 z(i1+49)=temp(7)
 z(i2+3)=temp(7)
 z(i1+50)=temp(8)
 z(i2+27)=temp(8)
 z(i1+51)=temp(9)
 z(i2+51)=temp(9)
 GO TO 100
 
! 1 X 1 PARTITION
 
 10 z(i1+1)=temp(1)
 z(i2+1)=temp(1)
 GO TO 100
 
 
! DO THE SAME IN DOUBLE PRECISION
 
 20 IF (ndof == 1) GO TO 30
 
 dz(i1+ 1)=dtemp(1)
 dz(i2+ 1)=dtemp(1)
 dz(i1+ 2)=dtemp(2)
 dz(i2+25)=dtemp(2)
 dz(i1+ 3)=dtemp(3)
 dz(i2+49)=dtemp(3)
 dz(i1+25)=dtemp(4)
 dz(i2+ 2)=dtemp(4)
 dz(i1+26)=dtemp(5)
 dz(i2+26)=dtemp(5)
 dz(i1+27)=dtemp(6)
 dz(i2+50)=dtemp(6)
 dz(i1+49)=dtemp(7)
 dz(i2+ 3)=dtemp(7)
 dz(i1+50)=dtemp(8)
 dz(i2+27)=dtemp(8)
 dz(i1+51)=dtemp(9)
 dz(i2+51)=dtemp(9)
 GO TO 100
 
 30 dz(i1+1)=dtemp(1)
 dz(i2+1)=dtemp(1)
 
 100 RETURN
END SUBROUTINE insert
