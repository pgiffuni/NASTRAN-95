SUBROUTINE ampc1(INPUT,output,ncol,z,mcb)
     
!     THE PURPOSE OF THIS ROUTINE IS TO COPY NCOL COLUMNS FROM INPUT
!      TO OUTPUT VIA UNPACK AND PACK.
 
!     THE PACK COMMON BLOCKS HAVE BEEN INITIALIZED OUTSIDE THE ROUTINE
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: output
 INTEGER, INTENT(IN)                      :: ncol
 INTEGER, INTENT(IN OUT)                  :: z(1)
 INTEGER, INTENT(IN OUT)                  :: mcb(7)
 
 
 COMMON /packx/it1,it2,ii,nn,incr
 
!-----------------------------------------------------------------------
 
 DO  i=1,ncol
   CALL unpack(*20,INPUT,z)
   CALL pack(z,output,mcb)
   CYCLE
   
!     NULL COLUMN
   
   20 CALL bldpk(it1,it2,output,0,0)
   CALL bldpkn(output,0,mcb)
 END DO
 RETURN
END SUBROUTINE ampc1
