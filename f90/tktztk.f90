SUBROUTINE tktztk(tk,z,nz,l,m,n)
     
!     THIS ROUTINE  PERFORMS A COORDINATE TRANSFORMATION ON THE
!     SYMMETRIC HALF OF A 3 BY 3 MATRIX
 
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: tk(3,3)
 DOUBLE PRECISION, INTENT(IN OUT)         :: z(1)
 INTEGER, INTENT(IN OUT)                  :: nz
 INTEGER, INTENT(IN OUT)                  :: l
 INTEGER, INTENT(IN OUT)                  :: m
 INTEGER, INTENT(IN OUT)                  :: n
 
 tk(1,1)=z(nz  )*z(l  )+z(nz+3)*z(l+1)+z(nz+6)*z(l+2)
 tk(2,1)=z(nz+1)*z(l  )+z(nz+4)*z(l+1)+z(nz+7)*z(l+2)
 tk(3,1)=z(nz+2)*z(l  )+z(nz+5)*z(l+1)+z(nz+8)*z(l+2)
 tk(1,2)=z(nz  )*z(l+1)+z(nz+3)*z(m  )+z(nz+6)*z(m+1)
 tk(2,2)=z(nz+1)*z(l+1)+z(nz+4)*z(m  )+z(nz+7)*z(m+1)
 tk(3,2)=z(nz+2)*z(l+1)+z(nz+5)*z(m  )+z(nz+8)*z(m+1)
 tk(1,3)=z(nz  )*z(l+2)+z(nz+3)*z(m+1)+z(nz+6)*z(n  )
 tk(2,3)=z(nz+1)*z(l+2)+z(nz+4)*z(m+1)+z(nz+7)*z(n  )
 tk(3,3)=z(nz+2)*z(l+2)+z(nz+5)*z(m+1)+z(nz+8)*z(n  )
 z(l  )=z(nz  )*tk(1,1)+z(nz+3)*tk(1,2)+z(nz+6)*tk(1,3)
 z(l+1)=z(nz  )*tk(2,1)+z(nz+3)*tk(2,2)+z(nz+6)*tk(2,3)
 z(l+2)=z(nz  )*tk(3,1)+z(nz+3)*tk(3,2)+z(nz+6)*tk(3,3)
 z(m  )=z(nz+1)*tk(2,1)+z(nz+4)*tk(2,2)+z(nz+7)*tk(2,3)
 z(m+1)=z(nz+1)*tk(3,1)+z(nz+4)*tk(3,2)+z(nz+7)*tk(3,3)
 z(n  )=z(nz+2)*tk(3,1)+z(nz+5)*tk(3,2)+z(nz+8)*tk(3,3)
 RETURN
END SUBROUTINE tktztk
