SUBROUTINE xsave
!     THE PURPOSE OF THIS ROUTINE IS TO PERFORM THE FUNCTIONS ASSIGNED
!     TO THE SAVE DMAP INSTRUCTION.
 
 COMMON/xvps/ ivps(1)
 COMMON/BLANK/ ipar(1)
 COMMON /oscent/ ioscr(7)
!     GET NUMBER OF PARAMETERS FROM OSCAR
 n = ioscr(7)*2 + 6
 DO  i1 = 8,n,2
!     GET VPS POINTER AND POINTER TO VALUE IN BLANK COMMON.
   j = ioscr(i1)
   k = ioscr(i1+1)
!     GET LENGTH OF VALUE FROM VPS
   l = ivps(j-1)
!     TRANSFER VALUE FROM BLANK COMMON TO VPS
   DO  i2 = 1,l
     ivps(j) = ipar(k)
     j = j + 1
     k = k + 1
   END DO
 END DO
 RETURN
END SUBROUTINE xsave
