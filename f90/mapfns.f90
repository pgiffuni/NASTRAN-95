FUNCTION mapfns (i)
     
!     THIS FUNCTION PROVIDES ENTRIES FOR VARIOUS FUNCTIONS
!     ON THE VAX VERSION OF NASTRAN
!     (THIS ROUTINE WAS PREVIOUSLY CALLED 'VAXFNS')
 
 
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER :: AND, andf, complf, orf, rshift, xorf
 COMMON /machin/ m(3), lqro
 
 mapfns = 0
 RETURN
 
 ENTRY AND (i,j)
!     ==============
 AND = IAND(i,j)
 RETURN
 
 ENTRY andf (i,j)
!     ================
 andf = IAND(i,j)
 RETURN
 
 ENTRY complf (i)
!     ================
 complf = NOT(i)
 RETURN
 
 ENTRY locfx (i)
!     ===============
 k = lqro/1000
 locfx = loc(i)/k
 RETURN
 
 ENTRY lshift (i,j)
!     ==================
 lshift = ishft(i,j)
 RETURN
 
 ENTRY orf (i,j)
!     ===============
 orf = ior (i,j)
 RETURN
 
 ENTRY rshift (i,j)
!     ==================
 rshift = ishft(i,-j)
 RETURN
 
 ENTRY xorf (i,j)
!     ================
 xorf = IEOR (i,j)
 RETURN
 
END FUNCTION mapfns
