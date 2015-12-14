SUBROUTINE shstts (tab,uab,vab)
     
!     TO CREATE STRESS TENSOR TRANSFORMATION MATRICES FROM AN ORTHOGONAL
!     TRANSFORMATION FOR SHELL ELEMENTS.
 
!     INPUT :
!           TAB    - ORTHOGONAL INPLANE ROTATION TRANSFORMATION
!     OUTPUT:
!           UAB    - TENSOR TRANSFORMATION FOR NORMAL AND INPLANE SHEAR
!                    COMPONENTS
!           VAB    - TENSOR TRANSFORMATION FOR OUT-OF-PLANE SHEAR
 
!     USAGE:
!           THE INPUT IS ASSUMED TO BE ROW-LOADED.
!           OUTPUTS ARE CREATED ROW-LOADED.
!           DEFINING:
!           [S]      AS A 2-D STRESS VECTOR;
!           [E]      AS A 2-D STRAIN VECTOR;
!           [Q]      AS A 2-D SHEAR FORCE VECTOR;
!           [G]      AS A 2-D STRESS/FORCE-STRAIN RELATION; AND
!           [ALPHA]  AS A VECTOR OF THERMAL EXPANSION COEFFICIENTS,
 
!           THEN THE FOLLOWING RELATIONSHIPS ARE TRUE:
 
!                       T                        T
!           [S]  = [UAB] [S]         [G]  = [UAB] [G] [UAB]
!              A            B           A            B
 
!                       T
!           [Q]  = [VAB] [Q]
!              A            B
 
!           IF [TBA] IS INPUT, THE OUTPUT WILL BE:
 
!                -1                      -1       T
!           [UAB]  = [UBA],   AND   [VAB]  = [VAB] = [VBA]
 
!           WHICH MAY BE USED IN THE FOLLOWING:
 
!           [E]  = [UBA] [E]         [ALPHA]  = [UBA] [ALPHA]
!              A            B               A                B
 
!           [Q]  = [VBA][Q]
!              A           B
 
 
 
 
 REAL, INTENT(IN)                         :: tab(9)
 REAL, INTENT(OUT)                        :: uab(9)
 REAL, INTENT(OUT)                        :: vab(4)
 
 
 uab(1) = tab(1)*tab(1)
 uab(2) = tab(4)*tab(4)
 uab(3) = tab(1)*tab(4)
 uab(4) = tab(2)*tab(2)
 uab(5) = tab(5)*tab(5)
 uab(6) = tab(2)*tab(5)
 uab(7) = tab(1)*tab(2)*2.0
 uab(8) = tab(4)*tab(5)*2.0
 uab(9) = tab(1)*tab(5) + tab(2)*tab(4)
 
 vab(1) = tab(5)*tab(9)
 vab(2) = tab(2)*tab(9)
 vab(3) = tab(4)*tab(9)
 vab(4) = tab(1)*tab(9)
 
 RETURN
END SUBROUTINE shstts
