BLOCK DATA ifx3bd
!IFX3BD
!     THIS TABLE CONTAINS TWO WORDS PER ENTRY (CARD TYPE)
!     FIRST  WORD IS USED AS THE CONICAL SHELL FLAG, AND
!     SECOND WORD IS USED INTERNALLY TO STORE THE NUMBER OF WORDS TO
!     BE OUTPUT TO THE GINO OUTPUT FILE
 
 COMMON /ifpx3/ i1(100),i2(100),i3(100),i4(100),i5(100),  &
     i6(100),i7(100),i8( 40)
 DATA i1/ -1, 0,    -1, 0,     0, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,     1, 0,  &
     1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,     1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,     0, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0, -1, 0,    -1, 0/
 DATA i2/ 0, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,     0, 0,     0, 0,    -1, 0,    -1, 0,  &
     -1, 0,     0, 0,     1, 0,     1, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,    -1, 0,     0, 0,    -1, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,    -1, 0,     1, 0, 0, 0,     0, 0/
 DATA i3/ 0, 0,    -1, 0,     0, 0,    -1, 0,    -1, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,    -1, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,    -1, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0, 0, 0,     0, 0/
 DATA i4/ 0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,     0, 0,     0, 0,     0, 0,    -1, 0,    -1, 0,  &
     -1, 0,     0, 0,     0, 0,     1, 0,     0, 0,    -1, 0,  &
     0, 0,     0, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0, -1, 0,     0, 0/
 DATA i5/ -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     0, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,     0, 0,     0, 0,    -1, 0,    -1, 0, -1, 0,    -1, 0/
 DATA i6/ -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,     0, 0,     0, 0,  &
     0, 0,     0, 0,     0, 0,     0, 0,    -2, 0,    -2, 0,  &
     -2, 0,    -2, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0, -2, 0,    -2, 0/
 DATA i7/ -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,     0, 0,     0, 0,     0, 0,     0, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,     0, 0,    -1, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0,    -1, 0, -1, 0,    -1, 0/
 DATA i8/ 0, 0,    -1, 0,    -1, 0,    -1, 0,     0, 0,    -1, 0,  &
     -1, 0,    -1, 0,    -1, 0,    22* 0/
 
END
