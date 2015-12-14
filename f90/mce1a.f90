SUBROUTINE mce1a
     
!     MCE1A PARTITIONS RG INTO RM AND RN
 
 INTEGER :: uset   ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,rn    ,  &
     um     ,un    ,ug    ,a     ,a11   ,a21   ,a12   ,a22   ,  &
     rule   ,usetxx,z     ,rect  ,square
 COMMON /BLANK / uset  ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,  &
     rn    ,l     ,u     ,mcb(7)
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua    ,  &
     uf    ,us    ,un    ,ug
 COMMON /parmeg/ a(7)  ,a11(7),a21(7),a12(7),a22(7),n     ,rule
 COMMON /patx  / nz    ,nsub1 ,nsub2 ,nsub3 ,usetxx
 COMMON /zzzzzz/ z(1)
 DATA    rect  / 2 /   ,square / 1 / ,i / 1 /
 
!     GENERATE ROW PARTITIONING VECTOR
 
 nz = korsz(z)
 usetxx = uset
 CALL calcv (scr1,ug,un,um,z)
 
!     GENERATE NULL COLUMN PARTITIONING VECTOR
 
 z(i  ) = 0
 z(i+2) = nsub2
 z(i+7) = 1
 z(i+8) = 2
 z(i+9) =-16777215
 
!     INITIALIZE MATRIX CONTROL BLOCKS
 
 n    = nz
 rule = 0
 a(1) = rg
 CALL rdtrl (a)
 a11(1) = rn
 a11(2) = nsub1
 a11(3) = nsub2
 a11(4) = rect
 a11(5) = a(5)
 a12(1) = rm
 a12(2) = nsub2
 a12(3) = nsub2
 a12(4) = square
 a12(5) = a(5)
 mcb(1) = scr1
 CALL rdtrl (mcb)
 a21(1) = 0
 a22(1) = 0
 
!     PARTITION RG INTO RM AND RN
 
 CALL partn (mcb,z,z)
 
!     WRITE TRAILERS FOR RM AND RN
 
 CALL wrttrl (a12)
 CALL wrttrl (a11)
 RETURN
END SUBROUTINE mce1a
