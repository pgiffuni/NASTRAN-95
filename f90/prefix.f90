SUBROUTINE prefix (iprefx,NAME)
     
 
 INTEGER, INTENT(OUT)                     :: iprefx
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 EXTERNAL        lshift,rshift,orf
 INTEGER :: rshift,orf,rword
 COMMON /system/ junk(38),nbpc,nbpw,ncpw
 DATA    iblnk / 4H     /
 
 iblank = iblnk
 
!     THIS ROUTINE PREFIXES THE TWO WORD VARIABLE 'NAME' WITH THE SINGLE
!     CHARACTER PREFIX 'IPREFX'.
 
!     SET RIGHT HAND PORTION OF WORDS TO ZERO.
 
 lword  = lshift( rshift( NAME(1),nbpw-4*nbpc ) , nbpw-4*nbpc )
 rword  = lshift( rshift( NAME(2),nbpw-4*nbpc ) , nbpw-4*nbpc )
 iprefx = lshift( rshift( iprefx,nbpw-nbpc ) , nbpw-nbpc )
 iblank = rshift( lshift( iblank,4*nbpc ) , 4*nbpc )
 
!     MOVE RIGHT WORD ONE CHARACTER AND PREFIX WITH LAST CHARACTER
!     OF LEFT WORD.
 
 rword = orf( lshift( lword,3*nbpc ) , rshift( rword,nbpc ) )
 rword = lshift( rshift( rword  ,nbpw-4*nbpc ) , nbpw-4*nbpc )
 rword = orf( rword , iblank )
 
!     MOVE LEFT WORD ONE CHARACTER TO RIGHT AND PREFIX WITH INPUT
!     VALUE.
 
 lword = orf( iprefx , rshift( lword,nbpc))
 lword = lshift( rshift( lword  ,nbpw-4*nbpc ) , nbpw-4*nbpc )
 lword = orf( lword , iblank )
 
 NAME(1) = lword
 NAME(2) = rword
 RETURN
END SUBROUTINE prefix
