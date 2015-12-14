SUBROUTINE dssend ( FILE )
     
! DSSEND (Dataset Set to End) will position a file to the end
! to allow for closing a file for read and opening it for write
! append.  This eliminates having to read sequentially to the end
! of the file before closing for read.
 
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 
 NAME  = FILE
 CALL dsgefl
 
! GET LAST BLOCK NUMBER IN THIS FILE FROM FCB
 
 nblock = fcb( 6, ifilex )
 
! GET CURRENT BLOCK NUMBER IN THIS FILE FROM FCB
 
 icblk  = fcb( 4, ifilex )
 IF ( icblk == nblock ) GO TO 10
 CALL dbmmgr( 6 )
 10    CONTINUE
 indclr = ibase( indbas+4) + indbas - 1
 indcbp = indclr
 CALL dssdcb
 RETURN
END SUBROUTINE dssend
