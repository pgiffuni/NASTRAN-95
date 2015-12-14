SUBROUTINE dsgncl
     INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 idsn =  mdsfcb( 2, ifilex )
 CALL dsclos( idsn )
 mdsfcb( 1,idsn ) = IAND( mdsfcb( 1,idsn ), maskh1 )
 RETURN
END SUBROUTINE dsgncl
