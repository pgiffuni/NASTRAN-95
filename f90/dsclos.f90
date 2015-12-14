SUBROUTINE dsclos ( iunit )
     INCLUDE 'DSIOF.COM'
!      print *,' dsclos,iunit=',iunit
 CLOSE ( iunit )
 numcls = numcls + 1
 RETURN
END SUBROUTINE dsclos
