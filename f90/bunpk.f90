INTEGER FUNCTION bunpk (ig,i,j)
     
    !     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
    !     UNPACK INTEGER GRID NO. FROM IG TABLE.   SEE BPACK FOR PACKING
    !     USE APPROP. PORTION OF THIS ROUTINE FOR DIFFERENT TYPE OF MACHINE.
 
 
    INTEGER*2, INTENT(IN)                    :: ig(1)
    INTEGER, INTENT(IN)                      :: i
    INTEGER, INTENT(IN)                      :: j
 
 
    COMMON /bandb /  nbit,     dum3b(3), ipass,    nw,       dum1b, nbpw
    COMMON /bands /  dum4s(4), ii1,      dum5s(5), mask
 
    ipass=ipass+1
    loc  =j-1
    n1=ii1*loc+i
    bunpk=ig(n1)

    RETURN
END FUNCTION bunpk
