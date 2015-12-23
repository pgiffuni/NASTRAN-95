SUBROUTINE bpack (ig,i,j,l)
     
 IMPLICIT INTEGER (a-z) 
 INTEGER*2, INTENT(OUT)                   :: ig(1)
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(IN)                      :: j
 INTEGER, INTENT(IN)                      :: l
 !IMPLICIT INTEGER 
 
!DC   NEXT 2 LINES FOR CDC AND UNIVAC ONLY
!     EXTERNAL         ORF,      LSHIFT
!     INTEGER          IG(1)
 
!     NEXT LINE FOR IBM, VAX, AND MACHINES THAT HAVE INTEGER*2
 
 
 COMMON /bandb /  nbit,     dum3b(3), ipass,    nw,       dum1b, nbpw
 COMMON /bands /  dum4s(4), ii1,      dum5s(5), mask
 
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     PACK INTERNAL GRID NO. INTO IG TABLE.  SEE BUNPK FOR UNPACKING
!     TABLE IG IS PACKED COLUMN-WISE.
!     USE APPROP. PORTION OF THIS ROUTINE FOR DIFFERENT TYPE OF MACHINE.
!     IPASS=COUNTER ON NUMBER OF CALLS TO PACK/UNPACK
 
!     NOTE - THIS ROUTINE DOES NOT CHECK NOR ZERO OUT THE PACKING SLOT
!            BEFORE PACKING.
!            L IS ASSUMED TO BE A POSITIVE INTEGER, NBIT BITS OR LESS
 
 ipass=ipass+1
 loc  =j-1
 
!     ********************************************
!     UNIVAC AND CDC MACHINES
!     (IG SHOULD BE IN INTEGER*4 HERE)
!     ********************************************
 
!     N1 =II1*(LOC/NW)+I
!     N2 =MOD(LOC,NW)*NBIT+NBIT
!     LOC=ORF(IG(N1),LSHIFT(L,NBPW-N2))
!     IG(N1)=LOC
 
!     RETURN
 
!     ********************************************
!     IBM AND VAX MACHINES
!     (IG IS SET TO INTEGER*2 IN BPACK AND BUNPK, ELSEWHERE INTEGER*4)
!     INTEGER*2     IG(1)
!     ********************************************
 
 n1=ii1*loc+i
 ig(n1)=l
 
 RETURN
END SUBROUTINE bpack
