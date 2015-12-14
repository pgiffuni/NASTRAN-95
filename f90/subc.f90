SUBROUTINE subc
     
 COMPLEX :: gusamp,sbkde1,sbkde2,  &
     f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,am5tt,am6,  &
     sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,am5t,  &
     ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1, c2p,c2n,  &
     c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,blam2,  &
     ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t, c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,  &
     ft3tst,c1a,c2a,c3a,cexp1,cexp2,cexp3,cexp1a, cexp2a,cexp3a,const,  &
     bc,bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,  &
     clift,cmomt,c4a,cexp4,cexp5,cexp4a,cexp5a,  &
     pres1,pres2,pres3,pres4,qres4,fqa,fqb,t1,t2,t3,fq7
 DIMENSION       pres1(21),pres2(21),pres3(21),pres4(21),qres4(21),  &
     sbkde1(201),sbkde2(201), sumsv1(201),sumsv2(201),svkl1(201),svkl2(201),  &
     xlsv1(21),xlsv2(21),xlsv3(21),xlsv4(21)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 COMMON /blk3  / sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,a,b,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1,c2p,  &
     c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,  &
     blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,c1p,c1n,  &
     bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,bc2,bc3,  &
     bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,pres1,pres2,  &
     pres3,pres4,qres4,fqa,fqb,fq7
 COMMON /blk4  / i,r,y,a1,b1,c4,c5,gl,i6,i7,jl,nl,ri,rt,r5,sn,sp,  &
     xl,y1,amu,gam,idx,inx,nl2,rl1,rl2,rq1,rq2,xl1,  &
     alp1,alp2,gamn,gamp,iner,iout,redf,stag,step,  &
     amach,betnn,betnp,bkap1,xlsv1,xlsv2,xlsv3,xlsv4,  &
     alpamp,amoaxs,gusamp,disamp,pitaxs,pitcor
 
 am4tst = 0.0
 s1 = sps*del - sigma
 s2 = sps/(dstr**2)
 s3 = sns/dstr
 s4 = sps + sns
 t3 = CEXP(-ai*sigma)
 DO  i = 1,200
   r  = i
   gamp = 2.0*pi*r + s1
   gamn =-2.0*pi*r + s1
   c1p = (gamp/dstr) - scrk
   c2p = (gamp/dstr) + scrk
   alp = gamp*s2 - s3*CSQRT(c1p)*CSQRT(c2p)
   t1  = alp - del
   CALL akapm (alp,bkdel3)
   sbkde1(i+1) = bkdel3
   sum1 = CEXP(ai*(alp*sps-gamp))*(alp*sps-gamp)*bkdel3/((alp*dstr**2  &
       - gamp*sps)*t1)*(f6s*t1/(t1+gl) + f5s  &
       + b*ai/(bkdel1*bkap1)*(del-amu)/(alp-del))
   c1n  = (gamn/dstr) - scrk
   c2n  = (gamn/dstr) + scrk
   aln  = gamn*s2 - s3*CSQRT(c1n)*CSQRT(c2n)
   t2   = aln - del
   CALL akapm (aln,bkdel3)
   sbkde2(i+1) = bkdel3
   sum2 = CEXP(ai*(aln*sps-gamn))*(aln*sps-gamn)*bkdel3/((aln*dstr**2  &
       - gamn*sps)*t2)*(f6s*(t2)/(t2+gl) + f5s  &
       + b*ai/(bkdel1*bkap1)*(del-amu)/(t2))
   c1p  = CEXP(-ai*(t1)*sps)
   c2p  = CEXP(-ai*(t1)*sns)
   c1n  = CEXP(-ai*(t2)*sps)
   c2n  = CEXP(-ai*(t2)*sns)
   f4   = f4 + sum1*t3*ai/(t1)*(c1p-c2p) + sum2*t3*ai/(t2)*(c1n-c2n)
   am4  = am4 + sum1*t3*(ai*sps*c1p/(t1) - ai*sns*c2p/(t1) + 1.0/  &
       ((t1)**2)*(c1p-c2p)+ai*(2.0-sps)/(t1)*(c1p-c2p)) +  &
       sum2*t3*(ai*sps*c1n/(t2)-ai*sns*c2n/(t2) + 1.0/  &
       ((t2)**2)*(c1n-c2n) + ai*(2.0-sps)/(t2)*(c1n-c2n))
   i6   = i + 1
   temp = (sps-sns)/rl1
   c1a  =-ai*t1
   c2a  =-ai*t2
   c3a  = ai*del
   cexp1  = CEXP(c1a*sns)
   cexp2  = CEXP(c2a*sns)
   cexp3  = CEXP(c3a*sns)
   cexp1a = CEXP(c1a*temp)
   cexp2a = CEXP(c2a*temp)
   cexp3a = CEXP(c3a*temp)
   const  = fq7/(2.0*pi)
   temp2  = 2.0*pi*r/s4
   c4a    =-ai*s1
   cexp4  = CEXP(c4a*(2.0*sns/s4+0.5))
   cexp5  = CEXP(c4a*0.5)
   cexp4a = CEXP(c4a*temp/s4)
   cexp5a = CEXP(c4a*temp/(sps+sns))
   xl     = sns
   DO  jl = 1,nl
     pres4(jl) = pres4(jl) - t3*(sum1*cexp1 + sum2*cexp2  &
         + const*cexp3*(cexp4*SIN(temp2*(sns+xl))/r  &
         - cexp5*SIN(temp2*(sps+xl))/r))
     xl    = xl + temp
     cexp1 = cexp1*cexp1a
     cexp2 = cexp2*cexp2a
     cexp3 = cexp3*cexp3a
     cexp4 = cexp4*cexp4a
     cexp5 = cexp5*cexp5a
   END DO
   IF (cabs((am4-am4tst)/am4) < 0.0006) GO TO 75
   am4tst = am4
 END DO
 GO TO 9994
 75 CONTINUE
 temp  = (sps-sns)/rl1
 temp1 = 2.0*sns/s4 + 0.5
 temp2 = 0.5 - (sps+sns)/s4
 c1a   = ai*del
 c2a   =-ai*s1
 c3a   =-c2a
 cexp1 = CEXP(c1a*sns)
 cexp2 = CEXP(c2a*temp1)
 cexp3 = CEXP(c3a*temp2)
 cexp1a= CEXP(c1a*temp)
 cexp2a= CEXP(c2a*temp/s4)
 const = t3*fq7/2.0
 xl    = sns
 DO  jl = 1,nl
   pres4(jl) = pres4(jl) - const*cexp1*(cexp2*((sns+xl)/s4-0.5)  &
       - cexp3*((sps+xl)/s4-1.5))
   xl    = xl + temp
   cexp1 = cexp1*cexp1a
   cexp2 = cexp2*cexp2a
   cexp3 = cexp3*cexp2a
 END DO
 CALL subcc
 RETURN
 
 9994 WRITE  (ibbout,3015) ufm
 3015 FORMAT (a23,' - AMG MODULE -SUBROUTINE SUBC.  AM4 LOOP DID NOT ',  &
     'CONVERGE.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE subc
