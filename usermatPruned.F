c *deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c*************************************************************************
c  Modified from usermatbm.F by gal, July 1999; by Gabriel Tardy, July 2024.
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,  cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, ONE, TWO, SMALL, SQTWOTHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 3
     &                 )
c
c --- local variables
c
      EXTERNAL         vmove, vzero, vapb1, vamb1,get_ElmData
      DOUBLE PRECISION sigElp(mcomp), stiffness(mcomp,mcomp), 
     &                 workingStrain(3), wk2(3), diagS(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      INTEGER          i, j, k, ncompgt
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 cpleq, dpleq,   pleq,    twoG,
     &				   Et, Ec, nut, nuc, Gt, Gc, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep,
     &                 c1, c2, c3, c4, c5
      DOUBLE PRECISION pv(3)
      data pv/TWOTHIRD, TWO, TWO/
c*************************************************************************
c Get the two Young's moduli and the Poisson's ratio
      Et        = prop(1)
	  Ec        = prop(2)
      nut       = prop(3)
	  nuc       = prop(4)
	  
c Calculate the shear moduli in tension/compression
      Gt        = Et/(TWO*(ONE+nut))
	  Gc        = Ec/(TWO*(ONE+nuc))
	  
c Get existing stress on the element
      call get_ElmInfo('NCOMP', ncompgt)
      call vzero(sigi(1),ncompgt)
      call get_ElmData ('ISIG', elemId, kDomIntPt, ncompgt, sigi)

c Determine the elastic stiffness matrix
      call vzero(stiffness(1,1), ncomp * ncomp)
      if (stress(1) .GE. 0) then
         stiffness (1,1) = Et
         stiffness (2,2) = Gt
         stiffness (3,3) = Gt
      else
         stiffness (1,1) = Ec
         stiffness (2,2) = Gc
         stiffness (3,3) = Gc
      end if

c Equivalent to:
c {sxx} = [ Ei 0  0  ]{exx}
c {sxz} = [ 0  Gi 0  ]{exz}
c {sxy} = [ 0  0  Gi ]{exy}
c where i=tension if sxx>0; otherwise, compression
	  
c Determine total strain at this timestep (=Strain + dStrain) and place it in vector workingStrain
      call vmove(Strain(1), workingStrain(1), ncomp)
      call vapb1(workingStrain(1), dStrain(1), ncomp)

c Copy elastic moduli stiffness to material Jacobian matrix (not convinced this is necessary)
      call vmove(stiffness(1,1), dsdePl(1,1), ncomp * ncomp)
	  
c Use constitutive law sigma = [S]epsilon (+sigi?) to determine stress sigElp
      do i=1,ncomp
         sigElp(i) = stiffness(i,i) * workingStrain(i) + sigi(i)
      end do
	  
c Update stress vector
      call vmove(sigElp(1), stress(1), ncomp)

c Update plastic strains (not convinced this is necessary)
      call vapb1(epsPl(1),wk2(1),ncomp)

c Determine elastic strain energy
      sedEl = ZERO
      do i=1,ncomp
		  sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i))
      end do
      sedEl = sedEl * HALF

 990  CONTINUE
c
      return
      end