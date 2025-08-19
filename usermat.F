*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain,epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, pVolDer, hrmflg, var3, var4,
     &                   var5, var6, var7)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and 3D/1D beam.
c
c           A 3D material constitutive model can be used for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be used.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model, which is the same as TB, BISO,
c       for different stress states. 
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c       This routine calls four routines,
c       usermat3d.F, usermatps.F usermatbm.F and usermat1d.F, w.r.t.
c       the corresponding stress states.
c       Each routine can be also a usermat routine for the specific 
c       element.
c
c*************************************************************************
c Copyright ANSYS.  All Rights Reserved.
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c      hrmflg   (dp,sc,io)                flag to indicate harmonic analysis 
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      pVolDer  (dp,ar(3),o)              derivatives of volumetric potential wrt to J
c                                         pVolDer(1) = dU/dJ
c                                         pVolDer(2) = d^2U/dJ^2
c                                         pVolDer(3) = d^3U/dJ^3
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,ncompgt,nStatev,nProp,
     &                 i, j
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp), sigi(ncomp), ustatev(nStatev),
     &                 wStress (ncomp), pVolDer (3),
     &                 Strain  (ncomp), dStrain (ncomp  ), 
     &                 wStrain (ncomp), strTrace, lambda, G,
     &                 epsPl   (ncomp  ), prop    (nProp  ),
     &                 dsdePl(ncomp, ncomp),	 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2),
     &                 stiffness (ncomp, ncomp),
     &                 E, Et, Ec, nu, nut, nuc, Kt, Kc, mut, muc,
     &                 kuzBeta, E0, kuzAlphaE, nu0, kuzAlphaNu,
     &                 ZERO, HALF, THIRD, ONE, TWO, THREE, PI,
     &                 mu1, mu2, mu3, K, sedElPrev, 
     &                 ed1, ed2, ed3, J2, J3, eDev(ncomp), 
     &                 dEDev(ncomp), eMean, dEMean
      DOUBLE PRECISION hrmflg

      EXTERNAL         usermat_harm

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      data             var1/0.0d0/
      data             var2/0.0d0/
c
c*************************************************************************
c

c ***    harmonic analysis

      IF(nint(hrmflg)/=0)THEN
         return
       ENDIF
	   
c *** ustatev(1) is used to determine which dimension was used
c *** ustatev(2) is used to determine the tensile/compressive/zero behavior of the element 

      ZERO = 0.d0
      HALF = 0.5d0
      THIRD = 1.d0/3.d0
      ONE = 1.d0
      TWO = 2.d0
      THREE = 3.d0
      PI = 4.D0*DATAN(1.D0)
c May need to be atan only

      keycut   = 0
      cutFactor = 0.d0

c The input property vector is expected to contain the following values in the following positions:
c prop(1) = Young's modulus in tension
c prop(2) = Young's modulus in compression
c prop(3) = Poisson's ratio in tension
c prop(4) = Poisson's ratio in compression
c prop(5) = Bimodulus method; 0 (Liu and Peddieson Step Function) or 1 (Kuznetsov Arctangent) (I am desperately trying to implement the Latorre/Montans hyperelastic formulation but am not doing well)
c prop(6) = Kuznetsov Beta (Smoothness constant)

c The output state vector can be expected to contain the following values in the following positions:
c ustatev(1) = Dimensionality; 1 (1D), 2 (Plane Stress), 3 (3D)
c ustatev(2) = Tensile/Compressive/Zero; 1 (Tensile), 0 (Zero), -1 (Compressive)
c ustatev(3) = Trace of stress tensor
c ustatev(4) = Young's modulus at current integration point
c ustatev(5) = Poisson's ratio at current integration point
c ustatev(6) = 
	  
c     If the fifth property is set (equal to 0 (Liu and Peddieson Step Function) or 1 (Kuznetsov Arctangent)), use invariant method
      if (nProp .eq. 4 .or. prop(5) .eq. 0 .or. prop(5) .eq. 1) then

c     Get the two Young's moduli and the Poisson's ratio
      Et        = prop(1)
      Ec        = prop(2)
      nut       = prop(3)
      nuc       = prop(4)
	  
c     Get initial stress (fbc) on the element
      call get_ElmInfo('NCOMP', ncompgt)
      call vzero(sigi(1),ncompgt)
      call get_ElmData ('ISIG', elemId, kDomIntPt, ncompgt, sigi)
	  
c     Determine total strain at this timestep (=Strain + dStrain) and place it in vector wStrain
      call vmove(Strain(1), wStrain(1), ncomp)
      call vapb1(wStrain(1), dStrain(1), ncomp)

C c Determine trace of strain tensor...
C       if(ncomp .ge. 4) then
C c ...for 3D, plane strain, axisymmetric elements
C          strTrace = wStrain(1) + wStrain(2) + wStrain(3)
C       else if(nDirect .eq. 2 .and. ncomp .eq. 3) then
C c ...for plane stress elements
C          strTrace = wStrain(1) + wStrain(2)
C       else
C c ...for beam and truss elements
C          strTrace = wStrain(1)
C       end if
	  
c ----- Liu and Peddieson Step Function Method -----
      if (prop(5) .eq. 0) then
c           Determine trace of stress tensor...
            if(ncomp .ge. 4) then
c                 ...for 3D, plane strain, axisymmetric elements
                  strTrace = stress(1) + stress(2) + stress(3)
            else if(nDirect .eq. 2 .and. ncomp .eq. 3) then
c                 ...for plane stress elements
                  strTrace = stress(1) + stress(2)
            else
c                 ...for beam and truss elements
                  strTrace = stress(1)
            end if  
 
c           Toggle tensile/compressive behavior based on trace of strain
            if(strTrace .gt. 0.d0) then
                  E = Et
                  nu = nut
                  ustatev(2) = 1.d0
            else if(strTrace .eq. 0.d0) then
                  E = Et*Ec/(Et+Ec)
                  nu = (nut*Ec + nuc*Et)/(Et+Ec)
                  ustatev(2) = 0.d0
            else
                  E = Ec
                  nu = nuc
                  ustatev(2) = -1.d0
            end if
      else
c --- Kuznetsov Arctangent Method (1D ONLY!) ---
            strTrace = wStrain(1)

c           Set kuzBeta to the sixth property given (or 0.01 by default)
            if (nProp .ge. 6) then
                  kuzBeta = prop(6)
            else 
                  kuzBeta = 0.01d0
            end if

c           Calculate required E0/aE
            E0 = (Ec + Et)/TWO
            kuzAlphaE = (Ec - Et)/(Ec + Et)

c           Extension of idea: calculate required nu0/aNu
            nu0 = (nuc + nut)/TWO
            kuzAlphaNu = (nuc - nut)/(nuc + nut)

c           Calculate smoothed elastic modulus/poisson's ratio
            E = E0 * (ONE - TWO*kuzAlphaE/pi*DATAN(strTrace/kuzBeta))
            nu = nu0 * (ONE - TWO*kuzAlphaNu/pi*DATAN(strTrace/kuzBeta))
      end if
	  
c     Set the third state variable to the trace stress (for debugging purposes)
      ustatev(3) = strTrace	

c     Set the fourth and fifth state variables to E, nu (for debugging purposes)
      ustatev(4) = E  
      ustatev(5) = nu

c Establish elastic stiffness matrix
      if(ncomp .ge. 4) then
c ...for 3D, plane strain, axisymmetric elements
         lambda = E*nu/((ONE+nu)*(ONE-TWO*nu))
         G = E/(TWO + TWO*nu)
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
         stiffness(1,1) = lambda + TWO*G
         stiffness(2,2) = lambda + TWO*G
         stiffness(3,3) = lambda + TWO*G
         stiffness(4,4) = G
         stiffness(5,5) = G
         stiffness(6,6) = G
         stiffness(1,2) = lambda
         stiffness(1,3) = lambda
         stiffness(2,3) = lambda
         stiffness(2,1) = lambda
         stiffness(3,1) = lambda
         stiffness(3,2) = lambda
         ustatev(1) = 3.d0
      else if(nDirect .eq. 2 .and. ncomp .eq. 3) then
c ...for plane stress elements
c          11, 22, 12                for plane stress
         stiffness(1,1) = E/(ONE-nu*nu)
         stiffness(2,2) = E/(ONE-nu*nu)
         stiffness(3,3) = E/(TWO + TWO*nu)
         stiffness(1,2) = E*nu/(ONE-nu*nu)
         stiffness(2,1) = E*nu/(ONE-nu*nu)
         epsZZ = 0.d0
         ustatev(1) = 2.d0
      else
c ...for beam and truss elements
         stiffness(1,1) = E
         stiffness(2,2) = E/(TWO*(ONE+nu))
         stiffness(3,3) = E/(TWO*(ONE+nu))
         ustatev(1) = 1.d0
      end if

c Copy elastic moduli stiffness to material Jacobian matrix (not convinced this is necessary for fully elastic analysis)
      call vmove(stiffness(1,1), dsdePl(1,1), ncomp * ncomp)
	  
c Use constitutive law sigma = [S]epsilon (+sigi) to determine stress workingStress
c Having a weird issue with stress, and I'm going to see if changing to this stress will work
      call vzero(wStress(1), ncomp)
      call vmove(sigi(1), wStress(1), ncomp)
      do i=1,ncomp
	     do j=1,ncomp
c may need to be wStrain(i) to work
            wStress(i) = wStress(i) + stiffness(i,j) * wStrain(j)
         end do
      end do
	  
c Update stress vector
      call vmove(wStress(1), stress(1), ncomp)

c Force overwrite state variables (because I cannot get this to work for 1D beam elements????)
      call put_ElmData('SVAR', elemId, kDomIntPt, nStatev, ustatev)

c Update plastic strains (not convinced this is necessary)
c      call vapb1(epsPl(1),wk2(1),ncomp)

c Determine elastic strain energy
      sedEl = ZERO
      do i=1,ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i))
      end do
      sedEl = sedEl * HALF
	  
      return
      end if
      end
