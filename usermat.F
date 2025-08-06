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
C          call usermat_harm (
C      &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
C      &                   ldstep,isubst,keycut,
C      &                   nDirect,nShear,ncomp,nProp,
C      &                   Time,dTime,Temp,stress,dsdePl, tsstif,
C      &                   Strain,prop,coords)
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
            E = E0 * (ONE - TWO*kuzAlphaE/pi*ATAN(strTrace/kuzBeta))
            nu = nu0 * (ONE - TWO*kuzAlphaNu/pi*ATAN(strTrace/kuzBeta))
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

c Copy elastic moduli stiffness to material Jacobian matrix (not convinced this is necessary)
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

c Update plastic strains (not convinced this is necessary)
c      call vapb1(epsPl(1),wk2(1),ncomp)

c Determine elastic strain energy
      sedEl = ZERO
      do i=1,ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i))
      end do
      sedEl = sedEl * HALF
	  
      return




c *** Use Latorre and Montans Uncoupled Strain Formulation (DO NOT USE! BROKEN RIGHT NOW AND I'M NOT PAID FOR THIS ONE!)
      else if (nProp .ge. 5 .and. prop(5) .eq. 2) then
c *** Get properties as before
      Et        = prop(1)
      Ec        = prop(2)
      nut       = prop(3)
      nuc       = prop(4)

c *** Get the previous elastic strain energy
      sedElPrev = ustatev(6)
	  
c *** Calculate addends mut, muc, Kt, Kc (Latorre and Montans, 2020, p.7, eqn. 69-72)
      mut = Et/(ONE+nut)-HALF*Ec/(ONE+nuc)
      muc = Ec/(ONE+nuc)-HALF*Et/(ONE+nut)
	  
      Kt = THIRD*Et/(ONE-TWO*nut)
      Kc = THIRD*Ec/(ONE-TWO*nuc)
	  
c *** Get deviatoric and volumetric strain tensors
c *** note: 11, 22, 33, 12, 23, 13 is order for Strain tensor in APDL
      eMean   = THIRD*(Strain(1) + Strain(2) + Strain(3))

      call vmove(Strain(1), eDev(1), ncomp)
      eDev(1) = eDev(1) - eMean
      eDev(2) = eDev(2) - eMean
      eDev(3) = eDev(3) - eMean
	  
c *** Get deviatoric and volumetric strain differential tensors
c *** note: 11, 22, 33, 12, 23, 13 is order for Strain tensor in APDL
      dEMean   = THIRD*(dStrain(1) + dStrain(2) + dStrain(3))

      call vmove(dStrain(1), dEDev(1), ncomp)
      dEDev(1) = dEDev(1) - dEMean
      dEDev(2) = dEDev(2) - dEMean
      dEDev(3) = dEDev(3) - dEMean
	  
c *** Determine deviatoric strain invariants J2, J3 (J1 = 0, always) (if there is a mistake it is here, probably J3!)
      J2 = -HALF*THIRD*((eDev(1)-eDev(2))*(eDev(1)-eDev(2))
     &   + (eDev(2)-eDev(3))*(eDev(2)-eDev(3))
     &   + (eDev(3)-eDev(1))*(eDev(3)-eDev(1))
     &   + TWO*THREE*(eDev(4)*eDev(4) + eDev(5)*eDev(5)
     &   + eDev(6)*eDev(6)))
	 
      J3 = ((eDev(3) - eMean)*eDev(2) + eMean*eMean - eMean*eDev(3)
     &    - eDev(5)*eDev(5))*eDev(1) + (eMean*eMean - eMean*eDev(3)
     &	  - eDev(6)*eDev(6))*eDev(2) + (eMean*eMean 
     &    - eDev(4)*eDev(4))*eDev(3) + eDev(4)*eDev(4)*eMean
     &	  + TWO*eDev(5)*eDev(4)*eDev(6) - eMean*(eMean*eMean 
     &    - eDev(5)*eDev(5) - eDev(6)*eDev(6))
c *** Directly solve the eigenvalue problem assuming that all three roots must necessarily be real (using depressed cubic trigonometric formulation, J1=0)
      ed3 = (TWO*sqrt(-THREE*J2)*sin((5.d0*Pi)/6.d0 
     &    + ACOS((THREE*J3*sqrt(-THREE/J2))/(TWO*J2))/THREE))/THREE
      ed2 = (TWO*sqrt(-THREE*J2)*sin((3.d0*Pi)/2.d0
     &    + ACOS((THREE*J3*sqrt(-THREE/J2))/(TWO*J2))/THREE))/THREE
      ed1 = (TWO*sqrt(-THREE*J2)*sin((13.d0*Pi)/6.d0
     &    + ACOS((THREE*J3*sqrt(-THREE/J2))/(TWO*J2))/THREE))/THREE
c *** Determine deviatoric and volumetric strain switches mu1, mu2, mu3, K (Latorre and Montans, 2020, p.3, eqn. 12-13)
c
      if (ed1 .GE. 0) then
         mu1 = mut
      else
         mu1 = muc
      end if
	  
      if (ed2 .GE. 0) then
         mu2 = mut
      else
         mu2 = muc
      end if
	  
      if (ed3 .GE. 0) then
         mu3 = mut
      else
         mu3 = muc
      end if

      if (eMean .GE. 0) then
         K = Kt
      else
         K = Kc
      end if
      
      ustatev(1) = ed1
      ustatev(2) = ed2
      ustatev(3) = ed3
      ustatev(4) = eMean
	  
c *** Determine strain energy for uncoupled strain-driven bilinear elastic isotropic material (Latorre and Montans, 2020, p.3, eqn. 14)
c 
      sedEl = mu1*ed1*ed1 + mu2*ed2*ed2 + mu3*ed3*ed3 
     &      + HALF*K*eMean*eMean
	 
c *** get initial stress
      call get_ElmInfo('NCOMP', ncompgt)
      call vzero(sigi(1),ncompgt)
      call get_ElmData ('ISIG', elemId,kDomIntPt, ncompgt, sigi)
	  
c *** Determine stress from derivative of strain energy function (eqn. 31) (THIS IS PROBABLY ALSO MISTAKE!)
      stress(1) = (sedEl - sedElPrev)/(dStrain(1)) + sigi(1)
      stress(2) = (sedEl - sedElPrev)/(dStrain(2)) + sigi(2)
      stress(3) = (sedEl - sedElPrev)/(dStrain(3)) + sigi(3)
      stress(4) = (sedEl - sedElPrev)/(dStrain(4)) + sigi(4)
      stress(5) = (sedEl - sedElPrev)/(dStrain(5)) + sigi(5)
      stress(6) = (sedEl - sedElPrev)/(dStrain(6)) + sigi(6)
      
c *** Set the previous elastic strain energy
      ustatev(6) = sedEl

      return
      end if
      end

*deck,usermat_harm    USERDISTRIB  parallel                    jmgerken
      subroutine usermat_harm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nProp,
     &                   freq,dfreq,Temp,stress,jacobi,tsstif,
     &                   strain,prop,coords)
c*************************************************************************
c
c           User defined harmonic constitutive model
c
c*************************************************************************
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
c      nProp     (int,sc,l)               Number of material constants
c      Temp      (dp,sc,in)               temperature
c      freq      (dp,sc,in)               frequency
c      dfreq      (dp,sc,in)              frequency increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp,2),io)       real and imaginary stress
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c
c     output arguments
c     ================
c      jacobi   (dp,ar(ncomp,ncomp,2),io)   real and imaginary material jacobian matrix
c
c*************************************************************************
#include "impcom.inc"
c
c***************** Arguments *************************************
c
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nProp
      DOUBLE PRECISION freq,dfreq, Temp
      DOUBLE PRECISION 
     &                 stress  (ncomp,2),
     &                 jacobi  (ncomp,ncomp,2),
     &                 strain  (ncomp,2),
     &                 prop    (nProp),
     &                 coords  (3),
     &                 tsstif  (2)
c
c***************** Local *************************************
c

      INTEGER i, j, iprony ! loop counter
      INTEGER nprnsh,nprnvl ! number of prony terms
      INTEGER nvar  ! index
      DOUBLE PRECISION c1, c2, c3 ! temporary
      DOUBLE PRECISION alpha, tau  ! prony parameters
      DOUBLE PRECISION Omega  ! frequency
      DOUBLE PRECISION smG,lmG ! storage and loss shear modulus factors
      DOUBLE PRECISION smK,lmK ! storage and loss bulk modulus factors
      DOUBLE PRECISION young, posn, G, bulk  ! elastic constants

      DOUBLE PRECISION pi
      PARAMETER (pi=3.14159265358979324D0)

!-----------------------------------------------------------------------

      keycut=0
        
! Prony series storage and loss modulus factors 
      
      Omega = freq*2.0D0*pi ! convert frequency to rad/sec
      nprnsh = nint(prop(3)) ! number of shear prony terms
      nprnvl = nint(prop(4)) ! number of bulk prony terms
      
      ! shear
      smG = 1.D0  ! storage modulus factor
      lmG = 0.D0  ! loss modulus factor
      nvar = 4 ! index to shear prony parameters
      
      do iprony = 1, nprnsh
        nvar = nvar + 2
        alpha = prop(nvar-1)
        tau   = prop(nvar)
        c1 = Omega * tau
        c2 = c1**2
        smG = smG + (c2  / ( 1.D0 + c2 ) - 1.D0 ) * alpha
        lmG = lmG + c1  / ( 1.D0 + c2 ) * alpha
      end do
      lmG = lmG / (freq*2.0D0*pi)  ! scale loss modulus by 1/Omega
      
      ! bulk
      smK = 1.D0  ! storage modulus factor
      lmK = 0.D0  ! loss modulus factor
      
      do iprony = 1, nprnvl
        nvar = nvar + 2
        alpha = prop(nvar-1)
        tau   = prop(nvar)
        c1 = Omega * tau
        c2 = c1**2
        smK = smK + (c2  / ( 1.D0 + c2 ) - 1.D0 ) * alpha
        lmK = lmK + c1  / ( 1.D0 + c2 ) * alpha
      end do
      lmK = lmK / (freq*2.0D0*pi)  ! scale loss modulus by 1/Omega
      
! stiffness      
      
      jacobi(:,:,:) = 0.D0
      young = prop(1)
      posn = prop(2)
      G = young / (2.D0*(1.D0+posn))
      bulk = young/(3.D0*(1.D0-2.D0*posn))
      
      ! transverse shear stiffness = storage shear modulus
      tsstif(1) = G*smG
      tsstif(2) = G*smG

      IF(ncomp .GE. 4) THEN
      ! ***    3d, plane strain and axisymmetric
        
        ! storage modulus
        do i=1,3
          jacobi(i,i,1) = 2.D0*G*smG
        end do 
        do i=4,ncomp
          jacobi(i,i,1) = G*smG
        end do 
        
        do i=1,3
          do j=1,3
            jacobi(j,i,1) = jacobi(j,i,1) - 2.D0/3.D0*G*smG + bulk*smK
          end do
        end do 
        
        ! loss modulus
        do i=1,3
          jacobi(i,i,2) = 2.D0*G*lmG
        end do 
        do i=4,ncomp
          jacobi(i,i,2) = G*lmG
        end do 
        
        do i=1,3
          do j=1,3
            jacobi(j,i,2) = jacobi(j,i,2) - 2.D0/3.D0*G*lmG + bulk*lmK
          end do
        end do 

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
      ! ***    plane stress
        
        ! storage modulus
        c1 = bulk*smK
        c2 = G*smG
        posn = 0.5D0*(3.D0*c1-2.D0*c2)/(3.D0*c1+c2)
        young = 9.D0*c1*c2/(3.D0*c1+c2)
        c1 = 1.D0 - posn * posn
        c2 = young / c1
        c3 = posn * c2
        
        jacobi(1,1,1) = c2
        jacobi(1,2,1) = c3
        jacobi(2,1,1) = c3
        jacobi(2,2,1) = c2
        jacobi(3,3,1) = G*smG
        
        ! loss modulus
        c1 = bulk*lmK
        c2 = G*lmG
        posn = 0.5D0*(3.D0*c1-2.D0*c2)/(3.D0*c1+c2)
        young = 9.D0*c1*c2/(3.D0*c1+c2)
        c1 = 1.D0 - posn * posn
        c2 = young / c1
        c3 = posn * c2
        
        jacobi(1,1,2) = c2
        jacobi(1,2,2) = c3
        jacobi(2,1,2) = c3
        jacobi(2,2,2) = c2
        jacobi(3,3,2) = G*lmG

      ELSE IF(ncomp .EQ. 3) THEN
      ! ***    3d beam
        
        ! storage modulus
        c1 = bulk*smK
        c2 = G*smG
        jacobi(1,1,1) = 9*c1*c2/(3*c1+c2)
        jacobi(2,2,1) = G*smG
        jacobi(3,3,1) = G*smG
        
        ! loss modulus
        c1 = bulk*lmK
        c2 = G*lmG
        jacobi(1,1,2) = 9*c1*c2/(3*c1+c2)
        jacobi(2,2,2) = G*lmG
        jacobi(3,3,2) = G*lmG

      ELSE IF(ncomp .EQ. 1) THEN
      ! ***    1d beam

        ! storage modulus
        jacobi(1,1,1) = 2.D0/3.D0*young*smG + 1.D0/3.D0*young*smK
        
        ! loss modulus
        jacobi(1,1,2) = 2.D0/3.D0*young*lmG + 1.D0/3.D0*young*lmK
         
      END IF

! stress

      ! real
      stress(:,1) = 
     &    matmul (jacobi(1:ncomp,1:ncomp,1),strain(1:ncomp,1)) -
     &    matmul (jacobi(1:ncomp,1:ncomp,2),strain(1:ncomp,2))
     &    *(freq*2.0D0*pi)
      
      ! imaginary
      stress(:,2) = 
     &    matmul (jacobi(1:ncomp,1:ncomp,1),strain(1:ncomp,2)) +
     &    matmul (jacobi(1:ncomp,1:ncomp,2),strain(1:ncomp,1))
     &    *(freq*2.0D0*pi)
     
      return
      end
