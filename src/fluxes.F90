module fluxes
   integer, parameter :: WHICH_FLUX = 1
contains
subroutine inv_flux(ql,qr,dt,flux)

real(kind=8), intent( in) :: QL     (:) ! Input: conservative variables
real(kind=8), intent( in) :: QR     (:) ! Input: conservative variables
real(kind=8), intent( in) :: dt
!Output
real(kind=8), intent(out) :: flux  (:)        ! Output: numerical flux

if (WHICH_FLUX == 1) then
   call ausm(ql,qr,flux)
else if (WHICH_FLUX == 2) then
   call hll(ql,qr,flux)
else if (WHICH_FLUX == 3) then
   call central(ql,qr,flux)
else
   call lax_friedrich_i(ql,qr,dt,flux)
end if

end subroutine inv_flux

subroutine central(QL, QR, flux)
implicit none

real(kind=8), parameter :: gamma    = 1.4d0
real(kind=8), parameter :: one      = 1.0d0
real(kind=8), parameter :: half     = 0.5d0
!Input
real(kind=8), intent( in) :: QL     (:) ! Input: conservative variables
real(kind=8), intent( in) :: QR     (:) ! Input: conservative variables
!Output
real(kind=8), intent(out) :: flux  (:)        ! Output: numerical flux


!Local variables
!real(kind=8) :: nx, ny, nz                   ! Normal vector
real(kind=8) :: uL, uR                       ! Velocity components.
real(kind=8) :: rhoL, rhoR, pL, pR           ! Primitive variables.
real(kind=8) :: qnL, qnR                     ! Normal velocities
real(kind=8) :: aL, aR, HL, HR               ! Speed of sound, Total enthalpy
real(kind=8) :: fL(3), fR(3)                 ! Fluxes ad dissipation term

! Face normal vector (unit vector)
!  Left state
rhoL  = QL(1)
uL    = QL(2) / rhoL
qnL   = uL
pL    = (gamma-one)*( QL(3) - half*rhoL*(uL*uL) )
aL    = sqrt(gamma*pL/rhoL)
HL    = (QL(3) + pL ) / rhoL
fL(1) = rhoL*qnL
fL(2) = rhoL*qnL * uL + pL
fL(3) = rhoL*qnL * HL

!  Right state
rhoR  = QR(1)
uR    = QR(2) / rhoR
qnR   = uR
pR    = (gamma-one)*( QR(3) - half*rhoR*(uR*uR) )
aR    = sqrt(gamma*pR/rhoR)
HR    = (QR(3) + pR ) / rhoR
fR(1) = rhoR*qnR
fR(2) = rhoR*qnR * uR + pR
fR(3) = rhoR*qnR * HR
flux(:)  = HALF * ( fL + fR   )
end subroutine central


subroutine lax_friedrich_i(QL, QR,dt, flux)
! QL = rho * [1, u ,v , E]
implicit none

real(kind=8), parameter :: gamma    = 1.4d0
real(kind=8), parameter :: one      = 1.0d0
real(kind=8), parameter :: half     = 0.5d0
!Input
real(kind=8), intent( in) :: QL     (:) ! Input: conservative variables
real(kind=8), intent( in) :: QR     (:) ! Input: conservative variables
real(kind=8), intent( in) :: dt
!Output
real(kind=8), intent(out) :: flux  (:)        ! Output: numerical flux


!Local variables
!real(kind=8) :: nx, ny, nz                   ! Normal vector
real(kind=8) :: uL, uR                       ! Velocity components.
real(kind=8) :: rhoL, rhoR, pL, pR           ! Primitive variables.
real(kind=8) :: qnL, qnR                     ! Normal velocities
real(kind=8) :: aL, aR, HL, HR               ! Speed of sound, Total enthalpy
real(kind=8) :: fL(3), fR(3)                 ! Fluxes ad dissipation term

! Face normal vector (unit vector)
!  Left state
rhoL  = QL(1)
uL    = QL(2) / rhoL
qnL   = uL
pL    = (gamma-one)*( QL(3) - half*rhoL*(uL*uL) )
aL    = sqrt(gamma*pL/rhoL)
HL    = (QL(3) + pL ) / rhoL
fL(1) = rhoL*qnL
fL(2) = rhoL*qnL * uL + pL
fL(3) = rhoL*qnL * HL

!  Right state
rhoR  = QR(1)
uR    = QR(2) / rhoR
qnR   = uR
pR    = (gamma-one)*( QR(3) - half*rhoR*(uR*uR) )
aR    = sqrt(gamma*pR/rhoR)
HR    = (QR(3) + pR ) / rhoR
fR(1) = rhoR*qnR
fR(2) = rhoR*qnR * uR + pR
fR(3) = rhoR*qnR * HR
flux(:)  = HALF * ( fL + fR   &
                  - (QR-QL)/dt)
         
!flux(:)  = HALF * (fL + fR) &
!         + QUARTER * (aL + aR) &
!                   * (QL(:) + QL(:))
end subroutine lax_friedrich_i

!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Roe(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: Roe(3)       ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two, four
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR, HL, HR               ! Speeds of sound.
 real :: RT,rho,v,H,a                 ! Roe-averages
 real :: drho,du,dP,dV(3)
 real :: ws(3),Da, R(3,3)
 integer :: j, k

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0
      four = 4.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*v*v) )

!Differences in primitive variables.
   drho = rhoR - rhoL
     du =   vR - vL
     dP =   pR - pL

!Wave strength (Characteristic Variables).
   dV(1) =  half*(dP-rho*a*du)/(a*a)
   dV(2) = -( dP/(a*a) - drho )
   dV(3) =  half*(dP+rho*a*du)/(a*a)

!Absolute values of the wave speeds (Eigenvalues)
   ws(1) = abs(v-a)
   ws(2) = abs(v  )
   ws(3) = abs(v+a)

!Modified wave speeds for nonlinear fields (to remove expansion shocks).
!There are various ways to implement an entropy fix. This is just one
!example.
   Da = max(zero, four*((vR-aR)-(vL-aL)) )
   if (ws(1) < half*Da) ws(1) = ws(1)*ws(1)/Da + quarter*Da
   Da = max(zero, four*((vR+aR)-(vL+aL)) )
   if (ws(3) < half*Da) ws(3) = ws(3)*ws(3)/Da + quarter*Da

!Right eigenvectors
   R(1,1) = one
   R(2,1) = v - a
   R(3,1) = H - v*a

   R(1,2) = one
   R(2,2) = v
   R(3,2) = half*v*v

   R(1,3) = one
   R(2,3) = v + a
   R(3,3) = H + v*a

!Compute the average flux.
   !Roe = half*( physical_flux(uL) + physical_flux(uR) )

!Add the matrix dissipation term to complete the Roe flux.
  do j = 1, 3
   do k = 1, 3
    Roe(j) = Roe(j) - half*ws(k)*dV(k)*R(j,k) 
   end do
  end do

 end function Roe
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- AUSM Flux Function ---
!*
!* M.-S. Liou and C. J. Steffen, A New Flux Splitting Scheme, Journal of 
!* Computational Physics, 107, pp. 23-39, 1993.
!*
!* NB: This may require a low CFL number to get started.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 subroutine AUSM(uL,uR,flux)
 implicit none
 real(kind=8),intent(in) :: uL(:), uR(:)  !  Input (conservative variables rho*[1, v, E])
 real(kind=8),intent(out) :: flux(:)       ! Output (numerical flux across L and R states)
!Local constants
 real(kind=8) :: gamma                        ! Ratio of specific heat.
 real(kind=8) :: zero, quarter, half, one, two
!Local variables
 real(kind=8) :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real(kind=8) :: HL, HR                       ! Specific enthaply (per unit mass).
 real(kind=8) :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real(kind=8) :: Pp, Pm, Mp, Mm
 real(kind=8) :: Fp(3), Fm(3)                 ! F_plus and F_minus

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR
      HR = ( uR(3) + pR ) / rhoR

!Positive M and p in the LEFT cell.
 if (ML <= -one) then
   Mp = zero
   Pp = zero
 elseif (ML < one) then
   Mp = quarter*(ML+one)*(ML+one)
   Pp = quarter*PL*(one+ML)*(one+ML)*(two-ML) ! or use Pp = half*(one+ML)*pL
 else
   Mp = ML
   Pp = PL
 endif

!Negative M and p in the RIGHT cell.
 if   (MR <= -one) then
   Mm = MR
   Pm = PR;
 elseif (MR < one) then
   Mm = -quarter*(MR-one)*(MR-one)
   Pm =  quarter*pR*(one-MR)*(one-MR)*(two+MR) ! or use Pm = half*(one-MR)*pR
 else
   Mm = zero
   Pm = zero
 endif

!Positive Part of Flux evaluated in the left cell.
 Fp(1) = max(zero,Mp+Mm)*aL * rhoL
 Fp(2) = max(zero,Mp+Mm)*aL * rhoL*vL  + Pp
 Fp(3) = max(zero,Mp+Mm)*aL * rhoL*HL

!Negative Part of Flux evaluated in the right cell.
 Fm(1) = min(zero,Mp+Mm)*aR * rhoR
 Fm(2) = min(zero,Mp+Mm)*aR * rhoR*vR  + Pm
 Fm(3) = min(zero,Mp+Mm)*aR * rhoR*HR

!Compute the flux: Fp(uL)+Fm(uR).
    flux = Fp + Fm

 end subroutine AUSM
!-----------------------------------------------------------------------------
!*****************************************************************************
!* --- HLL Flux Function ---
!*
!* A. Harten, P. D. Lax, and B. van Leer,On Upstream Differencing and 
!* Godunov-Type Schemes for Hyperbolic Conservation Laws, SIAM Review,
!* 25(1), pp. 35-61, 1983.
!*
!* With wave speeds evaluated by Einfeldt's method:
!* B. Einfeldt, On Godunov-Type Methods for Gas Dynamics, SIAM Journal of
!* Numerical Analysis, 25(2), pp. 294-318, 1988.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 subroutine HLL(uL,uR,flux)
 implicit none
 real(kind=8),intent(in):: uL(:), uR(:)     !  Input (conservative variables rho*[1, v, E])
 real(kind=8),intent(out):: flux(:)           ! Output (numerical flux across L and R states)
!Local constants
 real(kind=8) :: gamma                        ! Ratio of specific heat.
 real(kind=8) :: zero, quarter, half, one
!Local variables
 real(kind=8) :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real(kind=8) :: aL, aR, HL, HR               ! Speeds of sound.
 real(kind=8) :: RT,rho,u,H,a
 real(kind=8) :: SRp,SLm, uma,upa
 real(kind=8) :: fl(3) , fr(3)

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!Evaluate the two wave speeds: Einfeldt.
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*u*u) )

   uma = u - a
   upa = u + a
   SLm = min(vL-aL, uma, zero)
   SRp = max(vR+aR, upa, zero)

!Compute the HLL flux.
   call physical_flux(uL,fl)
   call physical_flux(uR,fr)
   flux = (SRp*fl-SLm*fr + SLm*SRp*(uR-uL))/(SRp-SLm)

 end subroutine HLL
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- 1D physical Euler flux --- 
!*
!* This is called in LaxFriedrichs, Richtmyer, MacCormack, StegerWarming,
!* Godunov, Osher, Roe, HLL, HLLL.
!*
!* The vector f(U) in u_t + f(u)_x = 0, as a function of 
!* the conservative variables, u = [density, density*velocity, total energy].
!
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 subroutine physical_flux(u,flux)
 real(kind=8),intent(in) :: u(:)             !  Input (conservative variables [rho, rho*v, rho*E])
 real(kind=8),intent(out) :: flux(:) ! Output (physical flux of the Euler equations)
!Local variables
 real(kind=8) :: density, velocity, pressure, enthalpy, gamma

!Define and compute some quantities.
     gamma = 1.4
   density = u(1)
  velocity = u(2)/u(1)
  pressure = (gamma-1.0)*( u(3) - 0.5*density*velocity*velocity )
  enthalpy = u(3) + pressure

!Evaluate the physical flux (mass, momentum, and energy fluxes).
  flux(1) =           density * velocity
  flux(2) = (density*velocity)* velocity + pressure
  flux(3) =          enthalpy * velocity

 end subroutine physical_flux
!*****************************************************************************

end module fluxes
