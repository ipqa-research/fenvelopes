subroutine read2PcubicNC(nc, nin, nout)
   !! Subroutine to read a GPECIN-like file and set up the system

   ! Critical constants must be given in K and bar
   ! b will be in L/mol and ac in bar*(L/mol)**2
   ! -> EoS parameters will be calculated from the critical constants
   !    to assure thermodynamic consistency

   use constants
   use system, only: setup, SRK_factory, PR76_factory, PR78_factory, &
                     z, bij, kij_mod => kij
   use system, only: z, nmodel => thermo_model, &
                     tc, pc, dceos => dc, om => w, &
                     ac, b, del1, k, kij, &
                     ntdep => tdep, ncomb => mixing_rule, bij, kinf, tstar, lij

   implicit none
   common /NAMES/ fluid

   integer, parameter :: nco=64

   integer, intent(in) :: nc !! number of components
   integer, intent(in) :: nin, nout !! IO units

   character*18 fluid(nco)

   real(pr) :: vc(nc)

   integer :: i, j

   read (NIN, *) ncomb, NTDEP

   call setup(nc, nmodel, ntdep, ncomb)

   Tstar = 0.d0
   if (nmodel .eq. 1) then
      del1 = 1.0D0
      write (nout, *) ' Model: Soave-Redlich-Kwong (1972)'
   else
      del1 = 1.0D0 + sqrt(2.0)
      write (nout, *) ' Model: Peng-Robinson (1976)'
   end if

   write (nout, *) ' Fluid           Tc(K)       Pc(bar)  Vceos(L/mol)    W'
   do i = 1, nc
      read (NIN, '(A)') fluid(i)
      read (NIN, *) Tc(i), Pc(i), OM(i), Vc(i)
      dceos(i) = 1/Vc(i)
      write (nout, 1) fluid(i), Tc(i), Pc(i), Vc(i), OM(i)
      read (NIN, *) ac(i), b(i), k(i)

      Kij(i, i) = 0.0D0
      Lij(i, i) = 0.0D0
      if (i .gt. 1) then
         if (ncomb .lt. 2) then
            read (NIN, *) (Kij(j, i), j=1, i - 1)
            Kij(i, :i - 1) = Kij(:i - 1, i)
            if (NTDEP >= 1) read (NIN, *) (Tstar(j, i), j=1, i - 1)
            Tstar(i, :i - 1) = Tstar(:i - 1, i)
            if (NTDEP == 2) read (NIN, *) (Kinf(j, i), j=1, i - 1)
            Kinf(i, :i - 1) = Kinf(:i - 1, i)
            read (NIN, *) (lij(j, i), j=1, i - 1)
            lij(i, :i - 1) = lij(:i - 1, i)
         end if
      end if
   end do

   write (nout, *) 'Fluid     ac(bar*L2/mol2)  b(L/mol)    d1      k'
   do I = 1, NC
      write (nout, 1) fluid(i), ac(i), b(i), del1(i), k(i)
   end do

   write (NOUT, *)
   if (ncomb .lt. 2) then
      if (NTDEP .eq. 0) then
         write (NOUT, *) '  Kij MATRIX'
      else
         write (NOUT, *) '  K0ij MATRIX'
      end if
      do I = 1, NC
         write (NOUT, 6) FLUID(I), (Kij(j, i), j=1, i - 1)
      end do
      if (NTDEP .eq. 1) then
         write (NOUT, *)
         write (NOUT, *) '  T* MATRIX'
         do I = 1, NC
            write (NOUT, 6) FLUID(I), (Tstar(j, i), j=1, i - 1)
         end do
      end if
      write (NOUT, *)
      write (NOUT, *) '  LIJ MATRIX'
      do I = 1, NC
         write (NOUT, 6) FLUID(I), (Lij(j, i), j=1, i - 1)
      end do
   end if


   select case(nmodel)
      case (1)
         call SRK_factory(z, tc_in=tc, pc_in=pc, w_in=om)
      case (2)
         call PR76_factory(z, tc_in=tc, pc_in=pc, w_in=om)
      case (3)
         call PR78_factory(z, tc_in=tc, pc_in=pc, w_in=om)
   end select
   
   write (NOUT, *)
   write (NOUT, *) ' Combining rules:'
   if (ncomb .eq. 0) then
      write (NOUT, *) ' 0: Classical or van der Waals '
      do i = 1, nc
         do j = i, nc
            bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
            bij(j, i) = bij(i, j)
         end do
      end do
   end if

1  format(A18, F8.3, 5x, F7.3, 3x, F7.3, 3x, F7.3)
6  format(A18, 20F10.5)
7  format(9x, F7.4, 2x, F7.4)
8  format(9x, F7.2, 2x, F7.2)
end

subroutine HelmSRKPR(nc, ND, NT, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   use constants, only: pr, R
   use system, only: del1, mixing_rule

   implicit real(pr)(A - H, O - Z)

   real(pr) :: rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)
   real(pr) :: dBi(nc), dBij(nc, nc)
   real(pr) :: dDi(nc), dDij(nc, nc), dDiT(nc)

   TOTN = sum(rn)
   D1 = del1(1)
   D2 = (1 - D1)/(1 + D1)

   if (mixing_rule .lt. 2) then
      call Bnder(nc, rn, Bmix, dBi, dBij)
      call DandTnder(NT, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
   end if

   ! The f's and g's used here are for Ar, not F (reduced Ar)
   ! This requires to multiply by R all g, f and its derivatives as defined by Mollerup
   f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
   g = R*log(1 - Bmix/V)
   fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
   fB = -(f + V*fv)/Bmix
   gv = R*Bmix/(V*(V - Bmix))
   fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
   gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

   ! Reduced Helmholtz Energy and derivatives
   Ar = -TOTN*g*T - D*f
   ArV = -TOTN*gv*T - D*fv
   ArV2 = -TOTN*gv2*T - D*fv2

   AUX = R*T/(V - Bmix)
   FFB = TOTN*AUX - D*fB
   FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
   FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2
   do i = 1, nc
      Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i)
      ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i)
      if (ND .eq. 2) then
         do j = 1, i
            Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                         + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
            Arn2(j, i) = Arn2(i, j)
         end do
      end if
   end do
   ! TEMPERATURE DERIVATIVES
   if (NT .eq. 1) then
      ArT = -TOTN*g - dDdT*f
      ArTV = -TOTN*gv - dDdT*fV
      ArTT = -dDdT2*f
      do i = 1, nc
         ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i)
      end do
   end if
end subroutine HelmSRKPR