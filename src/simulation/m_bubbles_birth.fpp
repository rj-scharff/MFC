!>
!! @file
!! @brief Contains module m_bubbles_birth

!> @brief Nucleation birth source for one Euler-Euler bubble class.
!!
!! Phase 2 scaffold. The birth rate is identically zero here, so enabling the
!! feature must change nothing: every contribution below is a multiple of the
!! rate, and adding exact zero is exact. That identity is the test this stage
!! exists to pass, before any physics is allowed to depend on the plumbing.
!!
!! The rate and the newborn state are separate routines on purpose. Selecting a
!! nucleation model replaces the first; the equation-of-state work replaces the
!! second; neither requires touching the bookkeeping that applies them.
!!
!! Birth sources the number density, the first radius moment, and the first
!! wall-velocity moment. It does not source the void fraction, which
!! `s_comp_alpha_from_n` derives from those under `adv_n`.
#:include 'macros.fpp'

module m_bubbles_birth

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: f_bubble_birth_rate, f_blake_pressure, s_bubble_newborn_state

contains

    !> Birth rate per unit mixture volume per unit time.
    !! @param cell_number_density Bubble number density in the cell
    !! @param cell_pressure Liquid pressure in the cell
    !> Liquid pressure below which a nucleus of radius @p fR0 has no equilibrium. The polytropic wall pressure this solver already
    !! integrates is p_l(R) = pv + p_g0 (R0/R)**(3 gam) - 2 sigma/R, p_g0 = Ca + 2/(Web R0), which is the Blake relation. It is
    !! non-monotone, and its minimum is the tension past which the nucleus runs away. Stationarity gives R_c = (3 gam A
    !! Web/2)**(1/(3 gam - 1)), A = p_g0 R0**(3 gam), in closed form, so no iteration is needed on device.
    !!
    !! Without surface tension the relation is monotone and its infimum is pv
    !! itself, which is the threshold this routine then returns -- so the
    !! surface-tension-free case reproduces the earlier vapour-pressure gate
    !! exactly rather than approximately.
    !! @param fR0 Equilibrium radius of the nucleus
    function f_blake_pressure(fR0) result(critical_pressure)

        $:GPU_ROUTINE(function_name='f_blake_pressure', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: fR0
        real(wp)             :: gas_pressure, amount, critical_radius, exponent
        real(wp)             :: critical_pressure

        if (f_is_default(Web)) then
            critical_pressure = pv
            return
        end if

        exponent = 3._wp*gam
        gas_pressure = Ca + 2._wp/(Web*fR0)
        amount = gas_pressure*fR0**exponent
        critical_radius = (exponent*amount*Web/2._wp)**(1._wp/(exponent - 1._wp))
        critical_pressure = pv + amount/critical_radius**exponent - 2._wp/(critical_radius*Web)

    end function f_blake_pressure

    function f_bubble_birth_rate(cell_number_density, cell_pressure, fR0) result(birth_rate)

        $:GPU_ROUTINE(function_name='f_bubble_birth_rate', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: cell_number_density, cell_pressure, fR0
        real(wp)             :: birth_rate, remaining_sites

        ! A stable liquid does not nucleate, and a stabilised nucleus does not
        ! either until the tension exceeds what its own curvature can carry. The
        ! gate is therefore the Blake threshold of the nucleus this class would
        ! create, not the vapour pressure: the two differ by orders of magnitude
        ! at nanometre sizes, and the vapour-pressure form lets every newborn
        ! grow unconditionally regardless of how small it is.
        !
        ! An unset vapour pressure is the sentinel dflt_real, which is negative,
        ! so the comparison fails and birth stays silent. Silence is the safe
        ! direction: without a vapour pressure there is no way to tell whether
        ! the liquid is metastable, and nucleating anyway would be a guess.
        if (cell_pressure < f_blake_pressure(fR0)) then
            birth_rate = bubble_birth_rate
        else
            birth_rate = 0._wp
        end if

        ! A quenched population is finite. The sites are in the liquid before the
        ! shot and each is spent when it activates, so the operator that models
        ! one must run out. The count already spent needs no field of its own:
        ! the number density is exactly that count, birth only ever adds to it,
        ! and it is transported with the liquid that carries the sites.
        !
        ! What is left is limited, not clipped after the fact. Clipping the state
        ! would let the source ask for sites that do not exist and then quietly
        ! discard the excess, which is the same concealment the void-fraction
        ! guard exists to prevent.
        !
        ! Compression raises the number density without activating anything, so a
        ! parcel can carry more bubbles per unit volume than the uncompressed
        ! liquid had sites. That reads as an empty inventory and closes the gate,
        ! which is the right direction: an inventory per unit mass would be exact,
        ! and the difference is the density ratio, one percent at the tensions
        ! where this gate is open at all.
        if (.not. f_is_default(bubble_site_density)) then
            remaining_sites = max(bubble_site_density - cell_number_density, 0._wp)
            birth_rate = min(birth_rate, remaining_sites/dt)
        end if

    end function f_bubble_birth_rate

    !> Radius and wall velocity of a newborn bubble in class @p bin.
    !! @param bin Bubble class index
    !! @param newborn_radius Radius assigned to newborn bubbles
    !! @param newborn_velocity Wall velocity assigned to newborn bubbles
    subroutine s_bubble_newborn_state(bin, newborn_radius, newborn_velocity)

        $:GPU_ROUTINE(function_name='s_bubble_newborn_state', parallelism='[seq]', cray_inline=True)

        integer, intent(in)   :: bin
        real(wp), intent(out) :: newborn_radius, newborn_velocity

        ! A newborn enters its own class at that class's equilibrium radius, at
        ! rest. The critical radius and the nucleus wall state replace this once
        ! the equation of state can supply them; see the corridor binding record.
        newborn_radius = R0(bin)
        newborn_velocity = 0._wp

    end subroutine s_bubble_newborn_state

end module m_bubbles_birth
