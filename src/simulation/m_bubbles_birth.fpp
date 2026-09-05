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

    private; public :: f_bubble_birth_rate, s_bubble_newborn_state

contains

    !> Birth rate per unit mixture volume per unit time.
    !! @param cell_void_fraction Void fraction in the cell
    !! @param cell_pressure Liquid pressure in the cell
    function f_bubble_birth_rate(cell_void_fraction, cell_pressure) result(birth_rate)

        $:GPU_ROUTINE(function_name='f_bubble_birth_rate', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: cell_void_fraction, cell_pressure
        real(wp)             :: birth_rate

        ! Identically zero for the Phase 2 no-op identity test. A nucleation
        ! model replaces this body; the arguments are the state it will need,
        ! and are referenced so the interface is not silently wrong when it is.
        birth_rate = 0._wp*(cell_void_fraction + cell_pressure)

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
