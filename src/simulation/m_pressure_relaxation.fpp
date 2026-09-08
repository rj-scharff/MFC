!>
!! @file
!! @brief Contains module m_pressure_relaxation

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Pressure relaxation for the six-equation multi-component model via Newton--Raphson equilibration and volume-fraction
!! correction
module m_pressure_relaxation

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel, f_pressure, f_phase_internal_energy, &
        & f_sg_thermal

    implicit none

    ! Liquid and its vapour, following m_phase_change's convention for the reacting pair.
    integer, parameter :: lp = 1
    integer, parameter :: vp = 2

    ! Volume fraction of a fresh nucleus. A trigger, not a physical size: the relaxation below carries the cell to an
    ! equilibrium that moved in the sixth significant figure over a six-decade sweep of this value.
    real(wp), parameter :: nucleus_volume_seed = 1.e-6_wp

    private; public :: s_pressure_relaxation_procedure

contains

    !> The main pressure relaxation procedure
    subroutine s_pressure_relaxation_procedure(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: i, j, k, l

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        real(wp) :: rho, gamma, pi_inf, qv_mix

        ! Formed here, not one call deeper: CCE OpenACC accepts a num_fluids-sized array passed to a device routine from a
        ! parallel-loop body, and rejects the same call from inside another acc routine seq.
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, alpha_rho, alpha, rho, gamma, pi_inf, qv_mix]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (mpp_lim) call s_correct_volume_fractions(q_cons_vf, j, k, l)

                    if (spall_pressure < 0._wp) call s_nucleate_vapor(q_cons_vf, j, k, l)

                    if (s_needs_pressure_relaxation(q_cons_vf, j, k, l)) then
                        call s_equilibrate_pressure(q_cons_vf, j, k, l)
                    end if

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                        alpha(i) = q_cons_vf(eqn_idx%E + i)%sf(j, k, l)
                    end do

                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_mix, alpha, alpha_rho)

                    call s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pressure_relaxation_procedure

    !> Open a vapour nucleus in a cell whose liquid has been stretched past the spall threshold
    !!
    !! Only a seed is placed. The state that follows is not chosen here: the relaxation below carries the cell to its own
    !! equilibrium, which is the liquid springing back towards its release density and vacating the volume its elastic strain
    !! had occupied. That equilibrium is a property of the cell, not of the seed - across a six-decade sweep of the seeded
    !! mass the resulting void fraction moved in the sixth significant figure, and at every tension it matched the liquid's
    !! elastic strain to four figures.
    !!
    !! The seed must sit at a positive pressure. A vapour has no stiffness, so its floor in s_equilibrate_pressure is zero;
    !! seeded at the liquid's own tension it is clamped there, its isentrope reference degenerates to unity, and the volume
    !! constraint is never satisfied - the cell is left with sum(alpha) < 1. Vaporising the liquid already present in the seed
    !! volume, at constant mass and temperature, places it well above zero and introduces no further constant.
    subroutine s_nucleate_vapor(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: alpha_l, alpha_rho_l, rho_l, pres_l, temp_l, pres_v, mass_v

        alpha_l = q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l)
        alpha_rho_l = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l)

        ! Nothing to nucleate from, or a nucleus is already open.
        if (alpha_l > 1._wp - nucleus_volume_seed .and. alpha_rho_l > sgm_eps) then
            rho_l = alpha_rho_l/alpha_l
            pres_l = ((q_cons_vf(lp + eqn_idx%int_en%beg - 1)%sf(j, k, l) - alpha_rho_l*qvs(lp))/alpha_l - pi_infs(lp))/gammas(lp)

            if (pres_l <= spall_pressure) then
                temp_l = f_sg_thermal(pres_l, rho_l, isentrope_n(lp), isentrope_B(lp), cvs(lp))
                pres_v = (isentrope_n(vp) - 1._wp)*cvs(vp)*rho_l*temp_l - isentrope_B(vp)
                mass_v = nucleus_volume_seed*rho_l

                q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l) = alpha_l - nucleus_volume_seed
                q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = alpha_rho_l - mass_v
                q_cons_vf(vp + eqn_idx%adv%beg - 1)%sf(j, k, l) = nucleus_volume_seed
                q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = mass_v

                ! Phasic energies supply the relaxation's isentrope references only; the mixture energy is untouched, so the
                ! latent heat of the seeded mass is drawn from the cell rather than invented.
                q_cons_vf(lp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = f_phase_internal_energy(pres_l, &
                          & alpha_l - nucleus_volume_seed, alpha_rho_l - mass_v, gammas(lp), pi_infs(lp), qvs(lp))
                q_cons_vf(vp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = f_phase_internal_energy(pres_v, nucleus_volume_seed, mass_v, &
                          & gammas(vp), pi_infs(vp), qvs(vp))
            end if
        end if

    end subroutine s_nucleate_vapor

    !> Check if pressure relaxation is needed for this cell
    logical function s_needs_pressure_relaxation(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        integer                                             :: i

        s_needs_pressure_relaxation = .true.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > (1._wp - sgm_eps)) then
                s_needs_pressure_relaxation = .false.
            end if
        end do

    end function s_needs_pressure_relaxation

    !> Correct volume fractions to physical bounds
    subroutine s_correct_volume_fractions(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: sum_alpha
        integer                                                :: i

        sum_alpha = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if ((q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) < 0._wp) .or. (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) < 0._wp)) then
                q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = 0._wp
            end if
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > 1._wp) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 1._wp
            sum_alpha = sum_alpha + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)/sum_alpha
        end do

    end subroutine s_correct_volume_fractions

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: pres_relax, f_pres, df_pres
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: pres_K_init, rho_K_s
        #:else
            real(wp), dimension(num_fluids) :: pres_K_init, rho_K_s
        #:endif
        integer, parameter :: MAX_ITER = 50
        ! Pressure relaxation convergence tolerance
        real(wp), parameter :: TOLERANCE = 1.e-10_wp
        integer             :: iter, i

        ! Initialize pressures
        pres_relax = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            ! A phase needs mass here, not just volume: everything below is built on rho_K_s, which is its partial density
            ! over its volume fraction, and is then divided by. The leading edge of an opening void reaches exactly that
            ! state - the non-conservative volume fraction advects outwards while the partial density is clipped to zero -
            ! and it is a vacuum, not an error. Its volume is real and is carried through the constraint below unchanged.
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                & l) > sgm_eps) then
                ! Phasic internal energy carries the formation energy: alpha_rho_k*qv_k must be
                ! removed before inverting the stiffened-gas EOS, or a nonzero qv inflates the
                ! phasic pressure by rho_k*qv_k/gamma_k (this is what breaks the reactive burn).
                pres_K_init(i) = ((q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, &
                            & k, l)*qvs(i))/q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) - pi_infs(i))/gammas(i)
                if (pres_K_init(i) <= -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) pres_K_init(i) = -(1._wp - 1.e-8_wp) &
                    & *isentrope_B(i) + 1.e-8_wp
            else
                pres_K_init(i) = 0._wp
            end if
            pres_relax = pres_relax + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)*pres_K_init(i)
        end do

        ! Newton-Raphson iteration
        f_pres = 1.e-9_wp
        df_pres = 1.e9_wp
        $:GPU_LOOP(parallelism='[seq]')
        do iter = 0, MAX_ITER - 1
            if (abs(f_pres) > TOLERANCE) then
                pres_relax = pres_relax - f_pres/df_pres

                ! Enforce pressure bounds. Only a phase carrying mass bounds the pressure: a vacuum has neither stiffness
                ! nor matter, and letting its zero isentrope_B floor the cell holds a liquid under real tension at zero
                ! instead, which the volume constraint below then pays for.
                do i = 1, num_fluids
                    if (q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) > sgm_eps) then
                        if (pres_relax <= -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) pres_relax = -(1._wp - 1.e-8_wp) &
                            & *isentrope_B(i) + 1.e-8_wp
                    end if
                end do

                ! Newton-Raphson step
                f_pres = -1._wp
                df_pres = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps) then
                        if (q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) > sgm_eps) then
                            ! Isentropic relation: rho = rho0 * (p/p0)^(1/gamma), Saurel et al. JFM (2009)
                            rho_K_s(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                    & l)/max(q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l), &
                                    & sgm_eps)*((pres_relax + isentrope_B(i))/(pres_K_init(i) + isentrope_B(i))) &
                                    & **(1._wp/isentrope_n(i))
                            f_pres = f_pres + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
                            df_pres = df_pres - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                                          & l)/(isentrope_n(i)*rho_K_s(i)*(pres_relax + isentrope_B(i)))
                        else
                            ! A vacuum holds its volume however the pressure moves, so it enters the constraint as a
                            ! constant and the phases carrying mass share what is left. Deleting it instead would make the
                            ! liquid expand into it at fixed mass, which reads as tension the flow never applied.
                            f_pres = f_pres + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
                        end if
                    end if
                end do
            end if
        end do

        ! Update volume fractions
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                & l) > sgm_eps) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                & l)/rho_K_s(i)
        end do

    end subroutine s_equilibrate_pressure

    !> Correct internal energies using equilibrated pressure
    subroutine s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp), intent(in)                                   :: rho, gamma, pi_inf, qv_mix
        real(wp)                                               :: dyn_pres, pres_relax
        integer                                                :: i

        dyn_pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            dyn_pres = dyn_pres + 5.e-1_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do

        pres_relax = f_pressure(q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_pres, gamma, pi_inf, qv_mix)

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = f_phase_internal_energy(pres_relax, &
                      & q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l), q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l), &
                      & gammas(i), pi_infs(i), qvs(i))
        end do

    end subroutine s_correct_internal_energies

end module m_pressure_relaxation
