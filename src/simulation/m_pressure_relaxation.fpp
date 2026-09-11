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

                    ! Nucleation runs BEFORE the mechanical relaxation and has to. The seed places the vapour at a large
                    ! positive pressure so that the relaxation has a valid isentrope reference; s_correct_internal_energies
                    ! below would overwrite that with the mixture pressure, which for a cell still almost entirely liquid is
                    ! the tension itself. A vapour with pi_inf = 0 has no state there, so it clamps at zero and the volume
                    ! constraint is never satisfied. The seed must therefore be relaxed in the same visit that places it.
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

                    ! Mass and energy transfer runs AFTER the mechanical relaxation, at the equilibrated mixture pressure.
                    ! Evaporation reads the liquid's phasic pressure out of its internal energy and volume fraction, and
                    ! s_equilibrate_pressure updates the volume fractions without touching the energies - so evaluated any
                    ! earlier than the correction above it would read a post-relaxation alpha against a pre-relaxation
                    ! energy, which is not a state. Moving the mass then changes the mixture pressure through the formation
                    ! energies, so the phasic energies are slaved to the conserved total a second time.
                    if (spall_pressure < 0._wp) then
                        call s_evaporate_into_void(q_cons_vf, j, k, l)

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                            alpha(i) = q_cons_vf(eqn_idx%E + i)%sf(j, k, l)
                        end do

                        call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_mix, alpha, alpha_rho)

                        call s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)
                    end if
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

    !> Feed an open void from the liquid wall, as far as the cell's own energy allows
    !!
    !! A void that no mass enters is a vacuum, and a vacuum is what breaks a spall calculation: the liquid drains away as
    !! the plane opens, the mixture density follows it to zero, and the sound speed - a finite bulk modulus over that
    !! vanishing mass - diverges. Measured on the trajectory a plane actually takes, the mixture sound speed reaches
    !! 40000 m/s and the run stops on its Courant number.
    !!
    !! Evaporation is what fills it, and the mass stays behind when the liquid leaves. The target is the density the
    !! vapour would have in thermal equilibrium with the wall it is evaporating from. The cell rarely affords that: the
    !! latent heat of reaching it can exceed half the cell's internal energy, so the transfer is capped where the pressure
    !! would reach zero, and again by the liquid actually present. All three limits are properties of the cell, so this
    !! introduces no constant. It is self-limiting - once the vapour is at target the transfer is zero - and the mixture
    !! energy is untouched, so the latent heat is paid by the cell rather than invented.
    subroutine s_evaporate_into_void(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: alpha_l, alpha_v, alpha_rho_l, alpha_rho_v
        real(wp)                                               :: pres_l, temp_l, rho_v_eq, gamma_mix, latent, mass_e
        real(wp)                                               :: gibbs_diff

        alpha_l = q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l)
        alpha_v = q_cons_vf(vp + eqn_idx%adv%beg - 1)%sf(j, k, l)
        alpha_rho_l = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l)
        alpha_rho_v = q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)
        latent = qvs(vp) - qvs(lp)

        if (alpha_v > sgm_eps .and. alpha_l > sgm_eps .and. alpha_rho_l > sgm_eps .and. latent > 0._wp) then
            pres_l = ((q_cons_vf(lp + eqn_idx%int_en%beg - 1)%sf(j, k, l) - alpha_rho_l*qvs(lp))/alpha_l - pi_infs(lp))/gammas(lp)

            ! A liquid already in tension has nothing to spend on evaporating itself.
            if (pres_l > 0._wp) then
                temp_l = f_sg_thermal(pres_l, alpha_rho_l/alpha_l, isentrope_n(lp), isentrope_B(lp), cvs(lp))

                if (temp_l > 0._wp) then
                    ! Whether the liquid wants to evaporate at all, which is not a question the energy budget can answer.
                    ! This is m_phase_change's own Gibbs equality residual, which is g_liquid - g_vapour; it is positive
                    ! exactly where the liquid is superheated for the local pressure. Without it the transfer has no
                    ! thermodynamic stop: the flow re-pressurises the cell between calls, the energy cap admits another
                    ! transfer, and a void sitting happily above its saturation pressure is fed until its vapour reaches
                    ! the density of the liquid.
                    gibbs_diff = temp_l*((cvs(lp)*isentrope_n(lp) - cvs(vp)*isentrope_n(vp))*(1._wp - log(temp_l)) - (qvps(lp) &
                                         & - qvps(vp)) + cvs(lp)*(isentrope_n(lp) - 1._wp)*log(pres_l + isentrope_B(lp)) - cvs(vp) &
                                         & *(isentrope_n(vp) - 1._wp)*log(pres_l + isentrope_B(vp))) + qvs(lp) - qvs(vp)
                else
                    gibbs_diff = 0._wp
                end if

                rho_v_eq = (pres_l + isentrope_B(vp))/((isentrope_n(vp) - 1._wp)*cvs(vp)*max(temp_l, sgm_eps))
                gamma_mix = alpha_l*gammas(lp) + alpha_v*gammas(vp)

                mass_e = min(alpha_v*rho_v_eq - alpha_rho_v, pres_l*gamma_mix/latent, alpha_rho_l)

                if (mass_e > 0._wp .and. gibbs_diff > 0._wp) then
                    q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = alpha_rho_l - mass_e
                    q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = alpha_rho_v + mass_e
                    q_cons_vf(lp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = f_phase_internal_energy(pres_l, alpha_l, &
                              & alpha_rho_l - mass_e, gammas(lp), pi_infs(lp), qvs(lp))
                    q_cons_vf(vp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = f_phase_internal_energy(pres_l, alpha_v, &
                              & alpha_rho_v + mass_e, gammas(vp), pi_infs(vp), qvs(vp))
                end if
            end if
        end if

    end subroutine s_evaporate_into_void

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

        ! A cell whose every phase has been zeroed above leaves nothing to normalise against, and dividing by that sum
        ! turns each fraction into a NaN rather than leaving the cell empty. Reachable at a spall plane, where the last
        ! of the liquid leaves.
        if (sum_alpha > sgm_eps) then
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)/sum_alpha
            end do
        end if

    end subroutine s_correct_volume_fractions

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: pres_relax, f_pres, df_pres
        real(wp)                                               :: alpha_v_old, driving, d_alpha_max
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: pres_K_init, rho_K_s, alpha_eq
            logical, dimension(3)  :: is_vacuum
        #:else
            real(wp), dimension(num_fluids) :: pres_K_init, rho_K_s, alpha_eq
            logical, dimension(num_fluids)  :: is_vacuum
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
            is_vacuum(i) = .true.
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                & l) > sgm_eps) then
                ! Phasic internal energy carries the formation energy: alpha_rho_k*qv_k must be
                ! removed before inverting the stiffened-gas EOS, or a nonzero qv inflates the
                ! phasic pressure by rho_k*qv_k/gamma_k (this is what breaks the reactive burn).
                pres_K_init(i) = ((q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, &
                            & k, l)*qvs(i))/q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) - pi_infs(i))/gammas(i)
                ! A phase pinned at its own floor carries no usable isentrope reference: it is in a state its equation
                ! of state cannot represent, and the ratio below would be taken about a fiction. A vapour reaches that
                ! whenever it finds itself in a liquid still under tension, since with no stiffness its floor is zero.
                if (pres_K_init(i) > -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) then
                    is_vacuum(i) = .false.
                else if (vapor_saturation_floor > 0._wp .and. i == vp .and. pres_K_init(lp) < 0._wp) then
                    ! A vapour in a stretched liquid is not a vacuum. It sits near its saturation pressure while the
                    ! liquid carries the tension, and that difference is what opens a cavity. Declaring it a vacuum
                    ! instead excludes it from the volume fraction update at the end of this routine, so the cavity
                    ! cannot grow at all however a growth law is written. Held at the supplied saturation pressure it
                    ! keeps a usable isentrope reference and stays a participating phase.
                    !
                    ! The rescue is conditioned on the LIQUID being in tension, and must be. The volume fraction
                    ! advects at the speed of sound while the mass under it does not, so the leading edge of an
                    ! opening void is volume with almost no vapour in it - measured four decades below the saturated
                    ! density. The test that a phase carries mass is alpha_rho > sgm_eps, and sgm_eps is 1e-16, so it
                    ! does not exclude those cells. Rescuing them anchors an isentrope at the saturation pressure on a
                    ! density that is nowhere near it; the phases are then out of equilibrium by the whole of the
                    ! cell's tension, and because s_correct_internal_energies refreshes the reference on the way out,
                    ! the volume the solve hands to the vapour is never handed back. Every dip below zero pressure is
                    ! banked. That turns any cell holding a trace of transported alpha into a cavitation site whose
                    ! threshold is zero rather than spall_pressure, and since a unit of volume fraction is worth
                    ! (pi_l - pi_v)/Gamma - about 2 GPa for water - a cell that accumulates a fifth of one reads
                    ! hundreds of megapascals and drives a compression wave back into the liquid.
                    !
                    ! Requiring tension is not a guard bolted on: it is the condition under which a cavity grows at
                    ! all, and it is the same test the Rayleigh limiter below already applies to its driving pressure.
                    ! pres_K_init(lp) is available because lp = 1 and vp = 2 and this loop is sequential in i, so the
                    ! liquid has been visited first. A cell holding no liquid leaves it at zero, which is not tension,
                    ! and the rescue correctly does not fire there either.
                    pres_K_init(i) = vapor_saturation_floor
                    is_vacuum(i) = .false.
                else
                    pres_K_init(i) = -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp
                end if
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
                    if (.not. is_vacuum(i)) then
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
                        if (.not. is_vacuum(i)) then
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
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) > sgm_eps .and. .not. is_vacuum(i)) alpha_eq(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
        end do

        ! Give the void a finite expansion rate, if one was asked for. Without this the
        ! relaxation carries a nucleated cell to mechanical equilibrium inside one step, so
        ! the void opens as fast as the mesh allows and the model has no growth kinetics at
        ! all - which is what makes a threshold on pressure rate independent by construction.
        !
        ! For n_s sites per unit volume sharing a void alpha_v, each has radius
        ! R = (3 alpha_v/(4 pi n_s))**(1/3) and the interfacial area per unit volume is
        ! 3 alpha_v/R, so an interface moving at the Rayleigh speed sqrt(2 dp/(3 rho_l))
        ! opens volume at
        !     d(alpha_v)/dt = 3 (4 pi n_s/3)**(1/3) alpha_v**(2/3) sqrt(2 dp/(3 rho_l)).
        ! One parameter, no new field, and alpha_v**(1/3) grows linearly under it.
        !
        ! The driving pressure is the LIQUID's tension, not p_vapour - p_liquid: the seed
        ! places the vapour at a large positive pressure to give the isentrope a valid
        ! reference (s_nucleate_vapor), and that is a numerical device, not the pressure of a
        ! cavitation bubble. Using it would drive growth at 200 MPa instead of the tension.
        !
        ! Only growth is limited. A collapsing void is left to the equilibrium solve, which
        ! keeps this to the smallest change that buys the kinetics.
        if (nucleus_site_density > 0._wp .and. .not. is_vacuum(vp)) then
            alpha_v_old = q_cons_vf(vp + eqn_idx%adv%beg - 1)%sf(j, k, l)
            if (alpha_eq(vp) > alpha_v_old .and. alpha_v_old > sgm_eps) then
                driving = max(0._wp, -pres_K_init(lp))
                if (driving > 0._wp) then
                    ! The procedure runs once per Runge-Kutta stage, so the per-call
                    ! increment is the step divided by the number of stages: the operator is
                    ! a splitting rather than an RHS contribution, so this is as accurate as
                    ! the splitting allows and it makes the per-step total exactly dt.
                    d_alpha_max = 3._wp*(4._wp*pi*nucleus_site_density/3._wp)**(1._wp/3._wp)*alpha_v_old**(2._wp/3._wp) &
                                         & *sqrt(2._wp*driving/(3._wp*q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, &
                                         & l)/max(q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l), &
                                         & sgm_eps)))*dt/real(time_stepper, wp)
                    if (alpha_eq(vp) - alpha_v_old > d_alpha_max) then
                        ! Hold back the vapour and return the volume it did not take to the
                        ! liquid, so the fractions still close. The cell is then left short of
                        ! equilibrium, which is the point: its tension is only partly relieved.
                        alpha_eq(lp) = alpha_eq(lp) + (alpha_eq(vp) - alpha_v_old - d_alpha_max)
                        alpha_eq(vp) = alpha_v_old + d_alpha_max
                    end if
                end if
            end if
        end if

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) > sgm_eps .and. .not. is_vacuum(i)) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = alpha_eq(i)
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
