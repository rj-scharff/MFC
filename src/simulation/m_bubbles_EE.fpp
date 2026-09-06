!>
!! @file
!! @brief Contains module @ref m_bubbles_ee "m_bubbles_EE"

#:include 'macros.fpp'

!> @brief Computes ensemble-averaged (Euler--Euler) bubble source terms for radius, velocity, pressure, and mass transfer
module m_bubbles_EE

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_bubbles
    use m_bubbles_birth

    implicit none

    real(wp), allocatable, dimension(:,:,:)   :: bub_adv_src
    real(wp), allocatable, dimension(:,:,:)   :: bub_n_src
    real(wp), allocatable, dimension(:,:,:,:) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    $:GPU_DECLARE(create='[bub_adv_src, bub_n_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src]')

    type(scalar_field) :: divu  !< matrix for div(u)
    $:GPU_DECLARE(create='[divu]')

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    $:GPU_DECLARE(create='[rs, vs, ms, ps]')

    !> State of a cell whose adaptive sub-integration did not converge, so that the abort can say where it failed instead of only
    !! that it failed. Written from inside the parallel loop without synchronisation, so with more than one failing cell it holds an
    !! arbitrary one of them; that is enough to diagnose with, and nothing computes from it.
    real(wp), dimension(7) :: adap_dt_fail_state
    $:GPU_DECLARE(create='[adap_dt_fail_state]')

contains

    !> Initialize the Euler-Euler bubble module
    impure subroutine s_initialize_bubbles_EE_module

        integer :: l

        @:ALLOCATE(rs(1:nb))
        @:ALLOCATE(vs(1:nb))
        @:ALLOCATE(ps(1:nb))
        @:ALLOCATE(ms(1:nb))

        do l = 1, nb
            rs(l) = qbmm_idx%rs(l)
            vs(l) = qbmm_idx%vs(l)
            if (.not. polytropic) then
                ps(l) = qbmm_idx%ps(l)
                ms(l) = qbmm_idx%ms(l)
            else
                ps(l) = rs(l)
                ms(l) = rs(l)
            end if
        end do

        $:GPU_UPDATE(device='[rs, vs]')
        $:GPU_UPDATE(device='[ps, ms]')

        @:ALLOCATE(divu%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(divu)

        @:ALLOCATE(bub_adv_src(0:m, 0:n, 0:p))
        @:ALLOCATE(bub_n_src(0:m, 0:n, 0:p))
        @:ALLOCATE(bub_r_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_v_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_p_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_m_src(0:m, 0:n, 0:p, 1:nb))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

    end subroutine s_initialize_bubbles_EE_module

    !> Compute the bubble volume fraction alpha from the bubble number density
    subroutine s_comp_alpha_from_n(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp)                                               :: nR3bar, alf_max
        integer(wp)                                            :: i, j, k, l

        alf_max = 0._wp

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, nR3bar]', reduction='[[alf_max]]', reductionOp='[max]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    nR3bar = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        nR3bar = nR3bar + weight(i)*(q_cons_vf(rs(i))%sf(j, k, l))**3._wp
                    end do
                    q_cons_vf(eqn_idx%alf)%sf(j, k, l) = (4._wp*pi*nR3bar)/(3._wp*q_cons_vf(eqn_idx%n)%sf(j, k, l)**2._wp)
                    alf_max = max(alf_max, q_cons_vf(eqn_idx%alf)%sf(j, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The void fraction is derived here rather than transported, so nothing
        ! downstream bounds it: with the number density roughly fixed it follows
        ! R**3, and a runaway radius carries it past one. A cell holding more
        ! bubble than cell is impossible rather than merely outside the closure's
        ! range, so it is reported and stopped instead of clipped -- clipping
        ! would hide the runaway and let every later step inherit it. The check
        ! lives here, where alpha is formed, rather than beside the ICFL test,
        ! because a realizability violation must not depend on whether run-time
        ! diagnostics happen to be switched on.
        if (alf_max >= 1._wp) then
            print *, 'max void fraction', alf_max
            call s_mpi_abort('Void fraction reached one: a cell holds more bubble ' // 'volume than cell volume. Exiting.')
        end if

    end subroutine s_comp_alpha_from_n

    !> Compute the right-hand side for Euler-Euler bubble transport
    subroutine s_compute_bubbles_EE_rhs(idir, q_prim_vf, divu_in)

        integer, intent(in)                                 :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(inout)                   :: divu_in  !< matrix for div(u)
        integer                                             :: j, k, l

        if (idir == 1) then
            if (.not. qbmm) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu_in%sf(j, k, l) = 0._wp
                            divu_in%sf(j, k, l) = 5.e-1_wp/dx(j)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j + 1, k, &
                                       & l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j - 1, k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (idir == 2) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dy(k)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, &
                                   & k + 1, l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k - 1, l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (idir == 3) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dz(l)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, &
                                   & l + 1) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubbles_EE_rhs

    !> Compute the Euler-Euler bubble source terms
    impure subroutine s_compute_bubble_EE_source(q_cons_vf, q_prim_vf, rhs_vf, divu_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(in) :: divu_in  !< matrix for div(u)
        real(wp) :: rddot
        real(wp) :: pb_local, mv_local, vflux, pbdot
        real(wp) :: n_tait, B_tait
        real(wp) :: chi_vw_l, k_mw_l, rho_mw_l     !< Per-thread bubble-wall scratch (avoid module-scalar race)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: Rtmp, Vtmp
            real(wp), dimension(3) :: myalpha, myalpha_rho
        #:else
            real(wp), dimension(nb)         :: Rtmp, Vtmp
            real(wp), dimension(num_fluids) :: myalpha, myalpha_rho
        #:endif
        real(wp)           :: myR, myV, alf, myP, myRho, R2Vav, R3
        real(wp)           :: nbub                                  !< Bubble number density
        integer            :: i, j, k, l, q, ii                     !< Loop variables
        integer            :: adap_dt_stop_sum, adap_dt_stop        !< Fail-safe exit if max iteration count reached
        real(wp)           :: entry_R, entry_V                      !< Sub-integration entry state, diagnostic only
        character(len=250) :: fail_message
        real(wp)           :: dmMass_n, dmBeta_c, dmBeta_t, dmCson  !< Lagrange-only arguments, unused here
        real(wp)           :: birth_rate, newborn_radius, newborn_velocity

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    bub_adv_src(j, k, l) = 0._wp
                    bub_n_src(j, k, l) = 0._wp

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        bub_r_src(j, k, l, q) = 0._wp
                        bub_v_src(j, k, l, q) = 0._wp
                        bub_p_src(j, k, l, q) = 0._wp
                        bub_m_src(j, k, l, q) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        adap_dt_stop_sum = 0
        $:GPU_PARALLEL_LOOP(private='[j, k, l, Rtmp, Vtmp, myalpha_rho, myalpha, myR, myV, alf, myP, myRho, R2Vav, R3, nbub, &
                            & pb_local, mv_local, vflux, pbdot, rddot, n_tait, B_tait, adap_dt_stop, chi_vw_l, k_mw_l, rho_mw_l, &
                            & birth_rate, newborn_radius, newborn_velocity, entry_R, entry_V]', collapse=3, copy='[adap_dt_stop_sum]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (adv_n) then
                        nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                    else
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            Rtmp(q) = q_prim_vf(rs(q))%sf(j, k, l)
                            Vtmp(q) = q_prim_vf(vs(q))%sf(j, k, l)
                        end do

                        R3 = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            R3 = R3 + weight(q)*Rtmp(q)**3._wp
                        end do

                        nbub = (3._wp/(4._wp*pi))*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/R3
                    end if

                    if (.not. adap_dt) then
                        R2Vav = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            R2Vav = R2Vav + weight(q)*Rtmp(q)**2._wp*Vtmp(q)
                        end do

                        bub_adv_src(j, k, l) = 4._wp*pi*nbub*R2Vav
                    end if

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        $:GPU_LOOP(parallelism='[seq]')
                        do ii = 1, num_fluids
                            myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                            myalpha(ii) = q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)
                        end do

                        if (num_fluids == 1) then
                            myRho = myalpha_rho(1)
                            n_tait = gammas(1)
                            B_tait = pi_infs(1)/pi_fac
                        else
                            myRho = 0._wp
                            n_tait = 0._wp
                            B_tait = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do ii = 1, num_fluids
                                myRho = myRho + myalpha_rho(ii)
                                n_tait = n_tait + myalpha(ii)*gammas(ii)
                                B_tait = B_tait + myalpha(ii)*pi_infs(ii)/pi_fac
                            end do
                        end if

                        n_tait = 1._wp/n_tait + 1._wp  ! make this the usual little 'gamma'
                        B_tait = B_tait*(n_tait - 1)/n_tait  ! make this the usual pi_inf

                        myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                        alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                        myR = q_prim_vf(rs(q))%sf(j, k, l)
                        myV = q_prim_vf(vs(q))%sf(j, k, l)

                        ! Nucleation is the operator that creates a population, so it cannot be
                        ! conditioned on one already existing. The branch below skips bubble
                        ! *dynamics* where there are no bubbles, which is correct; birth is
                        ! evaluated here, outside it, and applied after it.
                        if (bubble_birth) then
                            birth_rate = f_bubble_birth_rate(alf, myP, R0(q))
                            call s_bubble_newborn_state(q, newborn_radius, newborn_velocity)
                        end if

                        if (alf < small_alf) then
                            bub_adv_src(j, k, l) = 0._wp
                            bub_n_src(j, k, l) = 0._wp
                            bub_r_src(j, k, l, q) = 0._wp
                            bub_v_src(j, k, l, q) = 0._wp
                            if (.not. polytropic) then
                                bub_p_src(j, k, l, q) = 0._wp
                                bub_m_src(j, k, l, q) = 0._wp
                            end if
                        else
                            if (.not. polytropic) then
                                pb_local = q_prim_vf(ps(q))%sf(j, k, l)
                                mv_local = q_prim_vf(ms(q))%sf(j, k, l)
                                call s_bwproperty(pb_local, q, chi_vw_l, k_mw_l, rho_mw_l)
                                call s_vflux(myR, myV, pb_local, mv_local, q, vflux, fchi_vw=chi_vw_l, frho_mw=rho_mw_l)
                                pbdot = f_bpres_dot(vflux, myR, myV, pb_local, mv_local, q, fk_mw=k_mw_l)
                                bub_p_src(j, k, l, q) = nbub*pbdot
                                bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                            else
                                pb_local = 0._wp; mv_local = 0._wp; vflux = 0._wp; pbdot = 0._wp
                            end if

                            adap_dt_stop = 0

                            ! Adaptive time stepping
                            if (adap_dt) then
                                ! f_advance_step updates myR and myV in place, so the entry state is
                                ! kept here: a failure is far more informative for where it started
                                ! than for where it gave up.
                                entry_R = myR
                                entry_V = myV

                                ! The class index and the vapour mass are real arguments here, not
                                ! placeholders: under polytropic = F the sub-integrator advances pb and
                                ! mv, and it needs the bin index to evaluate the wall transfer.
                                adap_dt_stop = f_advance_step(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, &
                                                              & bub_adv_src(j, k, l), divu_in%sf(j, k, l), q, mv_local, dmMass_n, &
                                                              & dmBeta_c, dmBeta_t, dmCson)

                                if (adap_dt_stop /= 0) then
                                    adap_dt_fail_state(1) = real(j, wp)
                                    adap_dt_fail_state(2) = entry_R
                                    adap_dt_fail_state(3) = entry_V
                                    adap_dt_fail_state(4) = myR
                                    adap_dt_fail_state(5) = myV
                                    adap_dt_fail_state(6) = alf
                                    adap_dt_fail_state(7) = myP
                                end if

                                q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                                q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV
                                if (.not. polytropic) then
                                    q_cons_vf(ps(q))%sf(j, k, l) = nbub*pb_local
                                    q_cons_vf(ms(q))%sf(j, k, l) = nbub*mv_local
                                end if
                            else
                                rddot = f_rddot(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, bub_adv_src(j, &
                                                & k, l), divu_in%sf(j, k, l), dmCson)
                                bub_v_src(j, k, l, q) = nbub*rddot
                                bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            adap_dt_stop_sum = adap_dt_stop_sum + adap_dt_stop
                        end if

                        ! Applied whether or not the cell had a population to begin with. The
                        ! adaptive path advances the class directly, so birth is an increment over
                        ! the same half step the sub-integrator covered; otherwise it enters as a
                        ! rate beside the other sources. Newborns source the number density and the
                        ! moments, never the void fraction, which adv_n derives from them.
                        if (bubble_birth) then
                            if (adap_dt) then
                                q_cons_vf(eqn_idx%n)%sf(j, k, l) = q_cons_vf(eqn_idx%n)%sf(j, k, l) + 5.e-1_wp*dt*birth_rate
                                q_cons_vf(rs(q))%sf(j, k, l) = q_cons_vf(rs(q))%sf(j, k, l) + 5.e-1_wp*dt*birth_rate*newborn_radius
                                q_cons_vf(vs(q))%sf(j, k, l) = q_cons_vf(vs(q))%sf(j, k, &
                                          & l) + 5.e-1_wp*dt*birth_rate*newborn_velocity
                            else
                                bub_n_src(j, k, l) = bub_n_src(j, k, l) + birth_rate
                                bub_r_src(j, k, l, q) = bub_r_src(j, k, l, q) + birth_rate*newborn_radius
                                bub_v_src(j, k, l, q) = bub_v_src(j, k, l, q) + birth_rate*newborn_velocity
                            end if
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (adap_dt .and. adap_dt_stop_sum > 0) then
            $:GPU_UPDATE(host='[adap_dt_fail_state]')
            write (fail_message, '(A,I0,6(A,ES13.6))') "Adaptive time stepping failed to converge. Cell i = ", &
                   & int(adap_dt_fail_state(1)), "; entry R = ", adap_dt_fail_state(2), ", entry Rdot = ", adap_dt_fail_state(3), &
                   & "; exit R = ", adap_dt_fail_state(4), ", exit Rdot = ", adap_dt_fail_state(5), "; alpha = ", &
                   & adap_dt_fail_state(6), ", p = ", adap_dt_fail_state(7)
            call s_mpi_abort(trim(fail_message))
        end if

        if (.not. adap_dt) then
            $:GPU_PARALLEL_LOOP(private='[i, k, l, q]', collapse=3)
            do l = 0, p
                do q = 0, n
                    do i = 0, m
                        rhs_vf(eqn_idx%alf)%sf(i, q, l) = rhs_vf(eqn_idx%alf)%sf(i, q, l) + bub_adv_src(i, q, l)
                        if (bubble_birth) rhs_vf(eqn_idx%n)%sf(i, q, l) = rhs_vf(eqn_idx%n)%sf(i, q, l) + bub_n_src(i, q, l)
                        if (num_fluids > 1) rhs_vf(eqn_idx%adv%beg)%sf(i, q, l) = rhs_vf(eqn_idx%adv%beg)%sf(i, q, &
                            & l) - bub_adv_src(i, q, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do k = 1, nb
                            rhs_vf(rs(k))%sf(i, q, l) = rhs_vf(rs(k))%sf(i, q, l) + bub_r_src(i, q, l, k)
                            rhs_vf(vs(k))%sf(i, q, l) = rhs_vf(vs(k))%sf(i, q, l) + bub_v_src(i, q, l, k)
                            if (polytropic .neqv. .true.) then
                                rhs_vf(ps(k))%sf(i, q, l) = rhs_vf(ps(k))%sf(i, q, l) + bub_p_src(i, q, l, k)
                                rhs_vf(ms(k))%sf(i, q, l) = rhs_vf(ms(k))%sf(i, q, l) + bub_m_src(i, q, l, k)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubble_EE_source

end module m_bubbles_EE
