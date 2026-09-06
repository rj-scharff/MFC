!>
!! @file
!! @brief Contains module m_bubbles_validity

!> @brief Reports how close the dispersed-phase closure is to its validity limits.
!!
!! This is a read-only observer. It stores no field, changes no solver state, and
!! stops nothing: it reports numbers so that a validity crossing is a
!! deterministic, recorded event rather than something inferred afterwards from a
!! failure.
!!
!! Three limits are watched, each for a reason that was measured rather than assumed:
!!
!!   1. **Radius collapse.** Rayleigh's solution for a nearly empty cavity gives a
!!      wall velocity growing as (R0/R)**(3/2), so a cavity driven below roughly a
!!      tenth of its reference radius reaches the liquid sound speed on its own.
!!      That is a finite-time singularity, and no integrator resolves one.
!!
!!   2. **Wall Mach number.** Keller-Miksis is a first-order expansion in Rdot/c.
!!      Its retardation factor 1 + Rdot/c vanishes at Rdot = -c and inverts beyond
!!      it, so the pressure difference that should decelerate a collapse
!!      accelerates it instead.
!!
!!   3. **Void fraction.** The ensemble closure is dilute. An independent 0-D
!!      reference measured its energy accounting departing from first-order
!!      behaviour in the void fraction at around one per cent.
!!
!! Both continuous extrema and counts at reference thresholds are written. The
!! extrema are the result; the thresholds only make the counts readable, and they
!! are recorded in the header so the log can be re-thresholded without rerunning.
#:include 'macros.fpp'

module m_bubbles_validity

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_variables_conversion
    use m_constants, only: sgm_eps

    implicit none

    private; public :: s_initialize_bubbles_validity_module, s_write_bubbles_validity, s_finalize_bubbles_validity_module

    integer, parameter :: validity_unit = 18  !< File unit for the validity log

    !> Reference thresholds. Nothing computes from these: they count cells for readability while the extrema carry the result.
    real(wp), parameter :: radius_ratio_floor = 1.e-1_wp  !< Below this, Rayleigh collapse is transonic
    real(wp), parameter :: wall_mach_ceiling = 5.e-1_wp  !< Half the singular wall velocity, where the expansion is visibly wrong
    real(wp), parameter :: dilute_void_ceiling = 1.e-2_wp  !< Where the 0-D ledger's closure deficit leaves first order
    real(wp), parameter :: cell_model_void_ceiling = 3.e-1_wp  !< Where sphericity and coalescence end the cell model
    real(wp), parameter :: eos_margin_floor = 1.e-1_wp  !< A tenth of the way from ambient to the stiffened-gas singularity

    !> The void-fraction limit is a property of the closure in use, not of the problem, so it is set at initialization rather than
    !! fixed here: the dilute derivation expires near 1e-2, while the cell-model correction carries to roughly 3e-1.
    real(wp) :: void_fraction_ceiling

    !> Whether each limit was ever crossed, and where it first happened. A run that touches a validity limit must say so in its own
    !! output rather than relying on somebody reading the log: a result that carries its own caveat cannot be quoted by accident.
    !! Order: collapse, wall Mach, void fraction, equation-of-state margin.
    logical  :: ever_crossed(4)
    integer  :: first_step(4)
    real(wp) :: first_time(4)
    real(wp) :: worst(4)

contains

    !> Open the log with a labelled header.
    impure subroutine s_initialize_bubbles_validity_module

        character(len=path_len + name_len) :: file_path

        ever_crossed = .false.
        first_step = -1
        first_time = 0._wp
        worst = [huge(1._wp), 0._wp, 0._wp, huge(1._wp)]

        ! The void-fraction limit belongs to the closure that is running. The
        ! dilute derivation expires near 1e-2; with the cell-model correction the
        ! binding limits become sphericity and coalescence, which is roughly 3e-1.
        void_fraction_ceiling = dilute_void_ceiling
        if (bubble_confinement) void_fraction_ceiling = cell_model_void_ceiling

        if (proc_rank == 0) then
            file_path = trim(case_dir) // '/bubbles_validity.dat'
            open (validity_unit, file=trim(file_path), form='formatted', status='replace')

            write (validity_unit, '(A)') '# Dispersed-phase closure validity, one row per time step.'
            write (validity_unit, '(A)') '# Read-only diagnostic: no solver state is changed and nothing is stopped.'
            write (validity_unit, '(A)') '# min_radius_ratio  smallest R/R0; a finite-time collapse drives this to zero.'
            write (validity_unit, '(A)') '# max_wall_mach     largest |Rdot| over the singular wall velocity, which is the'
            write (validity_unit, '(A)') '#                   sound speed alone only when the confinement correction is off.'
            write (validity_unit, '(A)') '# max_void_fraction largest alpha; the ceiling depends on the closure in use.'
            write (validity_unit, '(A)') '# min_eos_margin    smallest (p + pi_inf)/pi_inf; one at ambient, zero where the'
            write (validity_unit, '(A)') '#                   stiffened gas admits no liquid at all. For water the'
            write (validity_unit, '(A)') '#                   homogeneous cavitation limit sits near 0.6, so a run below'
            write (validity_unit, '(A)') '#                   that sustains a tension real water could not.'
            write (validity_unit, '(A)') '# The n_* columns count cells past a reference threshold, for readability only.'
            write (validity_unit, '(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') '# Thresholds: radius ratio ', radius_ratio_floor, &
                   & ', wall Mach ', wall_mach_ceiling, ', void fraction ', void_fraction_ceiling, ', EOS margin ', eos_margin_floor
            write (validity_unit, &
                   & '(A)') '# t_step time min_radius_ratio max_wall_mach max_void_fraction min_eos_margin ' &
                   & // 'n_collapsing n_transonic n_dense n_strained'
        end if

    end subroutine s_initialize_bubbles_validity_module

    !> Scan the interior and write one row.
    !! @param q_prim_vf Primitive variables
    !! @param q_cons_vf Conservative variables
    !! @param t_step Current time step
    impure subroutine s_write_bubbles_validity(q_prim_vf, q_cons_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf, q_cons_vf
        integer, intent(in)                                 :: t_step
        real(wp)                                            :: ratio_min_loc, mach_max_loc, void_max_loc, margin_min_loc
        real(wp)                                            :: ratio_min_glb, mach_max_glb, void_max_glb, margin_min_glb
        real(wp)                                            :: collapsing_loc, transonic_loc, dense_loc, strained_loc
        real(wp)                                            :: collapsing_glb, transonic_glb, dense_glb, strained_glb
        real(wp)                                            :: myR, myV, alf, myP, myRho, n_tait, B_tait, c_liquid, ratio, mach
        real(wp)                                            :: beta, singular, margin
        integer                                             :: j, k, l, q, ii

        ratio_min_loc = huge(1._wp)
        mach_max_loc = 0._wp
        void_max_loc = 0._wp
        margin_min_loc = huge(1._wp)
        collapsing_loc = 0._wp
        transonic_loc = 0._wp
        dense_loc = 0._wp
        strained_loc = 0._wp

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, q, ii, myR, myV, alf, myP, myRho, n_tait, B_tait, c_liquid, ratio, &
                            & mach, beta, singular, margin]', reduction='[[mach_max_loc, void_max_loc], [ratio_min_loc, &
                            & margin_min_loc], [collapsing_loc, transonic_loc, dense_loc, strained_loc]]', reductionOp='[max, &
                            & min, +]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)

                    ! Cells the solver itself skips carry stale state, so reporting
                    ! them would describe the gate rather than the physics.
                    if (alf >= small_alf) then
                        void_max_loc = max(void_max_loc, alf)
                        if (alf > void_fraction_ceiling) dense_loc = dense_loc + 1._wp

                        ! The mixture sound speed the wall equation itself uses, so the
                        ! reported Mach number is the one the model is expanded in.
                        myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                        if (num_fluids == 1) then
                            myRho = q_cons_vf(1)%sf(j, k, l)
                            n_tait = gammas(1)
                            B_tait = pi_infs(1)/pi_fac
                        else
                            myRho = 0._wp
                            n_tait = 0._wp
                            B_tait = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do ii = 1, num_fluids
                                myRho = myRho + q_cons_vf(ii)%sf(j, k, l)
                                n_tait = n_tait + q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)*gammas(ii)
                                B_tait = B_tait + q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)*pi_infs(ii)/pi_fac
                            end do
                        end if
                        n_tait = 1._wp/n_tait + 1._wp
                        B_tait = B_tait*(n_tait - 1._wp)/n_tait
                        c_liquid = sqrt(n_tait*max((myP + B_tait)/(myRho*(1._wp - alf)), sgm_eps))

                        ! How much tension the stiffened gas has left. It admits no
                        ! liquid below -B_tait, where the sound speed vanishes, so the
                        ! margin runs from one at ambient to zero at that wall. For
                        ! water the homogeneous cavitation limit near -120 MPa sits at
                        ! a margin of about 0.6, so a run below that is sustaining a
                        ! tension real water could not.
                        margin = (myP + B_tait)/B_tait
                        margin_min_loc = min(margin_min_loc, margin)
                        if (margin < eos_margin_floor) strained_loc = strained_loc + 1._wp

                        ! Confinement moves the wall equation's singularity from the
                        ! sound speed to (1 - beta) times it, so the Mach number is
                        ! reported against the velocity that is actually singular.
                        ! Otherwise the reported margin flatters the model exactly
                        ! where the void fraction is largest.
                        beta = 0._wp
                        if (bubble_confinement) beta = max(alf, 0._wp)**(1._wp/3._wp)
                        singular = max(1._wp - beta, sgm_eps)*c_liquid

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            myR = q_prim_vf(qbmm_idx%rs(q))%sf(j, k, l)
                            myV = q_prim_vf(qbmm_idx%vs(q))%sf(j, k, l)
                            ratio = myR/R0(q)
                            mach = abs(myV)/singular

                            ratio_min_loc = min(ratio_min_loc, ratio)
                            mach_max_loc = max(mach_max_loc, mach)
                            if (ratio < radius_ratio_floor) collapsing_loc = collapsing_loc + 1._wp
                            if (mach > wall_mach_ceiling) transonic_loc = transonic_loc + 1._wp
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (num_procs > 1) then
            call s_mpi_allreduce_min(ratio_min_loc, ratio_min_glb)
            call s_mpi_allreduce_max(mach_max_loc, mach_max_glb)
            call s_mpi_allreduce_max(void_max_loc, void_max_glb)
            call s_mpi_allreduce_min(margin_min_loc, margin_min_glb)
            call s_mpi_allreduce_sum(collapsing_loc, collapsing_glb)
            call s_mpi_allreduce_sum(transonic_loc, transonic_glb)
            call s_mpi_allreduce_sum(dense_loc, dense_glb)
            call s_mpi_allreduce_sum(strained_loc, strained_glb)
        else
            ratio_min_glb = ratio_min_loc
            mach_max_glb = mach_max_loc
            void_max_glb = void_max_loc
            margin_min_glb = margin_min_loc
            collapsing_glb = collapsing_loc
            transonic_glb = transonic_loc
            dense_glb = dense_loc
            strained_glb = strained_loc
        end if

        worst(1) = min(worst(1), ratio_min_glb)
        worst(2) = max(worst(2), mach_max_glb)
        worst(3) = max(worst(3), void_max_glb)
        worst(4) = min(worst(4), margin_min_glb)
        call s_note_crossing(1, ratio_min_glb < radius_ratio_floor, t_step)
        call s_note_crossing(2, mach_max_glb > wall_mach_ceiling, t_step)
        call s_note_crossing(3, void_max_glb > void_fraction_ceiling, t_step)
        call s_note_crossing(4, margin_min_glb < eos_margin_floor, t_step)

        if (proc_rank == 0) then
            write (validity_unit, '(I9,1X,5(ES24.16,1X),4(I12,1X))') t_step, mytime, ratio_min_glb, mach_max_glb, void_max_glb, &
                   & margin_min_glb, nint(collapsing_glb), nint(transonic_glb), nint(dense_glb), nint(strained_glb)
            flush (validity_unit)
        end if

    end subroutine s_write_bubbles_validity

    !> Record the first time a limit is crossed, and nothing thereafter.
    !! @param which Which limit: 1 collapse, 2 wall Mach, 3 void fraction
    !! @param crossed Whether the limit is crossed on this step
    !! @param t_step Current time step
    impure subroutine s_note_crossing(which, crossed, t_step)

        integer, intent(in) :: which
        logical, intent(in) :: crossed
        integer, intent(in) :: t_step

        if (crossed .and. (.not. ever_crossed(which))) then
            ever_crossed(which) = .true.
            first_step(which) = t_step
            first_time(which) = mytime
        end if

    end subroutine s_note_crossing

    !> Close the log and write the summary that marks the run.
    impure subroutine s_finalize_bubbles_validity_module

        character(len=15), parameter :: label(4) = [character(len=15)::'radius_collapse','wall_mach      ', 'void_fraction  ', &
                  & 'eos_margin     ']
        real(wp)                           :: threshold(4)
        character(len=path_len + name_len) :: summary_path
        integer                            :: i

        threshold = [radius_ratio_floor, wall_mach_ceiling, void_fraction_ceiling, eos_margin_floor]

        if (proc_rank == 0) then
            close (validity_unit)

            summary_path = trim(case_dir) // '/bubbles_validity_summary.dat'
            open (validity_unit, file=trim(summary_path), form='formatted', status='replace')
            write (validity_unit, '(A)') '# Did this run leave the dispersed-phase closure it is derived under?'
            write (validity_unit, '(A)') '# Written so a result carries its own caveat rather than relying on the log being read.'
            write (validity_unit, '(A)') '# limit  crossed  first_step  first_time  worst_value  threshold'
            do i = 1, 3
                write (validity_unit, '(A,2X,L1,2X,I9,2X,3(ES24.16,2X))') label(i), ever_crossed(i), first_step(i), &
                       & first_time(i), worst(i), threshold(i)
            end do
            if (any(ever_crossed)) then
                write (validity_unit, '(A)') 'VERDICT OUTSIDE_VALIDITY'
                write (validity_unit, '(A)') '# At least one limit was crossed. Results from this run are not'
                write (validity_unit, '(A)') '# predictions of the model and must not be quoted as such.'
                print '(A)', ' WARNING: the dispersed-phase closure left its validity limits during this run.'
                print '(A)', '          See bubbles_validity_summary.dat.'
            else
                write (validity_unit, '(A)') 'VERDICT WITHIN_VALIDITY'
            end if
            close (validity_unit)
        end if

    end subroutine s_finalize_bubbles_validity_module

end module m_bubbles_validity
