!>
!! @file
!! @brief Contains module m_conservation_monitor

!> @brief Reports domain-integrated conserved quantities each time step.
!!
!! This is a read-only observer. It stores no field, changes no solver state,
!! and classifies nothing: it reports numbers so that a later local source can
!! be judged against the unmodified solver's own behaviour.
!!
!! Each net integral is accompanied by the integral of the absolute value.
!! A small net change obtained from two large independently computed terms is
!! less trustworthy than one obtained from small terms, and the pair makes the
!! difference visible.
#:include 'macros.fpp'

module m_conservation_monitor

    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

    private; public :: s_initialize_conservation_monitor_module, s_write_conservation_monitor, &
        & s_finalize_conservation_monitor_module

    integer, parameter :: monitor_unit = 17  !< File unit for the monitor output

    !> Accumulators, laid out as [energy, per-fluid mass, per-fluid volume, momentum], with the net integrals first and the absolute
    !! integrals in the second half. Shaped (1, 2*n_quantity) because the MPI vector reduction takes rank-2 arguments.
    real(wp), allocatable, dimension(:,:) :: integral_loc, integral_glb
    $:GPU_DECLARE(create='[integral_loc, integral_glb]')
    integer :: n_quantity  !< Number of integrated quantities (half of the accumulator length)

contains

    !> Allocate the accumulators and open the output file with a labelled header.
    impure subroutine s_initialize_conservation_monitor_module

        character(len=path_len + name_len) :: file_path
        integer                            :: i

        n_quantity = 1 + 2*num_fluids + num_vels

        @:ALLOCATE(integral_loc(1:1, 1:2*n_quantity))
        @:ALLOCATE(integral_glb(1:1, 1:2*n_quantity))

        if (proc_rank == 0) then
            file_path = trim(case_dir) // '/conservation_monitor.dat'
            open (monitor_unit, file=trim(file_path), form='formatted', status='replace')

            write (monitor_unit, '(A)') '# Domain-integrated conserved quantities, one row per time step.'
            write (monitor_unit, '(A)') '# Read-only diagnostic: no solver state is changed and nothing is classified.'
            write (monitor_unit, '(A)') '# abs_* integrate the absolute value of the same integrand over the same volume.'
            if (cyl_coord) then
                write (monitor_unit, '(A)') '# Axisymmetric: the volume element carries the radial factor y_cc, without 2*pi.'
            end if

            write (monitor_unit, '(A)', advance='no') '# t_step time energy'
            do i = 1, num_fluids
                write (monitor_unit, '(A,I0)', advance='no') ' mass_', i
            end do
            do i = 1, num_fluids
                write (monitor_unit, '(A,I0)', advance='no') ' volume_', i
            end do
            do i = 1, num_vels
                write (monitor_unit, '(A,I0)', advance='no') ' momentum_', i
            end do
            write (monitor_unit, '(A)', advance='no') ' abs_energy'
            do i = 1, num_fluids
                write (monitor_unit, '(A,I0)', advance='no') ' abs_mass_', i
            end do
            do i = 1, num_fluids
                write (monitor_unit, '(A,I0)', advance='no') ' abs_volume_', i
            end do
            do i = 1, num_vels
                write (monitor_unit, '(A,I0)', advance='no') ' abs_momentum_', i
            end do
            write (monitor_unit, *)
        end if

    end subroutine s_initialize_conservation_monitor_module

    !> Integrate the conserved rows over the interior and write one row.
    !! @param q_cons_vf Conservative variables
    !! @param t_step Current time step
    impure subroutine s_write_conservation_monitor(q_cons_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: t_step
        real(wp)                                            :: cell_volume, contribution
        integer                                             :: i, j, k, l, slot

        integral_loc = 0._wp

        $:GPU_UPDATE(device='[integral_loc]')

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, slot, cell_volume, contribution]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    cell_volume = dx(j)
                    if (n > 0) cell_volume = cell_volume*dy(k)
                    if (p > 0) cell_volume = cell_volume*dz(l)
                    ! Axisymmetric and cylindrical grids weight by radius; the constant 2*pi is omitted and the header records that.
                    if (cyl_coord) cell_volume = cell_volume*y_cc(k)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, n_quantity
                        if (i == 1) then
                            slot = eqn_idx%E
                        else if (i <= 1 + num_fluids) then
                            slot = eqn_idx%cont%beg + (i - 2)
                        else if (i <= 1 + 2*num_fluids) then
                            slot = eqn_idx%adv%beg + (i - 2 - num_fluids)
                        else
                            slot = eqn_idx%mom%beg + (i - 2 - 2*num_fluids)
                        end if
                        contribution = real(q_cons_vf(slot)%sf(j, k, l), wp)*cell_volume
                        $:GPU_ATOMIC(atomic='update')
                        integral_loc(1, i) = integral_loc(1, i) + contribution
                        $:GPU_ATOMIC(atomic='update')
                        integral_loc(1, n_quantity + i) = integral_loc(1, n_quantity + i) + abs(contribution)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[integral_loc]')

        if (num_procs > 1) then
            call s_mpi_allreduce_vectors_sum(integral_loc, integral_glb, 1, 2*n_quantity)
        else
            integral_glb = integral_loc
        end if

        if (proc_rank == 0) then
            write (monitor_unit, '(I9,1X)', advance='no') t_step
            write (monitor_unit, '(ES24.16,1X)', advance='no') mytime
            do i = 1, 2*n_quantity
                write (monitor_unit, '(ES24.16,1X)', advance='no') integral_glb(1, i)
            end do
            write (monitor_unit, *)
            flush (monitor_unit)
        end if

    end subroutine s_write_conservation_monitor

    !> Close the output file and deallocate the accumulators.
    impure subroutine s_finalize_conservation_monitor_module

        if (proc_rank == 0) close (monitor_unit)

        @:DEALLOCATE(integral_loc)
        @:DEALLOCATE(integral_glb)

    end subroutine s_finalize_conservation_monitor_module

end module m_conservation_monitor
