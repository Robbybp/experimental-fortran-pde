        program main
        implicit none

        real(8) :: k1 = 2, k2 = 1, H1 = 10, H2 = 1, H3 = 3
        real(8) :: u_bc = 1, v_bc = 1, w_bc = 1
        real(8) :: u_ic = 0, v_ic = 0, w_ic = 1
        real(8) :: ht, hx
        real(8) :: nz_tol, tol_rel, tol_abs, resid_norm, b_norm
        
        real(8), allocatable :: u(:,:), v(:,:), w(:,:)
        real(8), allocatable :: rhs_matrix(:,:), rhs_vector(:)
        real(8), allocatable :: A_matrix(:,:), b_vector(:),&
                                x_estimate(:), ic_vector(:),&
                                identity(:,:), resid(:)
        integer(4) :: nt = 100, nx = 100, nvar_x, dimen
        ! nt, nx are the number of elements in the discretization
        ! there are +1 points for each variable, one of which
        ! is specified by a boundary condition

        integer(4) :: i, j, it, eq_index
        integer(4) :: itr_max, mr
        integer(4) :: nnz

        integer(4), allocatable :: ia(:), ja(:)
        real(8), allocatable :: nza(:)

        character(100) :: argument, filename, which_problem

        if (command_argument_count() .gt. 0) then
          call get_command_argument(1, argument)
          read (argument,*) which_problem
        endif  
        if (command_argument_count() .gt. 1) then
          call get_command_argument(2, argument)
          read (argument,*) nt
        endif
        if (command_argument_count() .gt. 2) then
          call get_command_argument(3, argument)
          read (argument,*) nx
        endif

        ht = 1/real(nt)
        hx = 1/real(nx)

        nvar_x = 3*nx
        allocate (rhs_matrix(nvar_x, nvar_x))
        allocate (rhs_vector(nvar_x))
        allocate (identity(nvar_x, nvar_x))

        do i=1,nvar_x
          identity(i,i) = 1
        enddo

        rhs_matrix = 0
        rhs_vector = 0
        ! initialize system at lower boundary (v only)
        rhs_matrix(1,1) = -1/hx - k2
        rhs_matrix(1,3) = 1/hx
        rhs_vector(1) = -k1*u_bc
        ! initialize system at upper boundary (u, w only)
        rhs_matrix(nvar_x-1,nvar_x-1) = 1/hx - k1
        rhs_matrix(nvar_x-1,nvar_x-4) = -1/hx
        rhs_vector(nvar_x-1) = -k2*v_bc
        rhs_matrix(nvar_x,nvar_x) = 1/hx - H3
        rhs_matrix(nvar_x,nvar_x-3) = -1/hx
        rhs_matrix(nvar_x,nvar_x-1) = k1*H1
        rhs_vector(nvar_x) = k2*v_bc*H2

        ! 
        rhs_matrix(2,2) = -k1
        rhs_matrix(2,3) = k2
        rhs_matrix(2,5) = 1/(2*hx)
        rhs_vector(2) = u_bc/(2*hx)

        rhs_matrix(3,3) = -k2
        rhs_matrix(3,2) = k1
        rhs_matrix(3,1) = -1/(2*hx)
        rhs_matrix(3,6) = 1/(2*hx)

        rhs_matrix(4,4) = -H3
        rhs_matrix(4,2) = k1*H1
        rhs_matrix(4,3) = -k2*H2
        rhs_matrix(4,7) = 1/(2*hx)
        rhs_vector(4) = w_bc/(2*hx)

        rhs_matrix(nvar_x-4,nvar_x-4) = -k1
        rhs_matrix(nvar_x-4,nvar_x-3) = k2
        rhs_matrix(nvar_x-4,nvar_x-7) = -1/(2*hx)
        rhs_matrix(nvar_x-4,nvar_x-1) = 1/(2*hx)

        rhs_matrix(nvar_x-3,nvar_x-3) = -k2
        rhs_matrix(nvar_x-3,nvar_x-4) = k1
        rhs_matrix(nvar_x-3,nvar_x-6) = -1/(2*hx)
        rhs_vector(nvar_x-3) = -v_bc/(2*hx)

        rhs_matrix(nvar_x-2,nvar_x-2) = -H3
        rhs_matrix(nvar_x-2,nvar_x-4) = k1*H1
        rhs_matrix(nvar_x-2,nvar_x-3) = k2*H2
        rhs_matrix(nvar_x-2,nvar_x) = 1/(2*hx)
        rhs_matrix(nvar_x-2,nvar_x-5) = -1/(2*hx)

        ! fill in rest (interior elements) of matrix
        do i=2,nx-2
          ! u
          eq_index = 1 + 3*i - 2
          rhs_matrix(eq_index, eq_index) = -k1
          rhs_matrix(eq_index, eq_index+1) = k2
          rhs_matrix(eq_index, eq_index+3) = 1/(2*hx)
          rhs_matrix(eq_index, eq_index-3) = -1/(2*hx)

          ! v
          eq_index = 1 + 3*i - 1
          rhs_matrix(eq_index, eq_index) = -k2
          rhs_matrix(eq_index, eq_index-1) = k1
          rhs_matrix(eq_index, eq_index+3) = 1/(2*hx)
          rhs_matrix(eq_index, eq_index-3) = -1/(2*hx)

          ! w
          eq_index = 1 + 3*i
          rhs_matrix(eq_index, eq_index) = -H3
          rhs_matrix(eq_index, eq_index-2) = k1*H1
          rhs_matrix(eq_index, eq_index-1) = -k2*H2
          rhs_matrix(eq_index, eq_index+3) = 1/(2*hx)
          rhs_matrix(eq_index, eq_index-3) = -1/(2*hx)
        enddo

        open(unit=20, file='rhs_matrix.txt',&
            status='replace', action='write')
        open(unit=21, file='rhs_vector.txt',&
            status='replace', action='write')
        write(20,*) rhs_matrix
        write(21,*) rhs_vector
        close(20)
        close(21)

        nz_tol = 1.0D-8
        tol_abs = 1.0D-8
        tol_rel = 1.0D-8
        itr_max = 20
        mr = 15 

        ! construct and solve steady state problem

        allocate(A_matrix(nvar_x, nvar_x))
        allocate(b_vector(nvar_x))
        allocate(x_estimate(nvar_x))
        allocate(resid(nvar_x))

        A_matrix = rhs_matrix
        b_vector = rhs_vector

        call get_nnz(A_matrix, nvar_x, nvar_x, nz_tol, nnz)

        allocate(ia(nvar_x+1))
        allocate(ja(nnz))
        allocate(nza(nnz))
        call dense2csr(A_matrix, nvar_x, nvar_x, ia, ja, nza, nz_tol)

        write(*,'(a)') 'Linear, steady state call:'
        call pmgmres_ilu_cr(nvar_x, nnz, ia, ja, nza, x_estimate, &
                            b_vector, itr_max, mr, tol_abs, tol_rel)

        if (trim(which_problem).eq.'steady') then
          stop
        endif 

        ! construct and solve dynamic problem

        allocate(ic_vector(nvar_x))
        
        ic_vector(1) = 0
        ic_vector(nvar_x) = 1
        ic_vector(nvar_x-1) = 0
        do i=1,nvar_x-2
          if (modulo(i,3) .eq. 2) then
            ic_vector(i) = 0
          endif
          if (modulo(i,3) .eq. 0) then
            ic_vector(i) = 0
          endif
          if (modulo(i,3) .eq. 1) then
            ic_vector(i) = 1
          endif
        enddo


        A_matrix = identity - ht*rhs_matrix
        
        do i=1,nt
          b_vector = ic_vector + rhs_vector

          call dense2csr(A_matrix, nvar_x, nvar_x, ia, ja, nza, nz_tol)

          write(*,*) ' '
          write(*,'(a,I0)') 'Linear, dynamic call for t', i
          write(*,*) ' '

          call pmgmres_ilu_cr(nvar_x, nnz, ia, ja, nza, x_estimate, &
                              b_vector, itr_max, mr, tol_abs, tol_rel)

          resid = matmul(A_matrix,x_estimate) - b_vector
          resid_norm = sqrt(dot_product(resid,resid))
          b_norm = sqrt(dot_product(b_vector,b_vector))

          if (resid_norm/b_norm.gt.tol_rel) then
            write(*,*) ' '
            write(*,'(a,I0)') 'Failed to solve dynamic linear problem &
                              within tolerance at step ', i
            write(*,'(a,E12.4)') 'Residual = ', resid_norm
            write(*,'(a)') 'Exiting integration loop'
            exit
          endif
          ic_vector = x_estimate
        enddo

        deallocate(identity)

        deallocate(resid)
        deallocate(ic_vector)
        deallocate(x_estimate)
        deallocate(b_vector)
        deallocate(A_matrix)

        deallocate (rhs_vector)
        deallocate (rhs_matrix)
        contains

        subroutine get_nnz(A,n,m,nz_tol,nnz)
          ! Get number of nonzeros of matrix A        
          integer(4), intent(in) :: n, m
          real(8), intent(in) :: A(:,:)
          real(8), intent(in) :: nz_tol
          integer(4), intent(out) :: nnz
          integer(4) :: ntmp

          nnz = 0
          do i=1,n
          do j=1,m
            if (abs(A(i,j)).gt.nz_tol) then
              ntmp = nnz
              nnz = ntmp + 1
            endif
          enddo
          enddo
        end subroutine

        subroutine dense2csr(A,n,m,ia,ja,nza,nz_tol)
          ! Convert dense-represented matrix A to a csr representation 
          integer(4), intent(in) :: n, m
          real(8), intent(in) :: A(:,:)
          real(8), intent(in) :: nz_tol
          integer(4), intent(out) :: ia(:), ja(:)
          real(8), intent(out) :: nza(:)
          integer(4) :: i, j, nnz, ntmp

          nnz = 0
          ia(1) = 1
          do i=1,n
          do j=1,m
            if (abs(A(i,j)).gt.nz_tol) then
              ntmp = nnz
              nnz = ntmp + 1
              ja(nnz) = j
              nza(nnz) = A(i,j)
            endif
          enddo
            ia(i+1) = nnz + 1
          enddo
        end subroutine

        end program
