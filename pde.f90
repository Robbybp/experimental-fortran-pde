      module util
      implicit none
      external pmgmres_ilu_cr

      contains

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

        subroutine get_nnz(A,n,m,nz_tol,nnz)
          ! Get number of nonzeros of dense-represented matrix A
          integer(4), intent(in) :: n, m
          real(8), intent(in) :: A(:,:)
          real(8), intent(in) :: nz_tol
          integer(4), intent(out) :: nnz
          integer(4) :: ntmp, i, j

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

        subroutine newton(x0, x, newton_tol, newton_iter,&
                          nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)
        real(8), intent(in) :: hx, k1, k2, &
                               H1, H2, H3, u_bc, v_bc, w_bc
        integer(4), intent(in) :: nx
        real(8), intent(in) :: x0(:), newton_tol
        integer(4), intent(in) :: newton_iter
        real(8), intent(out) :: x(:)
        real(8), dimension(size(x0), size(x0)) :: Jf
        real(8), dimension(size(x0)) :: f, x_new, b_vec
        real(8), dimension(size(x0)) :: search_d
        real(8) :: norm_f, norm_fprev
        real(8) :: alf
        integer(4) :: k, i, n, nnz
        real(8), allocatable :: nza(:)
        integer(4), allocatable :: ia(:), ja(:)
        write(*,'(a)') "Beginning Newton's method"
        n=size(x0)
        alf = 1
        x = x0

        call get_Jfx1(x, f, Jf,&
                      nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)
        call get_nnz(Jf,n,n,1.0D-8,nnz)
        allocate(nza(nnz))
        allocate(ja(nnz))
        allocate(ia(n+1))
        call dense2csr(Jf,n,n,ia,ja,nza,1.0D-8)
                      
        norm_fprev = sqrt(dot_product(f,f))
        do k=1,newton_iter
          write(*,'(a,i0)') 'Entering Newton iteration ', k
          b_vec = -f
          call pmgmres_ilu_cr(n, nnz, ia, ja, nza, search_d, &
                             b_vec, 20, 20, 1.0D-8, 1.0D-8)
          x_new = x + search_d
          call get_Jfx1(x_new, f, Jf,&
                        nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)
          call get_nnz(Jf,n,n,1.0D-8,nnz)
          call dense2csr(Jf,n,n,ia,ja,nza,1.0D-8)
          !do i=1,n
          !write(*,*) x(i)
          !enddo
          norm_f = sqrt(dot_product(f,f))
          write(*,'(a,E12.4)') 'Initial ||f(x)|| = ', norm_f
          alf = 1
          do i=1,20
          if (norm_f.ge.norm_fprev) then
            alf = alf/2
            x_new = x + search_d/alf
            call get_Jfx1(x_new, f, Jf,&
                        nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)
            norm_f = sqrt(dot_product(f,f))
            write(*,'(a,E12.4)') 'Cut alpha, ||f(x)|| = ', norm_f

            cycle
          else
            write(*,'(a,E8.1)') 'Accepted iterate with alpha = ', alf
            exit
          endif
          enddo

          if (norm_f.ge.norm_fprev) then
            write(*,'(a,i0)') 'Error: Line search failed at iterate ',k
            return
          endif
          norm_fprev = norm_f
          write(*,'(a,E12.4)') "||f(x)|| = ", norm_f

          ! calculate search direction
          ! update x
          ! recalculate Jfx
          ! convert to CSR
          if (norm_f.lt.newton_tol) then
             write(*,'(a)') "Newton's method terminated successfully"
             write(*,'(a,E12.4)') "||f(x)|| = ", norm_f
             return
           endif
        enddo
        write(*,'(a)') "Newton's method reached max iteration"
        write(*,'(a,E12.4)') "||f(x)|| = ", norm_f
        

        end subroutine newton

        subroutine get_Jfx1(x, f, Jf,&
                           nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)
        real(8), intent(in) :: x(:)
        real(8), intent(out), dimension(size(x),size(x)) :: Jf
        real(8), intent(out), dimension(size(x)) :: f
        real(8), intent(in) :: hx, k1, k2, &
                               H1, H2, H3, u_bc, v_bc, w_bc
        integer(4), intent(in) :: nx
        integer(4) :: i, eq, n
        n = size(x)
        
        f = 0
        Jf = 0

        f(1) = 1/hx*(x(3)-x(1)) - k2*x(1)**1 + k1*u_bc**2
        f(2) = 1/(2*hx)*(x(5)-u_bc) - k1*x(2)**2 + k2*x(3)**1
        f(3) = 1/(2*hx)*(x(6)-x(1)) + k1*x(2)**2 - k2*x(3)**1
        f(4) = 1/(2*hx)*(x(7)-w_bc) + H1*k1*x(2)**2 - &
                                      H2*k2*x(3)**1 - H3*x(4)
        do i=2,nx
          eq = 3*i - 1
          f(eq) = 1/(2*hx)*(x(eq+3)-x(eq-3)) - k1*x(eq)**2 + &
                                               k2*x(eq+1)**1
          eq = 3*i
          f(eq) = 1/(2*hx)*(x(eq+3)-x(eq-3)) + k1*x(eq-1)**2 - &
                                               k2*x(eq)**1
          eq = 3*i + 1
          f(eq) = 1/(2*hx)*(x(eq+3)-x(eq-3)) + H1*k1*x(eq-2)**2 -&
                                               H2*k2*x(eq-1)**1 -&
                                               H3*x(eq)
        enddo
        f(n-4) = 1/(2*hx)*(x(n-1)-x(n-7)) - k1*x(n-4)**2 + &
                                             k2*x(n-3)**1
        f(n-3) = 1/(2*hx)*(v_bc-x(n-6)) + k1*x(n-4)**2 - &
                                             k2*x(n-3)**1
        f(n-2) = 1/(2*hx)*(x(n)-x(n-5)) + H1*k1*x(n-4)**2 -&
                                             H2*k2*x(n-3)**1 -&
                                             H3*x(n-2)
        f(n-1) = 1/(2*hx)*(x(n-1)-x(n-4)) - k1*x(n-1)**2 + &
                                             k2*v_bc**1
        f(n) = 1/(2*hx)*(x(n)-x(n-2)) + H1*k1*x(n-1)**2 -&
                                             H2*k2*v_bc**1 -&
                                             H3*x(n)

        ! initialize system at lower boundary (v only)
        Jf(1,1) = -1/hx - k2*1
        Jf(1,3) = 1/hx
        ! initialize system at upper boundary (u, w only)
        Jf(n-1,n-1) = 1/hx - k1*2*x(n-1)**1
        Jf(n-1,n-4) = -1/hx 
        Jf(n,n) = 1/hx - H3
        Jf(n,n-2) = -1/hx
        Jf(n,n-1) = 2*k1*H1*x(n-1)**1

        ! 
        Jf(2,2) = -k1*2*x(2)**1
        Jf(2,3) = k2*1
        Jf(2,5) = 1/(2*hx)

        Jf(3,3) = -k2*1
        Jf(3,2) = k1*2*x(2)**1
        Jf(3,1) = -1/(2*hx)
        Jf(3,6) = 1/(2*hx)

        Jf(4,4) = -H3
        Jf(4,2) = k1*H1*2*x(2)**1
        Jf(4,3) = -k2*H2*1
        Jf(4,7) = 1/(2*hx)

        Jf(n-4,n-4) = -k1*2*x(n-4)**1
        Jf(n-4,n-3) = k2*1
        Jf(n-4,n-7) = -1/(2*hx)
        Jf(n-4,n-1) = 1/(2*hx)

        Jf(n-3,n-3) = -k2*1
        Jf(n-3,n-4) = k1*2*x(n-4)**1
        Jf(n-3,n-6) = -1/(2*hx)

        Jf(n-2,n-2) = -H3
        Jf(n-2,n-4) = k1*H1*2*x(n-4)**1
        Jf(n-2,n-3) = -k2*H2*1
        Jf(n-2,n) = 1/(2*hx)
        Jf(n-2,n-5) = -1/(2*hx)

        ! fill in rest (interior elements) of matrix
        do i=2,nx-2
          ! u
          eq = 1 + 3*i - 2
          Jf(eq, eq) = -k1*2*x(eq)**1
          Jf(eq, eq+1) = k2*1
          Jf(eq, eq+3) = 1/(2*hx)
          Jf(eq, eq-3) = -1/(2*hx)

          ! v
          eq = 1 + 3*i - 1
          Jf(eq, eq) = -k2*1
          Jf(eq, eq-1) = k1*2*x(eq-1)**1
          Jf(eq, eq+3) = 1/(2*hx)
          Jf(eq, eq-3) = -1/(2*hx)

          ! w
          eq = 1 + 3*i
          Jf(eq, eq) = -H3
          Jf(eq, eq-2) = k1*H1*2*x(eq-2)**1
          Jf(eq, eq-1) = -k2*H2*1
          Jf(eq, eq+3) = 1/(2*hx)
          Jf(eq, eq-3) = -1/(2*hx)
        enddo

        end subroutine get_Jfx1

      end module util

      program main
        use util
        implicit none

        real(8) :: k1 = 2, k2 = 1, H1 = 10, H2 = 1, H3 = 5
        real(8) :: u_bc = 1, v_bc = 1, w_bc = 1
        real(8) :: u_ic = 0, v_ic = 0, w_ic = 1
        real(8) :: ht, hx
        real(8) :: nz_tol, tol_rel, tol_abs, resid_norm, b_norm

        real(8) :: newton_tol
        integer(4) :: newton_iter
        real(8), allocatable :: init_guess(:)

        
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

        if (which_problem.eq.'nonlinear1') then
          allocate(x_estimate(nvar_x))
          allocate(init_guess(nvar_x))
          init_guess = 1
          newton_tol = 1.0D-8
          newton_iter = 20
          call newton(init_guess, x_estimate, newton_tol, newton_iter,&
                      nx, hx, k1, k2, H1, H2, H3, u_bc, v_bc, w_bc)

          stop
        endif

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


      end program
