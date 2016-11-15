! Solving Differential Equation ( Boundary Value Problem ) using Forward Difference Method

program main
    implicit none

    integer :: i, j, steps
    real :: xmin,xmax,ta,dx,h
    real,dimension(:),allocatable :: t, b
    real,dimension(:,:),allocatable :: coeff,invcoeff

    xmin = 0.0
    xmax = 10.0
    ta = 20.0
    h = 0.01

    print*, "Enter the value of dx :"
    read*,dx
    steps = ( xmax-xmin )/dx - 1

    allocate(t(steps))
    allocate(b(steps))
    allocate(coeff(steps,steps))
    allocate(invcoeff(steps,steps))

    ! Making everything 0 for avoiding funny values
    do i=1,steps
        do j=1,steps
            coeff(i,j) = 0.0
            invcoeff(i,j) = 0.0
        end do 
        t(i) = 0.0
        b(i) = 0.0
    end do

    ! Populating the Constants matrix
    do i = 1,steps
        b(i) = h*dx*dx*ta
    end do
    b(1) = b(1) + 40.0
    b(steps) = b(steps) + 200.0
    ! Calculating the Coeffecients matrix 
    do i=1,steps
        coeff(i,i) = 2.0 + h*dx*dx
        if(i .ne. 1) then
            coeff(i,i-1) = -1.0
        end if
        if (i .ne. steps) then
            coeff(i,i+1) = -1.0
        endif

    end do

    print*,"The tridiagonal matrix :"
    do i=1,steps
        write(*,*),(coeff(i,j),j=1,steps)
    end do
    print*, "The constants column vector :"
    do i=1,steps
        write(*,*) b(i)
    end do

    ! Finding the inverse of the matrix
    call inverse(coeff,invcoeff,steps)

    t = matmul(invcoeff,b)

    print*, "The solutions :"
    do i=1,steps
        write(*,*) t(i)
    end do
    
end program main


subroutine inverse(a,c,n)
    implicit none 
    integer n
    real a(n,n), c(n,n)
    real L(n,n), U(n,n), b(n), d(n), x(n)
    real coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
end subroutine inverse
