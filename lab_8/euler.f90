! Solving a differential equation (Initial value problem) using Euler method
! Error : O(h^2)

program main
    implicit none

    real::x,stepsize,func
    integer::i
! Plotting the original function    
    x = 0.0
    stepsize = 0.01
    open(unit = 5,file = "original.txt")
    do i=1,401
        write(5,*)x,func(x)
        x = x + stepsize
    end do
! Calling the function for estimation of the values
    call euler(0.0,4.0,0.5,1.0)
    call euler(0.0,4.0,0.25,1.0)
    call euler(0.0,4.0,0.1,1.0)
    call euler(0.0,4.0,0.05,1.0)

end program main

! Subroutine for estimation of function using euler method
subroutine euler(xmin,xmax,stepsize,yinit)
    implicit none
    integer::i,j,iter
    real::xmin,xmax,yinit,stepsize,diff
    real,dimension(:),allocatable::x,y
    character(20)::filename

    ! xmin = 0.0
    ! xmax = 4.0
    ! stepsize = 0.5
    ! yinit = 1.0
    iter = (xmax-xmin)/stepsize

    allocate(x(iter+1))
    allocate(y(iter+1))

    x(1) = xmin
    do i=2,iter+1
        x(i) = x(i-1) + stepsize
    end do

    y(1) = yinit
    do i=2,iter+1
        y(i) = y(i-1) + diff(x(i-1))*stepsize
    end do

    write(filename,1)stepsize
    1 format(f4.2,".txt")
    open(unit = 10, file = filename)
    do i=1,iter+1
        write(10,*)x(i),y(i)
    end do

end subroutine euler

! Function to calculate the differential of the function
function diff(x)
    implicit none
    real :: x,diff
    diff = -2*(x**3) + 12*(x**2) - 20*x + 8.5
    return
end function diff

! Function to calculate the function 
function func(x)
    implicit none
    real:: x,func
    func = -0.5*x**4 + 4*x**3 -10*x**2 + 8.5*x + 1
    return 
end function func
