! Solving differential equation ( Initial Value Problem ) using Heun's method
! Heun's method uses average slope of current point and next point
! Error : O(h^3)
program main
    implicit none

    ! Calling the function for estimation of the values
    call heun_euler(0.0,4.0,0.5,1.0)

end program main

! Subroutine for estimation of function using euler and heun's method
subroutine heun_euler(xmin,xmax,stepsize,yinit)
    implicit none
    integer::i,j,iter
    real::xmin,xmax,yinit,stepsize,diff,func
    real,dimension(:),allocatable::x,y,yimp
    character(20)::filename

    iter = (xmax-xmin)/stepsize

    allocate(x(iter+1))
    allocate(y(iter+1))
    allocate(yimp(iter+1))

    x(1) = xmin
    do i=2,iter+1
        x(i) = x(i-1) + stepsize
    end do

    ! Euler method 
    y(1) = yinit
    do i=2,iter+1
        y(i) = y(i-1) + diff(x(i-1))*stepsize
    end do

    ! Improved euler method (Heun's)
    yimp(1) = yinit
    do i=2,iter+1
        yimp(i) = yimp(i-1) + (diff(x(i-1)) + diff(x(i)))*stepsize/2.0
    end do

    write(filename,1)stepsize
    1 format(f4.2,".txt")
    open(unit = 10, file = filename)
    do i=1,iter+1
        write(10,*)x(i),func(x(i)),yimp(i),100*(yimp(i)-func(x(i)))/func(x(i)), y(i), 100*(y(i)-func(x(i)))/func(x(i))
    end do

end subroutine heun_euler

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
