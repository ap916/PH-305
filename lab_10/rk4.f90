! Solving a differential equation ( Initial value Problem ) using RK4 Method
! Error : O(h^4)

program main
    implicit none

    ! Calling the function for estimation of the values
    call euler_heun_rk4(0.0,4.0,0.5,4.0,6.0)

end program main

! Subroutine for estimation of function using Rk4, euler and heun's method
subroutine euler_heun_rk4(xmin,xmax,stepsize,yinit1,yinit2)
    implicit none
    integer::i,j,iter
    real::xmin,xmax,yinit1,yinit2,stepsize,diff1,diff2,func1,func2,y1,y2
    real k1,k2,k3,k4,m1,m2,m3,m4
    real,dimension(:,:),allocatable::y,yimp,yrk
    real,dimension(:),allocatable::x
    character(20)::filename

    iter = (xmax-xmin)/stepsize

    allocate(x(iter+1))
    allocate(y(2,iter+1))
    allocate(yimp(2,iter+1))
    allocate(yrk(2,iter+1))

    x(1) = xmin
    do i=2,iter+1
        x(i) = x(i-1) + stepsize
    end do

    ! Euler method 
    y(1,1) = yinit1
    y(2,1) = yinit2
    do i=2,iter+1
        y(1,i) = y(1,i-1) + diff1(y(1,i-1),y(2,i-1))*stepsize
        y(2,i) = y(2,i-1) + diff2(y(1,i-1),y(2,i-1))*stepsize
    end do

    ! Improved euler method (Heun's)
    yimp(1,1) = yinit1
    yimp(2,1) = yinit2
    do i=2,iter+1
        y1 = yimp(1,i-1) + diff1(yimp(1,i-1),yimp(2,i-1))*stepsize
        y2 = yimp(2,i-1) + diff2(yimp(1,i-1),yimp(2,i-1))*stepsize
        
        yimp(1,i) = yimp(1,i-1) + ( diff1(yimp(1,i-1),yimp(2,i-1)) + diff1(y1,y2) )*stepsize/2.0
        yimp(2,i) = yimp(2,i-1) + ( diff2(yimp(1,i-1),yimp(2,i-1)) + diff2(y1,y2) )*stepsize/2.0
    end do

    ! Rk4 method 
    yrk(1,1) = yinit1
    yrk(2,1) = yinit2
    do i = 2, iter+1
        k1 = diff1( yrk(1,i-1),yrk(2,i-1) )
        m1 = diff2( yrk(1,i-1),yrk(2,i-1) )
        k2 = diff1( yrk(1,i-1)+ stepsize*k1/2.0, yrk(2,i-1)+ stepsize*m1/2.0 )
        m2 = diff2( yrk(1,i-1)+ stepsize*k1/2.0, yrk(2,i-1)+ stepsize*m1/2.0 )
        k3 = diff1( yrk(1,i-1)+ stepsize*k2/2.0, yrk(2,i-1)+ stepsize*m2/2.0 )
        m3 = diff2( yrk(1,i-1)+ stepsize*k2/2.0, yrk(2,i-1)+ stepsize*m2/2.0 )
        k4 = diff1( yrk(1,i-1)+ stepsize*k3, yrk(2,i-1)+ stepsize*m3 )
        m4 = diff2( yrk(1,i-1)+ stepsize*k3, yrk(2,i-1)+ stepsize*m3 )

        yrk(1,i) = yrk(1,i-1) + stepsize*(k1+2*k2+2*k3+k4)/6.0
        yrk(2,i) = yrk(2,i-1) + stepsize*(m1+2*m2+2*m3+m4)/6.0
    end do

    ! Writing stuff to file with name as stepsize.txt
    write(filename,1)stepsize
    1 format(f4.2,".txt")
    open(unit = 10, file = filename)
    do i=1,iter+1
        write(10,11)x(i),yrk(1,i),yrk(2,i), yimp(1,i), yimp(2,i), y(1,i), y(2,i), func1(x(i)), func2(x(i))
        11 format(f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " ,f9.6" " )
    end do

end subroutine euler_heun_rk4

! Function to calculate the differential of the function
function diff1(y1,y2)
    implicit none
    real :: y1,y2,diff1
    diff1 = -0.5*y1
    return
end function diff1

function diff2(y1,y2)
    implicit none
    real :: y1,y2,diff2
    diff2 = 4.0 - 0.3*y2 - 0.1*y1
    return
end function diff2

! Function to calculate the function 
function func1(x)
    implicit none
    real:: x,func1
    func1 = exp(-0.5*x)*4.0
    return 
end function func1

function func2(x)
    implicit none
    real:: x,func2
    func2 = 12.0 - 2*exp(-0.5*x) - 4.0*exp(-0.3*x)
    return 
end function func2
