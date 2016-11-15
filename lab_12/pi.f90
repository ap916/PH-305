! Finding the value of pi using monte carlo simulation 
! and Linear Congruential random number generator

program main

    integer::iter,countCircle,a,b,c,x,y,sqc,temp
    iter = 10000
    countCircle = 0
    ! Seeding the values ( Usually Prime )
    a = 3361
    b = 5023
    c = 16033
    x = 9999
    y = 1123
    sqc = c*c

    open(unit = 10,file = "output.txt")
    do i=1,iter

        ! Checking whether the point is inside circle or not
        if(x*x + y*y .lt. c*c) then
            countCircle = countCircle + 1
        end if

        ! Writing to file
        write(10,*)i,4.0*countCircle/i

        ! Updating the values ( Finding the next random number )
        x = mod((a*x + b),c)
        y = mod((a*y + b), c)

    end do

    print*,"Value of Pi :", 4.0*countCircle/iter

    close(10)

end program main
