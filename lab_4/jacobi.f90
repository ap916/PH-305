! Solving system of linear equations using Jacobi Iteration

program jacobi
    implicit none
! Initialisation of variables
    real*8,parameter::maxerror = 0.00000001
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff,Lmat
    real*8,dimension(:),allocatable::Constants,Xvalues,Yvalues
    integer*8::i,j,k,iter
    real*8::summation,lamda,temp,error

    print*,"Enter the no of variables (size of coeff matrix) followed by Coeff and Constants:"
    read*,len

    allocate(Coeff(len,len))
    allocate(Constants(len))
    allocate(Xvalues(len))
    allocate(Yvalues(len))

! Reading values from file "input.txt"
    do i=1,len
        read (*,*)(Coeff(i,j),j=1,len)
    end do
    do i=1,len
        read (*,*)Constants(i)
    end do

! Printing the matrix and vector
    print *, "The matrix A and vector B :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""
    do i=1,len
        write (*,*)Constants(i)
    end do

    do i=1,len
        Xvalues(i) = 1.0
        Yvalues(i) = 1.0
    end do

! Jacobi method
    iter = 0
    do
        iter = iter + 1
        ! Printing values of X
        print *, "Solutions after iteration :", iter
        print *, "-------------------------------------"
        do i=1,len
            write (*,*)Xvalues(i)
        end do
        
        do i=1,len
            temp = 0.0
            do j=1,len
                if(j.ne.i) then 
                    temp = temp + Coeff(i,j)*Xvalues(j)
                end if
            end do
            Yvalues(i) = (Constants(i) - temp)/Coeff(i,i)
        end do
        ! error calculation 
        error = abs(Yvalues(1) - Xvalues(1))
        do i=1,len
            temp = abs(Yvalues(i) - Xvalues(i))
            Xvalues(i) = Yvalues(i)
            if(temp .ge. error) error = temp
        end do
        if (error .le. maxerror) exit


    end do

! Printing values of X
    print *, "Final Solutions :"
    print *, "-----------------"
    do i=1,len
        write (*,*)Xvalues(i)
    end do


end program jacobi
