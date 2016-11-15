! Finding largest Eigenvalue of a matrix using Power method
! max(abs(lamda(i)))

program power
    implicit none
! Initialisation of variables
    real*8,parameter::maxerror = 0.00000001
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff
    real*8,dimension(:),allocatable::Xvalues,Yvalues
    integer*8::i,j,k,iter
    real*8::temp,error,Largest

    print*,"Enter the no of variables (size of coeff matrix) followed by Coeff:"
    read*,len

    allocate(Coeff(len,len))
    allocate(Xvalues(len))
    allocate(Yvalues(len))

! Reading values from file "input.txt"
    do i=1,len
        read (*,*)(Coeff(i,j),j=1,len)
    end do


! Printing the matrix and vector
    print *, "The matrix A and vector B :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""

    do i=1,len
        Xvalues(i) = 1.0
        Yvalues(i) = 1.0
    end do

! Power method
    iter = 0
    Largest = -9e10
    do
        iter = iter + 1
        ! Printing values of X
        print *, "Eigenvector after iteration :", iter
        print *, "-------------------------------------"
        do i=1,len
            write (*,*)Xvalues(i)
        end do
        print *, "Eigenvalue :", Largest
        print *
        
        Yvalues = matmul(Coeff,Xvalues)
        
        Largest = Yvalues(1)
        do i=1,len
            if (abs(Yvalues(i)) > abs(Largest)) Largest = Yvalues(i)
        end do
        Yvalues = Yvalues/Largest
        
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
    print *, "Final Eigenvector :"
    print *, "-----------------"
    do i=1,len
        write (*,*)Xvalues(i)
    end do
! Printing the Dominant Eigenvalue
    print *, "Dominant Eigenvalue:", Largest

end program power
