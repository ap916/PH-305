! Finding Smallest Eigenvalue of a matrix using Inverse Power Method 
! Smallest i.e. smallest(abs(lamda(i)))

program powerinv
    implicit none
! Initialisation of variables
    real*8,parameter::maxerror = 0.00000001
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff,InverseCoeff
    real*8,dimension(:),allocatable::Xvalues,Yvalues
    integer*8::i,j,k,iter
    real*8::temp,error,Smallest

    print*,"Enter the no of variables (size of coeff matrix) followed by Coeff:"
    read*,len

    allocate(Coeff(len,len))
    allocate(InverseCoeff(len,len))
    allocate(Xvalues(len))
    allocate(Yvalues(len))

! Reading values from file "input.txt"
    do i=1,len
        read (*,*)(Coeff(i,j),j=1,len)
    end do


! Printing the matrix
    print *, "The matrix A  :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""

    do i=1,len
        Xvalues(i) = 1.0
        Yvalues(i) = 1.0
    end do

! Calculating the inverse of Coeff matrix
        call Inverse(Coeff,InverseCoeff,len)
        
    do i=1,len
        write (*,*)(InverseCoeff(i,j)*Coeff(i,j),j=1,len)
    end do
    print *, ""

! Power method
    iter = 0
    Smallest = -9e10
    do
        iter = iter + 1
        ! Printing values of X
        print *, "Eigenvector after iteration :", iter
        print *, "-------------------------------------"
        do i=1,len
            write (*,*)Xvalues(i)
        end do
        print *, "Smallest Eigenvalue :", 1.0/Smallest
        print *
        
        Yvalues = matmul(InverseCoeff,Xvalues)
        
        Smallest = Yvalues(1)
        do i=1,len
            if (abs(Yvalues(i)) > abs(Smallest)) Smallest = Yvalues(i)
        end do
        Yvalues = Yvalues/Smallest
        
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
    print *, "Smallest Eigenvalue:", 1.0/Smallest

end program powerinv


subroutine inverse(a,c,n)
implicit none 
integer*8 n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer*8 i, j, k

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
