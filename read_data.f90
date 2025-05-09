!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module read_data

    implicit none

!-----------------------------------------------------------------------------------------------------------------!
!                                                                                                                 !
!  Contains subroutines used to load in a file containing pulsar signal data from the provided collection.        !
!  Converts data in file into a usable format, then creates an array containing columns for index, time, and      !
!  signal amplitude.                                                                                              !
!                                                                                                                 !
!  Author: Rachel Wood                                                                                            !
!                                                                                                                 !
!-----------------------------------------------------------------------------------------------------------------!

contains

    subroutine load_file(input, A, n)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !   Subroutine for use loading in .txt files from the provided collection and converting to arrays            !
    !   for use in fft calculations. Input file format is limited to that of those files which are provided.      !
    !   Provided data, for use with this program, is modified in format from the PULSE@Parkes data.               !
    !                                                                                                             !        
    !              input : input file path of desired pulsar signal.                                              !
    !                  A : output array containing columns for index, time, and amplitude for each data point.    !
    !                  n : number of data points in file.                                                         !
    !                 io : unique integer unit assigned to the current file.                                      !
    !                  i : temporary index variable, for use in do loops.                                         !
    !               stat : stores the I/O status of the previous attempt, used to identify errors in file         !
    !                       reading.                                                                              !
    !              index : name of first column in array A; index.                                                !
    !                  t : name of second column in array A; time.                                                !
    !                  x : name of third column in array A; signal amplitude.                                     !
    !               line : temporary string variable used to store each line of data read from the input file.    !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        character(len=73), intent(IN) :: input
        real, allocatable, intent(OUT) :: A(:, :)
        integer, intent(OUT) :: n
        integer :: io, i, stat
        real :: index, t, x
        character(len=132) :: line

        ! Open file specified by "input". Assign to unique integer unit io
        open(newunit=io, file=input, status="old", action="read")

        ! Begin at zero data points
        n = 0

        ! Do loop which runs through input file and assigns n as the number of lines present
        do
            read(io, '(A)', iostat=stat) line
            if (stat /= 0) exit  ! Ends run if an error occured while reading the file
            n = n + 1
        end do

        ! Rewind the input file back to the first line
        rewind(io)
 
        ! Allocate array A with dimension 3 by n
        allocate(A(n, 3))

        ! Record the value of index, time, and signal amplitude at each data point to its respective column
        do i = 1, n
            read(io, *, iostat=stat) index, t, x
            if (stat /=0) exit
            A(i, 1) = index
            A(i, 2) = t
            A(i, 3) = x
        end do

        close(io)
    end subroutine load_file


end module read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!