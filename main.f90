!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main

!-----------------------------------------------------------------------------------------------------------------!
!                                                                                                                 !          
!  Main program for use with fft_module.f90 and read_data.f90.                                                    !
!                                                                                                                 !
!  Implements subroutines from both modules to read in an input data file (from the files provided) containing    !
!  pulsar signal data; Perform an fft on the data, and output a file containing information on the prominent      !
!  frequencies present in the signal.                                                                             !
!                                                                                                                 !
!  Author: Rachel Wood                                                                                            !
!                                                                                                                 !
!-----------------------------------------------------------------------------------------------------------------!
!                                                                                                                 !
!    Valid input files are: 'J0006+1834.txt', 'J0108-1431.txt', 'J0152-1637.txt', 'J0206-4028.txt',               !
!                           'J0452-1759.txt', 'J0946+0951.txt', 'J2145-0750.txt', 'J1832+0029.txt',               !
!                           'J1807-2715.txt', 'J1125-5825.txt'.                                                   !
!                                                                                                                 !
!                                                                                                                 !
!  *****REQUIRED INPUTS FROM USER*****                                                                            !
!                                                                                                                 !
!        - Data file ("data_file.txt") : From provided collection of files containing pulsar signal data.         !
!           User must specify file path.                                                                          !
!                                                                                                                 !
!  *****OPTIONAL INPUTS FROM USER*****                                                                            !
!                                                                                                                 !
!        - Number of top peaks desired ("num_peaks") : Number of desired top frequencies to be written to         !
!           output file. Default set to 10.                                                                       !
!                                                                                                                 !
!        - Threshold value ("threshold") : Minimum fraction of max amplitude desired to define a peak. Default    !
!           set to 0.1.                                                                                           !
!                                                                                                                 !
!        - Minimum distance ("min_distance") : Minimum distance between defined peaks. Default set to 50.         !
!                                                                                                                 !
!-----------------------------------------------------------------------------------------------------------------!
!                                                                                                                 !
!    Variables:                                                                                                   !
!                                                                                                                 !
!                  A : array for use storing values of index, time, and signal from the input file data.          !
!               t, x : arrays for use storing time and signal amplitude values for each data point.               !
!    freqs, spectrum : arrays for use storing intermediate frequency and spectrum values.                         !
!         peak_freqs : array for use storing frequency data for identified peaks.                                 !
!          peak_mags : array for use storing magnitude data for identified peaks.                                 !
!         freqs_half : array for use storing the positive values of frequency to avoid unnecesary redundancy      !
!                       which is otherwise present between positive and negative values.                          ! 
!          mags_half : array for use storing the positive values of magnitude to avoid unneccesary redundancy     !
!                       which is otherwise present between positive and negative values.                          !
!          top_freqs : array for use storing frequency data for the top identified peak candidates. number of     !
!                       components determined by num_peaks.                                                       !
!           top_mags : array for use storing magnitude data for the top identified peak candidates. number of     !
!                       components determined by num_peaks.                                                       !
!            x_cmplx : used to store converted complex x values which have been calulated from the real           !
!                       x values.                                                                                 !
!         fft_result : output result from the fft.                                                                !
!       peak_indices : array for use storing new indices for the candidate peaks.                                 !
!          num_peaks : user-input integer value for desired number of top peaks to be written to the              !
!                       output file.                                                                              !
!          N, half_N : original number of data points, and half that number, respectively.                        !
!                  i : temporary index for use in do loops.                                                       !
!                 io : unique integer unit assigned to the current file.                                          !
!                 dt : time step for input data.                                                                  !
!                 fs : sampling frequency of input data.                                                          !
!          threshold : parameter which determines how tall a peak must be relative to the maximum to              !
!                       be considered. can be changed by user as desired.                                         !
!       min_distance : parameter which sets minimum distance between peaks. can be changed by user as desired.    !
!           filename : file name for output text file containing information on final frequency candidates.       !
!                       can be changed by user as desired.                                                        !
!          tolerance : for use in identifying harmonics of the potential fundamental frequency. determines        !
!                       within what tolerance a frequency will be considered a multiple of the minimum.           !
!             labels : array for use storing information on each of the top peaks, indicating which are           !
!                       the fundamental frequency, which are harmonics of it, and which are neither.              !
!           min_freq : minimum frequency value among the values in top_freqs. identifies likely fundamental       !
!                       frequency.                                                                                !
!              ratio : numerical value, indicates how many times min_freq fits into the current frequency.        !
!      rounded_ratio : the ratio value rounded to the nearest integer or half integer.                            !
!          min_index : index of min_freq.                                                                         !
!        is_multiple : logical flag used to indicate whether a given frequency is a harmonic of the               !
!                       fundamental frequency.                                                                    !
!                                                                                                                 !
!-----------------------------------------------------------------------------------------------------------------!

    ! Module containing subroutines for use in reading the input file
    use read_data

    ! Module containing subroutines for use in completing ffts and sorting / filtering data
    use fft_module

    implicit none

    real, allocatable :: A(:,:), t(:), x(:), freqs(:), spectrum(:), peak_freqs(:), peak_mags(:)
    real, allocatable :: freqs_half(:), spectrum_half(:), top_freqs(:), top_mags(:)
    complex, allocatable :: x_cmplx(:), fft_result(:)
    integer, allocatable :: peak_indices(:)
    integer :: num_peaks, N, half_N, i, io
    real :: dt, fs, threshold, min_distance
    character(len=100) :: filename
    real(8) :: tolerance
    character(len=32), dimension(10) :: labels
    real(8) :: min_freq, ratio, rounded_ratio
    integer :: min_index
    logical :: is_multiple
    character(len=100) :: spectrum_filename
    integer :: io_spec, j

    ! *****REPLACE FILE PATH WITH THE PATH FOR THE DESIRED INPUT FILE ON YOUR MACHINE*****

    ! Load data from file (subroutine contained in module 'read_data')
    call load_file('/home/rachel-wood/h_PHYS4840_labs/Final Project/Data_Files/J0452-1759.txt', A, n)

    ! Allocate arrays of length n to store time and signal values
    allocate(t(n), x(n))
    t = A(:, 2)
    x = A(:, 3)

    ! Calculate time step of input data
    dt = t(2) - t(1)

    ! Calculate sampling frequency of input data
    fs = 1.0 / dt

    ! Round N down to next lower power of 2, to ensure no odd-numbered signal lengths
    N = 2**floor(log(real(n))/log(2.0))
    half_N = N / 2

    ! Allocate arrays of length N for the signal and fft results
    allocate(x_cmplx(N), fft_result(N))

    ! Convert the real-valued signal into complex numbers and store them in the x_cmplx array
    do i = 1, N
        x_cmplx(i) = cmplx(x(i), 0.0)
    end do

    ! Perform the Fast Fourier Transform
    call fft_recursive(x_cmplx, fft_result, N)

    ! Allocate frequency array
    allocate(freqs_half(half_N))

    ! Populate the array freqs_half with positive frequency values
    do i = 1, half_N
        freqs_half(i) = (i - 1) * fs / real(N)
    end do

    ! Allocate spectrum array
    allocate(spectrum_half(half_N))

    ! Populate the array spectrum_half with the magnitudes of the fft result for positive frequencies
    do i = 1, half_N
        spectrum_half(i) = abs(fft_result(i))
    end do

    ! *****DEFAULT threshold VALUE IS 0.1. ADJUST AS NEEDED*****
    ! *****DEFAULT min_distance VALUE IS 50. ADJUST AS NEEDED***** 

    ! Find and sort peaks by descending magnitude
    call find_peaks(spectrum_half, freqs_half, threshold, min_distance, peak_freqs, peak_mags)

    num_peaks = 10 ! *****INTEGER NUMBER OF TOP PEAKS DESIRED FOR OUTPUT. ADJUST AS NEEDED*****

    ! Allocate arrays for the following suborutines
    allocate(top_freqs(num_peaks), top_mags(num_peaks))

    call top_peaks(peak_freqs, peak_mags, top_freqs, top_mags, num_peaks)

    ! set tolerance for use in identifying harmonics of the fundamental frequency 
    tolerance = 1.0d-4   ! *****ADJUST AS NEEDED*****

    ! Identify the fundamental frequency and any harmonics of it present in the results
    call identify_harmonics(top_freqs, num_peaks, tolerance, labels)

    ! Write results to file
    filename = "fft_peaks_output_J0452-1759.txt" ! *****OUTPUT FILE NAME. ADJUST AS NEEDED*****
    open(unit=10, file=filename, status="replace", action="write", iostat=io)
    if (io /= 0) then
        print *, "Error opening file ", trim(filename)
        stop
    end if

    write(10, '(A)') "Potential Candidates for Pulsar Frequency"
    write(10, '(A)') "Frequency(Hz)     Magnitude     Harmonics"
    do i = 1, size(top_freqs)
        write(10, '(F12.6, 3X, F12.6, 3X, A)') top_freqs(i), top_mags(i), labels(i)
    end do

    close(10)
    print *, "Results written to ", trim(filename)
    
    ! Output frequency spectrum to file
   

    spectrum_filename = "fft_spectrum_output_J0452-1759.txt"  ! Change name if needed
    open(unit=20, file=spectrum_filename, status="replace", action="write", iostat=io_spec)
    if (io_spec /= 0) then
        print *, "Error opening file ", trim(spectrum_filename)
        stop
    end if

    write(20, '(A)') "Frequency(Hz)     Magnitude"
    do j = 1, size(freqs_half)
        write(20, '(F12.6, 3X, F12.6)') freqs_half(j), spectrum_half(j)
    end do

    close(20)
    print *, "Frequency spectrum written to ", trim(spectrum_filename)


end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!