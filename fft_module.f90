!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fft_module

    implicit none

!-----------------------------------------------------------------------------------------------------------------!
!                                                                                                                 !
!  Contains subroutines to perform a Fast Fourier Transform on data, as well as to sort and identify              !
!  frequency peaks and their magnitudes.                                                                          !
!                                                                                                                 !
!  Author: Rachel Wood                                                                                            !
!                                                                                                                 !
!-----------------------------------------------------------------------------------------------------------------!

contains

    recursive subroutine fft_recursive(x, result, N)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !    Subroutine to perform a recursive Fast Fourier Transform on a set of data, using the Cooley-Tukey        !
    !    algorithm for ffts.                                                                                      !
    !                                                                                                             !
    !                  x : input signal in time space. must have an even number of samples.                       !
    !             result : output after the fft is performed. Returns a complex number with real part             !
    !                       representing the magnitude of a particular frequency, and imaginary part              !
    !                       representing the phase of that specific frequency.                                    !
    !                  N : number of samples present in data file.                                                !
    !                  k : temporary index variable for use in do loop.                                           !
    !          even, odd : arrays to store the even and odd indexes (respectively) of the data.                   !
    !  even_fft, odd_fft : arrays to store the even and odd indexed (respectively) results of the fft.            !
    !            twiddle : array to store the "twiddle factor" corresponding to each output of the fft.           !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        complex, dimension(:), intent(IN) :: x
        complex, dimension(:), intent(OUT) :: result
        integer, intent(IN) :: N
        integer :: k
        complex, allocatable :: even(:), odd(:), even_fft(:), odd_fft(:), twiddle(:)

        ! If the signal length is 1, then the fft is just the point itself:
        if (N == 1) then
            result(1) = x(1)
            return
        ! End program if signal length is 1
        end if

        ! Ensure that the signal length is a power of 2
        if (mod(log(real(N)) / log(2.0), 1.0) /= 0) then
            print *, "Signal length must be a power of 2"
            stop
        end if

        ! Allocate arrays
        allocate(even(N/2), odd(N/2), even_fft(N/2), odd_fft(N/2), twiddle(N/2))

        ! Split signal into even and odd indices
        do k = 1, N/2
            even(k) = x(2*k - 1)
            odd(k) = x(2*k)
        end do

        ! Call the recursive fft subroutine within itself
        call fft_recursive(even, even_fft, N/2)
        call fft_recursive(odd, odd_fft, N/2)

        ! Perform the fft 
        do k = 1, N/2
            ! Compute the "twiddle factors"
            twiddle(k) = exp(-2.0 * 3.14159 * cmplx(0.0, 1.0) * real(k - 1) / real(N))
            result(k) = even_fft(k) + twiddle(k) * odd_fft(k)
            result(k + N/2) = even_fft(k) - twiddle(k) * odd_fft(k)
        end do

        ! Deallocate arrays
        deallocate(even, odd, twiddle, even_fft, odd_fft)

    end subroutine fft_recursive


    subroutine find_peaks(spectrum, frequencies, threshold, min_distance, peak_freqs, peak_mags)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !    Subroutine to find peaks in the frequency spectrum, following execution of the fft.                      !
    !                                                                                                             !
    !           spectrum : frequency spectrum as output by the fft subroutine above.                              !
    !        frequencies : array of frequency values.                                                             !
    !          threshold : parameter that determines how tall a peak must be (relative to the maximum) to         !
    !                       be considered. *CAN BE ADJUSTED AS NEEDED*                                            !
    !       min_distance : sets minimum distance between peaks. *CAN BE ADJUSTED AS NEEDED*                       !
    !         peak_freqs : frequencies of identified peaks, used in extract_and_sort_peaks subroutine.            !
    !          peak_mags : magnitudes of identified peaks, used in extract_and_sort_peaks subroutine.             !
    !         candidates : array used to store possible peak candidates, which meet the threshold and             !
    !                       minimum distance requirements.                                                        !
    !       peak_indices : array used to store peak indices.                                                      !
    !            max_val : value of maximum point in frequency spectrum.                                          !
    !            min_val : value of minimum considered magnitude for peak candidates. determined by               !
    !                       threshold (min_value = threshold * max_value).                                        !
    !               i, j : temporary index variables, for use in do loops.                                        !
    !          candidate : temporary value for each element of candidates array, for use in do loops.             !
    !           distance : temporary value for the distance between point for each candidate, for use in          !
    !                       do loops.                                                                             !
    !         peak_count : used to track number of peaks which meet the requirements to be candidates.            !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        real, dimension(:), intent(IN) :: spectrum, frequencies
        real, intent(IN) :: threshold, min_distance
        real, dimension(:), intent(OUT), allocatable :: peak_freqs, peak_mags

        integer, dimension(:), allocatable :: candidates, peak_indices
        real :: max_val, min_val
        integer :: i, j, candidate, distance, peak_count

        logical :: isolated

        ! Allocate arrays
        allocate(candidates(0), peak_indices(0))
        ! Determine maximum value within spectrum array
        max_val = maxval(spectrum)
        ! Calculate minimum value based on magnitude of maximum value
        min_val = threshold * max_val

        ! Loop to identify local peaks which are greater than the minimum value
        do i = 2, size(spectrum) - 1
            if (spectrum(i) > spectrum(i-1) .and. spectrum(i) > spectrum(i+1) .and. spectrum(i) > min_val) then
                call append_to_array(candidates, i)
            end if
        end do

        ! Set initial peak count to zero
        peak_count = 0

        ! Find peaks that are isolated enough from others, distance determined by min_distance
        do i = 1, size(candidates)
            candidate = candidates(i)
            isolated = .true.
            do j = 1, peak_count
                distance = abs(candidate - peak_indices(j))
                if (distance < min_distance) then
                    isolated = .false.
                    exit
                end if
            end do
            if (isolated) then
                peak_count = peak_count + 1
                call append_to_array(peak_indices, candidate)
            end if
        end do

        ! Subroutine defined below, used to sort the peaks and identify the largest ones
        call extract_and_sort_peaks(frequencies, spectrum, peak_indices, peak_freqs, peak_mags)

    end subroutine find_peaks



    subroutine extract_and_sort_peaks(frequencies, spectrum, indices, peak_freqs, peak_mags)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !    Subroutine to sort the identified peaks and restrict to peaks which are tall enough (determined          !
    !    by threshold) and distant enough (determined by min_distance) from other nearby peaks.                   !
    !                                                                                                             !
    !          indices : array to store the indices of peaks                                                      !
    !         spectrum : frequency spectrum as output by the fft subroutine above.                                !
    !      frequencies : array of frequency values.                                                               !
    !       peak_freqs : array used to store frequencies of identified peaks.                                     !
    !        peak_mags : array used to store magnitudes of identified peaks.                                      !
    !          i, j, k : temporary index variables, for use in do loops.                                          !
    !                n : length of indices array.                                                                 !
    !         temp_mag : temporary placeholder value for magnitude, for use when sorting.                         !
    !        temp_freq : temporary placeholder value for frequency, for use when sorting.                         !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        integer, dimension(:), intent(IN) :: indices
        real, dimension(:), intent(IN) :: spectrum, frequencies
        real, intent(OUT), allocatable :: peak_freqs(:), peak_mags(:)
        integer :: i, j, k, n
        real :: temp_mag, temp_freq

        n = size(indices)

        ! Allocate arrays
        allocate(peak_freqs(n), peak_mags(n))

        ! Extract peak data from full frequencies and spectrum arrays, isolating only significant peaks
        do i = 1, n
            peak_freqs(i) = frequencies(indices(i))
            peak_mags(i) = spectrum(indices(i))
        end do

        ! Sort frequencies by magnitude
        do j = 1, n - 1
            do k = j + 1, n
                if (peak_mags(j) < peak_mags(k)) then
                    temp_mag = peak_mags(j)
                    peak_mags(j) = peak_mags(k)
                    peak_mags(k) = temp_mag

                    temp_freq = peak_freqs(j)
                    peak_freqs(j) = peak_freqs(k)
                    peak_freqs(k) = temp_freq
                end if
            end do
        end do
    end subroutine extract_and_sort_peaks


    subroutine append_to_array(array, value)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !  Short subroutine for using appending values to an array.                                                   !
    !                                                                                                             !
    !              array : the input array to which we plan to append.                                            !
    !               temp : temporary array to store values.                                                       !
    !              value : specific value that is being appended to array.                                        !
    !                  n : integer size of array.                                                                 !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        integer, dimension(:), allocatable :: array, temp
        integer :: value, n

        n = size(array)

        allocate(temp(n+1))

        if (n > 0) temp(1:n) = array
        temp(n+1) = value

        deallocate(array)

        allocate(array(n+1))

        array = temp

        deallocate(temp)

    end subroutine append_to_array



    subroutine top_peaks(peak_freqs_in, peak_mags_in, peak_freqs_out, peak_mags_out, num_peaks)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !    Short subroutine to extract the top N peaks, previously sorted by magnitude from                         !
    !    extract_and_sort_peaks subroutine.                                                                       !
    !                                                                                                             ! 
    !      peak_freqs_in : input array containing values of all peak frequencies.                                 !
    !       peak_mags_in : input array containing values of all peak magnitudes.                                  !
    !     peak_freqs_out : output array containing frequency values of a specified number of peaks                !
    !                       (number determined by num_peaks).                                                     !
    !      peak_mags_out : output array containing magnitude values of a specified number of peaks                !
    !                       (number determined by num_peaks).                                                     !
    !          num_peaks : integer number of top peaks desired for the output. *CAN BE ADJUSTED AS NEEDED*        !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        real, dimension(:), intent(IN) :: peak_freqs_in, peak_mags_in
        real, allocatable, intent(OUT) :: peak_freqs_out(:), peak_mags_out(:)
        integer, intent(IN) :: num_peaks

        ! Return only top N peaks
        allocate(peak_freqs_out(num_peaks), peak_mags_out(num_peaks))

        peak_freqs_out = peak_freqs_in(1:num_peaks)
        peak_mags_out = peak_mags_in(1:num_peaks)

    end subroutine top_peaks


    subroutine identify_harmonics(top_freqs, num_peaks, tolerance, labels)

    !-------------------------------------------------------------------------------------------------------------!
    !                                                                                                             !
    !    Subroutine for use identifying the likely fundamental frequency detected from the pulsar signal          !
    !    data, and assesing which of the other top frequency values are just multiples (harmonics) of the         !
    !    fundamental frequency. Produces an additional array that contains information on which peak is likely    !
    !    to be the fundamental frequency and which other frequencies are harmonics.                               !
    !                                                                                                             !
    !          num_peaks : integer number of top peaks desired in output.                                         !
    !          top_freqs : array which contains the frequency values of the top peaks. number of components       !
    !                       determined by num_peaks.                                                              !
    !          tolerance : the tolerance value which determines within what margin a frequency will be            !
    !                       considered a multiple of the minimum frequency (fundamental frequency).               !
    !             labels : array containing harmonic info for each of the frequencies in top_freqs.               !
    !                       specifically labels the fundamental frequency and any harmonics present.              !
    !                                                                                                             !
    !-------------------------------------------------------------------------------------------------------------!

        implicit none

        integer, intent(IN) :: num_peaks
        real(4), intent(IN) :: top_freqs(num_peaks)
        real(8), intent(IN) :: tolerance
        character(len=32), intent(OUT) :: labels(num_peaks)

        real(8) :: min_freq, ratio, rounded_ratio
        integer :: i, min_index
        logical :: is_multiple

        ! Find minimum frequency (most likely to be the fundamental frequency)
        min_index = 1
        min_freq = top_freqs(1)
        do i = 2, num_peaks
            if (top_freqs(i) < min_freq) then
                min_freq = top_freqs(i)
                min_index = i
            end if
        end do

        ! Label each frequency either "fundamental frequency", "harmonic" or "n/a"
        do i = 1, num_peaks
            if (i == min_index) then
                labels(i) = "fundamental frequency"
            else
                ratio = top_freqs(i) / min_freq
                rounded_ratio = nint(ratio * 2.0d0) / 2.0d0
                is_multiple = abs(ratio - rounded_ratio) < tolerance

                if (is_multiple) then
                    write(labels(i), '(A, F4.1, A)') "harmonic (", rounded_ratio, "× f₀)"
                else
                    labels(i) = "n/a"
                end if
            end if
        end do
    end subroutine identify_harmonics


end module fft_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!