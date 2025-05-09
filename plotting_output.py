import numpy as np
import matplotlib.pyplot as plt

###################################################################################################################
#                                                                                                                 #
#    This program is intended to be used with main.f90 and associated files to visualize pulsar signal data.      #
#                                                                                                                 #
#    Author: Rachel Wood                                                                                          #
#                                                                                                                 #
###################################################################################################################

####################################
#                                  #
###   Plotting original signal   ###
#                                  #
####################################

# Data file which contains index, time, and signal information for each data point
og_data = np.loadtxt('Data_Files/J0452-1759.txt', usecols=(1, 2), unpack=True)

# Define time and signal arrays
time, signal = og_data

fig, ax = plt.subplots(2, 1, figsize=(10, 6))

# Plot of zoomed-in portion of signal data
ax[0].plot(time, signal, color='teal')
ax[0].set_xlim(0, 0.2) # Adjust as desired
ax[0].set_title('Zoomed in signal')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Signal Strength')
ax[0].grid(True)

# Plot of full signal data
ax[1].plot(time, signal, color='darkslategrey')
ax[1].set_xlim(0, np.max(time))
ax[1].set_title('Full signal')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Signal Strength')
ax[1].grid(True)

plt.tight_layout()
plt.show()



##############################################
#                                            #
###   Plotting signal in frequency space   ###
#                                            #
##############################################


# Load full frequency spectrum (from your previously output file)
freqs, mags = np.loadtxt("fft_spectrum_output_J0452-1759.txt", usecols=(0, 1), skiprows=1, unpack=True)

# Load peak data (from the specified file)
peak_freqs, peak_mags = np.loadtxt("fft_peaks_output_J0452-1759.txt", usecols=(0, 1), skiprows=2, unpack=True)

# Plot the full frequency spectrum
plt.figure(figsize=(10, 6))
plt.plot(freqs, mags, label="FFT Spectrum", color='maroon')

# Plot the identified peaks
plt.scatter(peak_freqs, peak_mags, color='darkorange', marker='o', label="Identified Peaks", zorder=5)

# Labels, title, and legend
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.title("FFT Frequency Complete Spectrum with Identified Peaks")
plt.legend()
plt.grid(True)

# Display the plot
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 1, figsize=(10, 12))

# Plot of zoomed-in portion of signal data
ax[0].plot(freqs, mags, color='maroon', label="FFT Spectrum")
ax[0].set_xlim(peak_freqs.min() - 0.5, peak_freqs.max() + 10) # Adjust as desired
ax[0].scatter(peak_freqs, peak_mags, marker='o', color='darkorange')
ax[0].set_title('Zoomed in signal')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Signal Strength')
ax[0].legend()
ax[0].grid(True)

# Plot of full signal data
ax[1].plot(freqs, mags, color='maroon', label="FFT Spectrum")
ax[1].set_xlim(np.min(freqs), np.max(freqs))
ax[1].scatter(peak_freqs, peak_mags, marker='o', color='darkorange')
ax[1].set_title('Full signal')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Signal Strength')
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.show()
