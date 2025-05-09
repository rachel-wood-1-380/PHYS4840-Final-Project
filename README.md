# PHYS4840-Final-Project
Repository containing programs which use Fourier Transforms to identify potential frequency candidates for various pulsars, given an input file containing signal data.
The primary purpose of this program is to analyze pulsar data from the provided data files, and create an output file containing information on the frequencies detected in the signal, in additon to plotting graphs to visualize the data.
## Installation
Required:
- Compiler for Fortran 90 (gfortran, G95, etc.).
- Python (only necessary if you want to create any plots or figures related to the data. also necessary for the use of additional data files other than the ones in the 'Data Files' folder within this repository.
- Numpy Python module.
- matplotlib Python module

To use this program, you will need to download the provided files: 

**Necessary:**
 - main.f90 : primary program file. this is where you will input file names and any specification preferences. 
 - read_data.f90 : module containing subroutines relevant to reading input files.
 - fft_module.f90 : module containing subroutines relevant to performing a Fast Fourier Transform on a set of data. also contains subroutines which identify and sort peaks in the frequency spectrum.
 - One or more of the data files: 'J0006+1834.txt', 'J0108-1431.txt', 'J0152-1637.txt', 'J0206-4028.txt', 'J0452-1759.txt', 'J0946+0951.txt', 'J2145-0750.txt', 'J1832+0029.txt', 'J1807-2715.txt', 'J1125-5825.txt'
 
**Optional:**
- file_conversion.py : program designed to take PULSE@Parkes data and create new data files that are in a format useable by read_data.f90 and main.f90. only necessary if you want to analyze pulsar data that is not available in the provided files.
- true_values.txt : text file containing info on each pulsar's actual rotation period / frequency.
- plotting_output.py : file containing a number of different plots that show the results.

## Usage
This program is designed to be as user-friendly as possible. There are a limited number of data files provided. This is because it was necessary to convert data files from other sources into a format that made running it through this program as simple as possible. It is possible to use this program to analyze additional pulsar data, but it may require additional steps or complications. 
**To analyze the signal from one of the included pulsar data files:**
1. Download read_data.f90, fft_module.f90, main.f90, and one or more of the .txt data files included in this repository. Save all of these to one folder for simplicity.
2. Open main.f90 in a compatible text editor.
3. Edit line 108 in main.f90 so that the file path reflects that of the data file which you want to run.
4. Edit line 172 in main.f90 so that the file path matches the name and destination desired for the output file.
5. Save main.f90.
6. In command line, navigate to folder where you have saved the relevant files.
7. Compile files: if using gfortran, enter the command: >$ gfortran -o pulsar_program.exe read_data.f90 fft_module.f90 main.f90       hit enter, then enter the command: >$ ./pulsar_program.exe      hit enter again.
8. Check that output files saved correctly.
If you did everything right, you should now have two files (with whatever names you specified for the output files) containing information on the results of the FFT. One file contains the full frequency spectrum, and the other contains info on the most prominent frequency peaks.
In addition to performing an FFT on the data, the program also identifies and sorts frequency peaks, ultimately writing info on the top ten highest-magnitude peaks to a data file. It also identifies the smallest frequency value present in the top ten (most likely to be the true rotational frequency of the pulsar) and labels any harmonic frequencies (multiples of the fundamental frequency). These frequencies aren't likely to be the true frequency of the pulsar.

**Additionally:**

10. Check to see if the program correctly identified the rotational frequency of the pulsar. Period and frequency information for each of the included pulsars is available in the true_values.txt file.
11. Use the plotting_output.py file to visualize the data. You will need to change the file path names in the program to reflect the corresponding ones on your machine before running it.


## Optional Changes
There are a few supplemental features that can be utilized if the user desires. These will be outlined below.
### Change "threshold" Value
The threshold variable is used to specify the minimum magnitude that will be considered as a peak in the frequency spectrum. This filters out noise and other unwanted fluctuations. The default value is 0.1. 
The minimum peak magnitude is found via a simple calculation.
Minimum Value = threshold * Maximum Value
The maximum value is equal to the tallest point in the frequency spectrum. The minimum value is some fraction of that - with default 0.1.

Increasing the threshold value will result in the program considering only peaks which are taller than the new value. This will result in less peaks meeting the requirements to be candidates.
Decreasing the threshold value will result in the program considering a larger number of peaks, with more points meeting the minimum requirements to be considered candidates. This may result in noise being read as a signal, and identified peaks being less accurate. 

### Change "min_distance" Value
The min_distance variable is used to set a minimum distance between peak candidates. This prevents one peak from being interpreted as two in the case of a slight dip in the spectrum. The default value is 50 Hz.

Increasing the min_distance value will result in less peaks being considered, with nearby peaks being excluded. 
Decreasing the min_distance will result in more peaks being considered, at the risk of peaks being considered more than once. 

### Run Additional PULSE@Parkes Signals
If the user desires to use this program to perform an fft and identify pulsar frequencies on pulsars that are not present in the provided set of data, it is not a complicated process. 
This program will work for any three-column data file, with first column index, second column time (in seconds) and third column signal magnitude.
To use any of the other PULSE@Parkes data available at https://research.csiro.au/pulseatparkes/data-analysis/database/
The file_conversion.py file will take any data file from PULSE@Parkes and convert it to a useable format for the main program. It is necessary, however, that some changes be made to the file_conversion.py program. 
**Necessary steps to analyze additional data files:**
1. Copy data over from database at https://research.csiro.au/pulseatparkes/data-analysis/database/  to a .txt file on your computer with header excluded, and file name format 'pulsar_name_old.txt', where pulsar_name is replaced by the J-name of the new pulsar (format JXXXX+XXXX / JXXXX-XXXX).
2. Ensure that file containing 'pulsar_name_old.txt' is in the same folder as the file_conversion.py program.
3. In the file_converion.py program, replace all file paths with the path to the current folder (the same one that you have the 'pulsar_name_old.txt' file in). These replacements should happen on lines 73, 77, 83, and 85.
4. Locate known rotational period (in seconds) of the new pulsar. Set that as P_11 in file_conversion.py.
5. Set h_11 = P_11 / 1024. (Replace 11 with additional numbers if adding multiple new pulsars.)
6. Run file_conversion.py. Ensure new data files are now present in the folder.
7. In main.f90, replace the file ath on line 108 with one appropriate to the data file you want to analyze.
8. In main.f90, replace the destination file path on line 172 with the desired destination path.
9. Compile and run read_data.f90, fft_module.f90, and main.f90 as usual.

## Run Additional Pulsar Signal Data (excluding that from PULSE@Parkes database)
The program can be used to analyze other pulsar signal data as well, but it may be more complicated, depending on the file format.
The data file should be a txt file with no header and three columns. The first column should contain indices of the data points, starting from 1. The second column should contain the time in seconds for each data point. The third column should contain signal magnitude. 
If it is in this exact format, the user should be able to replace the current file path on line 108 with the path for the new data file, make any desired changes to the destination location on line 172, and compile and run the program as usual.
If you do choose to use data from another source, keep in mind that this program was specifically designed for the data format of PULSE@Parkes. You may have to make adjustments to the threshold and min_distance variables in order to get a useful selection of candidate frequencies.


## Methods / Information
Contains information and mathematical methods for a number of the different concepts relevent to the project.
### Pulsars
Pulsars are rapidly rotating neautron stars. They tend to have very short rotational periods ($\lesssim$ 1 second). 
Pulsars emit beams of radio waves that appear as pulses occuring at regular intervals. 
Because of their remarkably stable rotational periods, pulsars can be used to measure time and perform high-precision experiments.

### Fast Fourier Transforms
The Fast Fourier Transoform (FFT) is a remarkably versatile tool and is used in virtually every branch of physics and astronomy. A Fourier Transform takes a given signal in timespace and transforms it into frequency space. This new signal can be used to identify peak frequencies present in the signal.The FFT reduces the number of computations needed compared to a DFT (Discrete Fourier Transform) performed on the same set of data.
This program specifically utilizes the Radix-2 Cooley-Tukey algorithm, which uses a "divide and conquer" strategy to break the DFT into two interleaved DFTs of size N/2.
Given an input signal $x_{n}$, the output frequency space signal $X_{k}$ is given by: 

![image](https://github.com/user-attachments/assets/b4feffad-1715-4eaf-9287-ce3252b01b91)


Which we can rewrite as:



![image](https://github.com/user-attachments/assets/2fa9ccb9-d16b-4d89-ac36-fd9e6e84165b)



Where $E_{k}$ represents the DFT of the even-indexed inputs, and $O_{k}$ represents the DFT of the odd-indexed inputs. 


### Harmonics
Harmonics, or Harmonic Frequencies, are frequencies that are integer multiples of a certain "fundamental frequency".


![image](https://github.com/user-attachments/assets/d130b540-789c-4bdc-a03a-2b708995430c)


Harmonics appear in FFTs when the input signal isn't a pure sine wave. Most of the top frequencies in the FFTs of the pulsar data are harmonics of the fundamental frequency, which, in this case, is the true rotational frequency of our pulsar. 

### PULSE@Parkes

"
About Pulse@Parkes

Through the PULSE@Parkes program, high school students take control of Murriyang our Parkes radio telescope to observe pulsars under the guidance of professional astronomers.

Pulsars, the rapidly spinning remnants of stars after supernova explosions, are studied by astronomers to test the fundamental laws of physics. They have allowed us to investigate the stability of atomic clocks on Earth, and might even help us to detect the gravitational waves predicted by Albert Einstein.

Over the course of a two hour session, high school students may discover a new pulsar, identify unusual ones or detect sudden glitches in their rotation. The data collected is added to a growing database of results and is used by astronomers for ongoing research. Sessions are held onsite in our headquarters at Marsfield, Sydney or online so that schools across Australia can participate.

The program is designed for students in the senior years of high school (Years 10 to 12). There is a simple application process for schools to be involved, and it's free to participate."

- from the PULSE@Parkes website




























## Acknowledgements
