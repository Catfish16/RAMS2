import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks


def get_sampling_freq(input_frequency):
    """ Processes input sampling frequency
        Parameters:
           input_frequency (str) : File to be loaded and processed

        Returns:
            int : Sampling frequency as integer
    """
    if input_frequency.isdigit():
        return int(input_frequency)
    else:
        print("Please, enter a valid number.")


# Function to get the start and end times input from the user
def get_times(input_times):
    """ Process input start and end times
        Parameters:
           input_times (str) : File to be loaded and processed

        Returns:
            list[float] : A list of start or end times
    """
    while True:
        try:
            times = list(map(float, input(input_times).strip().split(',')))
            return times
        except ValueError:
            print("Please enter valid numbers separated by commas.")


def load_signal(filename):
    """ Processes .mat file
        Parameters:
           filename (str) : File to be loaded and processed

        Returns:
            any : 1 channel signal from matlab file
    """
    mat_data = loadmat(f"signals/{filename}")
    print(f"Keys in the .mat file: {mat_data.keys()}")

    # Console input: key name
    signal_key = input("Please enter the key for the signal data: ")

    if signal_key in mat_data:
        signal = mat_data[signal_key][0]  # Adjusting indexing based on typical .mat file structure
        return signal
    else:
        raise KeyError(f"Key '{signal_key}' not found in the .mat file.")


def process_intervals(signal, start_times, end_times, fs):
    """ Processes signal in intervals between start and end time pairs
        Parameters:
           signal (any) : ECG signal to be processed
           start_times (list[float]): List of start times in second
           end_times (list[float]): List of end times in seconds
           fs (int): sampling frequency of the signal

        Returns:
            Prints results to console
    """
    # Variables for averaging several results
    hr_values = []
    sdnn_values = []
    sdsd_values = []

    for start, end in zip(start_times, end_times):
        start_x = int(start * fs)
        end_x = int(end * fs)
        ecg_signal_interval = signal[start_x:end_x]

        # Find peaks in the interval
        # Height: minimal amplitude to qualify as a peak
        # Distance: minimal distance between the peaks (=> assumption average hr ~60 BPM)
        peaks, _ = find_peaks(ecg_signal_interval, height=0.5, distance=fs / 4)

        # Calculate R-R peak intervals in seconds
        rr_intervals = np.diff(peaks) / fs

        # Calculate heart rate in bpm
        heart_rate = [round(60 / rr) for rr in rr_intervals]
        hr_values.extend(heart_rate)

        avg_heart_rate = int(np.average(heart_rate))
        min_heart_rate = np.min(heart_rate)
        max_heart_rate = np.max(heart_rate)

        print(f"\nProcessing interval from {start}s to {end}s (indices {start_x} to {end_x}):")
        # print(f"Heart rates: {heart_rate}")
        print(f"Average heart rate: {avg_heart_rate}")
        print(f"Minimal heart rate: {min_heart_rate}, maximal heart rate: {max_heart_rate}")

        if len(rr_intervals) > 1:
            # Calculate SDNN
            sdnn = np.std(rr_intervals)
            sdnn_values.append(sdnn)

            # Calculate successive differences and SDSD
            successive_diffs = np.diff(rr_intervals)
            sdsd = np.std(successive_diffs)
            sdsd_values.append(sdsd)

            print(f"SDNN: {sdnn:.4f} seconds")
            print(f"SDSD: {sdsd:.4f} seconds")
        else:
            print(f"Interval from {start}s to {end}s has insufficient peaks for heart rate variability (HRV) analysis.")

    # Calculate statistics for more than 7 values
    if len(start_times) > 7:
        if hr_values:
            hr_mean = np.mean(hr_values)
            hr_std = np.std(hr_values)
            hr_min = np.min(hr_values)
            hr_max = np.max(hr_values)

            print("\nHeart Rate Statistics:")
            print(f"Mean HR: {hr_mean:.2f} bpm")
            print(f"Standard Deviation HR: {hr_std:.2f} bpm")
            print(f"Minimum HR: {hr_min:.2f} bpm")
            print(f"Maximum HR: {hr_max:.2f} bpm")

        if sdnn_values:
            sdnn_mean = np.mean(sdnn_values)
            sdnn_std = np.std(sdnn_values)
            sdnn_min = np.min(sdnn_values)
            sdnn_max = np.max(sdnn_values)

            print("\nSDNN Statistics:")
            print(f"Mean SDNN: {sdnn_mean:.4f} seconds")
            print(f"Standard Deviation SDNN: {sdnn_std:.4f} seconds")
            print(f"Minimum SDNN: {sdnn_min:.4f} seconds")
            print(f"Maximum SDNN: {sdnn_max:.4f} seconds")

        if sdsd_values:
            sdsd_mean = np.mean(sdsd_values)
            sdsd_std = np.std(sdsd_values)
            sdsd_min = np.min(sdsd_values)
            sdsd_max = np.max(sdsd_values)

            print("\nSDSD Statistics:")
            print(f"Mean SDSD: {sdsd_mean:.4f} seconds")
            print(f"Standard Deviation SDSD: {sdsd_std:.4f} seconds")
            print(f"Minimum SDSD: {sdsd_min:.4f} seconds")
            print(f"Maximum SDSD: {sdsd_max:.4f} seconds")


# Main function
def main():
    # User manual:
    print("The following console application can be used analyse single channel ECG signal of a matlab (.mat) file."
          "\n \n The user will be requested for the following inputs:"
          "\n* File name: The .mat file must be located in the signals folder of this program."
          "\n* Key: variable in which the signal is stored in the .mat file"
          "\n* Sampling frequency of the signal"
          "\n* Start and end times: Intervals of the signal to be examined"
          "\n \n The application will return the following values:"
          "\n* Average, minimum and maximum heart rate"
          "\n* HRV (heart rate variability measures:"
          "\n  ~SDNN: Standard deviation of NN intervals. "
          "\n  ~SDSD: Standard deviation of successive differences of NN intervals."
          "\n* If more than 7 interval pairs were defined, the application will return statistics for these."
          "\n \n Please, note that the application will interrupt execution if it encounters problems. "
          "\n --------------------------------------------------------------------------------------------------"
          "\n")

    # Console input: file name
    mat_filename = input(f"Enter the name of the .mat file (include extension): ")

    try:
        # Load the signal from the .mat file
        signal = load_signal(mat_filename)

        fs = get_sampling_freq(input("Enter the sampling frequency of the signal in Herz (e.g. 200): "))

        signal_length_secs = len(signal) / fs

        # Console input: start and end times in seconds
        start_times = get_times("Enter start times (in seconds), separated by commas: ")
        end_times = get_times(f"Enter end times (in seconds), "
                              f"separated by commas (signal ends after {signal_length_secs} seconds): ")

        # Ensure start and end times lists are the same length
        if len(start_times) != len(end_times):
            print("Start times and end times must have the same number of elements.")
            return

        # Process each interval
        process_intervals(signal, start_times, end_times, fs)
    except FileNotFoundError:
        print(f"File {mat_filename} not found.")
    except KeyError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
