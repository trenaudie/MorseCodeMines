#%%
import numpy as np 
import matplotlib.pyplot as plt
import wave


""" 1. Read and Parse the WAV File
* Open the WAV file as a binary stream
* Parse the WAV header to extract sample rate, bit depth, number of channels
* Load the audio samples into memory as an array of amplitude values
 """
wavfile = wave.open("morse_code.wav", "r")
values = []
for i in range(0, wavfile.getnframes()):
    frame = wavfile.readframes(1)
    value = int.from_bytes(frame, byteorder='little', signed=True)
    values.append(value)
values = np.array(values)
plt.plot(values)

#%% 

fft_values_raw = np.fft.fft(values)
fft_values = np.abs(fft_values_raw)
values_rebuilt = np.fft.ifft(fft_values_raw)
plt.plot(values_rebuilt)
#%%
plt.plot(values)
# %%
# what is the sample rate of the signal?
sample_rate = wavfile.getframerate()
sample_rate

#%% 
# amplitude of the signal
def amplitude(values):
    return np.max(values) - np.min(values)
amplitude_ = np.max(values) - np.min(values)
amplitude_
threshold = 0.2 * amplitude_
#%%
sample_rate = wavfile.getframerate()
# fft with sample rate 
fft_values = np.fft.fft(values)
fft_values = np.abs(fft_values)
frequencies = np.fft.fftfreq(len(fft_values), 1/sample_rate)
plt.plot(frequencies, fft_values)
#%%
plt.plot(frequencies[:len(frequencies)//20], fft_values[:len(fft_values)//20])

#%%



def dft(arr:np.ndarray):
    """
    Compute the Discrete Fourier Transform (DFT) of a signal by hand.
    """
    from tqdm import tqdm
    N = len(arr)
    dft_values = np.zeros(N, dtype=complex)
    for k in tqdm(range(N)):
        for n in range(N):
            dft_values[k] += arr[n] * np.exp(-2j * np.pi * k * n / N)
    return dft_values


#%% 
non_null_index = 0 
while non_null_index < len(values) -1 : # keep moving the non_null_index until we find a non null value that satisfies the right range
    print(non_null_index)
    while values[non_null_index] == 0:
        non_null_index += 1
    # update the range end
    range_end = non_null_index + 1
    while values[range_end] != 0 or values[range_end +1] != 0:
        range_end += 1
        print(f'range: {range_end- non_null_index}')
    if range_end - non_null_index > 1000:
        print(f'breaking because minimum range is reached: {range_end - non_null_index}')
        break 
    if values[range_end] == 0 and values[range_end+1] == 0:
        print(f'found an empty signal at index {non_null_index}')
        non_null_index = range_end + 1
        continue

#%%
assert range_end - non_null_index > 999
values_sinusoidal = values[non_null_index+1 : range_end - 1][:3000]
dft_values = dft(values_sinusoidal)
frequencies = np.fft.fftfreq(len(dft_values), 1/sample_rate)
idx = np.argsort(np.abs(dft_values))[-10:]
print(frequencies[idx]) 
major_frequency = frequencies[idx[-1]]
print(f'top frequency: {major_frequency}')

#%%

#%% 
# top10 fft_values, frequencies pairs  
top10 = np.argsort(fft_values[:len(fft_values)//20])[-10:]
top10_freqs = frequencies[:len(frequencies)//20][top10]
top10_freqs
freq_found = 546
T_sinewave = 1 / freq_found
num_points_sinewave = T_sinewave * sample_rate
num_points_sinewave

#%%

import plotly.express as px
fig = px.line(x=frequencies, y=fft_values)
fig.show()
#%%
# how to normalize the amplitude of the signal using a window?
window_size = 1024
window = np.hanning(window_size)

# Create the Hanning window
window_size = 1024
window = np.hanning(window_size)

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
fig.suptitle('Hanning Window Analysis', fontsize=16)

# Plot 1: Hanning window in time domain
ax1.plot(window, 'b')
ax1.set_title('Hanning Window (Time Domain)')
ax1.set_xlabel('Sample Index')
ax1.set_ylabel('Amplitude')
ax1.grid(True, alpha=0.3)

# Plot 2: Histogram of window values
ax2.hist(window, bins=50, color='green', alpha=0.7)
ax2.set_title('Distribution of Hanning Window Values')
ax2.set_xlabel('Value')
ax2.set_ylabel('Frequency')
ax2.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

#%%

import numpy as np
import matplotlib.pyplot as plt

def windowed_normalization(values, window_size=1024, target_level=0.5, method='rms'):
    """
    Apply windowed normalization to an audio signal using 50% overlap.
    
    Parameters:
    - values: Input audio array
    - window_size: Size of the processing window
    - target_level: Target normalization level (0.0 to 1.0)
    - method: 'rms' or 'peak' normalization
    
    Returns:
    - normalized_values: Normalized audio array
    """
    # Create Hanning window
    window = np.hanning(window_size)
    
    # Pad the input signal to handle edge cases
    pad_size = window_size // 2
    padded_values = np.pad(values, (pad_size, pad_size), mode='edge')
    
    # Initialize output array
    normalized_values = np.zeros_like(values)
    
    # Create arrays to accumulate the processed samples and window weights
    output_accumulator = np.zeros(len(padded_values))
    window_accumulator = np.zeros(len(padded_values))
    
    # Step size for 50% overlap
    hop_size = window_size // 2
    
    # Process the signal in overlapping windows
    for i in range(0, len(values), hop_size):
        # Define the range for this window
        start = i
        end = start + window_size
        
        if end > len(padded_values):
            break
            
        # Extract the current frame
        frame = padded_values[start:end]
        
        # Apply window to the frame
        windowed_frame = frame * window
        
        # Compute normalization factor based on method
        if method == 'rms':
            # Root Mean Square level
            level = np.sqrt(np.mean(windowed_frame**2))
        else:  # peak
            # Peak level
            level = np.max(np.abs(windowed_frame))
        
        # Avoid division by zero
        if level < 1e-10:
            gain = 1.0
        else:
            gain = target_level / level
        
        # Apply gain to the frame
        normalized_frame = frame * gain
        
        # Apply window again for smooth overlap
        windowed_normalized_frame = normalized_frame * window
        
        # Accumulate the processed samples and window weights
        output_accumulator[start:end] += windowed_normalized_frame
        window_accumulator[start:end] += window
    
    # Divide by the accumulated window weights to normalize overlapping regions
    # Add small epsilon to avoid division by zero
    normalized_output = output_accumulator[pad_size:-pad_size] / (window_accumulator[pad_size:-pad_size] + 1e-10)
    
    return normalized_output
values[100:200] = 60000
plt.plot(values)
plt.show()
values_normalized = windowed_normalization(values, window_size=108, target_level=0.5, method='rms')
plt.plot(values_normalized)
plt.show()
#%%
amplitude_norm = amplitude(values_normalized)
threshold_norm = 0.2 * amplitude_norm
#%% 
values_normalized_abs = np.abs(values_normalized)

above_threshold = values_normalized_abs > threshold_norm
plt.imshow(above_threshold.reshape(1, -1)[:1,:100], cmap='gray', aspect='auto')
plt.title('Above Threshold')
plt.show()
# %%
above_threshold[:100]
# %%
# for each black pixel, if the black spotch is longer than x pixels, then its is a morse code for blank 
plt.plot(values[:100])
#%%
plt.plot(above_threshold[:100])

#%% 
# Determining the freq by hand
# it takes about 80 points to make a sine wave
wavelength = 80
wavelength_secs = wavelength / sample_rate
wavelength_secs
frequency = 1 / wavelength_secs
frequency # nice, it is 550 hz approx


# %%
"""
IT MAY BE IMPORTANT TO CALCULATE THE POINTS PER WAVELENGTH
why ? 
then we can put a threshold on the blank equal to half a wavelength
"""

#%% \
# must remove the end blanks 

ptr_end = 0
while above_threshold[-(ptr_end+1)] == 0:
    ptr_end += 1
    if ptr_end > len(above_threshold):
        break
above_threshold_noend = above_threshold[:-ptr_end]
#%% 
plt.imshow(above_threshold_noend.reshape(1, -1)[:1,:10_000], cmap='gray', aspect='auto')
#%% 
signal_or_blank = []
blank_counter = 0
blank_threshold = 80 * 20 # 20 wavelengths ie negligible, but safer than 2
for i in range(len(above_threshold_noend)):
    if above_threshold_noend[i] == 1:
        if blank_counter == 0:
            signal_or_blank.append(1)
        else: 
            blank_counter +=1 
            if blank_counter > blank_threshold:
                signal_or_blank.extend([0] * blank_counter)
            else:
                signal_or_blank.extend([1] * blank_counter)
        blank_counter = 0
    elif above_threshold_noend[i] == 0:
        blank_counter += 1
if blank_counter > blank_threshold:
    signal_or_blank.extend([0] * blank_counter)
else:
    signal_or_blank.extend([1] * blank_counter)

# %%
plt.imshow(np.array(signal_or_blank).reshape(1, -1)[:,:], cmap='gray', aspect='auto')

# %%
durations = []
block_values = []
ptr = 0
current_value = signal_or_blank[0]
num_conseq = 0
while ptr < len(signal_or_blank) -1 :
    ptr += 1
    if signal_or_blank[ptr]  == current_value:
        num_conseq += 1
    else:
        block_values.append(current_value)
        durations.append(num_conseq)
        num_conseq = 0
        current_value = signal_or_blank[ptr]
        num_conseq += 1
durations.append(num_conseq)
block_values.append(current_value)
durations
#%%
block_values
# %%
# 5200 points is how many seconds?
unitlen = min(durations)
unitlen_secs = 5200 / sample_rate
unitlen_secs
unitlen_samples_buffer = unitlen * 0.2
#%%
num_unit_lens = np.array(durations) / unitlen
num_unit_lens = np.round(num_unit_lens)
num_unit_lens_unique = np.unique(num_unit_lens)
print(f'number of unique unit lengths: {len(num_unit_lens_unique)}')
# should be 3, but can be 4

len_unit = 1 
len_letter_space_or_dash = num_unit_lens_unique[1]
len_word_space = num_unit_lens_unique[-1]
len_unit, len_letter_space_or_dash, len_word_space
# %%

symbols1 = ""
for i,block in enumerate(num_unit_lens):
    block_value = block_values[i]
    if block == len_unit:
        if block_value == 1:
            symbols1 += "."
        else:
            symbols1 += ""
    elif block == len_letter_space_or_dash:
        if block_value == 1:
            symbols1 += "-"
        else:
            symbols1 += " "
    elif block == len_word_space and block_value == 0:
        symbols1 += "/"
symbols1

# %%
symbols1 = '.... ../.. -.'
# %%
def morse_to_text(morse_string):
    # Morse code lookup table
    morse_dict = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E',
        '..-.': 'F', '--.': 'G', '....': 'H', '..': 'I', '.---': 'J',
        '-.-': 'K', '.-..': 'L', '--': 'M', '-.': 'N', '---': 'O',
        '.--.': 'P', '--.-': 'Q', '.-.': 'R', '...': 'S', '-': 'T',
        '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X', '-.--': 'Y',
        '--..': 'Z', '.----': '1', '..---': '2', '...--': '3', '....-': '4',
        '.....': '5', '-....': '6', '--...': '7', '---..': '8', '----.': '9',
        '-----': '0', '--..--': ',', '.-.-.-': '.', '..--..': '?',
        '-..-.': '/', '-.--.': '(', '-.--.-': ')', '.-...': '&',
        '---...': ':', '-.-.-.': ';', '-...-': '=', '.-.-.': '+',
        '-....-': '-', '..--.-': '_', '.-..-.': '"', '...-..-': '$',
        '.--.-.': '@', '-.-.--': '!'
    }
    
    # Split the input string into words (separated by '/')
    words = morse_string.strip().split('/')
    translated_text = []
    
    for word in words:
        # Split word into individual morse characters (separated by spaces)
        chars = word.strip().split()
        word_chars = []
        
        for char in chars:
            # Look up the character in the dictionary
            if char in morse_dict:
                word_chars.append(morse_dict[char])
            # Handle special case if the morse code isn't in our dictionary
            elif char == '':
                continue  # Skip empty characters
            else:
                word_chars.append('<' + char + '>')  # Mark unknown patterns
                
        # Join all characters in this word
        translated_text.append(''.join(word_chars))
    
    # Join all words with spaces
    return ' '.join(translated_text)

# Test with the given example
symbols1 = '.... ../.. -.'
result = morse_to_text(symbols1)
print(f"Morse: {symbols1}")
print(f"Text: {result}")
# %%
