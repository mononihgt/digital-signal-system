# main_exp2.py
import numpy as np
import pandas as pd
import sounddevice as sd
from psychopy import core, visual, event, gui, logging
import random
import os
import time
import json
from scipy import signal # Needed for filtering
from questplus.qp import QuestPlus
from questplus import _constants

# --- Helper Function ---
def find_closest_intensity(qp_instance, presented_intensity):
    """Finds the intensity in the Quest+ stim_domain closest to the presented one."""
    known_intensities = qp_instance.stim_domain['intensity']
    idx = np.argmin(np.abs(known_intensities - presented_intensity))
    closest_intensity = known_intensities[idx]
    return closest_intensity

# --- Experiment Configuration ---
# Participant Info GUI
exp_info = {'participant': 'test'}
dlg = gui.DlgFromDict(dictionary=exp_info, title='Pure Tone Threshold Experiment 2 (Bandwidth)')
if not dlg.OK:
    core.quit()
participant_id = exp_info['participant']

# Stimulus Parameters
FS = 16000  # Sampling rate (Hz) - Keep consistent with filter design
DURATION_NOISE = 1.0
DURATION_PURE = 0.100
TONE_ONSET_IN_NOISE = 0.5
FREQ_PURE = 1000  # Center frequency (Hz)
GAP_DURATION = 0.5

# Define the bandwidth parameters (A values in Hz)
A_VALUES_TO_TEST = [100, 300, 900] # Hz

# Filter Parameters
FILTER_ORDER = 64 # FIR filter order (n in the example script was order+1)

N_TRIALS_PRACTICE_PER_COND = 3 # Practice trials *per A value*
N_TRIALS_MAIN_PER_COND = 32  # Main trials *per A value*
N_TRIALS_PRACTICE = N_TRIALS_PRACTICE_PER_COND * len(A_VALUES_TO_TEST)
N_TRIALS_MAIN = N_TRIALS_MAIN_PER_COND * len(A_VALUES_TO_TEST)
BREAK_INTERVAL = 32 # Show break message every N total main trials

# Response Keys
KEY_INTERVAL_1 = '1'
KEY_INTERVAL_2 = '2'
KEY_QUIT = 'q'

# Data Saving
DATA_FOLDER = 'data_exp2' # Use a different folder for Exp 2 data
if not os.path.exists(DATA_FOLDER):
    os.makedirs(DATA_FOLDER)
filename_base = f"{participant_id}_bandwidth" # Updated base name
data_filename = os.path.join(DATA_FOLDER, f"{filename_base}_trials.csv")
results_filename = os.path.join(DATA_FOLDER, f"{filename_base}_results.txt")
# We will save QP states separately for each A value
qp_state_filenames = {
    A_val: os.path.join(DATA_FOLDER, f"{filename_base}_A{A_val}Hz_questplus_state.json")
    for A_val in A_VALUES_TO_TEST
}


# Check if data file exists to prevent overwriting
if os.path.exists(data_filename) or any(os.path.exists(f) for f in qp_state_filenames.values()):
    overwrite_dlg = gui.Dlg(title="File Exists")
    overwrite_dlg.addText(f"Data files for '{filename_base}' might already exist.")
    overwrite_dlg.addText("Do you want to overwrite?")
    overwrite_dlg.show()
    if not overwrite_dlg.OK:
        core.quit()
    else:
        print(f"Overwriting existing files for {filename_base}")

# --- Cosine Ramp Function ---
def apply_cosine_ramp(audio, sample_rate, fade_time=0.02):
    """
    cosine ramp for fade in/out
    """
    fade_length = int(fade_time * sample_rate)
    fade_in = (1 - np.cos(np.linspace(0, np.pi, fade_length))) / 2
    
    fade_out = (1 + np.cos(np.linspace(0, np.pi, fade_length))) / 2
    
    audio[:fade_length] *= fade_in
    audio[-fade_length:] *= fade_out
    
    return audio

# --- Quest+ Configuration (Common Settings) ---
# Define the domains for stimuli and parameters
stim_intensities_ratio = np.linspace(0.001, 1.6, 1200) # P/N ratio (relative to normalized noise)
thresholds_ratio = stim_intensities_ratio.copy()
slopes = np.linspace(1, 15, 10)
lower_asymptotes = [0.5] # Fixed for 2AFC
lapse_rates = np.linspace(0.0, 0.05, 5)

stim_domain = dict(intensity=stim_intensities_ratio)
param_domain = dict(threshold=thresholds_ratio,
                    slope=slopes,
                    lower_asymptote=lower_asymptotes,
                    lapse_rate=lapse_rates)
outcome_domain = dict(response=['Correct', 'Incorrect'])

# --- Initialize QuestPlus instances for each A value ---
qps = {}
for A_val in A_VALUES_TO_TEST:
    qps[A_val] = QuestPlus(
        stim_domain=stim_domain,
        param_domain=param_domain,
        outcome_domain=outcome_domain,
        func='weibull',
        stim_scale='linear',
        param_estimation_method='mean',
        stim_selection_method='min_entropy'
    )
    print(f"Quest+ Initialized for A = {A_val} Hz")


# --- PsychoPy Setup ---
win = visual.Window(
    size=[800, 600], fullscr=True, screen=0,
    winType='pyglet', allowGUI=True, allowStencil=False,
    monitor='testMonitor', color=[0.8, 0.8, 0.8], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    units='height')

# Instructions and Feedback Text (updated for bandwidth)
instructions_text = """
Welcome to Experiment 2!

In each trial, you will hear two sounds separated by a short silence.
Each sound contains FILTERED noise, and the bandwidth of the noise might change between trials.
One of the two sounds will ALSO contain a faint pure tone (1000 Hz).

Your task is to identify whether the PURE TONE was in the FIRST or SECOND sound interval.

Press '1' if you heard the tone in the FIRST interval.
Press '2' if you heard the tone in the SECOND interval.

Please listen carefully. The tone might be very quiet.
Guess if you are unsure.

There will be practice trials first.

Press any key to start the practice.
"""

practice_over_text = "Practice is over. Press any key to start the main experiment."
exp_start_text = "Press any key to begin the main experiment."
break_text = "Take a short break. Press any key to continue."
end_text = "Experiment finished! Thank you. Saving results..."

fixation = visual.TextStim(win, text='+', height=0.05, color='black')
message_stim = visual.TextStim(win, text='', height=0.04, color='black', wrapWidth=win.size[0]*0.8/win.size[1])
response_prompt_stim = visual.TextStim(win, text="Tone in: [1] First or [2] Second?", height=0.05, color='black', pos=(0, -0.2))

# --- Filter and Normalization Helper ---
def generate_filtered_normalized_noise(A_value, fs, n_samples, filter_order):
    """Generates, filters, and normalizes white noise to RMS=1."""
    center_freq = 1000 # Hz

    # Calculate filter cutoff frequencies
    f1 = center_freq - A_value
    f2 = center_freq + A_value

    # Ensure cutoffs are valid
    f1 = max(1, f1) # Avoid 0 Hz or negative frequencies
    f2 = min(fs/2 - 1, f2) # Avoid Nyquist frequency

    if f1 >= f2:
        # This shouldn't happen with the given A values, but good to check
        print(f"Warning: Invalid frequency range for A={A_value}. f1={f1}, f2={f2}. Using broad noise.")
        # Fallback: generate non-filtered normalized noise
        raw_noise = np.random.randn(n_samples)
        rms_raw = np.sqrt(np.mean(raw_noise**2))
        if rms_raw > 1e-9:
            return raw_noise / rms_raw
        else:
            return np.zeros(n_samples) # Return silence if RMS is zero


    # Normalize frequencies for the filter design function
    nyquist = fs / 2
    Wn = [f1 / nyquist, f2 / nyquist]

    # Design FIR bandpass filter
    # Note: Scipy's firwin uses numtaps = order + 1
    filter_coeffs = signal.firwin(filter_order + 1, Wn, pass_zero='bandpass',window='hamming')

    # Generate white noise
    raw_noise = np.random.randn(n_samples)

    # Apply filter
    filtered_noise = signal.lfilter(filter_coeffs, 1.0, raw_noise)

    # Calculate RMS of the filtered noise
    rms_filtered = np.sqrt(np.mean(filtered_noise**2))

    # Normalize to RMS = 1 (Power = 1)
    if rms_filtered > 1e-9: # Avoid division by zero
        normalized_noise = filtered_noise / rms_filtered
    else:
        normalized_noise = np.zeros(n_samples) # Should not happen often

    return normalized_noise

# --- Audio Generation Function (Modified for Filtering/Normalization) ---
def generate_trial_sound(A_value, pn_ratio, correct_interval, fs,
                         duration_noise, duration_pure, tone_onset_in_noise, freq_pure, gap_duration, filter_order):
    """Generates the 2AFC sound stimulus with filtered, normalized noise."""
    pn_ratio = max(0, pn_ratio)
    # Power of the tone relative to the *normalized* noise power (which is 1)
    pure_power = pn_ratio * 1.0
    amplitude_pure = np.sqrt(2 * pure_power)

    # Time vectors and sample counts
    t_pure = np.arange(0, duration_pure, 1/fs)
    n_samples_noise = int(duration_noise * fs)
    n_samples_pure = len(t_pure)
    n_samples_gap = int(gap_duration * fs)
    n_samples_onset = int(tone_onset_in_noise * fs)
    n_samples_before_tone = n_samples_onset
    n_samples_after_tone = max(0, n_samples_noise - n_samples_before_tone - n_samples_pure)

    # Adjust pure tone length if it exceeds remaining time
    if n_samples_before_tone + n_samples_pure > n_samples_noise:
         n_samples_pure = n_samples_noise - n_samples_before_tone
         pure_tone_signal = amplitude_pure * np.sin(2 * np.pi * freq_pure * t_pure[:n_samples_pure])
         n_samples_after_tone = 0
         print(f"Warning: Tone duration ({duration_pure}s) + onset ({tone_onset_in_noise}s) exceeds noise duration ({duration_noise}s). Truncating tone.")
    else:
         pure_tone_signal = amplitude_pure * np.sin(2 * np.pi * freq_pure * t_pure)

    # Generate two independent segments of filtered, normalized noise
    noise1 = generate_filtered_normalized_noise(A_value, fs, n_samples_noise, filter_order)
    noise2 = generate_filtered_normalized_noise(A_value, fs, n_samples_noise, filter_order)

    # Prepare padded tone
    padding_before = np.zeros(n_samples_before_tone)
    padding_after = np.zeros(n_samples_after_tone)

    # Apply cosine ramp to the pure tone signal
    pure_tone_signal_ramped = apply_cosine_ramp(pure_tone_signal, fs, fade_time=0.02)
    tone_padded = np.concatenate([padding_before, pure_tone_signal_ramped, padding_after])
    # Final length check
    if len(tone_padded) < n_samples_noise:
        tone_padded = np.concatenate([tone_padded, np.zeros(n_samples_noise - len(tone_padded))])
    elif len(tone_padded) > n_samples_noise:
        tone_padded = tone_padded[:n_samples_noise]

    # Combine tone with the correct noise interval
    if correct_interval == 1:
        interval1 = noise1 + tone_padded
        interval2 = noise2
    else:
        interval1 = noise1
        interval2 = noise2 + tone_padded

    # Create final signal
    silence_gap_signal = np.zeros(n_samples_gap)
    combined_signal = np.concatenate([interval1, silence_gap_signal, interval2])
    stereo_signal = np.vstack([combined_signal, combined_signal]).T
    return stereo_signal

# --- Create Trial List ---
def create_trial_list(n_trials_per_cond, a_values, block_type):
    """Creates a list of trials with randomly interleaved A values."""
    trial_conditions = []
    for A_val in a_values:
        trial_conditions.extend([A_val] * n_trials_per_cond)
    random.shuffle(trial_conditions) # Interleave the conditions randomly
    # Add trial number and block type
    trials = [{'trial_overall': i + 1, 'A_value': cond, 'block': block_type}
               for i, cond in enumerate(trial_conditions)]
    return trials


# --- Run Experiment Function (Modified for A Values) ---
def run_experiment_block(trial_list, is_practice):
    """Runs a block of trials with interleaved A value conditions."""
    trial_data_list = []
    global qps # Allow access to the global qps dictionary

    practice_high_pn_ratio = 1.2 # Use a relatively high ratio for initial practice
    practice_trials_done_per_cond = {level: 0 for level in A_VALUES_TO_TEST}
    initial_practice_phase_count = 1 # How many initial trials use fixed high ratio

    for trial_info_base in trial_list:
        trial_num_overall = trial_info_base['trial_overall']
        current_A_value = trial_info_base['A_value']
        current_qp = qps[current_A_value] # Select the QP instance for this A value

        # --- Break Time --- (Based on overall *main experiment* trial number)
        main_trial_index = trial_num_overall - N_TRIALS_PRACTICE # Index within main block
        if not is_practice and (main_trial_index % BREAK_INTERVAL == 1 and main_trial_index > 1):
            message_stim.text = break_text
            message_stim.draw()
            win.flip()
            # Use space to continue, escape to potentially abort
            keys = event.waitKeys(keyList=['space', 'escape'])
            if 'escape' in keys:
                 print("Escape pressed during break. Exiting.")
                 return None # Signal exit
            event.clearEvents() # Clear buffer after wait

        # --- Get Stimulus Intensity ---
        if is_practice and practice_trials_done_per_cond[current_A_value] < initial_practice_phase_count:
             presented_pn_ratio = practice_high_pn_ratio
             practice_trials_done_per_cond[current_A_value] += 1
        else:
             next_stim_dict = current_qp.next_stim
             presented_pn_ratio = next_stim_dict['intensity']
             presented_pn_ratio = max(0.0001, presented_pn_ratio) # Clamp

        # --- Prepare Trial Sound ---
        correct_interval = random.choice([1, 2])
        stimulus_sound = generate_trial_sound(
            A_value=current_A_value, # Pass A value
            pn_ratio=presented_pn_ratio,
            correct_interval=correct_interval,
            fs=FS,
            duration_noise=DURATION_NOISE,
            duration_pure=DURATION_PURE,
            tone_onset_in_noise=TONE_ONSET_IN_NOISE,
            freq_pure=FREQ_PURE,
            gap_duration=GAP_DURATION,
            filter_order=FILTER_ORDER # Pass filter order
        )

        # --- Present Fixation & Play Sound ---
        fixation.draw()
        win.flip()
        core.wait(0.2 + random.uniform(0.0, 0.2))

        play_start_time = core.getTime()
        sd.play(stimulus_sound, FS, blocking=True)
        play_end_time = core.getTime()

        # --- Get Response ---
        fixation.draw()
        response_prompt_stim.draw()
        win.flip()

        response_key = None
        rt = None
        response_timer = core.Clock()
        event.clearEvents()

        while response_key is None:
            keys = event.getKeys(keyList=[KEY_INTERVAL_1, KEY_INTERVAL_2, KEY_QUIT], timeStamped=response_timer)
            if keys:
                key_pressed, rt = keys[0]
                if key_pressed == KEY_QUIT:
                    print("Quit key pressed. Exiting.")
                    return None # Signal to exit
                elif key_pressed in [KEY_INTERVAL_1, KEY_INTERVAL_2]:
                    response_key = key_pressed
            core.wait(0.001)

        # --- Process Response & Update Quest+ ---
        response_interval = 1 if response_key == KEY_INTERVAL_1 else 2
        is_correct = (response_interval == correct_interval)
        outcome = 'Correct' if is_correct else 'Incorrect'

        # Provide feedback during practice
        if is_practice:
             feedback_text = "Correct!" if is_correct else "Incorrect."
             feedback_stim = visual.TextStim(win, text=feedback_text, height=0.05,
                                             color='green' if is_correct else 'red')
             feedback_stim.draw()
             win.flip()
             core.wait(0.75)

        # Find the closest intensity in the stim_domain for the update
        intensity_for_update = find_closest_intensity(current_qp, presented_pn_ratio)

        # Update the *correct* Quest+ instance
        stim_update = dict(intensity=intensity_for_update)
        outcome_update = dict(response=outcome)
        current_qp.update(stim=stim_update, outcome=outcome_update)

        # --- Store Trial Data ---
        trial_info = {
            'trial_overall': trial_num_overall,
            'A_value_Hz': current_A_value, # Record the A value used
            'block': trial_info_base['block'],
            'presented_pn_ratio': presented_pn_ratio,
            'qp_update_intensity': intensity_for_update,
            # Tone power is relative to normalized noise power (1.0)
            'presented_tone_power': presented_pn_ratio * 1.0,
            'correct_interval': correct_interval,
            'response_interval': response_interval,
            'response_key': response_key,
            'is_correct': is_correct,
            'rt': rt,
            'timestamp': time.time()
        }
        trial_data_list.append(trial_info)
        # print(f"Trial {trial_num_overall} (A={current_A_value}Hz, Block: {trial_info_base['block']}): "
        #       f"PresRatio={presented_pn_ratio:.4f}, UpdRatio={intensity_for_update:.4f} Correct={is_correct}, Resp={response_interval}, RT={rt:.3f}")

        # --- Inter-Trial Interval (ITI) ---
        fixation.draw()
        win.flip()
        core.wait(0.2 + random.uniform(0.0, 0.2))

    return trial_data_list


# --- Main Execution ---
try:
    all_trial_data = []

    # --- Instructions ---
    message_stim.text = instructions_text
    message_stim.draw()
    win.flip()
    event.waitKeys()

    # --- Practice Block ---
    print("\n--- Starting Practice ---")
    practice_trials = create_trial_list(N_TRIALS_PRACTICE_PER_COND, A_VALUES_TO_TEST, 'practice')
    practice_data = run_experiment_block(practice_trials, is_practice=True)
    if practice_data is None:
        raise KeyboardInterrupt("User quit during practice.")
    all_trial_data.extend(practice_data)

    # --- Transition to Main Experiment ---
    message_stim.text = practice_over_text
    message_stim.draw()
    win.flip()
    event.waitKeys()

    # --- Main Experiment Block ---
    print("\n--- Starting Main Experiment ---")
    main_trials = create_trial_list(N_TRIALS_MAIN_PER_COND, A_VALUES_TO_TEST, 'main')
    # Adjust trial_overall numbers to continue from practice
    start_trial_num = N_TRIALS_PRACTICE + 1
    for trial in main_trials:
        trial['trial_overall'] += start_trial_num - 1 # Renumber based on overall count

    main_data = run_experiment_block(main_trials, is_practice=False)
    if main_data is None:
        raise KeyboardInterrupt("User quit during main experiment.")
    all_trial_data.extend(main_data)

    # --- End Experiment ---
    message_stim.text = end_text
    message_stim.draw()
    win.flip()

    # --- Save Data ---
    df = pd.DataFrame(all_trial_data)
    df.to_csv(data_filename, index=False)
    print(f"\nTrial data saved to: {data_filename}")

    # --- Get Final Estimates & Save ---
    final_results_summary = []
    final_results_summary.append(f"Participant: {participant_id}\n")
    final_results_summary.append(f"Bandwidth Parameter A Values Tested (Hz): {A_VALUES_TO_TEST}\n")
    final_results_summary.append(f"Total Main Trials per Condition: {N_TRIALS_MAIN_PER_COND}\n")

    print("\n--- Final Parameter Estimates (Mean of Posterior) ---")
    for A_val, qp_instance in qps.items():
        # # Save Quest+ state for this condition
        # qp_state_json = qp_instance.to_json()
        # qp_filename = qp_state_filenames[A_val]
        # with open(qp_filename, 'w') as f:
        #     f.write(qp_state_json)
        # print(f"Quest+ state for A={A_val}Hz saved to: {qp_filename}")

        # Get estimates
        final_estimate = qp_instance.param_estimate
        final_threshold_ratio = final_estimate['threshold']
        final_threshold_slope = final_estimate['slope']
        final_threshold_lapse = final_estimate['lapse_rate']

        result_str = (f"\nA Value (Hz): {A_val}\n"
                      f"  Estimated Threshold (P_Tone / P_NormNoise Ratio): {final_threshold_ratio:.6f}\n"
                      f"  Estimated Slope: {final_threshold_slope:.3f}\n"
                      f"  Estimated Lapse Rate: {final_threshold_lapse:.4f}\n")
        print(result_str)
        final_results_summary.append(result_str)

    # Save final results summary
    with open(results_filename, 'w') as f:
        f.writelines(final_results_summary)
    print(f"Final results summary saved to: {results_filename}")

    core.wait(1.5)

except Exception as e:
    logging.error(f"An error occurred: {e}")
    print(f"ERROR: An error occurred: {e}")
    # Save any data collected so far
    if 'all_trial_data' in locals() and all_trial_data:
         df = pd.DataFrame(all_trial_data)
         error_filename = os.path.join(DATA_FOLDER, f"{filename_base}_trials_ERROR.csv")
         df.to_csv(error_filename, index=False)
         print(f"Partial data saved to: {error_filename}")
    import traceback
    traceback.print_exc()

finally:
    # --- Cleanup ---
    win.close()
    sd.stop()
    core.quit()