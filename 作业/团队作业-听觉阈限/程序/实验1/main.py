import numpy as np
import pandas as pd
import sounddevice as sd
from psychopy import core, visual, event, gui, logging
import random
import os
import time
import json # Needed for saving/loading QP state
from questplus.qp import QuestPlus
from questplus import _constants

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

# --- Helper Function ---
def find_closest_intensity(qp_instance, presented_intensity):
    """Finds the intensity in the Quest+ stim_domain closest to the presented one."""
    known_intensities = qp_instance.stim_domain['intensity']
    # Add a small epsilon to prevent potential floating point issues if presented_intensity is exactly 0
    # presented_intensity = max(presented_intensity, 1e-9)
    idx = np.argmin(np.abs(known_intensities - presented_intensity))
    closest_intensity = known_intensities[idx]
    # Optional: Warn if the difference is large, indicating a potential issue
    # if np.abs(closest_intensity - presented_intensity) > (known_intensities[1] - known_intensities[0]):
    #     print(f"Warning: Presented intensity {presented_intensity:.4f} mapped to {closest_intensity:.4f}")
    return closest_intensity

# --- Weibull Function Definition for Plotting ---
def weibull_psychometric_function(x, threshold, slope, lower_asymptote, lapse_rate):
    """
    Calculates the Weibull psychometric function value.
    P(x) = g + (1 - g - l) * (1 - exp(-(x / t)**s))
    """
    # Ensure inputs are numpy arrays for vectorized operations
    x = np.asarray(x)
    # Prevent division by zero or log(0) if threshold is near zero or x is zero
    # We clamp x/t to a small positive value if the base is zero or negative.
    # Also handle cases where threshold might be estimated as <= 0, although unlikely here.
    t_safe = max(threshold, 1e-9) # Prevent division by zero threshold
    base = x / t_safe
    # Calculate exponent safely
    try:
        # Use maximum to avoid issues with 0 or negative base for potentially fractional slopes
        exponent = np.power(np.maximum(base, 1e-12), slope)
    except ValueError: # Handle potential complex numbers if base is negative and slope is fractional
         # This case *should* be avoided by maximum(base, 1e-12), but as a fallback:
         print(f"Warning: ValueError during np.power calculation. Base: {base}, Slope: {slope}. Setting exponent to large value.")
         exponent = np.inf # Or some other indicator of failure/saturation

    prob = lower_asymptote + (1.0 - lower_asymptote - lapse_rate) * (1.0 - np.exp(-exponent))
    return prob

# --- Experiment Configuration ---
# Participant Info GUI
exp_info = {'participant': 'test'} # Removed noise_condition selection
dlg = gui.DlgFromDict(dictionary=exp_info, title='Pure Tone Threshold Experiment')
if not dlg.OK:
    core.quit()
participant_id = exp_info['participant']

# Stimulus Parameters
FS = 16000
DURATION_NOISE = 1.0
DURATION_PURE = 0.100
TONE_ONSET_IN_NOISE = 0.5
FREQ_PURE = 1000
GAP_DURATION = 0.5

# Define the two noise power levels
NOISE_POWERS_TO_TEST = [1.44, 0.64]

N_TRIALS_PRACTICE_PER_COND = 3 # Practice trials *per noise level*
N_TRIALS_MAIN_PER_COND = 32  # Main trials *per noise level*
N_TRIALS_PRACTICE = N_TRIALS_PRACTICE_PER_COND * len(NOISE_POWERS_TO_TEST)
N_TRIALS_MAIN = N_TRIALS_MAIN_PER_COND * len(NOISE_POWERS_TO_TEST)
BREAK_INTERVAL = 32 # Show break message every N total trials

# Response Keys
KEY_INTERVAL_1 = '1'
KEY_INTERVAL_2 = '2'
KEY_QUIT = 'q'

# Data Saving
DATA_FOLDER = 'data'
if not os.path.exists(DATA_FOLDER):
    os.makedirs(DATA_FOLDER)
filename_base = f"{participant_id}_mixed_noise" # Updated base name
data_filename = os.path.join(DATA_FOLDER, f"{filename_base}_trials.csv")
results_filename = os.path.join(DATA_FOLDER, f"{filename_base}_results.txt")
# We will save QP states separately for each noise level
qp_state_filenames = {
    noise_level: os.path.join(DATA_FOLDER, f"{filename_base}_noise{noise_level:.2f}_questplus_state.json")
    for noise_level in NOISE_POWERS_TO_TEST
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


# --- Quest+ Configuration (Common Settings) ---
# Define the domains for stimuli and parameters (same for both conditions initially)
stim_intensities_ratio = np.linspace(0.001, 0.3, 600) # P/N ratio
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

# --- Initialize QuestPlus instances for each noise level ---
qps = {}
for noise_level in NOISE_POWERS_TO_TEST:
    qps[noise_level] = QuestPlus(
        stim_domain=stim_domain,
        param_domain=param_domain,
        outcome_domain=outcome_domain,
        func='weibull',
        stim_scale='linear',
        param_estimation_method='mean', # Or 'mode'
        stim_selection_method='min_entropy'
    )
    print(f"Quest+ Initialized for Noise Power: {noise_level:.2f}")


# --- PsychoPy Setup ---
win = visual.Window(
    size=[800, 600], fullscr=True, screen=0,
    winType='pyglet', allowGUI=True, allowStencil=False,
    monitor='testMonitor', color=[0.8, 0.8, 0.8], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    units='height')

# Instructions and Feedback Text (slightly updated)
instructions_text = """
Welcome!

In each trial, you will hear two sounds separated by a short silence.
Each sound contains noise, and the noise level might change between trials.
One of the two sounds will ALSO contain a faint pure tone.

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

# --- Audio Generation Function (remains the same) ---
def generate_trial_sound(pn_ratio, noise_power, correct_interval, fs,
                         duration_noise, duration_pure, tone_onset_in_noise, freq_pure, gap_duration):
    """Generates the 2AFC sound stimulus for one trial."""
    pn_ratio = max(0, pn_ratio)
    pure_power = pn_ratio * noise_power
    amplitude_pure = np.sqrt(2 * pure_power)
    t_noise = np.arange(0, duration_noise, 1/fs)
    t_pure = np.arange(0, duration_pure, 1/fs)
    n_samples_noise = len(t_noise)
    n_samples_pure = len(t_pure)
    n_samples_gap = int(gap_duration * fs)
    n_samples_onset = int(tone_onset_in_noise * fs)
    n_samples_before_tone = n_samples_onset
    # Ensure we don't have negative samples after tone if onset+duration > noise duration
    n_samples_after_tone = max(0, n_samples_noise - n_samples_before_tone - n_samples_pure)
    # Adjust pure tone length if it exceeds remaining time
    if n_samples_before_tone + n_samples_pure > n_samples_noise:
         n_samples_pure = n_samples_noise - n_samples_before_tone
         pure_tone_signal = amplitude_pure * np.sin(2 * np.pi * freq_pure * t_pure[:n_samples_pure])
         n_samples_after_tone = 0 # No samples left after
         print(f"Warning: Tone duration ({duration_pure}s) + onset ({tone_onset_in_noise}s) exceeds noise duration ({duration_noise}s). Truncating tone.")
    else:
         pure_tone_signal = amplitude_pure * np.sin(2 * np.pi * freq_pure * t_pure)

    noise1 = np.sqrt(noise_power) * np.random.randn(n_samples_noise)
    noise2 = np.sqrt(noise_power) * np.random.randn(n_samples_noise)
    padding_before = np.zeros(n_samples_before_tone)
    padding_after = np.zeros(n_samples_after_tone)

    pure_tone_signal_ramped = apply_cosine_ramp(pure_tone_signal, fs, fade_time=0.02)
    tone_padded = np.concatenate([padding_before, pure_tone_signal_ramped, padding_after])
    # Final check for length consistency due to potential rounding
    if len(tone_padded) < n_samples_noise:
        tone_padded = np.concatenate([tone_padded, np.zeros(n_samples_noise - len(tone_padded))])
    elif len(tone_padded) > n_samples_noise:
        tone_padded = tone_padded[:n_samples_noise]

    if correct_interval == 1:
        interval1 = noise1 + tone_padded
        interval2 = noise2
    else:
        interval1 = noise1
        interval2 = noise2 + tone_padded

    silence_gap_signal = np.zeros(n_samples_gap)
    combined_signal = np.concatenate([interval1, silence_gap_signal, interval2])
    stereo_signal = np.vstack([combined_signal, combined_signal]).T
    return stereo_signal

# --- Create Trial List ---
def create_trial_list(n_trials_per_cond, noise_levels, block_type):
    """Creates a list of trials with randomly interleaved noise levels."""
    trial_conditions = []
    for noise_level in noise_levels:
        trial_conditions.extend([noise_level] * n_trials_per_cond)
    random.shuffle(trial_conditions) # Interleave the conditions randomly
    # Add trial number and block type
    trials = [{'trial_overall': i + 1, 'noise_level': cond, 'block': block_type}
               for i, cond in enumerate(trial_conditions)]
    return trials


# --- Run Experiment Function (Modified for Interleaving) ---
def run_experiment_block(trial_list, is_practice):
    """Runs a block of trials with interleaved conditions."""
    trial_data_list = []
    global qps # Allow access to the global qps dictionary

    # Fixed high P/N ratio for early practice trials
    practice_high_pn_ratio = 0.3
    practice_trials_done_per_cond = {level: 0 for level in NOISE_POWERS_TO_TEST}
    initial_practice_phase_count = 1 # How many initial trials use fixed high ratio

    for trial_info_base in trial_list:
        trial_num_overall = trial_info_base['trial_overall']
        current_noise_level = trial_info_base['noise_level']
        current_qp = qps[current_noise_level] # Select the QP instance for this noise level

        # --- Break Time --- (Based on overall trial number)
        if not is_practice and ((trial_num_overall-N_TRIALS_PRACTICE) % BREAK_INTERVAL == 1 and trial_num_overall != (1+N_TRIALS_PRACTICE)):
            message_stim.text = break_text
            message_stim.draw()
            win.flip()
            event.waitKeys(keyList=['space', 'escape']) # Use space to continue, escape to potentially abort
            if 'escape' in event.getKeys():
                 print("Escape pressed during break. Exiting.")
                 return None # Signal exit
            event.clearEvents() # Clear buffer after wait

        # --- Get Stimulus Intensity ---
        # Use fixed high intensity for the very first practice trial(s) per condition
        if is_practice and practice_trials_done_per_cond[current_noise_level] < initial_practice_phase_count:
             presented_pn_ratio = practice_high_pn_ratio
             practice_trials_done_per_cond[current_noise_level] += 1
        else:
             # Get intensity from the *correct* Quest+ instance
             next_stim_dict = current_qp.next_stim
             presented_pn_ratio = next_stim_dict['intensity']
             presented_pn_ratio = max(0.0001, presented_pn_ratio) # Clamp if needed


        # --- Prepare Trial ---
        correct_interval = random.choice([1, 2])
        stimulus_sound = generate_trial_sound(
            pn_ratio=presented_pn_ratio,
            noise_power=current_noise_level, # Use the noise level for this trial
            correct_interval=correct_interval,
            fs=FS,
            duration_noise=DURATION_NOISE,
            duration_pure=DURATION_PURE,
            tone_onset_in_noise=TONE_ONSET_IN_NOISE,
            freq_pure=FREQ_PURE,
            gap_duration=GAP_DURATION
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

        # *** IMPORTANT FIX FOR KEYERROR ***
        # Find the closest intensity in the stim_domain to the one actually presented
        intensity_for_update = find_closest_intensity(current_qp, presented_pn_ratio)

        # Update the *correct* Quest+ instance using the *closest known intensity*
        stim_update = dict(intensity=intensity_for_update)
        outcome_update = dict(response=outcome)
        current_qp.update(stim=stim_update, outcome=outcome_update)

        # --- Store Trial Data ---
        trial_info = {
            'trial_overall': trial_num_overall,
            'noise_level': current_noise_level, # Record the noise level used
            'block': trial_info_base['block'],
            'presented_pn_ratio': presented_pn_ratio,
            # Intensity used for QP update might differ slightly due to 'find_closest'
            'qp_update_intensity': intensity_for_update,
            'presented_tone_power': presented_pn_ratio * current_noise_level,
            'correct_interval': correct_interval,
            'response_interval': response_interval,
            'response_key': response_key,
            'is_correct': is_correct,
            'rt': rt,
            'timestamp': time.time()
        }
        trial_data_list.append(trial_info)
        print(f"Trial {trial_num_overall} (Noise: {current_noise_level:.2f}, Block: {trial_info_base['block']}): "
              f"PresRatio={presented_pn_ratio:.4f}, UpdRatio={intensity_for_update:.4f} Correct={is_correct}, Resp={response_interval}, RT={rt:.3f}")

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
    practice_trials = create_trial_list(N_TRIALS_PRACTICE_PER_COND, NOISE_POWERS_TO_TEST, 'practice')
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
    main_trials = create_trial_list(N_TRIALS_MAIN_PER_COND, NOISE_POWERS_TO_TEST, 'main')
    # Adjust trial_overall numbers to continue from practice
    start_trial_num = N_TRIALS_PRACTICE + 1
    for trial in main_trials:
        trial['trial_overall'] += start_trial_num -1 # Renumber based on overall count

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
    final_results_summary.append(f"Noise Levels Tested: {NOISE_POWERS_TO_TEST}\n")
    final_results_summary.append(f"Total Main Trials per Condition: {N_TRIALS_MAIN_PER_COND}\n")

    print("\n--- Final Parameter Estimates (Mean of Posterior) ---")
    for noise_level, qp_instance in qps.items():
        # # Save Quest+ state for this condition
        # qp_state_json = qp_instance.to_json()
        # qp_filename = qp_state_filenames[noise_level]
        # with open(qp_filename, 'w') as f:
        #     f.write(qp_state_json)
        # print(f"Quest+ state for noise={noise_level:.2f} saved to: {qp_filename}")

        # Get estimates
        final_estimate = qp_instance.param_estimate
        final_threshold_ratio = final_estimate['threshold']
        final_threshold_slope = final_estimate['slope']
        final_threshold_lapse = final_estimate['lapse_rate']

        result_str = (f"\nNoise Level: {noise_level:.2f}\n"
                      f"  Estimated Threshold (P/N Ratio): {final_threshold_ratio:.6f}\n"
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