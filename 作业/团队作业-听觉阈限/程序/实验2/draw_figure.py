# draw_figure_exp2.py

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines # For creating legend handles
import matplotlib.cm as cm        # For colormaps
import numpy as np
import glob # To find files
import re
import os
import shutil # For removing directory contents
from scipy.optimize import curve_fit



# --- Configuration for Experiment 2 (Bandwidth) ---

# !!! UPDATE THESE AS NEEDED FOR EXPERIMENT 2 !!!
# 1. Name of the independent variable column in CSV and descriptive name for plots
EXP2_CONDITION_VARIABLE_NAME = 'A_value_Hz' # Column name in CSV
EXP2_CONDITION_LABEL = '滤波半带宽 A (Hz)'    # Label for axes/legends, e.g., 'Frequency (Hz)', 'Duration (ms)'

# 2. Expected values of the independent variable
EXP2_EXPECTED_CONDITIONS = [100, 300, 900] # e.g., [1000, 4000], [100, 500], ['TypeA', 'TypeB']

# 3. File naming patterns for Experiment 2 data
EXP2_CSV_PATTERN = '*_bandwidth_trials.csv'
EXP2_TXT_PATTERN = '*_bandwidth_results.txt'
EXP2_FIGURE_PREFIX = 'EXP2_BW' # Prefix for saved figure filenames

# 4. Regex to find the condition identifier in the TXT file (adjust based on actual text)
# Example: looks for "A Value (Hz): 100" (case-insensitive)
EXP2_TXT_REGEX_CONDITION = r"A Value \(Hz\):\s*([\d\.]+)" # Adjust this regex!

# 5. Directory containing participant data
data_folder = 'data_exp2' # Assumes Exp 2 data is in 'data_exp2'
# Subdirectory for saving figures (will be cleared)
figure_folder = os.path.join(data_folder, f'figure_{EXP2_FIGURE_PREFIX}')
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = ['DejaVu Sans', 'SimHei']
axis_linewidth = 2
# --- End of Exp 2 Configuration ---


# --- Plotting Style Configuration ---
plt.rcParams['font.size'] = 18
# Ensure SimHei or another appropriate font is installed for Chinese characters if needed
plt.rcParams['font.family'] = ['DejaVu Sans', 'SimHei']
# ---

# --- Helper Functions ---
# Weibull function definition (same as before)
def weibull_2afc_fit_lapse(x, threshold, slope, lapse):
    """Weibull function for 2AFC (guess rate = 0.5), fits lapse."""
    g = 0.5
    x = np.maximum(x, 1e-9)
    threshold = np.maximum(threshold, 1e-9)
    slope = np.maximum(slope, 0.1)
    lapse = np.clip(lapse, 0, 0.2)
    weibull_term = np.exp(-(x / threshold)**slope)
    return g + (1.0 - g - lapse) * (1.0 - weibull_term)

# --- Prepare Figure Directory ---
if os.path.exists(figure_folder):
    print(f"Clearing contents of existing figure directory: {figure_folder}")
    for filename in os.listdir(figure_folder):
        file_path = os.path.join(figure_folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')
else:
    print(f"Creating figure directory: {figure_folder}")
    os.makedirs(figure_folder)


# --- Part 1: Generate Individual Psychometric Function Plots for Exp 2 ---

print(f"\n--- Generating Individual Psychometric Plots ({EXP2_FIGURE_PREFIX}) ---")
csv_file_pattern = os.path.join(data_folder, EXP2_CSV_PATTERN)
participant_csv_files = sorted(glob.glob(csv_file_pattern))
count = 0 # For sequential participant IDs in plots

if not participant_csv_files:
    print(f"Error: No participant CSV files found matching '{csv_file_pattern}'")
else:
    print(f"Found {len(participant_csv_files)} participant CSV files.")

    # Define colors using a matplotlib colormap based on expected conditions
    cmap = cm.get_cmap('viridis', len(EXP2_EXPECTED_CONDITIONS))
    # condition_colors = {level: cmap(i) for i, level in enumerate(sorted(EXP2_EXPECTED_CONDITIONS))}
    # Or define fixed colors if preferred:
    condition_colors = {100:'#5F97D2', 300:'#7BC9A3', 900:'#EF7A6D'}

    for csv_filename in participant_csv_files:
        # Extract ID from filename, assuming pattern like 'participantID_bandwidth_trials.csv'
        try:
            participant_id_from_file = os.path.basename(csv_filename).split('_bandwidth')[0]
        except IndexError:
            participant_id_from_file = f"Unknown_{count}" # Fallback ID
        participant_display_id = f"Participant {count + 1}" # Use sequential ID for titles/filenames
        count += 1
        print(f" Processing {participant_display_id} (File: {os.path.basename(csv_filename)})")


        # --- Load Data ---
        try:
            df = pd.read_csv(csv_filename)
        except Exception as e:
            print(f"  Error loading {csv_filename}: {e}. Skipping.")
            continue

        # --- Prepare Data ---
        intensity_col = 'qp_update_intensity' if 'qp_update_intensity' in df.columns else 'presented_pn_ratio'
        if EXP2_CONDITION_VARIABLE_NAME not in df.columns:
             print(f"  Error: Condition column '{EXP2_CONDITION_VARIABLE_NAME}' not found in {csv_filename}. Skipping.")
             continue
        if 'block' in df.columns:
            df_main = df[df['block'] == 'main'].copy()
            if df_main.empty and not df.empty:
                print("  Warning: No 'main' block trials found. Using all trials.")
                df_main = df.copy()
            elif df_main.empty:
                print("  Warning: No data found for fitting. Skipping plot.")
                continue
        else:
            df_main = df.copy()
        try:
            df_main['is_correct'] = pd.to_numeric(df_main['is_correct'], errors='coerce')
            df_main.dropna(subset=['is_correct'], inplace=True)
        except Exception as e:
            print(f"  Error cleaning 'is_correct' column: {e}. Skipping plot.")
            continue
        if df_main.empty:
            print("  Warning: No valid main trial data found for fitting. Skipping plot.")
            continue
        try:
            psych_data = df_main.groupby([EXP2_CONDITION_VARIABLE_NAME, intensity_col])['is_correct']\
                                .agg(['mean', 'count']).reset_index()\
                                .rename(columns={'mean': 'proportion_correct', 'count': 'n_trials', intensity_col: 'pn_ratio'})
        except Exception as e:
            print(f"  Error grouping data by '{EXP2_CONDITION_VARIABLE_NAME}': {e}. Skipping plot.")
            continue

        # --- Create Plot for this Participant ---
        fig_ind, ax_ind = plt.subplots(figsize=(12, 9))
        participant_conditions = sorted(psych_data[EXP2_CONDITION_VARIABLE_NAME].unique())
        fitted_params = {}

        # --- Plotting Loop for Conditions ---
        for condition_value in participant_conditions:
            subset = psych_data[psych_data[EXP2_CONDITION_VARIABLE_NAME] == condition_value].sort_values('pn_ratio')
            color = condition_colors.get(condition_value, 'gray')
            if subset.empty:
                continue
            condition_label_str = f"{condition_value:g} Hz"

            # 1. Scatter plot of raw proportions
            point_sizes = np.sqrt(np.maximum(subset['n_trials'], 1)) * 15
            ax_ind.scatter(subset['pn_ratio'], subset['proportion_correct'],
                           s=point_sizes,
                           color=color, alpha=0.6, edgecolor='black', linewidth=0.5,
                           label=f'A = {condition_label_str} (Data)')

            # 2. Fit and plot Weibull
            if len(subset) >= 3:
                try:
                    # --- CORRECTED BOUNDS DEFINITION ---
                    # Initial guesses
                    p0_threshold_guess = subset['pn_ratio'].median()
                    if p0_threshold_guess <= 1e-6 or not np.isfinite(p0_threshold_guess):
                         p0_threshold_guess = 0.1
                    p0 = [p0_threshold_guess, 5.0, 0.02] # T, S, L guess

                    # Define lower bounds first
                    lower_bounds = [1e-6, 0.1, 0.0]

                    # Calculate a reasonable upper bound for threshold based on data
                    # Make sure it's greater than the lower bound
                    upper_threshold_bound = max(subset['pn_ratio']) * 2
                    if upper_threshold_bound <= lower_bounds[0]:
                        upper_threshold_bound = lower_bounds[0] * 100 # Fallback

                    # Define upper bounds
                    upper_bounds = [upper_threshold_bound, 50.0, 0.2]

                    # Now define the full bounds tuple
                    bounds = (lower_bounds, upper_bounds)
                    # --- END OF CORRECTION ---

                    weights = 1.0 / np.sqrt(np.maximum(subset['n_trials'].values, 1))

                    params, covariance = curve_fit(weibull_2afc_fit_lapse,
                                                   subset['pn_ratio'].values,
                                                   subset['proportion_correct'].values,
                                                   p0=p0, bounds=bounds, sigma=weights,
                                                   absolute_sigma=False, maxfev=10000, method='trf')

                    fitted_params[condition_value] = params

                    min_x = max(1e-6, subset['pn_ratio'].min() * 0.8)
                    max_x = subset['pn_ratio'].max() * 1.2
                    x_smooth = np.linspace(min_x, max_x, 200)
                    y_fit = weibull_2afc_fit_lapse(x_smooth, *params)
                    ax_ind.plot(x_smooth, y_fit, color=color, linestyle='-', linewidth=2,
                                label=f'A = {condition_label_str} (Fit)')
                except Exception as e:
                    # Print specific error during fitting
                    print(f"  ERROR during curve_fit for Condition A = {condition_value}: {e}")
                    # Optionally import traceback and print full traceback
                    # import traceback
                    # traceback.print_exc()
            else:
                print(f"  Skipping fit for Condition A = {condition_value}: Not enough data points.")

        # --- Final Touches on Individual Plot ---
        ax_ind.set_xlabel("信噪比 (P_Tone / P_NormNoise Ratio)")
        ax_ind.set_ylabel("正确率")
        ax_ind.set_title(f"心理物理函数: {participant_display_id}")
        ax_ind.set_ylim(0.45, 1.05)
        # ax_ind.set_xlim(right=0.06) # Removed fixed xlim, let it be dynamic based on data
        if not psych_data.empty:
             ax_ind.axhline(0.5, color='grey', linestyle=':', alpha=0.7, linewidth=1.5, label='Chance Level')
        handles, labels = ax_ind.get_legend_handles_labels()
        from collections import OrderedDict
        unique_labels = OrderedDict(zip(labels, handles))
        ax_ind.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize='medium', title="Condition")
        ax_ind.spines['top'].set_visible(False)
        ax_ind.spines['right'].set_visible(False)
        ax_ind.spines['left'].set_linewidth(1.5)
        ax_ind.spines['bottom'].set_linewidth(1.5)
        plt.tight_layout()
        individual_plot_filename = os.path.join(figure_folder, f"{EXP2_FIGURE_PREFIX}_{participant_display_id.replace(' ','_')}_psychometric.png")
        try:
             fig_ind.savefig(individual_plot_filename, dpi=300)
             print(f"  Saved: {individual_plot_filename}")
        except Exception as e:
             print(f"  Error saving figure {individual_plot_filename}: {e}")
        plt.close(fig_ind)


# --- Part 2: Generate Combined Threshold Summary Plot for Exp 2 ---

print(f"\n--- Generating Combined Threshold Summary Plot ({EXP2_FIGURE_PREFIX}) ---")
txt_file_pattern = os.path.join(data_folder, EXP2_TXT_PATTERN)
participant_txt_files = sorted(glob.glob(txt_file_pattern))
num_participants = len(participant_txt_files)

if not participant_txt_files:
    print(f"Error: No participant TXT results files found matching '{txt_file_pattern}'")
else:
    print(f"Found {num_participants} participant results files.")

    # --- Data Extraction from TXT ---
    all_threshold_data = []
    all_conditions_found = set() # Store unique conditions found (A values)
    for file_idx, txt_filename in enumerate(participant_txt_files):
        try:
            participant_id_from_file = os.path.basename(txt_filename).split('_bandwidth')[0]
        except IndexError:
            participant_id_from_file = f"Unknown_{file_idx}"
        participant_display_id = f"Participant {file_idx + 1}"
        current_condition_value = None
        try:
            with open(txt_filename, 'r') as f:
                for line in f:
                    match_condition = re.search(EXP2_TXT_REGEX_CONDITION, line, re.IGNORECASE)
                    if match_condition:
                        try:
                            current_condition_value = float(match_condition.group(1))
                            all_conditions_found.add(current_condition_value)
                        except ValueError:
                             print(f"  Warning: Could not parse condition value as float in {txt_filename}, line: {line.strip()}")
                             current_condition_value = None
                    if current_condition_value is not None:
                        match_thresh = re.search(r"Estimated Threshold.*?Ratio\):\s*([\d\.]+)", line, re.IGNORECASE)
                        if match_thresh:
                            threshold = float(match_thresh.group(1))
                            all_threshold_data.append({
                                'participant_id': participant_display_id,
                                'participant_idx': file_idx,
                                EXP2_CONDITION_VARIABLE_NAME: current_condition_value,
                                'threshold': threshold
                            })
                            current_condition_value = None
        except Exception as e:
            print(f"  Warning: Error parsing {txt_filename}: {e}. Skipping participant.")

    if not all_threshold_data:
        print("Error: No threshold data could be extracted. Cannot create summary plot.")
    else:
        threshold_df = pd.DataFrame(all_threshold_data)
        print(f"Extracted {len(threshold_df)} threshold data points for summary plot.")

        # --- Plotting Setup for Summary ---
        fig_sum, ax_sum = plt.subplots(figsize=(16, 12))
        sorted_conditions = sorted(list(all_conditions_found))
        if not sorted_conditions:
             print("Error: No conditions (A values) identified from data.")
             exit()
        cmap_sum = cm.get_cmap('viridis', len(sorted_conditions))
        summary_condition_colors = {level: cmap_sum(i) for i, level in enumerate(sorted_conditions)}
        min_alpha = 0.3
        max_alpha = 0.6
        if num_participants > 1:
            alpha_values = np.linspace(min_alpha, max_alpha, num_participants)
        else:
            alpha_values = np.array([max_alpha])

        # --- Plotting Individual Participant Thresholds ---
        grouped_by_participant = threshold_df.groupby('participant_id')
        for participant_id, group in grouped_by_participant:
            group_sorted = group.sort_values(EXP2_CONDITION_VARIABLE_NAME)
            participant_idx = group_sorted['participant_idx'].iloc[0]
            participant_alpha = alpha_values[participant_idx]
            if len(group_sorted) >= 1:
                x_coords = np.log10(group_sorted[EXP2_CONDITION_VARIABLE_NAME])
                ax_sum.plot(x_coords, group_sorted['threshold'],
                           marker='o', linestyle='-', linewidth=3.0, markersize=20,
                           color="#5F97D2", alpha=participant_alpha,
                           markeredgecolor='none', label='_nolegend_')

        # --- Plotting Average Data ---
        average_thresholds = threshold_df.groupby(EXP2_CONDITION_VARIABLE_NAME)['threshold'].mean().reset_index()
        average_thresholds = average_thresholds.sort_values(EXP2_CONDITION_VARIABLE_NAME)
        if not average_thresholds.empty:
            for _, row in average_thresholds.iterrows():
                condition_val = np.log10(row[EXP2_CONDITION_VARIABLE_NAME])
                condition_color = summary_condition_colors.get(condition_val, 'black')
                ax_sum.plot(condition_val, row['threshold'],
                           marker='^', linestyle='None', color="#EF7A6D", alpha=1.0,
                           markersize=18, markeredgecolor='none', linewidth=1,
                           zorder=5, label='_nolegend_')
        else:
            print("Warning: Could not calculate average thresholds.")

        # --- Final Touches on Summary Plot ---
        ax_sum.set_xlabel(EXP2_CONDITION_LABEL,fontsize = 24)
        ax_sum.set_ylabel("信噪比 (P/N Ratio)", fontsize=24)
        sorted_conditions = np.log10(sorted_conditions)
        # ax_sum.set_xscale('log')
        ax_sum.set_xticks(sorted_conditions)
        ax_sum.set_xticklabels([f"{val:g}" for val in np.power(10,sorted_conditions)])
        if len(sorted_conditions) > 1:
            x_range = max(sorted_conditions) - min(sorted_conditions)
            #  ax_sum.set_xlim(min(sorted_conditions) - 0.1 * x_range,
            #                  max(sorted_conditions) + 0.1 * x_range)
            x_pad = 0.2
            ax_sum.set_xlim(np.log10(100)-x_pad, np.log10(900)+x_pad)
        elif len(sorted_conditions) == 1:
             ax_sum.set_xlim(sorted_conditions[0] * 0.8, sorted_conditions[0] * 1.2)
        if not threshold_df.empty:
             min_thresh_overall = threshold_df['threshold'].min()
             max_thresh_overall = threshold_df['threshold'].max()
             y_pad = (max_thresh_overall - min_thresh_overall) * 0.1
             ax_sum.set_ylim(min_thresh_overall - y_pad, 1.22)
        sorted_conditions = sorted(list(all_conditions_found))
        

        # --- Create Custom Legend for Summary Plot---
        summary_legend_handles = []
        summary_legend_handles.append(mlines.Line2D([], [], color="#5F97D2", marker='o', linestyle='-',
                                            linewidth=2.0, alpha=np.mean(alpha_values),
                                            markersize=10, label='各被试阈值'))
        for condition_value in sorted_conditions:
             color = "#EF7A6D"
             condition_label_str = f"{condition_value:g} Hz"
             summary_legend_handles.append(mlines.Line2D([], [], color=color, marker='^', linestyle='None',
                                                markersize=10, markeredgecolor='black',
                                                label=f'A = {condition_label_str}'))
        # ax_sum.legend(handles=summary_legend_handles, loc=(0.1,0.1), fontsize=18)
        ax_sum.tick_params(axis="both",which = "major",width = axis_linewidth, length = 6,pad = 8)
        ax_sum.spines['top'].set_visible(False)
        ax_sum.spines['right'].set_visible(False)
        ax_sum.spines['left'].set_linewidth(axis_linewidth)
        ax_sum.spines['bottom'].set_linewidth(axis_linewidth)

        # add significance text
        line_margin = 0.01
        x_1_1 = np.log10(100)+line_margin
        x_1_2 = np.log10(300)-line_margin
        x_2_1 = np.log10(300)+line_margin
        x_2_2 = np.log10(900)-line_margin
        y_line_1 = 0.95
        y_line_2 = 1.04
        text_height = 0.02
        ax_sum.hlines(y_line_1,x_1_1,x_1_2,color="k", linewidth = axis_linewidth)
        ax_sum.hlines(y_line_1,x_2_1,x_2_2,color="k", linewidth = axis_linewidth)
        ax_sum.hlines(y_line_2,x_1_1,x_2_2,color="k", linewidth = axis_linewidth)
        ax_sum.text(s = '$p = .183$',x = (x_1_1+x_1_2)/2,y = y_line_1+text_height,ha='center',fontsize = 24)
        ax_sum.text(s = '$p = .059$',x = (x_2_1+x_2_2)/2,y = y_line_1+text_height,ha='center',fontsize = 24)
        ax_sum.text(s = '$**$',x = (x_1_1+x_2_2)/2,y = y_line_2+text_height,ha='center',fontsize = 24)


        # --- Save Summary Plot & Data ---
        plt.tight_layout()
        combined_plot_filename = os.path.join(figure_folder, f"{EXP2_FIGURE_PREFIX}_ALL_participants_threshold_summary.png")
        csv_summary_filename = os.path.join(figure_folder, f"{EXP2_FIGURE_PREFIX}_ALL_participants_threshold_summary.csv")
        csv_average_filename = os.path.join(figure_folder, f"{EXP2_FIGURE_PREFIX}_ALL_participants_average_thresholds.csv")
        try:
            fig_sum.savefig(combined_plot_filename, dpi=300)
            print(f"\nCombined threshold summary plot saved to {combined_plot_filename}")
            threshold_df.to_csv(csv_summary_filename, index=False)
            average_thresholds.to_csv(csv_average_filename, index=False)
            print(f"Summary data saved to {csv_summary_filename}")
            print(f"Average data saved to {csv_average_filename}")
        except Exception as e:
            print(f"Error saving summary figure/data {combined_plot_filename}: {e}")
        plt.close(fig_sum)

print("\nProcessing finished.")