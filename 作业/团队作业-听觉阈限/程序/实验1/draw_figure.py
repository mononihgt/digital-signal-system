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

# --- Configuration ---
# Directory containing participant data (both CSV and TXT)
data_folder = 'data'
# Subdirectory for saving figures
figure_folder = os.path.join(data_folder, 'figure')
count = 0
# Define expected noise levels (update if different)
EXPECTED_NOISE_LEVELS = [0.64, 1.0]
# Define plotting style
# plt.style.use('seaborn-v0_8-whitegrid') # Use a matplotlib style
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = ['DejaVu Sans', 'SimHei']
axis_linewidth = 2
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


# --- Part 1: Generate Individual Psychometric Function Plots ---

print("\n--- Generating Individual Psychometric Plots ---")
csv_file_pattern = os.path.join(data_folder, '*_mixed_noise_trials.csv')
participant_csv_files = sorted(glob.glob(csv_file_pattern))

if not participant_csv_files:
    print(f"Error: No participant CSV files found matching '{csv_file_pattern}'")
else:
    print(f"Found {len(participant_csv_files)} participant CSV files.")

    # Define colors using a matplotlib colormap
    cmap = cm.get_cmap('viridis', len(EXPECTED_NOISE_LEVELS))
    # noise_colors = {level: cmap(i) for i, level in enumerate(sorted(EXPECTED_NOISE_LEVELS))}
    noise_colors = {0.64:'#5F97D2', 1.44:'#EF7A6D'} # Fixed colors as requested

    for csv_filename in participant_csv_files:
        participant_id = os.path.basename(csv_filename).split('_mixed_noise')[0]
        print(f" Processing Participant: {participant_id}")

        participant_id = "participant " + str(count+1)
        count += 1
        # --- Load Data ---
        try:
            df = pd.read_csv(csv_filename)
        except Exception as e:
            print(f"  Error loading {csv_filename}: {e}. Skipping.")
            continue

        # --- Prepare Data ---
        intensity_col = 'qp_update_intensity' if 'qp_update_intensity' in df.columns else 'presented_pn_ratio'
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
            psych_data = df_main.groupby(['noise_level', intensity_col])['is_correct']\
                                .agg(['mean', 'count']).reset_index()\
                                .rename(columns={'mean': 'proportion_correct', 'count': 'n_trials', intensity_col: 'pn_ratio'})
        except Exception as e:
            print(f"  Error grouping data: {e}. Skipping plot.")
            continue

        # --- Create Plot for this Participant ---
        fig_ind, ax_ind = plt.subplots(figsize=(12, 9)) # Use plt.subplots

        participant_noise_levels = sorted(psych_data['noise_level'].unique())
        fitted_params = {}

        # --- Plotting Loop for Noise Levels ---
        for noise_level in participant_noise_levels:
            subset = psych_data[psych_data['noise_level'] == noise_level].sort_values('pn_ratio')
            color = noise_colors.get(noise_level, 'gray') # Assign color

            if subset.empty:
                continue

            # 1. Scatter plot of raw proportions
            point_sizes = np.sqrt(np.maximum(subset['n_trials'], 1)) * 15
            ax_ind.scatter(subset['pn_ratio'], subset['proportion_correct'],
                           s=point_sizes,
                           color=color, alpha=0.6, edgecolor='black', linewidth=0.5,
                           label=f'Noise {noise_level:.2f} (Data)')

            # 2. Fit and plot Weibull
            if len(subset) >= 3: # Need enough points to fit
                try:
                    p0 = [subset['pn_ratio'].median(), 5.0, 0.02]
                    bounds = ([1e-6, 0.1, 0.0], [max(subset['pn_ratio'])*2, 50.0, 0.2])
                    weights = 1.0 / np.sqrt(np.maximum(subset['n_trials'].values, 1))

                    params, covariance = curve_fit(weibull_2afc_fit_lapse,
                                                   subset['pn_ratio'].values,
                                                   subset['proportion_correct'].values,
                                                   p0=p0, bounds=bounds, sigma=weights,
                                                   absolute_sigma=False, maxfev=10000, method='trf')

                    fitted_params[noise_level] = params
                    # print(f"  Fit Noise {noise_level:.2f}: T={params[0]:.3f}, S={params[1]:.1f}, L={params[2]:.3f}") # Optional

                    min_x = max(1e-6, subset['pn_ratio'].min() * 0.8)
                    max_x = subset['pn_ratio'].max() * 1.2
                    x_smooth = np.linspace(min_x, max_x, 200)
                    y_fit = weibull_2afc_fit_lapse(x_smooth, *params)
                    ax_ind.plot(x_smooth, y_fit, color=color, linestyle='-', linewidth=2,
                                label=f'Noise {noise_level:.2f} (Fit)')
                except Exception as e:
                    print(f"  Could not fit curve for Noise Level {noise_level}: {e}")
            else:
                print(f"  Skipping fit for Noise Level {noise_level}: Not enough data points.")

        # --- Final Touches on Individual Plot ---
        ax_ind.set_xlabel("信噪比 (P/N Ratio)")
        ax_ind.set_ylabel("正确率")
        ax_ind.set_title(f"心理物理函数: {participant_id}")
        # ax_ind.set_xlabel("Pure Tone / Noise Power Ratio (P/N Ratio)")
        # ax_ind.set_ylabel("Proportion Correct")
        # ax_ind.set_title(f"Psychometric Function: Participant {participant_id}")
        ax_ind.set_ylim(0.45, 1.05)
        ax_ind.set_xlim(right=0.06)
        # Add chance line AFTER potentially adding data points
        if not psych_data.empty: # Only add if there was data
             ax_ind.axhline(0.5, color='grey', linestyle=':', alpha=0.7, linewidth=1.5, label='Chance Level')

        # Create legend for this plot
        handles, labels = ax_ind.get_legend_handles_labels()
        # Keep unique labels in order, filtering out potential duplicates if fit failed etc.
        from collections import OrderedDict
        unique_labels = OrderedDict(zip(labels, handles))
        ax_ind.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=18)
        ax_ind.spines['top'].set_visible(False)
        ax_ind.spines['right'].set_visible(False)
        ax_ind.spines['left'].set_linewidth(1.5)
        ax_ind.spines['bottom'].set_linewidth(1.5)
        plt.tight_layout()
        # Save the individual plot
        individual_plot_filename = os.path.join(figure_folder, f"{participant_id}_psychometric_functions.png")
        try:
             fig_ind.savefig(individual_plot_filename, dpi=300)
             print(f"  Saved: {individual_plot_filename}")
        except Exception as e:
             print(f"  Error saving figure {individual_plot_filename}: {e}")

        plt.close(fig_ind) # Close the figure to free memory


# --- Part 2: Generate Combined Threshold Summary Plot ---

print("\n--- Generating Combined Threshold Summary Plot ---")
txt_file_pattern = os.path.join(data_folder, '*_mixed_noise_results.txt')
participant_txt_files = sorted(glob.glob(txt_file_pattern))
num_participants = len(participant_txt_files)

if not participant_txt_files:
    print(f"Error: No participant TXT results files found matching '{txt_file_pattern}'")
else:
    print(f"Found {num_participants} participant results files.")

    # --- Data Extraction from TXT ---
    all_threshold_data = []
    all_noise_levels_found = set()

    for file_idx, txt_filename in enumerate(participant_txt_files):
        participant_id = os.path.basename(txt_filename).split('_mixed_noise')[0]
        current_noise_level = None
        try:
            with open(txt_filename, 'r') as f:
                for line in f:
                    match_noise = re.search(r"Noise Level:\s*([\d\.]+)", line, re.IGNORECASE)
                    if match_noise:
                        current_noise_level = float(match_noise.group(1))
                        all_noise_levels_found.add(current_noise_level)

                    if current_noise_level is not None:
                        match_thresh = re.search(r"Estimated Threshold.*?Ratio\):\s*([\d\.]+)", line, re.IGNORECASE)
                        if match_thresh:
                            threshold = float(match_thresh.group(1))
                            all_threshold_data.append({
                                'participant_id': participant_id,
                                'participant_idx': file_idx,
                                'noise_level': current_noise_level,
                                'threshold': threshold
                            })
                            current_noise_level = None
        except Exception as e:
            print(f"  Warning: Error parsing {txt_filename}: {e}. Skipping participant.")

    if not all_threshold_data:
        print("Error: No threshold data could be extracted. Cannot create summary plot.")
    else:
        threshold_df = pd.DataFrame(all_threshold_data)
        print(f"Extracted {len(threshold_df)} threshold data points for summary plot.")

        # --- Plotting Setup for Summary ---
        fig_sum, ax_sum = plt.subplots(figsize=(9, 12)) # Use plt.subplots

        # Define Colors using a matplotlib colormap
        sorted_noise_levels = sorted(list(all_noise_levels_found))
        cmap_sum = cm.get_cmap('viridis', len(sorted_noise_levels))
        summary_noise_colors = {level: cmap_sum(i) for i, level in enumerate(sorted_noise_levels)}

        # Define Alpha levels (0.5 to 0.8)
        min_alpha = 0.3
        max_alpha = 0.6
        if num_participants > 1:
            alpha_values = np.linspace(min_alpha, max_alpha, num_participants)
        else:
            alpha_values = np.array([max_alpha]) # Handle single participant

        # --- Plotting Individual Participant Thresholds (Revised) ---
        grouped_by_participant = threshold_df.groupby('participant_id')

        for participant_id, group in grouped_by_participant:
            group_sorted = group.sort_values('noise_level')
            participant_idx = group_sorted['participant_idx'].iloc[0]
            participant_alpha = alpha_values[participant_idx]

            # *** Use ONE plot call per participant for line + markers ***
            if len(group_sorted) >= 1: # Check if there's data for the participant
                ax_sum.plot(group_sorted['noise_level'], group_sorted['threshold'],
                           marker='o',          # Marker style
                           linestyle='-',       # Line style
                           linewidth=3.0,       # Line width
                           markersize=20,       # Marker size
                           color="#5F97D2",     # Fixed color as requested
                           alpha=participant_alpha, # Apply alpha to both line and marker
                           markeredgecolor='none',# No edge color on marker
                           label='_nolegend_') # Hide from automatic legend

        # --- Plotting Average Data ---
        average_thresholds = threshold_df.groupby('noise_level')['threshold'].mean().reset_index()

        if not average_thresholds.empty:
            for _, row in average_thresholds.iterrows():
                noise_color = summary_noise_colors.get(row['noise_level'], 'black')
                ax_sum.plot(row['noise_level'], row['threshold'],
                           marker='^',       # Triangle marker
                           linestyle='-.',
                           color="#EF7A6D",
                           alpha=1.0,        # Full opacity
                           markersize=18,    # Larger size
                           markeredgecolor='none',
                           linewidth=3,
                           zorder=5,         # Plot on top
                           label='_nolegend_')
        else:
            print("Warning: Could not calculate average thresholds.")

        # --- Final Touches on Summary Plot ---
        ax_sum.set_xlabel("噪音功率",fontsize = 24,labelpad=8)
        ax_sum.set_ylabel("信噪比 (P/N Ratio)",fontsize = 24, labelpad=8)
        # ax_sum.set_title(f"Individual and Average Thresholds (N={num_participants} Participants)")

        # Set x-axis ticks and labels
        middle_x = np.mean(sorted_noise_levels)
        ax_sum.set_xlim([middle_x-0.7, middle_x+0.7])
        ax_sum.set_xticks(sorted_noise_levels)
        ax_sum.set_yticks(np.linspace(0.01, 0.05, num=5,endpoint=True))
        ax_sum.set_xticklabels([f"{nl:.2f}" for nl in sorted_noise_levels])
        ax_sum.tick_params(axis="both",which = "major",width = axis_linewidth, length = 6,pad = 8)
        ax_sum.spines['top'].set_visible(False)
        ax_sum.spines['right'].set_visible(False)
        ax_sum.spines['left'].set_linewidth(axis_linewidth)
        ax_sum.spines['bottom'].set_linewidth(axis_linewidth)

        #set significance
        y_line = 0.045
        text_height = 0.001
        ax_sum.hlines(y_line, 0.64, 1.44, color="k", linewidth = axis_linewidth)
        ax_sum.text(middle_x,y_line+text_height,s="$ns$",fontsize = 24, ha = 'center')

        # Adjust y-limits
        if not threshold_df.empty:
             min_thresh = threshold_df['threshold'].min()
             max_thresh = threshold_df['threshold'].max()
             ax_sum.set_ylim(0.003, 0.055) # Increased upper padding slightly

        # --- Create Custom Legend for Summary Plot---
        summary_legend_handles = []
        # Handles for Noise Level Colors (using 'o' marker)
        for noise_level in sorted_noise_levels:
            color = summary_noise_colors.get(noise_level, 'gray')
            summary_legend_handles.append(mlines.Line2D([], [], color=color, marker='o', linestyle='None',
                                                markersize=8, label=f'信噪比 {noise_level:.2f}'))
        # Handle for Average Marker ('^')
        summary_legend_handles.append(mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                                            markersize=9, label='平均阈值'))
        # Handle for Individual Lines
        summary_legend_handles.append(mlines.Line2D([], [], color='grey', marker=None, linestyle='-',
                                            linewidth=1.0, alpha=0.7, label='各被试阈值'))

        # ax_sum.legend(handles=summary_legend_handles, loc=(0.6,0.1),fontsize = 18)

        # --- Save Summary Plot ---
        plt.tight_layout()
        combined_plot_filename = os.path.join(figure_folder, f"ALL_participants_threshold_summary.png")
        try:
            fig_sum.savefig(combined_plot_filename, dpi=300)
            print(f"\nCombined threshold summary plot saved to {combined_plot_filename}")
            threshold_df.to_csv(os.path.join(figure_folder, "ALL_participants_threshold_summary.csv"), index=False)
            average_thresholds.to_csv(os.path.join(figure_folder, "ALL_participants_average_thresholds.csv"), index=False)
        except Exception as e:
            print(f"Error saving summary figure {combined_plot_filename}: {e}")

        plt.close(fig_sum) # Close the summary figure

print("\nProcessing finished.")