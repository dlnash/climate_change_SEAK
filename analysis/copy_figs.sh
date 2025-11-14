#!/bin/bash
######################################################################
# Filename:    copy_figs.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Copies finalized figures to the 'final_figs' folder and 
#              Overleaf manuscript directory, runs git pull/push for Overleaf.
######################################################################

set -euo pipefail

# --- Input directories ---
maindir="../figs"                            # Main figure folder
finaldir="${maindir}/final_figs"             # Folder for final figures
overleafdir="/cw3e/mead/projects/cwp140/repos/SEAK_climate_change_manuscript"  # Overleaf repo

# --- Figure name lists ---
fig_names=(
    "topographic_map"
    "clim/ONDJFM_95th_percentile_clim"
    "clim/ONDJFM_95th_percentile_clim_hexbin"
    "cfsr_ros_landslide_strict"
    "ros_case"
    "ros_strict/ONDJFM/ONDJFM_ROS_FREQ_DIFF"
    "ros_strict/ONDJFM/ONDJFM_ros_frequency_clim_hexbin"
    "ros_strict/ONDJFM/ONDJFM_ROS_INTENSITY_DIFF"
    "ros_strict/ONDJFM/ONDJFM_ros_intensity_clim_hexbin"
)

fig_nums=(1 2 3 4 5 6 7 8 9)

# --- Ensure directories exist ---
mkdir -p "$finaldir"
mkdir -p "$overleafdir"

# --- Pull latest from Overleaf git repo ---
echo "🔄 Pulling latest changes from Overleaf..."
if [[ -d "$overleafdir/.git" ]]; then
    git -C "$overleafdir" pull --rebase
else
    echo "⚠️  Warning: $overleafdir is not a git repository. Skipping git pull."
fi

# --- Copy figures ---
echo "📂 Copying figures to $finaldir and $overleafdir..."
for i in "${!fig_names[@]}"; do
    src="${maindir}/${fig_names[$i]}.png"
    base="fig${fig_nums[$i]}.png"

    if [[ -f "$src" ]]; then
        cp -v "$src" "${finaldir}/${base}"
        cp -v "$src" "${overleafdir}/${base}"
    else
        echo "⚠️  Missing source file: $src"
    fi
done

# --- Zip all final figs ---
echo "🗜️  Zipping final figures..."
cd "$finaldir"
zip -r figs.zip fig*.png > /dev/null
echo "✅ Created $(realpath figs.zip)"

# --- Commit and push to Overleaf ---
echo "⬆️  Committing and pushing updates to Overleaf..."
if [[ -d "$overleafdir/.git" ]]; then
    cd "$overleafdir"
    git add fig*.png || true

    # Only commit if there are staged changes
    if ! git diff --cached --quiet; then
        git commit -m "Update final figures on $(date '+%Y-%m-%d %H:%M')"
        git push
        echo "✅ Figures committed and pushed to Overleaf."
    else
        echo "ℹ️  No figure changes to commit."
    fi
else
    echo "⚠️  Warning: $overleafdir is not a git repository. Skipping git push."
fi
