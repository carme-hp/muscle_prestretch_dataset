#!/bin/bash

# Ask for dataset id
read -p "Enter dataset id: " id

# Define parameter sweeps
prestretch_forces=(0.0 1.0 2.0)
frequencies=(0.1 0.01)

# Go to build directory once
cd simulate_muscle/build_release || exit 1

for prestretch_force in "${prestretch_forces[@]}"; do
  for specific_states_call_frequency in "${frequencies[@]}"; do

    echo "Running simulation:"
    echo "  prestretch_force = $prestretch_force"
    echo "  frequency        = $specific_states_call_frequency"

    ./muscle \
      ../settings_contraction_with_prestretch.py \
      muscle.py \
      "$id" \
      "$prestretch_force" \
      "$specific_states_call_frequency"

    out_dir="out/muscle_${id}_prestretch_${prestretch_force}_frequency_${specific_states_call_frequency}"
    mv "$out_dir" ../../results

  done
done

cd ../../
