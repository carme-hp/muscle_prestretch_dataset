#!/bin/bash

# Ask for user input
read -p "Enter dataset id: " id
read -p "Enter prestretch force: " prestretch_force
read -p "Enter specific states call frequency: " specific_states_call_frequency

# Run OpenDiHu simulation
cd simulate_muscle/build_release
./muscle ../settings_contraction_with_prestretch.py muscle.py "$id" "$prestretch_force" "$specific_states_call_frequency"

out_dir="out/muscle_${id}_prestretch_${prestretch_force}_frequency_${specific_states_call_frequency}"
mv "$out_dir" ../../results
cd ../../

