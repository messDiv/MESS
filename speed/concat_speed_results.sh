##!/bin/bash
for community_assembly_model in "neutral" "filtering" "competition" "pairwise_competition"
do
    for m in 0.01 0.001
    do
        for s in 0.01 0.001
        do
            for J in 1000 5000
            do
                statfile="speed_test_${community_assembly_model}_J${J}_m${m}_s${s}.txt"
                cat statfile >> "speed_test_summary_${community_assembly_model}_J_${J}.txt"
                
            done
        done
    done
done

