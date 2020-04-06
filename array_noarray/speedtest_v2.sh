##!/bin/bash

create_param_array(){
    echo "------- MESS params file (v.0.1.1)----------------------------------------------
array_${community_assembly_model}_speed               ## [0] [simulation_name]: The name of this simulation scenario
./array_${community_assembly_model}      ## [1] [project_dir]: Where to save files
200                  ## [2] [generations]: Duration of simulations. Values/ranges Int for generations, or float [0-1] for lambda.
${community_assembly_model}          ## [3] [community_assembly_model]: Model of Community Assembly: neutral, filtering, competition
point_mutation       ## [4] [speciation_model]: Type of speciation process: none, point_mutation, protracted, random_fission
2.2e-08              ## [5] [mutation_rate]: Mutation rate scaled per base per generation
7000             ## [6] [alpha]: Abundance/Ne scaling factor
570                  ## [7] [sequence_length]: Length in bases of the sequence to simulate
------- Metacommunity params: --------------------------------------------------
100                  ## [0] [S_m]: Number of species in the regional pool
750000               ## [1] [J_m]: Total # of individuals in the regional pool
2                    ## [2] [speciation_rate]: Speciation rate of metacommunity
0.7                  ## [3] [death_proportion]: Proportion of speciation rate to be extinction rate
2                    ## [4] [trait_rate_meta]: Trait evolution rate parameter for metacommunity
0.3                    ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
------- LocalCommunity params: island1------------------------------------------
island1              ## [0] [name]: Local community name
1000                 ## [1] [J]: Number of individuals in the local community
0.005                ## [2] [m]: Migration rate into local community
0.001                ## [3] [speciation_prob]: Probability of speciation per timestep in local community" > $filename
}

create_param_noarray(){
    echo "------- MESS params file (v.0.1.1)----------------------------------------------
noarray_${community_assembly_model}_speed               ## [0] [simulation_name]: The name of this simulation scenario
./noarray_${community_assembly_model}      ## [1] [project_dir]: Where to save files
200                  ## [2] [generations]: Duration of simulations. Values/ranges Int for generations, or float [0-1] for lambda.
${community_assembly_model}          ## [3] [community_assembly_model]: Model of Community Assembly: neutral, filtering, competition
point_mutation       ## [4] [speciation_model]: Type of speciation process: none, point_mutation, protracted, random_fission
2.2e-08              ## [5] [mutation_rate]: Mutation rate scaled per base per generation
7000             ## [6] [alpha]: Abundance/Ne scaling factor
570                  ## [7] [sequence_length]: Length in bases of the sequence to simulate
------- Metacommunity params: --------------------------------------------------
100                  ## [0] [S_m]: Number of species in the regional pool
750000               ## [1] [J_m]: Total # of individuals in the regional pool
2                    ## [2] [speciation_rate]: Speciation rate of metacommunity
0.7                  ## [3] [death_proportion]: Proportion of speciation rate to be extinction rate
2                    ## [4] [trait_rate_meta]: Trait evolution rate parameter for metacommunity
0.3                    ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
------- LocalCommunity params: island1------------------------------------------
island1              ## [0] [name]: Local community name
1000                 ## [1] [J]: Number of individuals in the local community
0.005                ## [2] [m]: Migration rate into local community
0.001                ## [3] [speciation_prob]: Probability of speciation per timestep in local community" > $filename
}

run_simulations(){
    MESS -p $filename -s 50 -c 6
}



#for community_assembly_model in "neutral" "filtering" "competition" "pairwise_competition"
for community_assembly_model in "pairwise_competition"
do
    filename="params-speed-${community_assembly_model}-array.txt"
    echo $filename
    create_param_array
    run_simulations
    cp $filename ./array_${community_assembly_model}
done


cp MESS/LocalCommunity.py ../LocalCommunity_save_auto.py
cp MESS/Metacommunity.py ../Metacommunity_save_auto.py
cp MESS/Region.py ../Region_save_auto.py

cp ../LocalCommunity_noarray.py MESS/LocalCommunity.py
cp ../Region_noarray.py MESS/Region.py
cp ../MetaCommunity_noarray.py MESS/Metacommunity.py


#for community_assembly_model in "neutral" "filtering" "competition" "pairwise_competition"
for community_assembly_model in "pairwise_competition"
do
    filename="params-speed-${community_assembly_model}-noarray.txt"
    echo $filename
    create_param_noarray
    run_simulations
    cp $filename ./noarray_${community_assembly_model}
done

cp array_speed_test_* noarray_speed_test_* speed/

create_param_noarray_old(){
    echo "------- MESS params file (v.0.1.1)----------------------------------------------
${community_assembly_model}_${m}_${s}_${J}_${alpha}                 ## [0] [simulation_name]: The name of this simulation scenario
./noarray_${community_assembly_model}_${m}_${s}_${J}_${alpha}       ## [1] [project_dir]: Where to save files
200                    ## [2] [generations]: Duration of simulations. Values/ranges Int for generations, or float [0-1] for lambda.
${community_assembly_model}          ## [3] [community_assembly_model]: Model of Community Assembly: neutral, filtering, competition
point_mutation       ## [4] [speciation_model]: Type of speciation process: none, point_mutation, protracted, random_fission
2.2e-08              ## [5] [mutation_rate]: Mutation rate scaled per base per generation
${alpha}                 ## [6] [alpha]: Abundance/Ne scaling factor
570                  ## [7] [sequence_length]: Length in bases of the sequence to simulate
------- Metacommunity params: --------------------------------------------------
100                  ## [0] [S_m]: Number of species in the regional pool
750000               ## [1] [J_m]: Total # of individuals in the regional pool
2                    ## [2] [speciation_rate]: Speciation rate of metacommunity
0.7                  ## [3] [death_proportion]: Proportion of speciation rate to be extinction rate
2                    ## [4] [trait_rate_meta]: Trait evolution rate parameter for metacommunity
1                    ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
------- LocalCommunity params: island1------------------------------------------
island1              ## [0] [name]: Local community name
${J}                 ## [1] [J]: Number of individuals in the local community
${m}                 ## [2] [m]: Migration rate into local community
${s}                    ## [3] [speciation_prob]: Probability of speciation per timestep in local community" > $filename
}

run_simulations_old(){
    MESS -p "$filename" -s 20 -c 8
}


for community_assembly_model in "pairwise_competition"
do
    for m in 1e-2 1e-3
    do
        for s in 1e-2 1e-3
        do
            for J in 1000
            do
                for alpha in 2000 7000
                do
                    filename="params-speed-${community_assembly_model}-${m}-${s}-${alpha}-noarray.txt"
                    echo $filename
                    create_param_noarray_old
                    run_simulations_old
                    cp $filename ./noarray_${community_assembly_model}_${m}_${s}_${J}_${alpha}
                done
            done
        done
    done
done



cp ../LocalCommunity_save_auto.py MESS/LocalCommunity.py
cp ../Region_save_auto.py MESS/Region.py
cp ../Metacommunity_save_auto.py MESS/Metacommunity.py

