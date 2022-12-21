#!/bin/bash

##SBATCH -C haswell
##SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem 10000
#SBATCH --time=10:00:00
#SBATCH --output=batch_jobs/myjob.o%j
##SBATCH --mail-user=srinirag@illinois.edu
##SBATCH --mail-type=ALL

##modules
##module load Python/2.7.9-GCC-4.9.2-bare
##module load xcb-proto/1.11-intel-2016.u3-Python-2.7.9
##module use /home/sri/modulefiles/
##module load anaconda
##module load python

# Your script content goes here...
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for10years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for10years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for1years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for1years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for2years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for2years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for3years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for3years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for4years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for4years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for5years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for5years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for6years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for6years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for7years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for7years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for8years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for8years.npy
python3 get_fisher_withdelensing.py -ilc_fname ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for9years.npy -opfname results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs//s4deepv3r025/s4deepv3r025_ilc_galaxy0_27-39-93-145-225-278_TT-EE_for9years.npy
