#/bin/sh
nohup srun -A fungcat --partition gpu_q -t 3-1 --mem-per-cpu 48GB --gres=gpu:1  --nodelist pegsndgpu00$2 ./Plinko.py predict -u root -p 1 /scratch/fungcat/shreyas/coding_domain/one_vs_$3/model_dir/$1/ /scratch/fungcat/jsporter/Data/Data_Ensemble_CDs/$1/ &> /scratch/fungcat/jsporter/Prediction_Out/ensemble_cds.$1.predict.out &
