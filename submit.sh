#guppy
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_%a.out -a 1-1000 `pwd`/runGuppy.sh 0 >> submit.out
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_%a.out -a 1-1000 `pwd`/runGuppy.sh 1000 >> submit.out
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_%a.out -a 1-1000 `pwd`/runGuppy.sh 2000 >> submit.out
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_%a.out -a 1-1000 `pwd`/runGuppy.sh 3000 >> submit.out
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_%a.out -a 1-308 `pwd`/runGuppy.sh 4000 >> submit.out
#sbatch -D `pwd` --time=48:00:00 --partition=gpu --gres=gpu:p100:1 --cpus-per-task=4 --mem=60g -o `pwd`/%A_guppy.out `pwd`/runGuppy.sh >> submit.out

#fast5mod
#sbatch -D `pwd` --time=48:00:00  --cpus-per-task=4 --mem=60g -o `pwd`/%A_fast5mod.out `pwd`/runFast5mod.sh >> submit_fast5mod_2.out

#call
#sbatch -D `pwd` --time=48:00:00  --cpus-per-task=4 --mem=60g -o `pwd`/%A_call.out `pwd`/runCallmeth.sh >> submit_call.out

#fast5mod and call
#sbatch -D `pwd` --time=48:00:00  --cpus-per-task=4 --mem=60g -o `pwd`/%A_guppy_fast5mod.out `pwd`/runGuppyFast5mod.sh >> submit_guppy_fast5mod.out

#fast5_sub
sbatch -D `pwd` --time=70:00:00  --cpus-per-task=4 --mem=60g -o `pwd`/%A_sub.out `pwd`/fast5_sub_3.sh >> submit_sub.out

#isolate read names
#sbatch -D `pwd` --time=48:00:00  --cpus-per-task=4 --mem=60g -o `pwd`/%A_IDs.out `pwd`/pull_read_names.sh >> submit_IDs.out

