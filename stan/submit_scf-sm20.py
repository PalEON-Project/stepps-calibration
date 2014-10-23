
import subprocess

runs = [ ('cal_12taxa_lower_comp_v0.1_bigC_c1',
          './calibration_bigC.exe \
<<<<<<< HEAD
          sample num_warmup=250 num_samples=10000 save_warmup=1\
=======
          sample num_warmup=75 num_samples=200 save_warmup=1\
>>>>>>> aa5d5ffa0f6872cd07a84afde5e0603c05ebbb46
          data file=../r/dump/cal_data_12taxa_lower_comp_v0.1.dump \
          output file=../output/12taxa_lower_comp_v0.1_bigC_c1.csv\
          random seed=42'),
         ('cal_12taxa_upper_comp_v0.1_bigC_c2',
          './calibration_bigC.exe \
<<<<<<< HEAD
          sample num_warmup=250 num_samples=10000 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_upper_comp_v0.1.dump \
          output file=../output/12taxa_upper_comp_v0.1_bigC_c1.csv\
          random seed=56'),
         ('cal_12taxa_mid_comp_v0.1_vary_psi_bigC_c1',
          './calibration_vary_psi_bigC.exe \
          sample num_warmup=250 num_samples=10000 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/12taxa_mid_comp_v0.1_vary_psi_bigC_c1.csv\
          random seed=121')
=======
          sample num_warmup=75 num_samples=200 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_upper_comp_v0.1.dump \
          output file=../output/12taxa_upper_comp_v0.1_bigC_c1.csv\
          random seed=42')
>>>>>>> aa5d5ffa0f6872cd07a84afde5e0603c05ebbb46
]

# ('cal_12taxa_mid_comp_v0.1_chain2',
#           './calibration_v7.exe \
#           sample num_warmup=250 num_samples=10000 save_warmup=1\
#           data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
#           output file=../output/12taxa_mid_comp_v0.1_chain2.csv\
#           random seed=21'),
# )

qsub = """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com

cd $HOME/Documents/paleon/stepps2/calibration/stan
export OMP_NUM_THREADS={threads}
srun {command}
"""

dry_run = False

for name, command in runs:
#    sub = qsub.format(queue="low.q", walltime="672:00:00", command=run, threads=1)
    sub = qsub.format(command=command, threads=1, name=name)
    with open(name + ".sh", 'w') as f:
        f.write(sub)
    print "submitting:", name
    if not dry_run:
        subprocess.check_call(['sbatch', name + '.sh'])
    else:
        print sub