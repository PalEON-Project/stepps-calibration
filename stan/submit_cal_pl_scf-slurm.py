
import subprocess

runs = [ ('cal_pl',
          './cal_pl.exe \
          sample num_warmup=250 num_samples=500 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/cal_pl_v0.3.csv\
          random seed=42'),
         ('cal_pl_Kgamma_EPs_v0.3',
          './cal_pl_Kgamma_EPs.exe \
          sample num_warmup=250 num_samples=500 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/cal_pl_Kgamma_EPs_v0.3.csv\
          random seed=42'),
         ('cal_pl_Ka_Kb_Kgamma_EPs_v0.3',
          './cal_pl_Ka_Kb_Kgamma_EPs.exe \
          sample num_warmup=250 num_samples=500 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/cal_pl_Ka_Kb_Kgamma_EPs_v0.3.csv\
          random seed=42'),
         ('cal_pl_Ka_Kgamma_EPs_v0.3',
          './cal_pl_Ka_Kgamma_EPs.exe \
          sample num_warmup=250 num_samples=500 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/cal_pl_Ka_Kgamma_EPs_v0.3.csv\
          random seed=42'),
         ('cal_pl_Ka_EPs_v0.3',
          './cal_pl_Ka_EPs.exe \
          sample num_warmup=250 num_samples=500 save_warmup=1\
          data file=../r/dump/cal_data_12taxa_mid_comp_v0.1.dump \
          output file=../output/cal_pl_Ka_EPs_v0.3.csv\
          random seed=42')
]

qsub = """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com

cd $HOME/Documents/projects/stepps-calibration/stan
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
