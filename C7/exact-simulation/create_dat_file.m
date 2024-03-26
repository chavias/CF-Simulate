dirname = pwd;      % on scratch
dirname = strsplit(dirname,'/');
dirname = dirname{end};

load(strcat('/cluster/scratch/chavezm/',dirname,'/OptimData.mat'))
make_dat_file(seq,strcat('/cluster/scratch/chavezm/',dirname,'/sequence'));
make_dat_file(rep_seq(seq,10),strcat('/cluster/scratch/chavezm/',dirname,'/sequence10'));

