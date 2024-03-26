dirname = pwd;				 % on scratch
dirname = strsplit(dirname,'/');	 % 
dirname = dirname{end};			 %  


opt = load(strcat('/cluster/scratch/chavezm/',dirname,'/OptimData.mat'));
seq = seq_n_phases_n_amplitudes_noise(opt.xmin,1,0);

% detune in 10 Hz
for i=1:10
    detune = i*10;
    seq_detune(seq,detune,strcat('/cluster/scratch/chavezm/',dirname,'/sequence_detune_',mat2str(detune)));
end
% detune in 100 Hz
for i=2:10
    detune = i*100;
    seq_detune(seq,detune,strcat('/cluster/scratch/chavezm/',dirname,'/sequence_detune_',mat2str(detune)));
end
% detune in kHz
for i=2:5
    detune = i*1000;
    seq_detune(seq,detune,strcat('/cluster/scratch/chavezm/',dirname,'/sequence_detune_',mat2str(detune)));
end



