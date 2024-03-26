function [] = make_dat_file(seq, name)
% Creates a .dat file where the first column is the amplitude in Hz and the
% second the phase in degree
    name = strcat(string(name),'.dat');
    fileID = fopen(name,'w');
    fprintf(fileID,'%d \t %d \n',[ seq.w_rf./(2*pi); seq.phase]);
end 