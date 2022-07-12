function [new_sigprof,hkm_fr,hkm_nr] = gate_prof(sigprof,start_bin_fr,start_bin_nr,hkm)
hkm_fr = hkm(start_bin_fr:end);
new_sigprof.an287fr=sigprof.an287(start_bin_fr:end,:);
new_sigprof.pc287fr=sigprof.pc287(start_bin_fr:end,:);
new_sigprof.an299fr=sigprof.an299(start_bin_fr:end,:);
new_sigprof.pc299fr=sigprof.pc299(start_bin_fr:end,:);

hkm_nr = hkm(start_bin_nr:end);
new_sigprof.an287nr=sigprof.an287nr(start_bin_nr:end,:);
new_sigprof.pc287nr=sigprof.pc287nr(start_bin_nr:end,:);
new_sigprof.an299nr=sigprof.an299nr(start_bin_nr:end,:);
new_sigprof.pc299nr=sigprof.pc299nr(start_bin_nr:end,:);

