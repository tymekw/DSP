function [zz,pp,ggain] = bilinearTZ(z,p,gain,fs)
% Bilinear transformation: H(s) (analog filter) --> H(z) digital filter
% zeros, poles, gain (z,p,gain) --> zeros, poles, gain (zz,pp,ggain)

pp = []; zz = []; ggain = gain;
for  k=1:length(z)   % transforming zeros
     zz = [ zz (2*fs+z(k))/(2*fs-z(k)) ];
     ggain = ggain*(2*fs-z(k));
end
for  k=1:length(p)   % transforming poles
     pp = [ pp (2*fs+p(k))/(2*fs-p(k)) ];
     ggain = ggain/(2*fs-p(k));
end
if (length(p)>length(z)) zz = [ zz -1*ones(1,length(p)-length(z)) ]; end
if (length(p)<length(z)) pp = [ pp -1*ones(1,length(z)-length(p)) ]; end
