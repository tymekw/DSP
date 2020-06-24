function XdB = skaladB(X)
% skalowanie intensywnoœci pikseli obrazu w decybelach
XdB = log10(abs(X)+1); maxXdB = max(max(XdB)); minXdB = min(min(XdB));
XdB = (XdB-minXdB)/(maxXdB-minXdB)*255;