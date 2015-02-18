function mel = fhz2mel(hz)
%fhz2mel Conversion between mel/linear scales using O'Shaugnessy's formula

mel = 2595*log10(1+hz/700);

end

