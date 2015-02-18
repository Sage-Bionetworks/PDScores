function hz = fmel2hz( mel )
%fmel2hz Conversion between mel/linear scales using O'Shaugnessy's formula

hz = 700*(10.^(mel/2595)-1);

end

