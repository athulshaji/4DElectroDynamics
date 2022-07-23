function VOLTAGE = GetVoltageAcrossDiode(I,I0,eta,Vth)
    VOLTAGE = eta*Vth*log((I/I0)+1);
end