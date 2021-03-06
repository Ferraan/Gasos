function cp = cp (T0, Tf,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4)
if(Tf<1000)  
    cp=(Tf*a0+a1*Tf^2/2+a2*Tf.^3/3+a3*Tf.^4/4+a4*Tf.^5/5-(T0*a0+a1*T0^2/2+a2*T0.^3/3+a3*T0.^4/4+a4*T0.^5/5))/(Tf-T0);
elseif(Tf>=1000)
    cp=((1000*a0+a1*1000^2/2+a2*1000.^3/3+a3*1000.^4/4+a4*1000.^5/5-(T0*a0+a1*T0^2/2+a2*T0.^3/3+a3*T0.^4/4+a4*T0.^5/5))+(Tf*b0+b1*Tf^2/2+b2*Tf.^3/3+b3*Tf.^4/4+b4*Tf.^5/5-(1000*b0+b1*1000^2/2+b2*1000.^3/3+b3*1000.^4/4+b4*1000.^5/5)))/(Tf-T0);
end
  cp=cp*8.31447;
end
