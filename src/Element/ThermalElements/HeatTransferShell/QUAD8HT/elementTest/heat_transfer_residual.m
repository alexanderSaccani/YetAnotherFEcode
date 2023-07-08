function [Td] = heat_transfer_residual(t,T,Kc,invMc,qext)

q = qext(t);

Td = invMc*(q-Kc*T);

end