syms b1 b2 mu1 mu2 mu3 d1 d2 d3 d4 d5 d6 a1 a2 a3 g1 g2 cs ca vb hs ha hb es ea eb

hsdt=b1*cs-mu1*hs-d1*hs*ha-d2*hs*hb;
hadt=b1*ca-mu2*ha-d3*ha*hs-d4*ha*hb; 
hbdt=b2*vb-mu3*hb-d5*hb*hs-d6*hb*ha; 
csdt=a1*es-b1*cs-g1*cs*ca; 
cadt=a2*ea-b1*ca-g2*ca*cs; 
vbdt=a3*eb-b2*vb; 
esdt=mu1*hs-a1*es; 
eadt=mu2*ha-a2*ea; 
ebdt=mu3*hb-a3*eb;

hseq=solve(hsdt==0,hs);
haeq=solve(hadt==0,ha);
hbeq=solve(hbdt==0,hb);
cseq=solve(csdt==0,cs);
caeq=solve(cadt==0,ca);
vbeq=solve(vbdt==0,vb);
eseq=solve(esdt==0,es);
eaeq=solve(eadt==0,ea);
ebeq=solve(ebdt==0,eb);

vector_prueba=[hsdt,hadt,hbdt,csdt,cadt,vbdt,esdt,eadt,ebdt];
variables_solve=[hs,ha,hb,cs,ca,vb,es,ea,eb];
sol=solve(vector_prueba==0,variables_solve);

eq1=[-((d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*d5*g1*g2*(d1*d4*d5 + d2*d3*d6)),((d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*d6*g1*g2*(d1*d4*d5 + d2*d3*d6)),-((d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*g1*g2*(d1*d4*d5 + d2*d3*d6)),- (b1^2*d1*d4*d5 + b1^2*d2*d3*d6 + d2*d6*g2*mu1*mu2 + d4*d5*g1*mu1*mu2)/(b1*(d1*d4*d5*g2 + d2*d3*d6*g2)) - (mu2*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*b1*d6*g2^2*(d1*d4*d5 + d2*d3*d6)),(mu1*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*b1*d5*g1^2*(d1*d4*d5 + d2*d3*d6)) - (b1^2*d1*d4*d5 + b1^2*d2*d3*d6 + d2*d6*g2*mu1*mu2 + d4*d5*g1*mu1*mu2)/(b1*(d1*d4*d5*g1 + d2*d3*d6*g1)),-(mu3*(d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*b2*g1*g2*(d1*d4*d5 + d2*d3*d6)),-(mu1*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*a1*d5*g1*g2*(d1*d4*d5 + d2*d3*d6)),(mu2*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*a2*d6*g1*g2*(d1*d4*d5 + d2*d3*d6)),-(mu3*(d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) + d5*g1*mu2 - d6*g2*mu1))/(2*a3*g1*g2*(d1*d4*d5 + d2*d3*d6))];
eq2=[((d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*d5*g1*g2*(d1*d4*d5 + d2*d3*d6)),-((d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*d6*g1*g2*(d1*d4*d5 + d2*d3*d6)),  ((d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*g1*g2*(d1*d4*d5 + d2*d3*d6)),   (mu2*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*b1*d6*g2^2*(d1*d4*d5 + d2*d3*d6)) - (b1^2*d1*d4*d5 + b1^2*d2*d3*d6 + d2*d6*g2*mu1*mu2 + d4*d5*g1*mu1*mu2)/(b1*(d1*d4*d5*g2 + d2*d3*d6*g2)), - (b1^2*d1*d4*d5 + b1^2*d2*d3*d6 + d2*d6*g2*mu1*mu2 + d4*d5*g1*mu1*mu2)/(b1*(d1*d4*d5*g1 + d2*d3*d6*g1)) - (mu1*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*b1*d5*g1^2*(d1*d4*d5 + d2*d3*d6)),  (mu3*(d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*b2*g1*g2*(d1*d4*d5 + d2*d3*d6)),(mu1*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*a1*d5*g1*g2*(d1*d4*d5 + d2*d3*d6)), -(mu2*(d2*d6*g2 + d4*d5*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*a2*d6*g1*g2*(d1*d4*d5 + d2*d3*d6)),  (mu3*(d1*g2 - d3*g1)*(((4*d1*d4*b1^2*d5^2*d6*g1*g2 + 4*d2*d3*b1^2*d5*d6^2*g1*g2 + d4*d5^3*g1^3*mu2^2 + 2*d4*d5^2*d6*g1^2*g2*mu1*mu2 + d2*d5^2*d6*g1^2*g2*mu2^2 + d4*d5*d6^2*g1*g2^2*mu1^2 + 2*d2*d5*d6^2*g1*g2^2*mu1*mu2 + d2*d6^3*g2^3*mu1^2)/(d2*d6*g2 + d4*d5*g1))^(1/2) - d5*g1*mu2 + d6*g2*mu1))/(2*a3*g1*g2*(d1*d4*d5 + d2*d3*d6))];
eq3=[0,0,0,0,0,0,0,0,0];
eq4=[-(d1*b1^2*mu2 + g1*mu1*mu2^2)/(d3*g1*mu1*mu2),0,0,-(d1*b1^2 + g1*mu1*mu2)/(b1*d3*g1),0,0,-(d1*b1^2*mu1*mu2 + g1*mu1^2*mu2^2)/(a1*d3*g1*mu1*mu2),0,0];

J=jacobian(vector_prueba,variables_solve);
eqeva1=subs(J,variables_solve,eq1);
eqeva2=subs(J,variables_solve,eq2);
eqeva3=subs(J,variables_solve,eq3);
eqeva4=subs(J,variables_solve,eq4);

eigen1=eig(eqeva1);
eigen2=eig(eqeva2);
eigen3=eig(eqeva3);
eigen4=eig(eqeva4);

dete1=det(eqeva1);
dete2=det(eqeva2);
dete3=det(eqeva3);
dete4=det(eqeva4);

traza1=trace(eqeva1);
traza2=trace(eqeva2);
traza3=trace(eqeva3);
traza4=trace(eqeva4);