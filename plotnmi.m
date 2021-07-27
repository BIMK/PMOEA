
x=0:0.05:0.5;
% x=0:8;
plot(x,a0,'-R>',x,b0,'-bd',x,c0,'-ks',x,d0,'-pg');h=legend('MOEA/D-Net','MEME-Net','GN','Infomod',3);
xlabel('External degree of a node');
ylabel('NMI');