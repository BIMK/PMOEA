%打印算法迭代过程
%generation_size: 迭代次数

function plotGA(generation_size)
global fitness_value;
global best_Q;
x = 1:1:generation_size;
y = best_Q;
plot(x,y)