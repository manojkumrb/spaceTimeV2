

function tick=ML_ytick(a,b,step)

if nargin==2; step=round(100*(b-a)/10)/100; end
y1=0:-step:a;
y2=0:step:b;
tick=sort([y1(2:end) y2]);