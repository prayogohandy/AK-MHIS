clc;clear;close all;
gfunc=@limitstate;
tdist={'norm',0,1 
    'norm',0,1 
    };

boolPlot=true;
AK_MHIS(gfunc,tdist,boolPlot);

function gx = limitstate(x)
    c=1.2;
    gx=c-1/20.*(x(:,1).^2+4).*(x(:,2)-1)+sin(5/2*x(:,1));
end