% Plot
bound=6;inter=0.01;
[xs1, xs2] = meshgrid(-bound:inter:bound);
gtrue=gfunc(transDist([xs1(:),xs2(:)],tdist));
gkrig=predictor([xs1(:),xs2(:)],dmodel);
gtrue = reshape(gtrue, size(xs1));
gkrig = reshape(gkrig, size(xs1));

hold on
h1=scatter(xt(ysign,1),xt(ysign,2),5,'.','g');
h2=scatter(xt(~ysign,1),xt(~ysign,2),5,'.','r');
scatter(xd(1:N1,1),xd(1:N1,2),40,'o','b');
scatter(xd(N1+1:end,1),xd(N1+1:end,2),70,'s','b');
contour(xs1, xs2, gtrue,'levellist',0,'LineStyle','-', 'Color','b');
contour(xs1, xs2, gkrig,'levellist',0,'LineStyle','--', 'Color','k');
for i=1:thisLayer-1
    circle(3+(i-1)*0.2);
end
hold off
axis equal
axis([-bound,bound,-bound,bound])
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';

legend('Initial DoE','Added DoE','True Function','Kriging Prediction','Hyperball Expansion','Location','southwest');

function h = circle(r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);
h = plot(xunit, yunit,'LineStyle','-.','Color',[0.3 0.3 0.3]);
end