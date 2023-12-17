function plot_graph(x,y,A,n)
figure(n);
hold on
plot(x,y,'o','Markersize',15,'MarkerEdgeColor',[0.5,0,0],'MarkerFaceColor',[1,0,0]);
ind = find(A == 1);
[I,J] = ind2sub(size(A),ind);
for k = 1 : length(ind)
    plot([x(I(k)),x(J(k))],[y(I(k)),y(J(k))],'linewidth',4,'Color',[0,0,0.5]);
end
daspect([1,1,1])
axis off
end