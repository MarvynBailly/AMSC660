numLoops = 10;

figure;
hold on;

legend_entries = cell(1, numLoops);

for i = 1:numLoops
    n = i*100;
    rayleigh(n,10e-12);
    legend_entries{i} = ['n = ', num2str(n)];
end


legend(legend_entries);
hold off;