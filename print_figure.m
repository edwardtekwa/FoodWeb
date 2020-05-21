function [] = print_figure( outfile, file_format )
disp('A simple test to illustrate generating figures in a batch mode.');
x = 0:.1:1;
A = exp(x);
plot(x,A);
print(outfile,file_format);
end
