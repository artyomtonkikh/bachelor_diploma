function [] = change_parameters(A, n_Norton, c, kappa, B, m, sigma0, tN, N)


f1=fopen('parameters_set.txt','wt');
pstring=num2str(A,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(n_Norton,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(c,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(kappa,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(B,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(m,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(sigma0,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(tN,9);
fprintf(f1,'%s ',pstring);
pstring=num2str(N,9);
fprintf(f1,'%s ',pstring);
fclose(f1);
