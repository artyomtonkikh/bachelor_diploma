stresses=[150 170 190 230];
strains=[0.004 0.006 0.01];
ncycles=[9 80 200];
for n=1:3
    [strain,stress,cycles,stress_amplitude,creep_rate] = simulation_CF_explicit(230, strains(n), ncycles(n));
    if(strains(n)==0.01)
        semilogx(cycles,stress_amplitude,'-^');
    end
    if(strains(n)==0.006)
        semilogx(cycles,stress_amplitude,'-s');
    end
    if(strains(n)==0.004)
        semilogx(cycles,stress_amplitude,'-o');
    end
    %semilogx(cycles,stress_amplitude);
    xlabel('N cycles');
    ylabel('stress amplitude, MPa');
    hold on;
end
figure;
for n=1:4
    [strain,stress,cycles,stress_amplitude,creep_rate] = simulation_CF_explicit(stresses(n), 0.004, 110);
    %loglog(cycles,creep_rate);
    if(stresses(n)==150)
        loglog(cycles,creep_rate,'-d');
    end
    if(stresses(n)==170)
        loglog(cycles,creep_rate,'-^');
    end
    if(stresses(n)==190)
        loglog(cycles,creep_rate,'-s');
    end
    if(stresses(n)==230)
        loglog(cycles,creep_rate,'-o');
    end
    xlabel('N cycles');
    ylabel('creep rate, 1/h');
    hold on;
end
% [strain,stress,cycles,stress_amplitude,creep_rate] = simulation_creep_fatigue(stress_creep, strain_max, number_cycles);
% figure;
% plot(strain,stress);
% xlabel('eps');
% ylabel('stress, MPa');
% grid on;
% figure;
% semilogx(cycles,stress_amplitude);
% xlabel('N cycles');
% ylabel('stress amplitude, MPa');
% grid on;
% figure;
% loglog(cycles,creep_rate);
% xlabel('N cycles');
% ylabel('creep rate, 1/h');
% grid on;