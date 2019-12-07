true_file='../output/enkftemp/true_field.nc'
ncload(true_file,'x','y','nv','zeta')
tmp=[x y zeta([33 45 57],:)'];
save('true_initial_day1_day2_zeta.dat','tmp','-ascii')


nsize=20
for i=1:nsize
    if (i<10)
    ncload(['../output/enkftemp/' 'enkf_output_0' num2str(i) '.nc'],'zeta')
    else
    ncload(['../output/enkftemp/' 'enkf_output_' num2str(i) '.nc'],'zeta')	    
    end
    zeta=zeta(1:24,:);
    if i==1 ; zeta_ave=zeta*0; end
    zeta_ave=zeta/nsize+zeta_ave;
end
tmp=[x y zeta_ave([11 23],:)'];
save('analysis_day1_day2_zeta.dat','tmp','-ascii')

for i=1:nsize
    if (i<10)
    ncload(['../Input/' 'enkf_initial_0' num2str(i) '.nc'],'zeta')
    else
    ncload(['../Input/' 'enkf_initial_' num2str(i) '.nc'],'zeta')
    end

    if i==1 ; zeta_ave=zeta*0; end
    zeta_ave=zeta/nsize+zeta_ave;
end
tmp=[x y zeta_ave'];
save('initial_zeta.dat','tmp','-ascii')
