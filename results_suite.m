% Jan 16, 2023 
%'You can't get away from yourself by moving from one place to another'
%
%% set up directories and output

raw_res_dir = './NS_tf_ISAFSAF/';

raw_data = 'NS_ISAF_2700000_centreline.csv';
numbersteps = 2410000;

output_dir = './NS_tf_ISAFSAF';
overwrite = true;

if(~exist(output_dir,'dir') && overwrite == true)
   mkdir(output_dir) 
end

raw_data = readtable(strcat(raw_res_dir,raw_data));

configuration = 'centreline'; % centreline, radial @ HAB, integrated values


%% isolate key variables in context of each flame configuration

pbe_mode = '2PBE'; % 1PBE, 2PBE

if(strcmp(configuration,'centreline')) 
    z = raw_data.Points_0; % HAM (m)  
    [v_m,psds,v] = gen_PSD_at_HABs(raw_data,z,pbe_mode);
    T = raw_data.T;
    soot_vol_fraction = raw_data.X1; % first moment
    total_number_density = raw_data.X0; % zeroth moment
    total_primary_particle_ND = raw_data.X01; % zeroth moment in terms of n_p
    Nuc_fv = raw_data.Nuc_fv;
    Growth_fv = raw_data.Growth_fv;
    Condense_fv = raw_data.Conden_fv;
    Ox_O2_fv = raw_data.Oxidation_O2_fv;
    Ox_OH_fv = raw_data.Oxidation_OH_fv;
    
end
%% plot

if(strcmp(configuration,'centreline'))
   exp_res_dir = './santoro_centreline_exp_results';
   
   h = zeros(4,1);
   
   exp_soot_volfrac = dlmread(strcat(exp_res_dir,'/central_line_fv_exp'));% note z in cm 
   exp_soot_volfrac(:,1) = exp_soot_volfrac(:,1)/100;% cm --> m
  
   h(1) = figure('Name','Soot Volume Fraction along centreline');
   
   semilogy(z,soot_vol_fraction*10^6,'k-');
   hold on
   semilogy(exp_soot_volfrac(:,1),exp_soot_volfrac(:,2),'or');
        xlabel('HAB (m)');
        ylabel('f_v,  Soot Volume Fraction (ppm)');
        axis([0.01,0.08,10^-3,10^2])
    hold off
    
   T_exp = dlmread(strcat(exp_res_dir,'/central_line_T_exp'));
   T_exp(:,1) = T_exp(:,1)/100; %cm --> m
   
   h(2) = figure('Name','Temperature along centreline');
   hold on
   plot(z,T,'k-');
   plot(T_exp(:,1),T_exp(:,2),'ro');
        xlabel('HAB (m)');
        ylabel('T (K)');
        axis([0,0.1,200,1800]);
        
   hold off
        
   h(3) = figure('Name','Rates contribution to soot volume fraction along centreline');
   hold on
   plot(z,Nuc_fv,'k-');
   plot(z,Growth_fv,'k:');
   plot(z,Condense_fv,'k--');
   plot(z,Ox_O2_fv,'r-');
   plot(z,Ox_OH_fv,'r:');
        legend('Nuc','Growth','Cond','Ox_{O_2}','Ox_{OH}');
        xlabel('HAB (m)');
        ylabel('R need to check units');
        axis([0,0.1,-10^-4,10^-4]);
   hold off
   
   h(4) = figure('Name', 'PSD at HAB 15 mm');
   i_hab = HAB2row_no(0.015,z);
   loglog(v_m,psds(i_hab,:,2),'k-');
   
 
    
end
%% error calc, helps with convergence studies
disp('T %error: ')
T_error_calc(T_exp,z,T)

disp('fv %error: ')
fv_error_calc(exp_soot_volfrac, z, soot_vol_fraction*10^6)


%% save figures
s = true;

if(s) 
   if(strcmp(configuration,'centreline'))
       
       savefig(h,strcat(raw_res_dir,strcat('centreline_sil',num2str(numbersteps))));
   end
    
end
%%

if not(isfolder(raw_res_dir))
    mkdir(raw_res_dir)
end

disp(strcat(raw_res_dir,'centreline'));

save(strcat(raw_res_dir,'centreline'));
%%

%VF

exp_soot_volfrac = dlmread(strcat(exp_res_dir,'/central_line_fv_exp'));% note z in cm 
exp_soot_volfrac(:,1) = exp_soot_volfrac(:,1)/100;% cm --> m
% semilogy(exp_soot_volfrac(:,1),exp_soot_volfrac(:,2),'or');
hold on

tol = 0.1;

sp = spap2(20,5,z, soot_vol_fraction*10^6);
fnplt(sp)

axis([0.01,0.080,10^-3,10^2])























