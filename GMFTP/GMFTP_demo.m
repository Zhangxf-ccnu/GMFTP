% Choose a data set to test GMFTP.
data_set = 'Collins';
% data_set = 'Gavin';
% data_set = 'Krogan_core';
% data_set = 'Krogan_extended';
% data_set = 'DIP';
% data_set = 'BIOGRID';

repeat_times = 100;

tau = 0.2;
K =1000;
lambda = 4;
T = 400;
rho = 1e-6;
switch data_set
    case 'Collins'
        
        % Test GMFTP using the Collins network ('.\data\Collins_PPI.txt') and the
        % total GO annotations ('.\data\Collins_Functional_profile.txt').  The detect complexes
        % will be written into file  'Collins_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\Collins_PPI.txt', '.\data\Collins_Functional_profile.txt', 'Collins_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
        
    case 'Gavin'
        
        % Test GMFTP using the Gavin network ('.\data\Gavin_PPI.txt') and the
        % total GO annotations ('.\data\Gavin_Functional_profile.txt').  The detect complexes
        % will be written into file  'Gavin_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\Gavin_PPI.txt', '.\data\Gavin_Functional_profile.txt', 'Gavin_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
        
    case 'Krogan_core'
        
        % Test GMFTP using the Krogan_core network ('.\data\Krogan_core_PPI.txt') and the
        % total GO annotations ('.\data\Krogan_core_Functional_profile.txt').  The detect complexes
        % will be written into file  'Krogan_core_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\Krogan_core_PPI.txt', '.\data\Krogan_core_Functional_profile.txt', 'Krogan_core_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
        
    case 'Krogan_extended'
        
        % Test GMFTP using the Krogan_extended network ('.\data\Krogan_extended_PPI.txt') and the
        % total GO annotations ('.\data\Krogan_extended_Functional_profile.txt').  The detect complexes
        % will be written into file  'Krogan_extended_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\Krogan_extended_PPI.txt', '.\data\Krogan_extended_Functional_profile.txt', 'Krogan_extended_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
            
    case 'DIP'
        
        % Test GMFTP using the DIP network ('.\data\DIP_PPI.txt') and the
        % total GO annotations ('.\data\DIP_Functional_profile.txt').  The detect complexes
        % will be written into file 'DIP_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\DIP_PPI.txt', '.\data\DIP_Functional_profile.txt', 'DIP_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
        
    case 'BIOGRID'
        
        % Test GMFTP using the BIOGRID network ('.\data\BIOGRID_PPI.txt') and the
        % total GO annotations ('.\data\BIOGRID_Functional_profile.txt').  The detect complexes
        % will be written into file  'BIOGRID_complex_result.txt' in current folder.
        theta_star = GMFTP_main('.\data\BIOGRID_PPI.txt', '.\data\BIOGRID_Functional_profile.txt', 'BIOGRID_complex_result.txt', lambda, K, repeat_times, T, rho, tau);
        
end