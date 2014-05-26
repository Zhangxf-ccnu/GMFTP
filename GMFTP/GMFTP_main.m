function theta_star = GMFTP_main(PPI_profie, Functional_profile, output_file_name, lambda, K, repeat_times, T, rho, tau)
% GMFTP_main is the main function of our model. It reads data (PPI network and functional profile) from the text
% files (function: data_read) and detects protein complexes using the model
% described in the paper (function: GMFTP). It also writes the detected
% complex into file 'output_file_name' where each line is a detected
% complex of tab-separated labels.


% Inputs:
%   PPI_profie: the input file name of the PPI network, where each line
%   contains two proteins defining the interaction between them.
%   Example: YAL064W-B  YPR126C

%   Functional_profile: the input file name of the functional profile, where each line
%   contains a protein and a function defining the the association between them.
%   Example: YNL079C    GO:0000001
%   If the the functional information of proteins is not available, please
%   input an empty file which dose not include any characters.

%   output_file_name: the file into which GMFTP writes the detected
%   complexes. Each line is a detected complex of tab-separated labels.
%   Example: YNL071W YER178W YGR193C YBR221C

%   lambda: rate parameter of exponential distribution. The default value
%   is 4.
%   K: maximum number of possible protein complexes. The default value is
%   1000.
%   repeat_times: the number of times that we repeat the entire calculation
%   to avoid local minimum. The default value is 100.
%   T: the number of iterations limited in GMFTP. The default value is 400.
%   rho: the tolerance threshold of the stop criterion.  The default
%   value is 1e-6.
%	tau: the threshold used to obtain protein complexes from estimator of theta. The default
%   value is 0.2.


% Outputs:
%   theta_star: the protein-complex membership indication matrix.



if nargin < 9
    tau = 0.2;
end

if nargin < 8
    rho = 1e-6;
end

if nargin < 7
    T = 400;
end

if nargin < 6
    repeat_times = 100;
end

if nargin < 5
    K = 1000;
end


if nargin < 4
    lambda = 4;
end

if nargin < 3
    output_file_name = 'complex_result.txt';
end


if nargin < 2
    error('You need input the PPI and functional profile');
end

% Read data from the text file and write it into matlab format
Data_set = data_read(PPI_profie, Functional_profile, K);

% Detect protein complex using GMFTP.
theta_star = GMFTP(Data_set.PPI, Data_set.Function_profile, lambda, K, repeat_times, T, rho, tau);

% Write the result of detected complexes into file 'output_file_name'.
fid = fopen(output_file_name,'w');
for k = 1:size(theta_star,2)
    member_indices = find(theta_star(:,k));
    for i = 1:length(member_indices);
        fprintf(fid, '%s\t', cell2mat(Data_set.Protein(member_indices(i))) );
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf(['The detected complex have been written into file ', output_file_name])
fprintf('\n')







function Data_set = data_read(PPI_profie, Functional_profile, K)
fprintf('Reading data...')
fprintf('\n')
fid_ppi=fopen(PPI_profie);
temp_PPI=textscan(fid_ppi,'%s%s%*[^\n]','delimiter','\t');
fclose(fid_ppi);

Data_set.Protein = union(temp_PPI{1},temp_PPI{2});
[~,Locb_1] = ismember(temp_PPI{1}, Data_set.Protein);
[~,Locb_2] = ismember(temp_PPI{2}, Data_set.Protein);
Data_set.PPI = sparse(Locb_1,Locb_2,1,length(Data_set.Protein),length(Data_set.Protein));
Data_set.PPI = Data_set.PPI + Data_set.PPI';


fid_function=fopen(Functional_profile);
temp_function=textscan(fid_function,'%s%s%*[^\n]','delimiter','\t');
fclose(fid_function);

if isempty(temp_function{1})&& isempty(temp_function{2})
    Data_set.Function_profile = sparse(length(Data_set.Protein), K);
else  
    Data_set.Function = unique(temp_function{2});
    [~,Locb_1] = ismember(temp_function{1}, Data_set.Protein);
    [~,Locb_2] = ismember(temp_function{2}, Data_set.Function);
    Data_set.Function_profile = sparse(Locb_1,Locb_2,1,length(Data_set.Protein),length(Data_set.Function));
end





