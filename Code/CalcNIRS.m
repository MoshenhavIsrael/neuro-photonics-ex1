function [ dHbR , dHbO, fig ] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile )
%% CalcNIRS - calculate and plot HbR HbO
% Input:
% dataFile - .mat file with intensity data.
%    DS.Lambda : two wavelengths (in nm)
%    t : time vector
%    d : intensity data of 20 channels
%        20 first rows-> first wavelength, 20 last rows->second wavelength
% SDS - Sourse-Detector Separation distance in cm
% tissueType - one of the rows in DPFperTissueFile (for example 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' )
% plotChannelIdx - vector with numbers in the range of [1-20] indicating channels to plot. If empty - none is plotted. (default = [])
%
% extinctionCoefficientsFile - .csv file with the following columns : wavelength, Water, HbO2, HHb, FatSoybean
%        default = '.\ExtinctionCoefficientsData.csv' (if not passed or empty)
% DPFperTissueFile - .txt file with two columns: Tissue and DPF (Tissue is tissue type, corresponding with tissueType input variable)
%        measured at 807nm
%        default = '.\DPFperTissue.txt' (if not passed or empty)
% relDPFfile - relative DPF according to wavelength
%        default = '.\RelativeDPFCoefficients.csv' (if not passed or empty)
%
% Output :
%    dHbR - HbR concentration change for all channels (nx20) where n is time vector length
%    dHbO - HbO concentration change for all channels (nx20) where n is time vector length
%    fig - handle to figure. Empty if plotChannelIdx==[].

%% Initialization

% define default values
arguments
    dataFile 
    SDS 
    tissueType 
    plotChannelIdx =[]
    extinctionCoefficientsFile = '.\ExtinctionCoefficientsData.csv'
    DPFperTissueFile = '.\DPFperTissue.txt'
    relDPFfile = '.\RelativeDPFCoefficients.csv'
end

% Check params validity
% For each validation display the mistake if exist

valid = zeros([6,1]); % 6 validations
% valid 1 - variables of type string
% Define "isText" function for valid(1)
function tf = isTextScalar(x)
    tf = ischar(x) || (isstring(x) && isscalar(x));
end
valid(1) = isTextScalar(dataFile) && isTextScalar(tissueType) &&isTextScalar(extinctionCoefficientsFile) &&isTextScalar(DPFperTissueFile) && isTextScalar(relDPFfile);
if valid(1) == 0
    fprintf("At least one of the the following parameters set incorrectly: \n - dataFile \n - tissueType \n - extinctionCoefficientsFile \n - DPFperTissueFile \n - relDPFfile\n")
    fprintf("All those params must be strings\n\n")
end
% valid 2 - variables which are numeric
valid(2) = isnumeric(SDS) && isnumeric(plotChannelIdx);
if valid(2) == 0
    fprintf("At least one of the the following parameters set incorrectly: \n - SDS \n - plotChannelIdx\n")
    fprintf("All those params must be numeric\n\n")
end
% valid 3 - SDS is single value with reasonable value in cm
valid(3) = isscalar(SDS) && SDS>0.1 && SDS<15; % SDS must be scalar in units of [cm]
if valid(3) == 0
    fprintf("SDS must be a scalar in units of [cm] (and with reasonable value)\n\n")
end
% valid 1 - plotChannelIdx is 0-20 whole numbers
valid(4) = length(plotChannelIdx)<=20 && all(floor(plotChannelIdx)==plotChannelIdx);
if valid(4) == 0
    fprintf("The variable plotChannelIdx must contain 0-20 whole numbers\n\n")
end
% valid 5 - files existance
valid(5) = exist(dataFile, "file") && exist(extinctionCoefficientsFile, "file") && exist(DPFperTissueFile, "file") && exist(relDPFfile, "file");
if valid(5) == 0
    fprintf("At least one of the the following files are missing: \n - dataFile \n - extinctionCoefficientsFile \n - DPFperTissueFile \n - relDPFfile\n")
    fprintf("Check path or spelling\n\n")
end
% valid 6 - correctness of files extinsions 
valid(6) = 1;
filesNames = {dataFile, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile};
filesTypes = {'.mat', '.csv', '.txt', '.csv'};
for i=1:4
    [~, ~, ext] = fileparts(filesNames{i});
    valid(6) = valid(6) && string(ext) == string(filesTypes{i});
end
if valid(6) == 0
    fprintf("At least one of the the following files don't match their expected suffix: \n - dataFile (mat) \n - extinctionCoefficientsFile (csv) \n - DPFperTissueFile (txt) \n - relDPFfile (csv) \n")
    fprintf("Check path or spelling\n\n")
end
% Exit function if one of params doesn't valid
if sum(logical(valid))<length(valid)
    return
end

% load data file and verify its content is the expected content
% Exit function if there is a mistake and display the mistake
load(dataFile, 't', 'd', 'SD');
if ~(exist('t') && exist('d') && exist('SD'))
    disp("The data file is not appropriate")
    disp("File must contain the variables: t, d and SD")
    return
end
if ~(isfield(SD, 'Lambda'))
    disp("The struct 'SD' must have the field 'Lambda'")
    return
end

%% Main calculation

% Using this equation: 
%       OD(t) = log_10(I(0)/I(t)) = L_eff[epsilon_HbR/ln(10)*dHbR + epsilon_HbO/ln(10)*dHbO]
% Of course that L_eff = SDS * DPF
% Note that epsilons and DPF depend on wavelength
% For system of 2 wavelengths derive 2 equations we will get the form of:

%       [OD_w1/Leff_w1    =   epsilonsMat * [dHbR
%        OD_w2/Leff_w2]                      dHbO]

% Assumptions:
% - the coeffisients in the extinctionCoefficientsFile are in form of
%       epsilon/ln(10)
% - The DPF in DPFperTissueFile is suit for 807nm and can be corrected by
%       the relDPFs by: DPF(w) = DPF(807)*relDPF(w)/relDPF(807)
% - Each wavelength contributes one equation, so we have 2 wavelengths for
%       2 equations and then we can find the 2 vars (dHbO, dHbR)
% - Each point on the time axis can be solved by itself

% Initial variables 
dHbR = zeros([length(t), 20]);
dHbO = zeros([length(t), 20]);
w1 = SD.Lambda(1);
w2 = SD.Lambda(2);

% get the DPF compatible to the tissue. The dictionary of DPFperTissue is
% relevant for 807nm
DPFperTissue = readtable(DPFperTissueFile);
for idx = 1:length(DPFperTissue.Tissue)
    if strcmp(DPFperTissue.Tissue{idx}, tissueType)
        break;
    end
end
DPF_807 = DPFperTissue.DPF(idx);

% Correct the DPF to the real wavelengths by the relation DPF table. Assume
% that diviving by relDPF(807) gives the DPF of reference wavelength and
% then multiply the reference DPF with relDPF(w) gives DPF(w)
relDPFtable = readtable(relDPFfile);
relDPF_807 = table2array(relDPFtable(relDPFtable.wavelength==807, "relDPFcoeff"));
relDPF_w1 = table2array(relDPFtable(relDPFtable.wavelength==w1, "relDPFcoeff"));
relDPF_w2 = table2array(relDPFtable(relDPFtable.wavelength==w2, "relDPFcoeff"));
DPFw1 = DPF_807 * relDPF_w1 / relDPF_807;
DPFw2 = DPF_807 * relDPF_w2 / relDPF_807;
Leff_w1 = SDS*DPFw1;
Leff_w2 = SDS*DPFw2;

% Extract coeffisients of HbR (called also 'HHb') and HbO (called 'HbO2')
% for the 2 wavelengths and save them as matrix
% Assume coeffisient are in form of epsilon/ln(10)
extCoeffTable = readtable(extinctionCoefficientsFile);
extCoeff_w1 = extCoeffTable(extCoeffTable.wavelength==w1, ["HHb", "HbO2"]);
extCoeff_w2 = extCoeffTable(extCoeffTable.wavelength==w2, ["HHb", "HbO2"]);
extCoeffMat = table2array([extCoeff_w1; extCoeff_w2]);

% Solve the Linear Equation System for all the 20 channels
% For each channel we have two equations with 2 vars, and we solve the system by linear algebra
for ch = 1:20
    % Extract signals of channel from the 2 wavelength
    I_t = d(:, [ch, ch+20]);
    % Transpose so that each row will display signal of one wavelength
    I_t = I_t';     
    I_0 = I_t(:, 1);
    % Calc the Optical Densities over time (OD(t))
    dOD_t = log10(I_0./I_t);
    % Divide by L_eff
    B = dOD_t./[Leff_w1; Leff_w2];
    % Now iterate the time axis and solve each time point separately
    for i = 1:length(t)
        bVector = B(:, i);
        % Solve the linear system and save dHbO and dHbR
        %   The linear system is:
        %           extCoeffMat * x = bVector
        %   So the sulotion is:
        x = extCoeffMat \ bVector;
        dHbR(i, ch) = x(1);
        dHbO(i, ch) = x(2);
    end
end

%% Figure with results
% Initial
fig = figure;
nChannels = length(plotChannelIdx);
nRows = floor(sqrt(nChannels));
nCols = ceil(nChannels/nRows);
if length(plotChannelIdx) > 1
    sgtitle('Changes in Hemoglobin concentrations over time');
end
% Subplots loop
for i = 1:length(plotChannelIdx)
    ch = plotChannelIdx(i);
    subplot(nRows, nCols, i);
    hold;
    plot(t, dHbR(:, ch), 'b'); 
    plot(t, dHbO(:, ch), 'r'); 
    title(['channel ' num2str(ch)]);
    if length(plotChannelIdx) < 2
        title(['Changes in Hemoglobin concentrations over time for channel ' num2str(ch)]);
    end
    legend(["dHbR", "dHbO"]);
    xlabel("time (sec)");
    ylabel("Concentration change (mol/L)");
end

end
