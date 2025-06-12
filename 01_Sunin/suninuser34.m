% Sunin parameter file, see "sunin_readme.doc" for details
% Last modification 5/2/06, add trial by trial filtering
% Last modification 3/24/06, add analysis based on ivf files. 
%% file location 
clear;
system='v';             % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'P:\';     % Data disk name
% datafolder = 'HRXu_analysis\';   % Data folder name on data disk, results will be saved in 'expresult'
% expname = '00_data\L15\120104_L15\'; % Exp folder name (in both data folder and result folder)
% runname = 'Run09_BRtraining_f10_2size\';      % Run foler name (in both data folder and result folder)
% runname = 'Run01_G8\';      % Run foler name (in both data folder and result folder)
datafolder = 'HRXu_analysis\00_data\L1\';     % Data disk name
% % datafolder = '00_data\';   % Data folder name on data disk, results will be saved in 'expresult'
expname = '091218_L1\'; % Exp folder name (in both data folder and result folder)
% runname = 'run03_G8\';      % Run foler name (in both data folder and result folder)
% datafolder = 'expt\';     % Data disk name
% expname = '141119_L30\'; % Exp folder name (in both data folder and result folder)
% runname = 'Run01_G8\';      % Run foler name (in both data folder and result folder)
% datafolder = 'expt\';     % Data disk name
% expname = '141208_L45\'; % Exp folder name (in both data folder and result folder)
% runname = 'Run05_LHBar\E03\';      % Run foler name (in both data folder and result folder)
runname = 'run02_RGLum8\';      % Run foler name (in both data folder and result folder)

filename={
% 'Angle_E01B000.BLK'
};
%% grating subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim={};
% 
% G8
% Smap(1, :)={'HV', [3 7], [1 5]};
% Smap(2, :)={'AO', [4 8], [2 6]};
% Smap(3, :)={'H0', [3 7], [9]};
% Smap(4, :)={'V0', [1 5], [9]};
% Smap(5, :)={'A0', [4 8], [9]};
% Smap(6, :)={'O0', [2 6], [9]};% 
% % 
% % % Smap(1, :)={'H_V', [1 3 5], [2 4 6]};
% % % Smap(2, :)={'L_R', [3 4], [5 6]};
% % % Smap(3, :)={'D_S', [1 2], [3:6]};
% % 
% % % Smap(1, :)={'H_All', [3 7], [1:9]};
% % % Smap(2, :)={'V_All', [1 5], [1:9]};
% % % Smap(3, :)={'A_All', [4 8], [1:9]};
% % % Smap(4, :)={'O_All', [2 6], [1:9]};
% % % Smap(5, :)={'H_V', [3 7], [1 5]};
% % % Smap(6, :)={'V_H', [1 5], [3 7]};

% G8 OD_calculation_110831
%  Smap(1, :)={'AO0', [2 4 6 8], [9]};

% G4
% Smap(1, :)={'HV', [2 4], [1 3]};
% Smap(2, :)={'All0', [1:4], [5]};
% Smap(3, :)={'H0', [2 4], [5]};
% Smap(4, :)={'V0', [1 3], [5]};
% Smap(5, :)={'LU_RD', [2 3], [1 4]};
% % Smap(6, :)={'O0', [2 6], [9]};% 

% G16


% Condition name
% stim(1)={'VR'};  
% stim(3)={'OR'};
% stim(5)={'HU'};
% stim(7)={'AL'};
% stim(9)={'VL'};  
% stim(11)={'OL'};
% stim(13)={'HD'};
% stim(15)={'AR'};

% Subtraction maps
% Smap(1, :)={'HV', [5 13], [1 9]};
% Smap(2, :)={'AO', [7 15], [3 11]};
% Smap(1, :)={'HV', [3 7], [1 5]};
% Smap(2, :)={'AO', [4 8], [2 6]};
% Smap(1, :)={'H0', [3 7], [9]};
% Smap(2, :)={'V0', [1 5], [9]};
% Smap(3, :)={'A0', [4 8], [9]};
% Smap(4, :)={'O0', [2 6], [9]};


% Smap(1, :)={'H_All', [5 13], [1:17]};
% Smap(2, :)={'V_All', [1 9], [1:17]};
% Smap(3, :)={'A_All', [7 15], [1:17]};
% Smap(4, :)={'O_All', [3 11], [1:17]};

%% RGLum

% % RGLum
% from Xiao Chao
% % Condition name
% 
% % RGLumSFTF
% stim(1)={'LCA'};   % Name stim condistions here
% stim(2)={'LCO'};
% stim(3)={'RCA'};
% stim(4)={'RCO'};
% stim(5)={'LLA'};   % Name stim condistions here
% stim(6)={'LLO'};
% stim(7)={'RLA'};
% stim(8)={'RLO'};
% stim(9)={'Blank'};
% % stim={};       % leave this line if you don't want to define condition names
% 
% %           name    A  -  B    (List of subtraction maps)
Smap(1, :)={'RGLum', [1 2 3 4], [5 6 7 8]};
% Smap(2, :)={'RGlumA', [1 3], [5 7]};
% Smap(3, :)={'RGlumO', [2 4], [6 8]};
% Smap(4, :)={'AO', [1 3 5 7], [2 4 6 8]};   
% Smap(5, :)={'AOrg', [1 3], [2 4]};   
% Smap(6, :)={'AOlum', [5 7], [6 8]};
% Smap(7, :)={'RGo', [1 2 3 4], [9]};
% Smap(8, :)={'RGd', [1 2 3 4], []};
% Smap(9, :)={'Lumo', [5 6 7 8], [9]};
% Smap(10, :)={'Lumd', [5 6 7 8], []};
% Smap(11, :)={'allo', [1 2 3 4 5 6 7 8], [9]};
% Smap(12, :)={'alld', [1 2 3 4 5 6 7 8], []};
Smap(2, :)={'L-R', [1 2 5 6], [3 4 7 8]};
% Smap(14, :)={'L-Rrg', [1 2], [3 4]};
% Smap(15, :)={'L-Rlu', [5 6], [7 8]};

% From Chen Ming
% Condition name
% stim(1)={'LRgA'};   % Name stim condistions here
% stim(2)={'LRgO'};
% stim(3)={'RRgA'};
% stim(4)={'RRgO'};
% stim(5)={'LLumA'};   
% stim(6)={'LLumO'};
% stim(7)={'RLumA'};
% stim(8)={'RLumO'};
% stim(9)={'blank'};
% % 
% % % Subtraction maps
% Smap(1, :)={'RGlum', [1 2 5 6], [3 4 7 8]};% define subtraction maps (name, A, B)
% Smap(2, :)={'AO', [1 3 5 7], [2 4 6 8]};    
% Smap(3, :)={'RgAO', [1 3], [2 4]};    
% Smap(4, :)={'LumAO', [5 7], [6 8]};    
% Smap(5, :)={'ARgLum', [1 5], [3 7]};    
% Smap(6, :)={'ORgLum', [2 6], [4 8]}; 
% Smap(7, :)={'Allo', [1:8], [9]};
% Smap(8, :)={'Rgo', [1 2 5 6], [9] };
% Smap(9, :)={'Lumo', [3 4 7 8], [9]};
% Smap(10, :)={'Ao', [1 3 5 7], [9]};
% Smap(11, :)={'Oo', [2 4 6 8], [9]};
% Smap(12, :)={'Rgd', [1 2 5 6], []};
% Smap(13, :)={'Lumd', [3 4 7 8], []};
% Smap(14, :)={'Alld', [1:8], []};
% Smap(15, :)={'blankd', [9], []};
% Smap(16, :)={'L_R', [5:8], [1:4]};

% % % RGLum4
% stim(1)={'RgA'};   % Name stim condistions here
% stim(2)={'RgO'};
% stim(3)={'LumA'};
% stim(4)={'LumO'};
% stim(5)={'blank'};   
% % Subtraction maps
% Smap(1, :)={'RGlum', [1 2], [3 4]};% define subtraction maps (name, A, B)
% Smap(2, :)={'AO', [1 3], [2 4]};    
% Smap(3, :)={'RG-blank', [1 2], [5]};    
% Smap(4, :)={'Lum-blank', [3 4], [5]};    
% Smap(5, :)={'All-blank', [1:4], [5]};


%% RD 8

% % RD dir8 (dir 1-8,  blank 9)
% From Chen Ming
% % Condition name
% 
% stim(1)={'D1'};     
% stim(2)={'D2'};     
% stim(3)={'D3'};     
% stim(4)={'D4'};     
% stim(5)={'D5'};     
% stim(6)={'D6'};     
% stim(7)={'D7'};     
% stim(8)={'D8'};     
% stim(9)={'blank'};   
% 
% Smap(1, :)={'axisHV', [1 5], [3 7]};
% Smap(2, :)={'axisAO', [2 6], [4 8]};
% Smap(3, :)={'axisHA-VO', [1 2 5 6], [3 4 7 8]};
% Smap(4, :)={'axisHO-VA', [1 8 4 5], [2 3 6 7]};
% Smap(5, :)={'RL', [2 1 8], [4 5 6]};
% Smap(6, :)={'UD', [2 3 4], [6 7 8]};
% Smap(7, :)={'1-5', [1], [5]};
% Smap(8, :)={'2-6', [2], [6]};
% Smap(9, :)={'3-7', [3], [7]};
% Smap(10, :)={'4-8', [4], [8]};
% Smap(11, :)={'1o', [1], [9]};
% Smap(12, :)={'2o', [2], [9]};
% Smap(13, :)={'3o', [3], [9]};
% Smap(14, :)={'4o', [4], [9]};
% Smap(15, :)={'5o', [5], [9]};
% Smap(16, :)={'6o', [6], [9]};
% Smap(17, :)={'7o', [7], [9]};
% Smap(18, :)={'8o', [8], [9]};
% Smap(19, :)={'axisHo', [1 5], [9]};
% Smap(20, :)={'axisAo', [2 6], [9]};
% Smap(21, :)={'axisVo', [3 7], [9]};
% Smap(22, :)={'axisOo', [4 8], [9]};
% Smap(23, :)={'Ro', [2 1 8], [9]};
% Smap(24, :)={'Uo', [2 3 4], [9]};
% Smap(25, :)={'Lo', [4 5 6], [9]};
% Smap(26, :)={'Do', [6 7 8], [9]};
% Smap(27, :)={'1c', [1], [1:8]};
% Smap(28, :)={'2c', [2], [1:8]};
% Smap(29, :)={'3c', [3], [1:8]};
% Smap(30, :)={'4c', [4], [1:8]};
% Smap(31, :)={'5c', [5], [1:8]};
% Smap(32, :)={'6c', [6], [1:8]};
% Smap(33, :)={'7c', [7], [1:8]};
% Smap(34, :)={'8c', [8], [1:8]};
% Smap(35, :)={'1o', [1], [9]};
% Smap(36, :)={'2o', [2], [9]};
% Smap(37, :)={'3o', [3], [9]};
% Smap(38, :)={'4o', [4], [9]};
% Smap(39, :)={'5o', [5], [9]};
% Smap(40, :)={'6o', [6], [9]};
% Smap(41, :)={'7o', [7], [9]};
% Smap(42, :)={'8o', [8], [9]};
% Smap(43, :)={'allo', [1:8], [9]};
% Smap(44, :)={'DIr',[1 2 3 4],[5 6 7 8]};
% Smap(45, :)={'blanko0',[9],[]};
% Smap(46, :)={'1357-2468',[1 3 5 7],[2 4 6 8]};

%% epilepsy and MBC

% epilepsy
% Smap(1, :)={'1', [1],[]};
% Smap(2, :)={'2', [2],[]};
% Smap(3, :)={'3', [3],[]};
% Smap(4, :)={'4', [4],[]};
% Smap(5, :)={'5', [5],[]};
% Smap(6, :)={'6', [6],[]};% 
% Smap(7, :)={'7', [7],[]};% 
% Smap(8, :)={'8', [8],[]};% 

% MBClong
% Smap(1, :)={'OD_all', [1 2 4 6], [3 5]};
% Smap(2, :)={'OD_MBC_all', [4 6], [3 5]};
% Smap(3, :)={'OD_MBC_1', [4], [3]};
% Smap(4, :)={'OD_MBC_2', [6], [5]};
% Smap(5, :)={'H_V_MBC_all', [4 5], [3 6]};

% contour 16 
% Smap(1, :)={'OD', [1:2:15], [2:2:16]};% define subtraction maps (name, A, B)

% Bar5
% Smap(1, :)={'1', [1],[]};
% Smap(2, :)={'2', [2],[]};
% Smap(3, :)={'3', [3],[]};
% Smap(4, :)={'4', [4],[]};
% Smap(5, :)={'5', [5],[]};
% Smap(6, :)={'1-6', [1],[6]};% 
% Smap(7, :)={'2-6', [2],[6]};% 
% Smap(8, :)={'3-6', [3],[6]};% 
% Smap(9, :)={'4-6', [4], [6]};
% Smap(10, :)={'5-6', [5], [6]};
% Smap(11, :)={'1-2', [1], [2]};
% Smap(12, :)={'2-3', [2], [3]};
% Smap(13, :)={'3-4', [3], [4]};
% Smap(14, :)={'4-5', [4], [5]};

%% BR training MBC

% % BR_training_MBC_1207_L21
% %
% Smap(1, :)={'L-R', [1 3 5 7], [2 4 6 8]};
% Smap(2, :)={'A-O', [1 2 5 6], [3 4 7 8]};
% Smap(3, :)={'LargeL-R', [1 3], [2 4]};
% Smap(4, :)={'LargeA-O', [1 2], [3 4]};
% Smap(5, :)={'SmallL-R', [5 7], [6 8]};
% Smap(6, :)={'SmallA-O', [5 6], [7 8]};
% Smap(7, :)={'OD1-2', [5], [8]};
% Smap(8, :)={'OD3-4', [7], [6]};
% Smap(9, :)={'Large-Small', [1 2], [5 6]};
% Smap(10, :)={'1', [1], []};
% Smap(11, :)={'2', [2], []};
% Smap(12, :)={'3', [3], []};
% Smap(13, :)={'4', [4], []};
% Smap(14, :)={'5', [5], []};
% Smap(15, :)={'6', [6], []};
% Smap(16, :)={'7', [7], []};
% Smap(17, :)={'8', [8], []};
% Smap(18, :)={'1-9', [1], [9]};
% Smap(19, :)={'2-9', [2], [9]};
% Smap(20, :)={'3-9', [3], [9]};
% Smap(21, :)={'4-9', [4], [9]};
% Smap(22, :)={'5-9', [5], [9]};
% Smap(23, :)={'6-9', [6], [9]};
% Smap(24, :)={'7-9', [7], [9]};
% Smap(25, :)={'8-9', [8], [9]};
% Smap(26, :)={'1-all', [1], [1:8]};
% Smap(27, :)={'2-all', [2], [1:8]};
% Smap(28, :)={'3-all', [3], [1:8]};
% Smap(29, :)={'4-all', [4], [1:8]};
% Smap(30, :)={'5-all', [5], [1:8]};
% Smap(31, :)={'6-all', [6], [1:8]};
% Smap(32, :)={'7-all', [7], [1:8]};
% Smap(33, :)={'8-all', [8], [1:8]};
% Smap(34, :)={'all-9', [1:8], [9]};
% Smap(35, :)={'SmallR-L', [6 8], [5 7]};
%% BR training

% % BR_training_short1
% % 
% Smap(1, :)={'L-R', [1 4], [2 3]};
% Smap(2, :)={'LA-RO', [1], [2]};
% Smap(3, :)={'LO_RA', [4], [3]};
% Smap(4, :)={'LA-RA', [1], [3]};
% Smap(5, :)={'LO-RO', [4], [2]};
% Smap(6, :)={'G8A-O', [5 6], [7 8]};
% Smap(7, :)={'dir135-315', [5], [6]};
% Smap(8, :)={'dir45-225', [7], [8]};
% Smap(9, :)={'ODA-O', [1 3], [2 4]};
% Smap(10, :)={'LA-LO', [1], [4]};
% Smap(11, :)={'RA_RO', [2], [3]};

% training short 2
% 
% Smap(1, :)={'OD_Ori-RD', [1 4], [2 3]};
% Smap(2, :)={'OD_RD', [5], [6]};
% Smap(3, :)={'OD_plaid', [7], [8]};
% Smap(4, :)={'Ori_A-O', [1 3], [2 4]};
% Smap(5, :)={'plaidngrating_L', [7], [1 4]};
% Smap(6, :)={'plaidngrating_R', [8], [2 3]};
% Smap(7, :)={'plaidngrating_LR', [9], [1:4]};
% Smap(8, :)={'LA-RO', [1], [2]};
% Smap(9, :)={'LA-RA', [1], [3]};
% Smap(10, :)={'LO-RA', [4], [3]};
% Smap(11, :)={'LO-RO', [4], [2]};
% Smap(12, :)={'LA-LO', [1], [4]};
% Smap(13, :)={'RA-RO', [3], [2]};
% 

% BR_training_L21
% % 
% Smap(1, :)={'L-R', [1 3], [2 4]};
% Smap(2, :)={'LA-RO', [1], [2]};
% Smap(3, :)={'LO_RA', [3], [4]};
% Smap(4, :)={'LA-RA', [1], [4]};
% Smap(5, :)={'LO-RO', [3], [2]};
% Smap(6, :)={'G8A-O', [5], [6]};
% Smap(7, :)={'ODA-O', [1 4], [2 3]};
% Smap(8, :)={'LA-LO', [1], [3]};
% Smap(9, :)={'RA_RO', [4], [2]};
% Smap(10, :)={'A_O', [1 4 5], [2 3 6]};
%  Smap(1, :)={'R-L', [2 4], [1 3]};

% BR_training_color_L21
%
% Smap(1, :)={'L-R', [1 3], [2 4]};
% Smap(2, :)={'red-green', [1 4], [2 3]};
% Smap(3, :)={'LredA-greenO', [1], [3]};
% Smap(4, :)={'RredA-greenO', [4], [2]};
% Smap(5, :)={'OD1-4', [1], [4]};
% Smap(6, :)={'OD3-2', [3], [2]};
% Smap(7, :)={'OD1-2', [1], [2]};
% Smap(8, :)={'OD4-3', [4], [3]};
% % Smap(9, :)={'RA_RO', [4], [2]};

% Smap(1, :)={'L', [1], []};
% % BR_training_color6_L22
% %
% Smap(1, :)={'L-R', [1 3], [2 4]};
% Smap(2, :)={'red-green', [1 4 5], [2 3 6]};
% Smap(3, :)={'LredA-greenO', [1], [3]};
% Smap(4, :)={'RredA-greenO', [4], [2]};
% Smap(5, :)={'OD1-4', [1], [4]};
% Smap(6, :)={'OD3-2', [3], [2]};
% Smap(7, :)={'OD1-2', [1], [2]};
% Smap(8, :)={'OD4-3', [4], [3]};
% Smap(9, :)={'DredA-DgreenO', [5], [6]};
% Smap(10, :)={'1', [1], []};
% Smap(11, :)={'2', [2], []};
% Smap(12, :)={'3', [3], []};
% Smap(13, :)={'4', [4], []};
% Smap(14, :)={'5', [5], []};
% Smap(15, :)={'6', [6], []};
% Smap(16, :)={'7', [7], []};
% Smap(17, :)={'1-7', [1], [7]};
% Smap(18, :)={'2-7', [2], [7]};
% Smap(19, :)={'3-7', [3], [7]};
% Smap(20, :)={'4-7', [4], [7]};
% Smap(21, :)={'5-7', [5], [7]};
% Smap(22, :)={'6-7', [6], [7]};
% Smap(23, :)={'1-all', [1], [1:6]};
% Smap(24, :)={'2-all', [2], [1:6]};
% Smap(25, :)={'3-all', [3], [1:6]};
% Smap(26, :)={'4-all', [4], [1:6]};
% Smap(27, :)={'5-all', [5], [1:6]};
% Smap(28, :)={'6-all', [6], [1:6]};
% Smap(29, :)={'all-7', [1:6], [7]};
% Smap(30, :)={'R-L', [2 4], [1 3]};

% % BR_training_with_RD
% %
% Smap(1, :)={'L-R', [1 3 5], [2 4 6]};
% Smap(2, :)={'red-green', [1 4], [2 3]};
% Smap(3, :)={'LredA-greenO', [1], [3]};
% Smap(4, :)={'RredA-greenO', [4], [2]};
% Smap(5, :)={'OD1-4', [1], [4]};
% Smap(6, :)={'OD3-2', [3], [2]};
% Smap(7, :)={'OD1-2', [1], [2]};
% Smap(8, :)={'OD4-3', [4], [3]};
% Smap(9, :)={'RD_L-R', [5], [6]};
% Smap(10, :)={'1', [1], []};
% Smap(11, :)={'2', [2], []};
% Smap(12, :)={'3', [3], []};
% Smap(13, :)={'4', [4], []};
% Smap(14, :)={'5', [5], []};
% Smap(15, :)={'6', [6], []};
% Smap(16, :)={'7', [7], []};
% Smap(17, :)={'1-7', [1], [7]};
% Smap(18, :)={'2-7', [2], [7]};
% Smap(19, :)={'3-7', [3], [7]};
% Smap(20, :)={'4-7', [4], [7]};
% Smap(21, :)={'5-7', [5], [7]};
% Smap(22, :)={'6-7', [6], [7]};
% Smap(23, :)={'1-all', [1], [1:6]};
% Smap(24, :)={'2-all', [2], [1:6]};
% Smap(25, :)={'3-all', [3], [1:6]};
% Smap(26, :)={'4-all', [4], [1:6]};
% Smap(27, :)={'5-all', [5], [1:6]};
% Smap(28, :)={'6-all', [6], [1:6]};
% Smap(29, :)={'all-7', [1:6], [7]};
% Smap(30, :)={'gratingL-R', [1 3], [2 4]};
% Smap(31, :)={'R-L', [2 4 6], [1 3 5]};

% 120104_L15
% BR_training_6+6 easy subtraction
% Smap(1, :)={'L-R_s', [1 3], [2 4]};
% Smap(2, :)={'L-R_b', [7 9], [8 10]};
% Smap(3, :)={'A-O_s', [1 4 5], [2 3 6]};
% Smap(4, :)={'A-O_b', [7 10 11], [8 9 12]};

% 
% % BR_training_color4_from_balancer
% %
% Smap(1, :)={'L-R', [1 3], [2 4]};
% Smap(2, :)={'red-green', [1 4], [2 3]};
% Smap(3, :)={'LredA-greenO', [1], [3]};
% Smap(4, :)={'RredA-greenO', [4], [2]};
% Smap(5, :)={'OD1-4', [1], [4]};
% Smap(6, :)={'OD3-2', [3], [2]};
% Smap(7, :)={'OD1-2', [1], [2]};
% Smap(8, :)={'OD4-3', [4], [3]};

% OD mask from BR_training_color4

% Smap(1, :)={'L-R', [1], [2]};
% Smap(2, :)={'R-L', [2], [1]};
% Smap(3, :)={'L+D-blank', [1], [5]};
% Smap(4, :)={'R+D-blank', [2], [5]};
% Smap(5, :)={'L+D', [1], []};
% Smap(6, :)={'R+D', [2], []};
% Smap(7, :)={'OD1-2', [1], [2]};
% Smap(8, :)={'OD4-3', [4], [3]};




%% BR color balancer

% BR_color_balancer_L20
%
% Smap(1, :)={runname(end-11:end-9), [1], [2]};
% Smap(2, :)={runname(end-7:end-5), [3], [4]};
% Smap(3, :)={runname(end-3:end-1), [5], [6]};
% % 
% Smap(1, :)={'1-2', [1], [2]};
% Smap(2, :)={'1-4', [1], [4]};
% Smap(3, :)={'1-6', [1], [6]};
% Smap(4, :)={'3-2', [3], [2]};
% Smap(5, :)={'3-4', [3], [4]};
% Smap(6, :)={'3-6', [3], [6]};
% Smap(7, :)={'5-2', [5], [2]};
% Smap(8, :)={'5-4', [5], [4]};
% Smap(9, :)={'5-6', [5], [6]};
% Smap(10, :)={'1-3', [1], [3]};
% Smap(11, :)={'1-5', [1], [5]};
% Smap(12, :)={'3-5', [3], [5]};
% Smap(13, :)={'2-4', [2], [4]};
% Smap(14, :)={'2-6', [2], [6]};
% Smap(15, :)={'4-6', [4], [6]};
% Smap(16, :)={'1-0', [1], [7]};
% Smap(17, :)={'2-0', [2], [7]};
% Smap(18, :)={'3-0', [3], [7]};
% Smap(19, :)={'4-0', [4], [7]};
% Smap(20, :)={'5-0', [5], [7]};
% Smap(21, :)={'6-0', [6], [7]};

%% BR contrast

% BR_contrast
% 
% Smap(1, :)={'L-R', [2 3 7 8], [4 5 9 10]};
% Smap(2, :)={'L-R_strong', [3 8], [5 10]};
% Smap(3, :)={'L-R_mild', [2 7], [4 9]};
% Smap(4, :)={'L-R_equal', [1], [6]};
% Smap(5, :)={'A-O', [2 3 9 10], [4 5 7 8]};
% Smap(6, :)={'A-O_strong', [3 10], [5 8]};
% Smap(7, :)={'A-O_mild', [2 9], [4 7]};
% Smap(8, :)={'A-O_equal', [1], [6]};
% Smap(9, :)={'LA-RA_strong', [3], [10]};
% Smap(10, :)={'LO-RO_strong', [8], [5]};
% Smap(11, :)={'LA-RO_strong', [3], [5]};
% Smap(12, :)={'LO-RA_strong', [8], [10]};
% Smap(13, :)={'LA-RO_mild', [2], [4]};
% Smap(14, :)={'LO-RA_mild', [7], [9]};
% Smap(15, :)={'LA-RA_mild', [2], [9]};
% Smap(16, :)={'LO-RO_mild', [7], [4]};
% Smap(17, :)={'LA-LO_strong', [3], [8]};
% Smap(18, :)={'RA-RO_strong', [10], [5]};
% Smap(19, :)={'LA-LO_mild', [2], [7]};
% Smap(20, :)={'RA-RO_mild', [9], [4]};
% Smap(21, :)={'3-0', [3], [11]};
% Smap(22, :)={'5-0', [5], [11]};
% Smap(23, :)={'8-0', [8], [11]};
% Smap(24, :)={'10-0', [10], [11]};
% Smap(25, :)={'1-0', [1], [11]};
% Smap(26, :)={'6-0', [6], [11]};




%% OD 
% OD8
% Condition name
% stim(1)={'LH'};   % Name stim condistions here
% stim(2)={'LO'};
% stim(3)={'LV'};
% stim(4)={'LA'};
% stim(5)={'RH'};
% stim(6)={'RO'};
% stim(7)={'RV'};
% stim(8)={'RA'};
% stim(9)={'Blank'};
% % % Subtraction maps
% Smap(1, :)={'LRall', [1 2 3 4], [5 6 7 8]};% define subtraction maps (name, A, B)
% Smap(2, :)={'LRh', [1], [5]};    
% Smap(3, :)={'LRa', [2], [6]};    
% Smap(4, :)={'LRv', [3], [7]};    
% Smap(5, :)={'LRo', [4], [8]};    
% Smap(6, :)={'HV', [1 5], [3 7]};
% Smap(7, :)={'HVle', [1], [3]};
% Smap(8, :)={'HVre', [5], [7]};
% Smap(9, :)={'AO', [2 6], [4 8]};
% Smap(10, :)={'AOle', [2], [4]};
% Smap(11, :)={'AOre', [6], [8]};
% Smap(12, :)={'HA-VO', [1 2 5 6], [3 4 7 8]};
% Smap(13, :)={'HO-VA', [1 4 5 8], [2 3 6 7]};
% Smap(14, :)={'LEo', [1 2 3 4], [9]};
% Smap(15, :)={'LEd', [1 2 3 4], []};
% Smap(16, :)={'REo', [5 6 7 8], [9]};
% Smap(17, :)={'REd', [5 6 7 8], []};
% Smap(18, :)={'allo', [1 2 3 4 5 6 7 8], [9]};
% Smap(19, :)={'alld', [1 2 3 4 5 6 7 8], []};
% Smap(20, :)={'H+Vo', [1 3 5 7], [9]};
% Smap(21, :)={'H+Vd', [1 3 5 7], []};
% Smap(22, :)={'A+Oo', [2 4 6 8], [9]};
% Smap(23, :)={'A+Od', [2 4 6 8], []};
% Smap(24, :)={'Hd', [1 5], []};
% Smap(25, :)={'Ho', [1 5], [9]};
% Smap(26, :)={'Ad', [2 6], []};
% Smap(27, :)={'Ao', [2 6], [9]};
% Smap(28, :)={'Vd', [3 7], []};
% Smap(29, :)={'Vo', [3 7], [9]};
% Smap(30, :)={'Od', [4 8], []};
% Smap(31, :)={'Oo', [4 8], [9]};

% Smap(1, :)={'LH-All', [1], [1:9]};
% Smap(2, :)={'LO-All', [2], [1:9]};
% Smap(3, :)={'LV-All', [3], [1:9]};
% Smap(4, :)={'LA-All', [4], [1:9]};
% Smap(5, :)={'RH-All', [5], [1:9]};
% Smap(6, :)={'RO-All', [6], [1:9]};
% Smap(7, :)={'RV-All', [7], [1:9]};
% Smap(8, :)={'RA-All', [8], [1:9]};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OD12 %
% Condition name
% stim(1)={'LH'};   % Name stim condistions here
% stim(2)={'LO'};
% stim(3)={'LV'};
% stim(4)={'LA'};
% stim(5)={'RH'};
% stim(6)={'RO'};
% stim(7)={'RV'};
% stim(8)={'RA'};
% stim(9)={'BH'};
% stim(10)={'BO'};
% stim(11)={'BV'};
% stim(12)={'BA'};
% stim(13)={'Blank'};
% % Subtraction maps
% Smap(1, :)={'LHV', [1], [3]};% define subtraction maps (name, A, B)
% Smap(2, :)={'RHV', [5], [7]};    
% Smap(3, :)={'BHV', [9], [11]};    
% Smap(4, :)={'LRHV', [1 5], [3 7]};    
% Smap(5, :)={'HV', [1 5 9], [3 7 11]};    
% Smap(6, :)={'BiMo', [9:12], [1:8]};
% Smap(7, :)={'BiL', [9:12], [1:4]};
% Smap(8, :)={'BiR', [9:12], [5:8]};
% Smap(9, :)={'HV', [1 5 9], [3 7 11]};
% Smap(10, :)={'HVl', [1], [3]};
% Smap(11, :)={'HVr', [5], [7]};
% Smap(12, :)={'HVb', [9], [11]};
% Smap(13, :)={'AO', [2 6 10], [4 8 12]};
% Smap(14, :)={'AOl', [2], [4]};
% Smap(15, :)={'AOr', [6], [8]};
% Smap(16, :)={'AOb', [10], [12]};
% Smap(17, :)={'HA-VO', [1 2 5 6 9 10], [3 4 7 8 11 12]};
% Smap(18, :)={'HO-VA', [1 4 5 8 9 12], [2 3 6 7 10 11]};
% Smap(19, :)={'LEo', [1 2 3 4], [13]};
% Smap(20, :)={'LEd', [1 2 3 4], []};
% Smap(21, :)={'REo', [5 6 7 8], [13]};
% Smap(22, :)={'REd', [5 6 7 8], []};
% Smap(23, :)={'allo', [1 2 3 4 5 6 7 8 9 10 11 12], [13]};
% Smap(24, :)={'alld', [1 2 3 4 5 6 7 8 9 10 11 12], []};
% Smap(25, :)={'H+Vo', [1 3 5 7 9 11], [13]};
% Smap(26, :)={'H+Vd', [1 3 5 7 9 11], []};
% Smap(27, :)={'A+Oo', [2 4 6 8 10 12], [13]};
% Smap(28, :)={'A+Od', [2 4 6 8 10 12], []};
% Smap(29, :)={'Hd', [1 5 9], []};
% Smap(30, :)={'Ho', [1 5 9], [13]};
% Smap(31, :)={'Ad', [2 6 10], []};
% Smap(32, :)={'Ao', [2 6 10], [13]};
% Smap(33, :)={'Vd', [3 7 11], []};
% Smap(34, :)={'Vo', [3 7 11], [13]};
% Smap(35, :)={'Od', [4 8 12], []};
% Smap(36, :)={'Oo', [4 8 12], [13]};
% Smap(37, :)={'HV-AOa', [1 3 5 7 9 11], [2 4 6 8 10 12]};
% Smap(38, :)={'HV-AOm', [1 3 5 7], [2 4 6 8]};
% Smap(39, :)={'HV-AOb', [9 11], [10 12]};
% Smap(40, :)={'L-R', [1:4], [5:8]};
% Smap(41, :)={'HV-AOb', [9 11], [10 12]};
% Smap(42, :)={'HV-AOb', [9 11], [10 12]};
% Smap(43, :)={'HV-AOb', [9 11], [10 12]};
% Smap(44, :)={'R-L', [5:8], [1:4]};
% Smap(45, :)={'L-All', [1:4], [1:13]};
% Smap(46, :)={'R-All', [5:8], [1:13]};

% %% plot eye
% Smap(1, :)={'1', [1], [3]};% define subtraction maps (name, A, B)
% Smap(2, :)={'2', [1], [4]};    
%% plot eye
% Smap(1, :)={'1', [1], []};% define subtraction maps (name, A, B)
% Smap(2, :)={'2', [2], []};    
% Smap(3, :)={'3', [3], []};    
% Smap(4, :)={'4', [4], []};    
% Smap(5, :)={'5', [5], []};    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datasource=1;       % 1: from original OI block files, 2: from previously saved single-condition 'ivf' files.
fframe = [1:2];     % first frame range e.g. [], [1], [1, 2] or [1:3]
sumrange=[9:16];    % sum frame range e.g. [5 6 7 8 9 10 11 12] or [5:12]
% fframe = [1:4];     % first frame range e.g. [], [1], [1, 2] or [1:3]
% sumrange=[22:39];    % sum frame range e.g. [5 6 7 8 9 10 11 12] or [5:12]
operation=2;        % specify type of pixel value: 0: raw DC, 1: dR, 2:dR/R, 3:dR/Rb, 4:dR/R0end
clipmethod = 0;     % 0: no clipping; 1: clipping to +-SD (value);  2: clipping to +-SD (value) with a mask (default.bmp);  3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
clipvalue  = 1.5;     % How many SD to be used for clipping, usually =1 (range is plus and minus 1SD on both sides of median)	% clipvalue = [10, 10; 20, 20; 0.8, 0];   %for method 3 only (window cliping), represents [x1, y1; x2, y2; sd, nouse]
flagmasking=0;      % 0: no masking in average value calculation, '1' for map loading, '2' for coordinates loading
maskname=[          % Names of the masks (need to be same length), these names will be added '.bmp' for loading mask map or add '_x.txt'/'_y.txt' for loading coordinates.
    ];
domainradius=[      % Can specify different size for different domains (e.g. V2 orientation domains are larger than V1's), default 5 pixels (leave empty)
];
flagmap=3;          % 0: do not generate any maps, 1: generate only block-averaged Smaps, 2: generate both block-avg Smaps and block-avg single-condition maps, 3 and 4: more single-block maps, see readme
flagsp=0;           % '0': no superpixel (ROI) time course output, '1': for mask loading, '2' for coordinates loading (for '1' and '2' need change 'flagmasking' to the same value as 'flagsp')
superpixel=[0];		% This variable need to be removed, set it same as 'flagsp' for now.
FrameRate=4;        % Hz, usually 4Hz for VDAQ and 7Hz for RedShirt, not critical for analysis, just make output sp file nicer
% ext = 'a6';       % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)
% ext = 'gs';       % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)
ext = '';       % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)

%---------------- More Advanced Options (default: 0 or []) -------------%
flagrandom = 0;     % 1: if is awake & randomnized data (requires '_stimseq.txt' at block file folder. Usually 0
flagperform= 0;     % 1: if is awake data (requires '_performance.txt' at block file folder), usually 0
flaggoodstim = 0;   %align% 1: will looking for 'goodstim.txt' which contains '1's for good stim conditions, and used these for average map & quantification
saveaccum=0;        % 1: save the temporary accumulate maps, usually =0
accummapnum=[];     % Smap number to be saved accumulatly, no single condition map will be saved this way to save space.
tempsaveblocks = 0; % how often save the temp analysis maps, eg. 5 means save results every 5 blocks. % set tempsaveblocks=999 to avoid save tempmaps
blockselect =[];    % select blocks for processing, eg. [ 1 3 4 5], leave empty for all block being selected ([]).
flagsaveivf=1;      % 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
flagsaveblk=0;      % 1: save averaged data into a block file (under testing). 0: not save
flagquantify=0;     % 1: quantify map value based on masks (fun1-fun3.txt output), 2: detailed quantify (blk by blk, slow), 0: no quantification
flagtrialfilter=0;	% 0: no trial by trial filtering; 1: filtering subtraction map for each trial before summation.
flaghpfilter=0;		% for average maps only (takes time) 0: no high-pass filter, 1: highpass filter based on data
hpfmethod='ribot';	% filter kernel type, see 'filtermethod' below in vector section, default 'fastmean'
hpkernelsize=2;	    % diameter in pixels. 50-200 depends on baseline-noise size. slow for large kernel. For 'ribot' filter, it's order of polynomial (usually 2 or 3)
flaglpfilter=0;		% for average maps only 0: no low-pass filter, 1: lowpass filter based on data
lpfmethod='gaussian'; % filter kernel type, see 'filtermethod' below in vector section, default 'gaussian'
lpkernelsize=5;		% diameter in pixels. 5-10 depends on high frequency noise size. 

% Alignment, do not change if not do shift/alignment
flagalign=0;       %align% 0: no shift/align;  1: shift/align all frames with asigned frame; 2: shift/align within stim frames, 11: load shift parameters (shiftinput.txt in run folder);
shiftmethod = 41;   % 1: coorleation method, 2: min difference method, 3:masked, 4:normxcorr2 (fast)
shiftrange = [-9 9 9 9];  % how many pixels to shift to search for best fit, towards left, right, up, down
shiftstep1 = 3;    % Shift by 'step1' first, then use finer shift 'step2' to find better alignment around the result from first shift.
shiftstep2 = 0;    % ==0 if no fine shift is needed
firstframe = 'els_E00B050.blk';    % read the very first frame for shifting (some case the first block is not in file array), or you want to use a frame in the middle. only for flag 1,2

%---------------- Vector analysis --------------------------------%
flagvector=0;   % 0: no vector analysis, 1: vector map from raw data, 2: vector map from previous analyzed single condition maps (i.e. only if you have run vector analysis before, and you have single condition map in 'vector\' folder).
% specify single condition maps used for vector analysis:
vect(1,:)={'H', [1 5 9], []};     % e.g. 'Horizontal orientation' has [1 5 9]: left, right, both eye horizontal. the 2nd [] is usually empty.
vect(2,:)={'A', [2 6 10], []};
vect(3,:)={'V', [3 7 11], []};
vect(4,:)={'O', [4 8 12], []};
lowpass=0;     % 15? low-pass filter kernel size (diameter in pixels, typical: ~7s-20 for size 300-500)
highpass=0;	% 80? high-pass filter kenel size; (typical: ~80 for 504x504)
filtermethod='gaussian';    % choose from 'fastmean', 'slowmean', 'gaussian'
                            % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
                            % 'slowmean' 'disk-like' but better at edge
                            % 'gaussian' Gaussian filter with half sd
                            % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
                              
%---------------- Do not change following lines ------------------%
sunincore34;
clear all
return;
