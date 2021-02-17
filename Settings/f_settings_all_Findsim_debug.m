function [stg] = f_settings_all_Findsim()

%% Import

% True or false to decide whether to run import functions
% (Import)
 stg.import = true;

% Name of the folder where everything related to the model is stored
% (Folder Model)
stg.folder_model = "Model_Findsim_debug";

% Name of the excel file with the sbtab
% (SBtab excel name)
stg.sbtab_excel_name = "SBTAB_Findsim.xlsx";

% Name of the model
% (Name)
stg.name = "Findsim";

% Name of the default model compartment
% (Compartment name)
stg.cname = "Cell";

% Name of the sbtab saved in .mat format
% (SBtab name)
stg.sbtab_name = "sbtab_" + stg.name;

%% Analysis

% String with the analysis to be run, the options are "diag",
% "opt", "SA"
% and can be combined as for example "RS,diag", to not run any analysis set
% stg.analysis to equal to ""
% (Analysis)

stg.analysis = "diag";

% Experiments to run
% stg.ms.exprun = [1,3,4];
stg.exprun = [1];

% Choice between 0,1,2 and 3 to change either and how to apply log10 to the
% scores (check documentation)
% (Use logarithm)
stg.useLog = 4;

% True or false to decide whether to use multicore everywhere it is available
% (Optimization Multicore)
stg.optmc = false;

% Choice of ramdom seed
% (Ramdom seed)
stg.rseed = 2;

% True or false to decide whether to save results
% (Save results)
stg.save_results = false;

% True or false to decide whether to run detailed simulation for plots
% (Save results)
stg.simdetail = true;

%% Simulation

% Maximum time for each individual function to run in seconds
% (Maximum time)
stg.maxt = 30;

% Equilibration time
% (Equilibration time)
stg.eqt  = 5000;

% True or false to decide whether to do Dimensional Analysis
% (Dimensional Analysis)
stg.dimenanal = true;

% True or false to decide whether to do Unit conversion
% (Unit conversion)
stg.UnitConversion = true;

% Value of Relative tolerance
% (Relative tolerance)
stg.reltol = 1.0E-4;

% Value of Absolute tolerance
% (Absolute tolerance)
stg.abstol = 1.0E-4;

% Time units for simulation
% (Simulation time)
stg.simtime = "second";

% True or false to decide whether to run sbioaccelerate (after changing this value
% you need to run "clear functions" to see an effect)
% (sbioaccelerate)
stg.sbioacc = true;

% Max step size in the simulation (if empty matlab decides whats best)
% (Maximum step)
stg.maxstep = [];

% Max step size in the equilibration (if empty matlab decides whats best)
% (Maximum step)
stg.maxstepeq = [];

% Max step size in the detailed plots (if empty matlab decides whats best)
% (Maximum step)
stg.maxstepdetail = [10];

%% Model

% Number of parameters to optimize
% (Parameter number)
% stg.parnum = 15;

% stg.ms.parnum = 5;

% stg.parnum = 160;
% test = [2.77814998	-0.301029996	0.903081857	0.936252283	0.103974669	0.544390543	0	-1	0.079184548	-1	0.301029996	-0.698970004	0	1.698970004	-0.920815452	-1	-0.22185002	-1	1.255268186	0.301029996	2.99999873	-1	0.477119984	0.602059991	-0.397940009	3.77814998	-1	-0.769551079	3.477119984	0	0.397940009	-0.823908741	4.477119984	0	3.999968151	-0.301029996	0	-4	-1	0	3.623249327	-0.602059991	-2.698970004	-3.48148606	-0.698970004	2.698952272	-1	1.397944226	-1.774690718	-3	1.397944226	-1.774690718	2.99999873	0	5.255268186	1	4.079184548	1	-1.301029996	-1.15490196	2.977624806	1.204119983	0.602059991	3.752552377	1.204119983	0.602059991	3.752552377	1.204119983	0.602059991	3.104734896	1.334453751	0.73239376	3.415310209	1.645029065	1.042969073	3.928643636	2.158362492	1.556302501	4.150492386	2.380211242	1.77815125	4.451522382	2.681241237	2.079181246	3.375801668	1.602059991	1	4.655642364	2.283301229	1.681241237	3.565465734	1.903089987	1.301029996	3.264435738	1.602059991	1	4.264435738	1.602059991	1	4.183916142	-0.22184875	-0.823908741	4.183916142	-0.22184875	-0.823908741	3.948847023	0.079181246	-0.522878745	3.948847023	0.079181246	-0.522878745	0.976437817	-1.096910013	-1.698970004	2.271066317	-1.096910013	-1.698970004	4.66888375	1.602059991	1	3.271066317	-0.096910013	-0.698970004	3.451522382	-0.096910013	-0.698970004	3.053582373	-0.096910013	-0.698970004	3.271066317	-0.096910013	-0.698970004	1.146128036	-1.012237759	1.755874856	-1.703520935	3.901628047	0.301029996	-0.301029996	3.43373423	1.204119983	0.602059991	3.43373423	1.204119983	0.602059991	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.256828796	1.380211242	0.77815125];
% test = [-0.22185002	-0.301029996	-2.096918143	0.936252283	0.103974669	0.544390543	0	-1	-2.920815452	-1	0.301029996	-0.698970004	0	1.698970004	-3.920815452	-1	-3.22185002	-1	-1.744731814	0.301029996	-1.27031E-06	-1	-2.522880016	0.602059991	-0.397940009	0.77814998	-1	-0.769551079	0.477119984	0	0.397940009	-0.823908741	1.477119984	0	0.999968151	-0.301029996	0	-4	-1	0	0.623249327	-0.602059991	-2.698970004	-3.48148606	-0.698970004	-0.301047728	-1	-1.602055774	-1.774690718	-3	-1.602055774	-1.774690718	-1.27031E-06	0	2.255268186	1	1.079184548	1	-1.301029996	-1.15490196	-0.022375194	1.204119983	0.602059991	0.752552377	1.204119983	0.602059991	0.752552377	1.204119983	0.602059991	0.104734896	1.334453751	0.73239376	0.415310209	1.645029065	1.042969073	0.928643636	2.158362492	1.556302501	1.150492386	2.380211242	1.77815125	1.451522382	2.681241237	2.079181246	0.375801668	1.602059991	1	1.655642364	2.283301229	1.681241237	0.565465734	1.903089987	1.301029996	0.264435738	1.602059991	1	1.264435738	1.602059991	1	1.183916142	-0.22184875	-0.823908741	1.183916142	-0.22184875	-0.823908741	0.948847023	0.079181246	-0.522878745	0.948847023	0.079181246	-0.522878745	-2.023562183	-1.096910013	-1.698970004	-0.728933683	-1.096910013	-1.698970004	1.66888375	1.602059991	1	0.271066317	-0.096910013	-0.698970004	0.451522382	-0.096910013	-0.698970004	0.053582373	-0.096910013	-0.698970004	0.271066317	-0.096910013	-0.698970004	1.146128036	1.987762241	1.755874856	1.296479065	0.901628047	0.301029996	-0.301029996	0.43373423	1.204119983	0.602059991	0.43373423	1.204119983	0.602059991	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.256828796	1.380211242	0.77815125];

%Correct parameters
stg.parnum = 29;

test = [0.999968151	-0.301029996	0	-4	-1	0.264435738	1.602059991	1	1.183916142	-0.22184875	-0.823908741	1.183916142	-0.22184875	-0.823908741	0.948847023	0.079181246	-0.522878745	0.948847023	0.079181246	-0.522878745	-2.023562183	-1.096910013	-1.698970004	-0.728933683	-1.096910013	-1.698970004	1.66888375	1.602059991	1];

% test = [3.565465734	1.903089987	1.301029996	3.264435738	1.602059991	1	4.264435738	1.602059991	1	4.183916142	-0.22184875	-0.823908741 3.43373423	1.204119983	0.602059991];

% Array with the lower bound of all parameters
% (Lower bound)
stg.lb = test-1;
% stg.lb = zeros(1,stg.parnum)-5;

% Array with the upper bound of all parameters
% (Upper bound)
stg.ub = test+1;
% stg.ub = zeros(1,stg.parnum)+5;

%% Diagnostics

% if ispc
%     load("Model\" + stg.folder_model +"\Results\bestx.mat",'bestx');
%     load("Model\" + stg.folder_model +"\Results\padiag.mat",'pa');
% else
%     load("Model/" + stg.folder_model +"/Results/bestx.mat",'bestx');
%     load("Model/" + stg.folder_model +"/Results/padiag.mat",'pa');
% end

% Choice of what parameters in the array to test, the indices correspond to
% the parameters in the model and the numbers correspond to the parameters
% in the optimization array, usually not all parameters are optimized so
% there needs to be a match between one and the other.
% (Parameters to test)
% stg.partest(:,1) = [1,zeros(1,159)];

% stg.partest(:,1) = [zeros(1,90),1,2,3,4,5,6,7,8,9,10,11,12,10,11,12,zeros(1,34),13,14,15,13,14,15,zeros(1,15)];
% stg.partest(:,1) = [1:160];

%Correct parameters
stg.partest(:,1) = [zeros(1,34),1,2,3,4,5,zeros(1,54),6,7,8,zeros(1,3),9:29,zeros(1,40)];

% stg.ms.partest(:,1) = [0  ,0  ,0  ,0  ,0  ,0  ,0 ,0 ,1 ,2,...
%                        3  ,4  ,5];
% (Parameter array to test)

stg.pat = [1:1];

% All the parameter arrays, in this case there is only one
% (Parameter arrays)

%  stg.pa(5,:) = [4.34356751226841,0.923670780241194,0.673875847564401,3.45149572954397,2.50559833769720,0.00415615022018036,5.24486891989970,1.52326843941571,1.99064830038686,3.88062066075096,0.124508238614293,0.173288877821028,3.60553925673816,1.38243382569349,1.39831064384363];
% 
%  stg.pa(1,:) = [3.565465734	1.903089987	1.301029996	3.264435738	1.602059991	1	4.264435738	1.602059991	1	4.183916142	-0.22184875	-0.823908741 3.43373423	1.204119983	0.602059991];
% %  stg.pa(2,:) = [1.47051661397124,-4.58771378327779,4.88472131798021,-4.47919212582039,-0.793785681798660,1.71742167512146,3.51108966232211,4.25446620406408,3.31352552011681,3.98680623624496,4.32884031430711,4.07169863187831,4.70715213220524,4.30810875901718,2.42131569665967];
% %  stg.pa(3,:) = [-5,4.12135585802449,4.99146911624014,-2.54227308884864,2.54704499352413,-2.02597328941843,1.85757330256128,-3.04266952712744,4.20592980948317,4.65444168720118,-4.98129815674133,4.88305072306658,4.22788671025649,5,4.99997935458053];
%  stg.pa(2,:) = [3.00782327965918,4.44103049270336,4.98689658505079,-4.46679183533548,4.83522766774498,1.79750232381264,3.57310348907238,-3.59344694743980,4.80528800453270,3.36780261465713,2.60490617774707,4.94602680659227,2.63969992012040,-3.92136714482799,2.35542805692282];
% %  stg.pa(5,:) = [-4.26279353805345,2.88643200396286,4.41496882294233,0.123672432876000,4.09629756048136,4.34966183555361,2.65916568571896,1.36526071837012,2.20182274204792,3.79352258710029,-3.46305900574538,-0.0994467635636818,3.14330662044877,3.64289011525887,4.33806057175983];
% %  stg.pa(6,:) = [2.69719440794330,0.183113486119264,4.80716764448831,-4.07115143573228,1.37733833852969,-0.128928118114379,3.81317398116445,-1.56434352721003,0.546332576501559,4.32332566297067,-3.50440094711640,1.45818898376434,3.59073083362675,4.03069341930590,4.22839598948340];
%  stg.pa(3,:) = [4.39041799990433,-0.485184857967128,-4.72079121923554,4.97203264419559,4.18950214248611,-3.73493921705331,3.54203779163694,-4.18713079085033,-4.90627884689042,4.99999785495303,-1.78844283387361,1.76402762860888,4.70494274198218,-0.483157240386130,-0.911361474223092];
%  stg.pa(4,:) = [1.11389406834506,-0.722033566206175,3.17293462462387,1.40396873098917,-4.02336305666310,-4.24983854952258,3.87893750075415,-4.26923209185099,-5,5,-4.82474735566891,3.18705402011508,4.34754000362115,-3.44834931956091,-0.835477432766747];
%  stg.pa(6,:) = [3.30733573938338,1.12838826776291,2.16219714471150,2.57377638803188,1.24491683896527,0.925384094311856,5.25047673117878,1.90811367913012,1.17552933268312,4.49106982875976,0.232185344435488,-0.513195559242825,3.57449699250194,1.12639725433913,1.15400443014967];
%  stg.pa(7,:) = [3.61149652874897,1.24861468886348,0.942914395672525,3.59208166553196,1.28106748928935,0,4.33696003211232,0.728229078256353,1.42451976653212,3.96610121450282,-1.21570312448490,0.172054293858328,4.43373423000000,2.02721036487820,1.01005126553146];
% 
%  
%  stg.pa(2,:) = [4.29773286146038,2.49693911369044,1.73196980711270,4.08282603953207,2.43462074402764,1.43217333100689,4.75965682330789,0.797908727278588,0.920556213396455,4.20628288522142,0.0235857240088855,0.0395309828353139,3.60345957798613,0.791998627014539,0.984708301857106];
 
% All parameters
% for n = 1:10
% stg.pa(1,:) = [2.77814998	-0.301029996	0.903081857	0.936252283	0.103974669	0.544390543	0	-1	0.079184548	-1	0.301029996	-0.698970004	0	1.698970004	-0.920815452	-1	-0.22185002	-1	1.255268186	0.301029996	2.99999873	-1	0.477119984	0.602059991	-0.397940009	3.77814998	-1	-0.769551079	3.477119984	0	0.397940009	-0.823908741	4.477119984	0	3.999968151	-0.301029996	0	-4	-1	0	3.623249327	-0.602059991	-2.698970004	-3.48148606	-0.698970004	2.698952272	-1	1.397944226	-1.774690718	-3	1.397944226	-1.774690718	2.99999873	0	5.255268186	1	4.079184548	1	-1.301029996	-1.15490196	2.977624806	1.204119983	0.602059991	3.752552377	1.204119983	0.602059991	3.752552377	1.204119983	0.602059991	3.104734896	1.334453751	0.73239376	3.415310209	1.645029065	1.042969073	3.928643636	2.158362492	1.556302501	4.150492386	2.380211242	1.77815125	4.451522382	2.681241237	2.079181246	3.375801668	1.602059991	1	4.655642364	2.283301229	1.681241237	3.565465734	1.903089987	1.301029996	3.264435738	1.602059991	1	4.264435738	1.602059991	1	4.183916142	-0.22184875	-0.823908741	4.183916142	-0.22184875	-0.823908741	3.948847023	0.079181246	-0.522878745	3.948847023	0.079181246	-0.522878745	0.976437817	-1.096910013	-1.698970004	2.271066317	-1.096910013	-1.698970004	4.66888375	1.602059991	1	3.271066317	-0.096910013	-0.698970004	3.451522382	-0.096910013	-0.698970004	3.053582373	-0.096910013	-0.698970004	3.271066317	-0.096910013	-0.698970004	1.146128036	-1.012237759	1.755874856	-1.703520935	3.901628047	0.301029996	-0.301029996	3.43373423	1.204119983	0.602059991	3.43373423	1.204119983	0.602059991	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.271066317	1.397940009	0.77815125	3.256828796	1.380211242	0.77815125];
% stg.pa(1,:) = [-0.22185002	-0.301029996	-2.096918143	0.936252283	0.103974669	0.544390543	0	-1	-2.920815452	-1	0.301029996	-0.698970004	0	1.698970004	-3.920815452	-1	-3.22185002	-1	-1.744731814	0.301029996	-1.27031E-06	-1	-2.522880016	0.602059991	-0.397940009	0.77814998	-1	-0.769551079	0.477119984	0	0.397940009	-0.823908741	1.477119984	0	0.999968151	-0.301029996	0	-4	-1	0	0.623249327	-0.602059991	-2.698970004	-3.48148606	-0.698970004	-0.301047728	-1	-1.602055774	-1.774690718	-3	-1.602055774	-1.774690718	-1.27031E-06	0	2.255268186	1	1.079184548	1	-1.301029996	-1.15490196	-0.022375194	1.204119983	0.602059991	0.752552377	1.204119983	0.602059991	0.752552377	1.204119983	0.602059991	0.104734896	1.334453751	0.73239376	0.415310209	1.645029065	1.042969073	0.928643636	2.158362492	1.556302501	1.150492386	2.380211242	1.77815125	1.451522382	2.681241237	2.079181246	0.375801668	1.602059991	1	1.655642364	2.283301229	1.681241237	0.565465734	1.903089987	1.301029996	0.264435738	1.602059991	1	1.264435738	1.602059991	1	1.183916142	-0.22184875	-0.823908741	1.183916142	-0.22184875	-0.823908741	0.948847023	0.079181246	-0.522878745	0.948847023	0.079181246	-0.522878745	-2.023562183	-1.096910013	-1.698970004	-0.728933683	-1.096910013	-1.698970004	1.66888375	1.602059991	1	0.271066317	-0.096910013	-0.698970004	0.451522382	-0.096910013	-0.698970004	0.053582373	-0.096910013	-0.698970004	0.271066317	-0.096910013	-0.698970004	1.146128036	1.987762241	1.755874856	1.296479065	0.901628047	0.301029996	-0.301029996	0.43373423	1.204119983	0.602059991	0.43373423	1.204119983	0.602059991	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.271066317	1.397940009	0.77815125	0.256828796	1.380211242	0.77815125];
% stg.pa(2,:) = zeros(1,160);
% end

%  stg.pa(2,:) = [-0.615611714776611,-0.440851888542055,-1.78919414099660,1.92484642774317,-0.482136416143407,1.12931989531952,-0.0931336115651770,-1.07439135759926,-3.16921794580828,-1.43094316112731,0.578421665760323,-0.688514048750643,-0.976666508309151,2.27752510776872,-4.00681478444751,-0.546971851883420,-3.67220773204434,-0.354443577426494,-0.856740126446678,0.703237851602270,0.257369325193082,-1.74260819226670,-2.82646137257834,0.811846190603446,-0.570072196872103,1.71240386388230,-0.136870962661423,-0.489340484048084,1.05278411311692,0.494591799266926,-0.392868498164796,-1.38450146137070,1.56434201909131,0.625291171821348,0.810538056910869,-0.374539082602307,0.641085521233614,-4.47103854274751,-0.595743322074374,-0.0673893143122104,1.14059281971451,-0.719054710754626,-1.79870140787098,-3.73625270696668,0.263488498304777,-1.21056373210913,-0.441466275837582,-0.898796521982342,-1.34488758096527,-2.92590585860333,-0.718154609842875,-2.44407313706029,-0.459449177239094,0.156021402043514,2.29587281905002,0.978805863801138,1.53335979967517,0.204758933177235,-2.05124161966972,-0.523070770902631,0.448485046316039,0.920082587281629,0.968094622190901,1.05269028248203,0.467609107715365,-0.0344780529364931,0.0605466237066875,0.276127737514672,1.26811230755881,-0.355857355252218,1.32534847633656,0.872970930431388,-0.267152763783748,1.87308978946544,1.02875000498964,0.889670759973508,1.76564515274086,0.791993672678531,2.02516990523358,1.80665597874737,2.39199511760452,1.77780049302988,2.44805667446608,2.90537061315281,1.19223301102351,2.58531989820272,0.435872998688251,2.05704395957564,1.84900505262796,1.65711041784618,-0.328785627280028,1.80606888452392,1.29329903577513,-0.0925053808844482,1.10095026500251,1.76199125229054,2.14895893416018,1.16789813003138,0.0283835862909192,1.08725167247844,0.520769083690245,-1.20716851977582,1.07316754192843,-0.587809599664108,-0.391642652153959,0.0616169100154886,0.0624864160453047,0.120388588313395,1.68629499076718,0.717957583690764,-0.249926724893577,-2.79686820276623,-1.99222663922933,-1.60339251427485,-0.0209847199658022,-1.45385891428268,-2.11315165721224,1.73035401642550,2.41412599838649,1.73450472300871,0.103976220591857,-0.437546336531047,-0.377902132296508,-0.445220882503009,-0.781827185901649,-0.246106563558373,-0.247954773968512,0.825513591319837,-1.01049858880713,0.632357831897745,-0.0360403750218010,-1.23901948071283,0.514474031371236,2.49208003434809,1.97309826838352,1.54468062110941,0.930452725408674,1.12890154278119,0.0748626510063124,-0.0234516444742313,1.34733722600082,1.02580969535366,-0.526410202336678,1.26689201357337,-0.0873708566898416,-0.130147515879578,0.518527103576569,0.419189033783347,0.296869562040892,0.869074710982389,0.745922745953776,1.22101157093359,1.75654677372417,0.288120422384609,0.409878064705834,0.776603352711404,0.920817002062236,0.435188191009384,0.795342939094698,0.0482802078382720];
%  stg.pa(3,:) = [-0.137433954776403,-0.379855704698571,-2.00347968440452,0.851903580736005,0.0738757562155914,0.567274438466867,0.0420184036801745,-1.06589771126239,-2.82290305890639,-0.953897660113348,0.300340046707986,-0.761928865518452,0.0467844395284122,1.67691990522323,-3.82641899922688,-0.901966081972425,-3.14460944655693,-1.03592683879670,-1.83636834348679,0.326696378376652,0.0773859904897652,-1.08935927448051,-2.42857626835596,0.670605322047080,-0.432662637439656,0.766964071204939,-0.935521934840381,-0.671478785285526,0.549418950828975,-0.0214280598672633,0.364538296361401,-0.789125757519286,1.39745244029031,0.0438310788075769,1.06224750498794,-0.239452892190053,0.0548893260304967,-4.03826530263785,-0.985683844191041,0.0216668286159466,0.541913889413742,-0.502397263218142,-2.62218172349892,-3.56883553590389,-0.731205414029603,-0.369540952509520,-1.05939762763884,-1.51934163344028,-1.84590872106140,-3.08945923139655,-1.56677847153309,-1.79950815539864,-0.0994426786635791,0.0754605014713411,2.30485734405107,1.03246448912090,1.11703952058070,1.00211131644506,-1.33770366946127,-1.25078688217600,0.0374499435484306,1.21842164859605,0.624792689333885,0.703782738228541,1.14645758074707,0.511978058421103,0.702017097346163,1.15429058116269,0.588404175783262,0.0696607913445402,1.30699395633463,0.798488804784570,0.390388668764562,1.73882243609375,1.11540866538910,0.921259040868999,2.06079119044862,1.62342542003123,1.05223952166759,2.41259323442115,1.72397865077855,1.40103167990683,2.75854704640765,2.09289199633417,0.409974670284776,1.50498862671709,1.07248231207133,1.68741155594635,2.37469006703211,1.71854496103408,0.495566803893766,1.98678560202645,1.23041634559732,0.256482677416593,1.56279444474813,1.05474956681318,1.33808042073330,1.57839827268450,1.06584063634798,1.12968091468964,-0.180790640286052,-0.817138889602041,1.28151344737042,-0.262894276291246,-0.760253050153393,1.02125689716146,0.132859351242294,-0.528203595349150,1.00100065619790,0.0294420639216752,-0.495281533911711,-1.95651868155989,-1.01758410152233,-1.75995554526641,-0.687520155700325,-1.07013554013379,-1.77638329578294,1.59654872047749,1.61275388872239,1.05541316569339,0.237082512741124,-0.00407290410219036,-0.747065838880799,0.356101274717097,-0.0218639348724035,-0.599652781986926,0.0145980706337156,-0.170444532282978,-0.662665228375095,0.362477399727756,-0.152334997140538,-0.599566619347430,1.10975357746244,1.97211661141326,1.77319811788985,1.28457496391765,0.896614137527185,0.383243142844524,-0.327931828355556,0.521050187675847,1.15405052462589,0.668792607053215,0.391562085012437,1.15089666406769,0.520362408287462,0.266874259754844,1.38503937521907,0.797817382826592,0.353012767401907,1.30456912578396,0.873052701152993,0.237746171049801,1.44633148269558,0.792024767310219,0.268089080494583,1.43612408099596,0.814173823271639,0.260343752478030,1.30951587608804,0.698146750791006];
%  stg.pa(2,:) = [-0.145847147493724,-0.325757493755420,-2.03711329157073,0.903725434887016,0.0313547132028129,0.627364800846056,0.0420184036801745,-1.06574929473518,-2.98516849453917,-0.912295595537812,0.226911886085611,-0.641313699085434,0.0631328055134493,1.71212648202147,-3.95911447793407,-0.998295417544014,-3.15046444769907,-0.954200171388445,-1.75771573225610,0.236802483533286,0.0933053300653260,-0.997448613067898,-2.55953482731855,0.573807439223409,-0.356306032715192,0.752267330910851,-0.950449350887427,-0.677057787784168,0.387969103197220,-0.00609038836765310,0.392700989712696,-0.871628270048271,1.45069266284766,-0.00264888229880354,1.06258156039029,-0.392614245396734,0.0424728990854217,-3.93937209518299,-0.978364032517805,0.0711818700611538,0.530575588775806,-0.524200851810862,-2.61677239355039,-3.56695777988349,-0.633776625993322,-0.229408450267644,-1.06522260659939,-1.51934163344028,-1.85263073020106,-3.08906018968114,-1.67718515726085,-1.73794827309108,-0.0994426786635791,0.0891924490744112,2.16330325647132,1.09611330509458,1.10287560396843,0.915918476520377,-1.37263215663024,-1.21416049697017,0.0542303621949206,1.10822262984398,0.614491234329984,0.654047818908599,1.18193510661690,0.662417742199739,0.661869811049604,1.28910453064209,0.520784844391489,0.172997317686211,1.35425875232030,0.798488804784570,0.332003521547971,1.54805520994682,0.990915566457522,1.00742375721438,2.11187128455326,1.63545663334787,1.17523511624797,2.31091618796671,1.79762391915873,1.40103167990683,2.76668806867903,2.01118616709502,0.384179558821359,1.62592999953418,0.948587832157236,1.74491710168379,2.25982905224032,1.64303821938483,0.526539384263914,1.91619981956490,1.28281463557738,0.180449465912545,1.70113911594233,0.909720489873201,1.33545967197855,1.51891389994473,1.03544109053537,1.09008504155366,-0.138570829345955,-0.746334869869912,1.27889247318777,-0.260924395149570,-0.727370751805314,0.953699574738096,0.132859351242294,-0.468773403718022,1.04268523212638,-0.0171661272797053,-0.441144231198177,-1.95651868155989,-1.06842861464701,-1.71823732739696,-0.813265281338447,-0.999501340215661,-1.79070999763157,1.65515918640998,1.55161506760884,1.02472270902239,0.221447524354430,-0.189778920154925,-0.760347024938212,0.375070341285341,-0.0545357218739391,-0.627485121225427,-0.0248786764385085,-0.109311087713826,-0.662665228375095,0.366655228813116,-0.113597486372111,-0.599566619347430,1.10975357746244,2.00004651800760,1.85199326363170,1.29373255244059,0.935981326997760,0.353253410687620,-0.304706990504391,0.519124396212809,1.11942420429611,0.689223555300896,0.361069173990820,1.23721787923943,0.596274245003963,0.246076979244390,1.44067100181017,0.691284462970664,0.361096709793939,1.33341695309512,0.859627262663874,0.349292211152455,1.47053830357278,0.792024767310219,0.341640593755530,1.33558780396452,0.861795373984810,0.242816260822079,1.48004211496444,0.830987671048061];
%  stg.pa(3,:) = [0.723120447529783,-1.29632999093680,-2.80594542345659,1.76874574893514,-0.641318933602066,1.28445697006349,-0.545858761290583,-0.820099490045370,-3.55153740768877,-0.290640669674559,0.911720379778285,0.235943442894696,0.0537149535357528,1.04108409187857,-4.05624913832172,-0.535762514295776,-2.29414744142958,-0.853400616758330,-1.08836741750606,-0.244132026704977,0.900501180892137,-1.45042270079751,-2.95174121739637,0.811846190603446,-0.232563289500472,0.727194232031983,-1.82012081816488,-1.46673065993841,0.177382901688732,-0.637831606243431,-0.196552396531337,-0.594962347046488,1.02586183595856,0.625291171821348,1.07383405327978,0.145208254714125,0.835515843971186,-4.88455669425219,-1.47329415686279,-0.0596246894674168,1.46775821908367,-0.226862497844894,-1.94212621887284,-4.43843916271625,0.0668364142306250,0.0635829712877463,-1.34359528417939,-1.03504689149402,-2.32328182308288,-3.99746383086376,-0.718154609842875,-1.85004819382055,-0.940737534136699,0.322299404512153,2.51249564835970,1.54642358935293,0.748748497367774,0.202308007199383,-0.813214721508565,-0.505697847791020,0.136103018776271,0.392359485857693,-0.261847137987496,1.21829231728727,0.946277369160434,-0.349514122951306,0.941014703371235,1.74240350911355,0.0749390941055482,0.0474157887497289,1.72966119805630,0.798790097081877,-0.391806963160578,1.47422698297400,1.77730953174528,0.472460022494735,1.79626104163678,1.76154822595706,2.11481944867603,1.60386749699329,2.09671927196415,2.18153795212379,1.98124210269467,2.59950638228367,0.529769760751684,0.997569937241572,0.170888190849422,1.27304510591278,2.72461955073118,1.89571679400593,1.08318249328790,1.83498432771944,1.49215466732560,0.717820921919219,1.35786718962320,0.498973009219860,2.22189969207139,0.967115278511412,1.56366349836234,1.15879235044878,0.739412872470185,-0.379252516326899,1.39710499554684,-0.500261861797075,-0.164403031826502,1.36571772895928,0.384816080897481,0.443368309179214,0.892826069174720,0.523713308295309,-0.258994913593027,-2.37591209167806,-1.96726565132662,-2.24164356659429,0.227246621781905,-1.25038732827267,-2.43022904105667,2.04432644072557,1.53660190590971,1.85362210720020,1.02236352472954,-0.185866092781083,-1.18820054477819,0.0731868061975125,-0.951873713186329,-0.456253345179585,-0.908198546758043,0.625018659024384,0.0654114014988600,0.407817677841490,-0.0730548452829631,-1.04996596124731,1.38930540277745,2.88570372150518,2.15851795049664,0.742379428929238,0.955332189111798,0.849042099166973,0.00399782804121895,1.41934376899105,1.08620330443138,-0.298361312740128,-0.291141439748803,2.07946814849392,0.573208188842549,1.03635042634365,1.64678466499404,-0.158903300275916,0.243311256913606,1.86951104973674,0.848671353445475,-0.187292452045927,0.745040745252404,0.577712444033781,0.605688390364702,1.18262625782613,1.67577877945018,0.505496067761115,1.23083638199876,0.281113376551916];

%Correct parameters

stg.pa(1,:) = [0.999968151	-0.301029996	0	-4	-1	0.264435738	1.602059991	1	1.183916142	-0.22184875	-0.823908741	1.183916142	-0.22184875	-0.823908741	0.948847023	0.079181246	-0.522878745	0.948847023	0.079181246	-0.522878745	-2.023562183	-1.096910013	-1.698970004	-0.728933683	-1.096910013	-1.698970004	1.66888375	1.602059991	1];
% stg.pa(2,:) = [1.20327577993523,0.313004819657628,-0.847960779005249,-3.72858066814786,-0.550448716714675,-0.611592770793515,2.31886573537624,0.481456747980214,0.365768657633490,-0.683415872577579,0.153923390182066,1.06063326577892,-1.06384580589771,0.159271375843753,1.44599697691288,0.864651828387445,-0.307382131999887,1.62276312499277,0.242936206435653,0.359885193960285,-2.26116772023767,-1.89600181858380,-1.46727820457674,-1.54897665986418,-0.681400435613773,-2.53706951179994,1.02910373944186,1.78194651365667,1.73016958751427];
% stg.pa(3,:) = [1.50588410182241,-0.0935763452208480,-0.580527199032361,-3.99491581003785,-0.000656107723683963,-0.243271025299648,2.06203111071450,1.69163898074757,0.485084894460439,0.715559192270335,-0.0620581188561402,1.85109334045899,-0.489539341319553,-0.821244626833460,1.49970672966907,-0.414279111634506,0.399438708993961,1.12891576951045,0.331066029914223,-0.407987391885445,-1.38800206698805,-0.790371374222949,-2.19644292044499,-1.10291642667294,-1.13858162786427,-2.41106669689700,1.27846472440081,0.835896249725143,1.31868321088228];
% stg.pa(4,:) = [1.26293814585163	0.390552340008227	-0.811489291937400	-4.78506405242233	-0.512411529483115	-0.554133604035905	1.52467928680644	0.185130623371452	1.11539718496376	0.544223394598075	0.133004635163269	1.52352755229485	0.582122720498769	0.0966006028844693	1.04473661531405	0.701202615391585	0.460199613671783	1.69379190578062	-0.309882796765651	-0.346013231200285	-1.86043932353880	-2.04170739195718	-2.39222612053496	-1.02135885675066	-1.03220270086897	-2.43751951704707	1.63022071900023	1.07946346038527	0.466978377963250];
% stg.pa(5,:) = [1.26293814585163	0.390552340008227	-0.811489291937400	-4.26685138406449	-0.762411529483115	-0.554133604035905	1.52467928680644	0.185130623371452	1.11539718496376	0.544223394598075	0.133004635163269	2.14736234432327	0.582122720498769	-0.111761175702841	1.16973661531405	0.701202615391585	0.460199613671783	1.69379190578062	-0.309882796765651	-0.142942457885960	-1.29614717422326	-1.53521800986492	-2.39222612053496	-1.02135885675066	-0.521161784370842	-2.52385880566385	1.63022071900023	1.07946346038527	0.883289303349880];
stg.pa(2,:) = [0.889708873794741,0.496937639527883,-0.776274726085585,-4.94981951161300,-0.873605121366917,-0.692910689994836,2.58533301066836,0.0984976998997857,0.812521118636573,0.483855438647211,0.140891437672981,1.13520734032854,-0.403088207531380,0.125193578841168,1.19861893169260,-0.179677696622224,0.429927631703068,1.84125952132435,0.328283981123265,0.449391633443089,-2.08207472518182,-0.706830784093234,-0.764803757430592,-1.61740747552252,-1.11347331849189,-2.53346724627044,1.13066398461756,0.841491835677703,0.897905803786056];
stg.pa(3,:) = [1.44824918743424,0.327411032873311,-1,-3.98956611640901,-1.29445953277027,-0.502343599207032,2.09991689730418,0.223431039212322,0.484185038322967,-0.177439325752746,0.173834992315229,1.06840226722783,-0.552427286595109,0.166073198809973,1.38526346100751,0.468727451229625,0.471179799226031,1.72763797908063,0.184780807823937,0.321244171880714,-2.08979766234699,-2.08618589462227,-1.70180270786196,-1.19101127140512,-0.212820936243880,-2.61088682429332,2.11755042715234,1.40118218335712,1.30562783165817];


 % Best parameter array found so far for the model
% (Best parameter array)
stg.bestpa = stg.pa(1,:);

%% Plots

% True or false to decide whether to plot results
% (Plots)
stg.plot = true;

% True or false to decide whether to use long names in the title of the outputs
% plots in f_plot_outputs.m
% (Plot outputs long names)
stg.plotoln = true;

%% Sensitivity analysis

% Number of samples to use in SA
% (Sensitivity analysis number of samples)
stg.sansamples = 21600;

% True or false to decide whether to subtract the mean before calculating SI and
% SIT
% (Sensitivity analysis subtract mean)
stg.sasubmean = true;

% Choose the way you want to obtain the samples of the parameters for 
% performing the SA;
% 0 Log uniform distribution truncated at the parameter bounds
% 1 Log normal distribution with mu as the best value for a parameter and
% sigma as stg.sasamplesigma truncated at the parameter bounds
% 2 same as 1 without truncation
% 3 Log normal distribution centered at the mean of the parameter bounds and
% sigma as stg.sasamplesigma truncated at the parameter bounds
% 4 same as 3 without truncation.
% (Sensitivity analysis sampling mode)
stg.sasamplemode = 2;

% Sigma for creating the normal distribution of parameters to perform
% sensitivity analysis
% (Sensitivity analysis sampling sigma)
stg.sasamplesigma = 0.1;

%% Profile Likelihood

% Parameter(optimization array) that is being worked on in a specific
% iteration of PL (if -1 no parameter is being worked in PL)
% (Profile Likelihood Index)
stg.PLind = -1;

% Which parameters to do PL on, it should be all parameters but can also be
% a subset for testing purposes
% (Profile Likelihood parameters to Test)
stg.pltest = (1:15);

% True or false to decide whether to do plots after calculating PL
% (Profile Likelihood Plots)
stg.plplot = true;

% True or false to decide whether to run simulated annealing
% (Profile Likelihood Simulated Annealing)
stg.plsa = true;

% 0 or 1 to decide whether to run fmincon
% (Profile Likelihood FMincon)
stg.plfm = false;

%% Optimization

%  Time for the optimization in seconds (fmincon does not respect this
% time!!)
% (Optimization time)
stg.optt = 60*30;

% Population size for the algorithms that use populations
% (Population size)
stg.popsize = 3600;

% True or false to decide whether to display Plots (Plots doesn't work if using
% multicore)
% (Optimization plots)
stg.optplots = true;

% True or false to decide whether to run fmincon (no gradient so this doesn't work
% very well, no max time!!)
stg.fmincon = false;

% True or false to decide whether to run simulated annealing
% (Simulated annealing)
stg.sa = false;

% True or false to decide whether to run Pattern search
% (Pattern search)
stg.psearch = false;

% True or false to decide whether to run Genetic algorithm
% (Genetic algorithm)
stg.ga = true;

% True or false to decide whether to run Particle swarm
% (Particle swarm)
stg.pswarm = true;

% True or false to decide whether to run Surrogate optimization
% (Surrogate optimization)
stg.sopt = false;
end