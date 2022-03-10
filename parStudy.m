classdef parStudy < handle
    %PARSTUDY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Properties relevant to parameter study
        parStudyRct           % structure: contains the stationary_oneDimRct objects containing simulation results
        numVar                % double: number of varied parameters
        parVars               % cell array: cell i contains varied parameter (string), cell i+1 contains range over which parameter is varied (line vector)
        caseMode              % double: decides how cases for parStudy are created
        parVarArray           % cell array: contains all possible combinations of parameter values from their ranges
        parVarArrayLine       % cell array: parVarArray transformed to line vector of cells
        diaryPath             % string: path to diary txt-file
        caseStat              % logical array: contains index which case resulted in an error
        casesSuccess          % row vector: contains numbers of successful cases
        numSuccess            % double: number of successful simulations
        tableParFix           % cell array: contains name value pairs of parameters being constant for the output tables
        tableParVar           % cell array: contains the parameters which are being written into the table for each case 
        tablePS               % cell array: contains tables with variables chosen by user
        
        sDateTime             % string: current date and time in this format: 'yyyy-mm-dd_HH-MM-SS'
        sDate                 % string: current date in this format: 'yyyy-mm-dd'
        computationTime       % double: computation time for the entire program in s
        
    end
    
    methods
        %% Constructor method of the class
        function obj = parStudy(varargin)
            %PARSTUDY ...
            
            % Create input parser
            parser = inputParser;

            addParameter (parser,'numVar',[]);
            addParameter (parser,'parVars',[]);
            addParameter (parser,'timeStamp',[]);
            addParameter (parser,'dateStamp',[]);
            addParameter (parser,'caseMode',[]);

            % Parse for input
            parse (parser,varargin{:});

            % Results
            obj.numVar = parser.Results.numVar;
            obj.parVars = parser.Results.parVars;
            obj.sDateTime = parser.Results.timeStamp;
            obj.sDate = parser.Results.dateStamp;
            obj.caseMode = parser.Results.caseMode;
            
        end
        
        %% cases
        function cases(obj)
            % PARVARARRAY returns cell array, that contains all possible
            % combinations of the variables values in parVars
            
            % iX contains the number of values of each variable
            iX = zeros(1,obj.numVar);
            % for every other cell in parVars containing values
            for m = 2:2:length(obj.parVars)
                iX(m/2) = length(obj.parVars{m});
                if iX(m/2) == 0
                    iX(m/2) = 1;
                end
            end
            
            % if parXRan has less elements than largest parXRan, fill with NaN
            for m = 2:2:length(obj.parVars)
                if length(obj.parVars{m}) < max(iX)
                    obj.parVars{m}(end+1:max(iX)) = NaN(1,max(iX)-length(obj.parVars{m}));
                end
            end
            
            if obj.caseMode == 1
                % Create n dimensional array for all possible combinations of values of
                % all the variables
                % Go through parVarArray by linear indexing. Indexes are converted
                % to subscripts, eg. 1 -> 1,1; 2 -> 2,1; 3 -> 1,2; 4 -> 2,2
                % Then [par1Ran(1) par2Ran(1)]; [par1Ran(2) par2Ran(1)]... is being
                % written into the parVarArray

                % ind maps output cell array by lin. index (see sub2ind documentation)
                ind = 1;
                % preallocate output
                if length(iX) == 1
                    obj.parVarArray = cell(iX,1);
                else
                    obj.parVarArray = cell(iX);
                end
                % preset variable that keeps while loop running until its set to 1
                start = 1;
                while start == 1
                    % sub maps output cell array by subscript
                    sub = cell(1,obj.numVar);
                    [sub{:}] = ind2sub(size(obj.parVarArray),ind);
                    sub = [sub{:}];
                    % preallocate temp
                    temp = zeros(1,obj.numVar);
                    % assemble combination of variables by using ind
                    for n = 2:2:length(obj.parVars)
                        temp(1,n/2) = obj.parVars{n}(sub(n/2));
                    end
                    % write temporary variable to output variable
                    obj.parVarArray(ind) = {temp};
                    clear temp
                    % increase index by one
                    ind = ind+1;
                    % if ind is larger than the maximum amount of values in the array,
                    % stop while loop by setting start to zero
                    if ind > prod(iX)
                        start = 0;
                    end
                end
            elseif obj.caseMode == 2
                % check if all variables have the same number of values
                numVal = zeros(obj.numVar,1);
                for i = 2:2:length(obj.parVars)
                    % number of values for each variable
                    numVal(i/2) = length(obj.parVars{i});
                end
                if ~all(numVal==numVal(1))
                    error('Variables do not have the same number of values!');
                end
                % Create 1 dimensional array with the same length as the
                % parXRan vectors
                % number of values i (the same for all variables)
                for i = 1:numVal(1)
                    % in parVars every other cell contains values
                    for j = 2:2:length(obj.parVars)
                        % write values into parVarArray which are then used
                        % to create the cases in setupPS
                        obj.parVarArray{i}(j/2) = obj.parVars{j}(i);
                    end
                end
            end
        end
        
        %% setupPS
        function setupPS(obj,varargin)
            %SETUPOBJECTS ...
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'speciesNames', {'H2';'CO';'N2';'CO2';'H2O'; ...
            'CH4';   'C2H6';  'C3H8';  'C4H10'; 'C5H12'; 'C6H14'; 'C7H16'; 'C8H18'; 'C9H20';  'C10H22'; ...
            'C11H24';'C12H26';'C13H28';'C14H30';'C15H32';'C16H34';'C17H36';'C18H38';'C19H40'; 'C20H42'; ...
            'C21H44';'C22H46';'C23H48';'C24H50';'C25H52';'C26H54';'C27H56';'C28H58';'C29H60'; 'C30H62'; ...
            'C31H64';'C32H66';'C33H68';'C34H70';'C35H72';'C36H74';'C37H76';'C38H78';'C39H80'; 'C40H82'; ...
            'C41H84';'C42H86';'C43H88';'C44H90';'C45H92';'C46H94';'C47H96';'C48H98';'C49H100';'C50H102'});
            addParameter (parser,'numSpeciesInlet', 55);
            addParameter (parser,'numSpeciesGasPhase', 25);
            addParameter (parser,'ReportStatistics', 'on');
            addParameter (parser,'rctEnthalpy', 160000);
            addParameter (parser,'wallTemperature', 273.15);
            addParameter (parser,'meshFactor',1);
%             addParameter (parser,'meshHeight',0.009);
%             addParameter (parser,'meshLength',0.03144);
            addParameter (parser,'elementSize',0.01);
            addParameter (parser,'resTol',1e-5);
            addParameter (parser,'WHSV',6.7);
            addParameter (parser,'ratioInlet',2.1);
            addParameter (parser,'TInlet',523.15);
            addParameter (parser,'pInlet',20);
            addParameter (parser,'kRct',600);
            addParameter (parser,'kCool',1e5);
            addParameter (parser,'alphaRct',600);
            addParameter (parser,'alphaCool',1e5);
%             addParameter (parser,'areaCrossRct',63.617e-6);
            addParameter (parser,'lengthRct',0.03144);
            addParameter (parser,'wRct', 8);
            addParameter (parser,'dRct', 1e-3);
            addParameter (parser,'volRedRct', 8);
            addParameter (parser,'rhoCat', 8);
%             addParameter (parser,'timeStamp',[]);
%             addParameter (parser,'dateStamp',[]);
            addParameter (parser,'TConst',523.15);
            addParameter (parser,'TCool',373.15);
            addParameter (parser,'thicknessWall',0.001);
            addParameter (parser,'rateFactor',1);
            addParameter (parser,'aRctFactor',1);
            addParameter (parser,'aWallFactor',1);
            addParameter (parser,'DFactor',1);
            addParameter (parser,'constMIXTURESwitch', 1);
            addParameter (parser,'TRctSwitch', 1);
            addParameter (parser,'kModeSwitch', 3);
            addParameter (parser,'TWallSwitch', 1);
            
            % Parse for input
            parse (parser,varargin{:});
            
            speciesNames = parser.Results.speciesNames;
            numSpeciesInlet = parser.Results.numSpeciesInlet;
            numSpeciesGasPhase = parser.Results.numSpeciesGasPhase;
            reportStats = parser.Results.ReportStatistics;
            DeltaRctH = parser.Results.rctEnthalpy';
            TWall = parser.Results.wallTemperature';
            meshFactor = parser.Results.meshFactor;
%             meshHeight = parser.Results.meshHeight;
%             meshLength = parser.Results.meshLength;
            elementSize = parser.Results.elementSize;
            resTol = parser.Results.resTol;
            WHSV = parser.Results.WHSV;
            ratioInlet = parser.Results.ratioInlet;
            TInlet = parser.Results.TInlet;
            TConst = parser.Results.TConst;
            TCool = parser.Results.TCool;
            pInlet = parser.Results.pInlet;
            kRct = parser.Results.kRct;
            kCool = parser.Results.kCool;
            alphaCool = parser.Results.alphaCool;
            alphaRct = parser.Results.alphaRct;
%             areaCrossRct = parser.Results.areaCrossRct;
            lengthRct = parser.Results.lengthRct;
            wRct = parser.Results.wRct;
            dRct = parser.Results.dRct;
            volRedRct = parser.Results.volRedRct;
            rhoCat = parser.Results.rhoCat;
%             sDateTime = parser.Results.timeStamp;
%             sDate = parser.Results.dateStamp;
            thicknessWall = parser.Results.thicknessWall;
            rateFactor = parser.Results.rateFactor;
            aRctFactor = parser.Results.aRctFactor;
            aWallFactor = parser.Results.aWallFactor;
            DFactor = parser.Results.rateFactor;
            constMIXTURESwitch = parser.Results.constMIXTURESwitch;
            TRctSwitch = parser.Results.TRctSwitch;
            kModeSwitch = parser.Results.kModeSwitch;
            TWallSwitch = parser.Results.TWallSwitch;
            
            % reshape parVarArray into linevector
            obj.parVarArrayLine = reshape(obj.parVarArray,1,[]);
            
            i1 = length(obj.parVarArrayLine);
            for i = 1:i1
                    % current date and time in proper format as string
                    sDate2 = datestr(datetime('now'),'yyyy-mm-dd');
                    sDateTime2 = datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS');
                    
                    % overwrite the array of values of the varied parameter
                    % with the single value of the i-th case
                    for j = 2:2:length(obj.parVars)
                        eval([obj.parVars{j-1},'=obj.parVarArrayLine{i}(j/2);']);
                    end

                    % Create instances of stationary_oneDimRct
                    obj.parStudyRct{i} = stationary_oneDimRct('speciesNames', speciesNames, ...
                    'numSpeciesInlet', numSpeciesInlet,...
                    'numSpeciesGasPhase', numSpeciesGasPhase,...
                    'meshFactor', meshFactor,...
                    'elementSize', elementSize,...
                    'resTol', resTol,...
                    'rctEnthalpy', DeltaRctH,...
                    'rhoCat', rhoCat,...
                    'WHSV', WHSV,...
                    'ratioInlet', ratioInlet,...
                    'inletTemperature', TInlet, ...
                    'TConst', TConst,...
                    'TCool', TCool,...
                    'inletPressure', pInlet,...
                    'wallTemperature', TWall,...
                    'kRct', kRct,...
                    'kCool', kCool,...
                    'alphaRct', alphaRct,...
                    'alphaCool', alphaCool,...
                    'rateFactor', rateFactor,...
                    'aRctFactor', aRctFactor,...
                    'aWallFactor', aWallFactor,...
                    'DFactor', DFactor,...
                    'ReportStatistics', reportStats,...
                    'lengthRct', lengthRct, ...
                    'wRct',wRct,...
                    'dRct',dRct,...
                    'volRedRct',volRedRct,...
                    'timeStamp', sDateTime2,...
                    'dateStamp', sDate2,...
                    'constMIXTURESwitch', constMIXTURESwitch,...
                    'TRctSwitch', TRctSwitch, ...
                    'kModeSwitch', kModeSwitch,...
                    'TWallSwitch', TWallSwitch, ...
                    'thicknessWall', thicknessWall );
%                     'areaCrossRct', areaCrossRct, ...
%                     'meshHeight', meshHeight,...
%                     'meshLength', meshLength,...
                        
                    % Initialize PDE Models
                    obj.parStudyRct{i}.initModel
                    obj.parStudyRct{i}.domainMesh
                    obj.parStudyRct{i}.setCoeff
                    obj.parStudyRct{i}.setBC
                    obj.parStudyRct{i}.setIC
            end
            
        end
        
        %% solvePS
        function obj = solvePS(obj)
            %SOLVEPS starts simulation of each case in parStudyRct. Results
            %need to be saved temporarily in variable due to parfor
            %specifications
            
            tic
            % create parallel processing pool, only if no pool exists, with
            % number workers = number of physical cores
            if isempty(gcp('nocreate'))
                parpool(maxNumCompThreads);
            end
            
            % preallocate future objects
            futureTemp(1:length(obj.parVarArrayLine)) = parallel.FevalFuture;
            % execute solve method via parfeval in parallel for all cases
            for i = 1:length(obj.parVarArrayLine)
                futureTemp(i) = parfeval(@solve, 1, obj.parStudyRct{i});
            end
            
            % wait for simulations to be finished
            wait(futureTemp);
            
            % preallocate resultTemp
            resultTemp(1,1:length(obj.parVarArrayLine)) = pde.StationaryResults;
            for i = 1:length(obj.parVarArrayLine)
                if ~isempty(futureTemp(i).Error)
                    disp('errör in täsk')
                else
                    % get output from future array
                    resultTemp(i) = fetchOutputs(futureTemp(i))';
                end
            end
            
            % write diary entries into corresponding case object
            for i = 1:length(obj.parVarArrayLine)
                % use splitline to create new cell when line break occurs
                obj.parStudyRct{i}.diaryOutput = splitlines(futureTemp(i).Diary);
            end
            
            % delete parallel pool
%             delete(gcp('nocreate'));

            % stop method if no simulation produced any results
            if isempty(resultTemp(:))
                error('No simulation produced any results.')
            end
            
            % Assign temp variable to resultPDE property
            for i = 1:length(obj.parVarArrayLine)
                obj.parStudyRct{i}.resultPDE =resultTemp(i);
                % if NodalSolution is empty, empty resultPDE property
                if isempty(obj.parStudyRct{i}.resultPDE.NodalSolution)
                    obj.parStudyRct{i}.resultPDE = [];
                end
            end
            
            toc
        end
        
        %% evaluatePS
        function evaluatePS(obj)
            %SOLVEPS starts simulation of each case in parStudyRct. Results
            %need to be saved temporarily in variable due to parfor
            %specifications
            
            for i = 1:length(obj.parVarArrayLine)
                % save first residual as property
                % first residual is saved in diaryOutput(2), as residual of
                % Iteration step 0
                % look with regex for symbols matching the residual
                % 1 or more digits - dot - 1 or more digits - e - plus or
                % minus - 1 or more digits
                % convert cell result to double
                if size(obj.parStudyRct{i}.diaryOutput,1) > 1
                    obj.parStudyRct{i}.firstResid = str2double(cell2mat(regexp(obj.parStudyRct{i}.diaryOutput{2},'\d+\.\d+[e][\+\-]\d+','match')));
                end
                % check if simulation was successful
                % Case 1: simulation produces an error, therefore resultPDE is
                % empty
                if isempty(obj.parStudyRct{i}.resultPDE)
                    obj.parStudyRct{i}.caseStatRct = 0;
                    obj.parStudyRct{i}.caseStatRctMsg = 'Simulation error, no results in resultPDE, probably step size error';
                elseif obj.parStudyRct{i}.firstResid < obj.parStudyRct{i}.resTol
                    obj.parStudyRct{i}.caseStatRct = 0;
                    obj.parStudyRct{i}.caseStatRctMsg = 'Simulation error, first residual < ResidualTolerance';
                elseif ~isempty(obj.parStudyRct{i}.resultPDE)
                    obj.parStudyRct{i}.evaluate
                    obj.parStudyRct{i}.caseStatRct = 1;
                    obj.parStudyRct{i}.caseStatRctMsg = [];
                end
            end
            % temporary comma seperated list of reactor objects
            tempCSL = [obj.parStudyRct{:}];
            
            obj.caseStat(1,find(~[tempCSL.caseStatRct])) = 0;
            obj.casesSuccess = find([tempCSL.caseStatRct]);
            obj.numSuccess = length(find([tempCSL.caseStatRct]));
        end    
        
        %% plotInterpPS
        function plotInterpPS(obj, varargin)
            % PLOTINTERP plots the interpolated solution values

            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', false);
            addParameter (parser,'norm_n',false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            norm_nSwitch = parser.Results.norm_n;
            
            % Plots
            fig1 = figure('name', 'Independent variables along reactor axis');
            figure(fig1);
            
            
            %%%%%%%%%%%%%%%%  
            % Mole flow of CO
            %%%%%%%%%%%%%%%%
            subplot(3,1,1)
            % create color for each successful simulation
            % Modify colormap to cut off bright yellow parts
            cmPre = parula(obj.numSuccess*1.1);
            colorMap = cmPre(1:obj.numSuccess,:);
            % preallocate, not useful for nOutCO and xq_i if lengths differ
            % between cases
            legendArray = cell(1,obj.numSuccess);          
            if norm_nSwitch
                for i = 1:obj.numSuccess
                    nOutCO(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.nOut(2,:)./obj.parStudyRct{obj.casesSuccess(i)}.nOut(2,1);
                    xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
                    hold on
                    plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,nOutCO(i,:),'Color',colorMap(i,:),'Linewidth',2);
                    legendArray{i} = (num2str(obj.parVarArrayLine{obj.casesSuccess(i)}));
%                     ylabel('$\frac{\dot{\mathsf{n}}_\mathsf{CO}}{\dot{\mathsf{n}}_\mathsf{CO,in}} \mathsf{ ~/~ -}$','fontsize',25,'Interpreter','latex');
                    ylabel('\fontsize{16}{16} $\frac{\dot{\mathsf{n}}_\mathsf{CO}}{\dot{\mathsf{n}}_\mathsf{CO,in}}$ \fontsize{10}{10}$\mathsf{ ~/~ -}$','Interpreter','latex');
                end
            else
                for i = 1:obj.numSuccess
                    nOutCO(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.nOut(2,:);
                    xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
                    hold on
                    plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,nOutCO(i,:),'Color',colorMap(i,:),'Linewidth',2);
                    legendArray{i} = (num2str(obj.parVarArrayLine{obj.casesSuccess(i)}));
                    ylabel('$\dot{\mathsf{n}}_\mathsf{CO} \mathsf{ ~/~ mol~s}^{\mathsf{-1}}$','Interpreter','latex');
                end
            end
            hold off
            box on
            title('(a)    Molar flows of CO')
            ylim([0 inf])
            xlim([0 1])
            xlabel('Dimensionless reactor length / -')
            legend(legendArray,'Location','northeastoutside')
            set(gca,'linewidth',1.5)

%             %%%%%%%%%%%%%%%%
%             % Mole flow of H2
%             %%%%%%%%%%%%%%%%
%             subplot(3,1,1)
%             % create color for each successful simulation
%             % Modify colormap to cut off bright yellow parts
%             cmPre = parula(obj.numSuccess*1.1);
%             colorMap = cmPre(1:obj.numSuccess,:);
%             % preallocate, not useful for nOutCO and xq_i if lengths differ
%             % between cases
%             legendArray = cell(1,obj.numSuccess);          
%             if norm_nSwitch
%                 for i = 1:obj.numSuccess
%                     nOutH2(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.nOut(1,:)./obj.parStudyRct{obj.casesSuccess(i)}.nOut(1,1);
%                     xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
%                     hold on
%                     plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,nOutH2(i,:),'Color',colorMap(i,:),'Linewidth',2);
%                     legendArray{i} = (num2str(obj.parVarArrayLine{obj.casesSuccess(i)}));
% %                     ylabel('$\frac{\dot{\mathsf{n}}_\mathsf{CO}}{\dot{\mathsf{n}}_\mathsf{CO,in}} \mathsf{ ~/~ -}$','fontsize',25,'Interpreter','latex');
%                     ylabel('\fontsize{16}{16} $\frac{\dot{\mathsf{n}}_\mathsf{H2}}{\dot{\mathsf{n}}_\mathsf{H2,in}}$ \fontsize{10}{10}$\mathsf{ ~/~ -}$','Interpreter','latex');
%                 end
%             else
%                 for i = 1:obj.numSuccess
%                     nOutH2(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.nOut(1,:);
%                     xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
%                     hold on
%                     plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,nOutH2(i,:),'Color',colorMap(i,:),'Linewidth',2);
%                     legendArray{i} = (num2str(obj.parVarArrayLine{obj.casesSuccess(i)}));
%                     ylabel('$\dot{\mathsf{n}}_\mathsf{H2} \mathsf{ ~/~ mol~s}^{\mathsf{-1}}$','Interpreter','latex');
%                 end
%             end
%             hold off
%             box on
%             title('(a)    Molar flows of H2')
%             ylim([0 inf])
%             xlim([0 1])
%             xlabel('Dimensionless reactor length / -')
%             legend(legendArray,'Location','northeastoutside')
%             set(gca,'linewidth',1.5)


            % TRct
            if obj.parStudyRct{1}.TRctSwitch
                subplot(3,1,2)
                for i = 1:obj.numSuccess
                    TRctOutPS(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.TRctOutMean;
                    xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
                    hold on
                    plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,TRctOutPS(i,:),'Color',colorMap(i,:),'Linewidth',2)
                end
                hold off
                box on
                xlim([0 1])
                title('(b)    Reactor temperatures')
                xlabel('Dimensionless reactor length / -')
                ylabel('T_{Rct} / K')
                % placeholder legend
                l2 = legend(legendArray,'Location','northeastoutside');
                % make it go away
                set(l2,'visible','off')
                set(gca,'linewidth',1.5)
            end
            
            % TWall
            if obj.parStudyRct{1}.TWallSwitch
                subplot(3,1,3)
                for i = 1:obj.numSuccess
                    TWallOutPS(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.TWallOutMean;
                    xq_i(i,:) = obj.parStudyRct{obj.casesSuccess(i)}.xq;
                    hold on
                    plot(xq_i(i,:)./obj.parStudyRct{obj.casesSuccess(i)}.lengthRct,TWallOutPS(i,:),'Color',colorMap(i,:),'Linewidth',2)
                end
                hold off
                box on
                xlim([0 1])
                title('(c)    Wall temperatures')
                xlabel('Dimensionless reactor length / -')
                ylabel('T_{Wall} / K')
                % placeholder legend
                l3 = legend(legendArray,'Location','northeastoutside');
                % make it go away
                set(l3,'visible','off')
                set(gca,'linewidth',1.5)
            end
            
            % put legend in front of all subplots
            h = get(gcf,'children');
            ind = find(isgraphics(h,'Legend'));
            set(gcf,'children',h([ind:end,1:ind-1]))
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_PS_Results.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end
            
        end
        
        %% Result presentation method: Hydrocarbon distribution
        function showResultFTPS(obj, varargin)
            %SHOWRESULTFT will display hydrocarbon distribution
            %   Mass fraction over carbon number
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            
            % Plot product mass fractions and diesel fraction
            fig2 = figure('name', 'FT');
            figure(fig2);
            cmPre = parula(obj.numSuccess*1.1);
            colorMap = cmPre(1:obj.numSuccess,:);
            hold on
            for i = 1:obj.numSuccess
                plot(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(:,end),'color',colorMap(i,:),'LineWidth',1.5);
                text(20, 0.035 + 0.035*(i/10), cat(2, 'w_{C10-C20} = ', num2str(obj.parStudyRct{obj.casesSuccess(i)}.w_Diesel(end))),'FontSize',16,'color',colorMap(i,:))
            end
            hold off
            box on
            xlabel('n_C / -','FontSize',20);
            ylabel('w_i / -');
            axis([0 30 0 0.1]);
            ax = gca;
            ax.FontSize = 16;
            title('Mass fraction of diesel in all hydrocarbons');
            set(gca,'linewidth',1.5)
            
            % patches for bar diagram
            for i = 1:obj.numSuccess
                if obj.parStudyRct{obj.casesSuccess(i)}.numSpeciesInlet < 10
                    disp('No mass balance for diesel fraction calculated.')
                elseif obj.parStudyRct{obj.casesSuccess(i)}.numSpeciesInlet < 26
                    xverts = [(10:10+length(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:end,end))); (11:10+(length(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:end,end))+1)); (11:10+(length(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:end,end))+1)); (10:10+length(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:end,end)))];
                    yverts = [zeros(1,1+length(obj.parStudyRct{i}.w_Hydrocarbons(11:end,end))); zeros(1,1+length(obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:end,end))); obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(10:end,end)'; obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(10:end,end)'];
                    patch(xverts,yverts,colorMap(i,:),'LineWidth',0.5);
                else
                    xverts = [(10:20); (11:21); (11:21); (10:20)];
                    yverts = [zeros(1,11); zeros(1,11); obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(11:21,end)'; obj.parStudyRct{obj.casesSuccess(i)}.w_Hydrocarbons(10:20,end)'];
                    patch(xverts,yverts,colorMap(i,:),'LineWidth',0.5);
                end
            end
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_PS_ResultsFT.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end
            
        end
        
        %% Result presentation method: Selectivities
        function showResultSDiesel(obj, varargin)
            %SHOWRESULTFT will display hydrocarbon distribution
            %   Mass fraction over carbon number
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            
            % Plot product mass fractions and diesel fraction
            fig2 = figure('name', 'Selectivity');
            figure(fig2);
            cmPre = parula(obj.numSuccess*1.1);
            colorMap = cmPre(1:obj.numSuccess,:);
            hold on
            for i = 1:obj.numSuccess
                plot(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z,'color',colorMap(i,:),'LineWidth',1.5);
                text(20, 0.035 + 0.035*(i/10), cat(2, 'S_{CO,C10-C20} = ', num2str(round(obj.parStudyRct{obj.casesSuccess(i)}.SCarbon_CO_Diesel,4))),'FontSize',16,'color',colorMap(i,:))
            end
            hold off
            box on
            xlabel('n_C / -','FontSize',20);
            ylabel('S_{CO,i} / -');
            axis([0 30 0 0.1]);
            ax = gca;
            ax.FontSize = 16;
            title('Selectivity of CO to each hydrocarbon');
            set(gca,'linewidth',1.5)
            
            % patches for bar diagram
            for i = 1:obj.numSuccess
                if obj.parStudyRct{obj.casesSuccess(i)}.numSpeciesInlet < 10
                    disp('No mass balance for diesel fraction calculated.')
                elseif obj.parStudyRct{obj.casesSuccess(i)}.numSpeciesInlet < 26
                    xverts = [(10:10+length(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:end,1))); (11:10+(length(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:end,1))+1)); (11:10+(length(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:end,1))+1)); (10:10+length(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:end,1)))];
                    yverts = [zeros(1,1+length(obj.parStudyRct{i}.S_Carbon_CO_Z(11:end,1))); zeros(1,1+length(obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:end,1))); obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(10:end,1)'; obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(10:end,1)'];
                    patch(xverts,yverts,colorMap(i,:),'LineWidth',0.5);
                else
                    xverts = [(10:20); (11:21); (11:21); (10:20)];
                    yverts = [zeros(1,11); zeros(1,11); obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(11:21,1)'; obj.parStudyRct{obj.casesSuccess(i)}.S_Carbon_CO_Z(10:20,1)'];
                    patch(xverts,yverts,colorMap(i,:),'LineWidth',0.5);
                end
            end
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_PS_ResultsSCO.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end
            
        end
        
        %% Write Results into Table
        function tableOutput(obj, varargin)
            %TABLEOUTPUT creates cell array of tables containing results of
            %parameter study cases as requested by user.
            %first column of each table contains the fixed variable, second
            %shows if simulation was successful
            %both input arguments need to be passed to method, otherwise an
            %error occurs
            %if value of the fixed variable isn't being used in any of the
            %cases, an empty table is returned
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', false);
            addParameter (parser,'tableParFix', []);
            addParameter (parser,'tableParVar', []);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            obj.tableParFix = parser.Results.tableParFix;
            obj.tableParVar = parser.Results.tableParVar;
            
            tableParFixVar = cell(1,length(obj.tableParFix)/2);
            tableParFixVal = zeros(1,length(obj.tableParFix)/2);
            for i = 1:(length(obj.tableParFix)/2)
                % Names of fixed variables, creating a table for each of them.
                tableParFixVar{i} = obj.tableParFix{i*2-1};
                % Values of fixed variables
                tableParFixVal(i) = str2double(obj.tableParFix{i*2});
            end
            
            % preallocate
            tableArr = zeros(length(obj.parStudyRct),length(obj.tableParVar));
            % number of tables generated i
            for i = 1:length(obj.tableParFix)/2
                % write fixed variable in beginning of tableParVar so that
                % first column contains fixed variable
                obj.tableParVar = horzcat(tableParFixVar{i},obj.tableParVar);
                % number of colums j
                for j = 1:length(obj.tableParVar)
                    % number of rows k
                    for k = 1:length(obj.parStudyRct)
                        % check if variable that should be written in tableis not empty && check if fixed value of table is equal to the corresponding value of the case && check if case was successfull
                        % if one of those isn't true, entry in table = NaN
                        if ~isempty(eval(['obj.parStudyRct{k}.',obj.tableParVar{j}])) && eval(['obj.parStudyRct{k}.',tableParFixVar{i}]) == tableParFixVal(i) && obj.parStudyRct{k}.caseStatRct
                            tableArr(k,j) = eval(['obj.parStudyRct{k}.',obj.tableParVar{j}]);
                        else
                            tableArr(k,j) = NaN;
                        end
                    end
                end
                % add column to show which row (simulation) resulted in an error
                tableArr = horzcat(tableArr(:,1),obj.caseStat(2,:)',obj.caseStat(1,:)',tableArr(:,2:end));
                % modify fixed variable name to prevent duplicate variable
                % name in table
                obj.tableParVar(1) = {[tableParFixVar{i},'_fix']};
                % add column name to tableParVar
                obj.tableParVar = horzcat(obj.tableParVar(:,1),{'case number'},{'success?'},obj.tableParVar(:,2:end));
                % delete lines in which fixed variable is NaN
                tableArr(isnan(tableArr(:,1)),:) = [];
                % write temporary array into cell array of tables
                obj.tablePS{i} = array2table(tableArr,'VariableNames',obj.tableParVar);
                
                % write current table into xlsx file
                if saveSwitch
                    filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_tablePS.xlsx');
                    if ~isfolder(strcat('results/',char(obj.sDate)))
                        mkdir('results/',obj.sDate)
                    end
%                     if i == 1
%                         writetable(obj.tablePS{i},filename,'Delimiter',';');
%                     elseif isfile(filename) && i >= 2
%                         tableCurr = readtable(filename);
%                         tableCurr = horzcat(tableCurr,obj.tablePS{i});
%                         writetable(tableCurr,filename,'Delimiter',';');
%                     end
                    
                    if isfile(filename) && i == 1
                        delete(filename)
                    end
                    writetable(obj.tablePS{i},filename,'Sheet',string(strcat(tableParFixVar(i),'=',num2str(tableParFixVal(i)))));
                end
                
                % reset tableParVar
                obj.tableParVar(:,1:3) = [];
                % reset tableArr to zero array
                clear tableArr tableCurr range
                tableArr = zeros(length(obj.parStudyRct),length(obj.tableParVar));
            end
        end
    end
end

