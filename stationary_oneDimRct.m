classdef stationary_oneDimRct < handle
    %STATIONARY_ONEDIMRCT is a class that is intended to solve system of coupled
    % PDEs. The system stands for a plug flow reactor. The spatial coordinate
    % refers to the direction of flow. Such a reactor can be
    % described with PDEs for the conservation of species (plural) and enthalpy
    % within the reaction zone, as well as the conservation of enthalpy within
    % a zone that corresponds with the reaction zone, e.g. the reactor wall and
    % even a cooling zone.
    % In the one dimensional case it is assumed, that there are no
    % gradients in the independent variables perpendicular to the spatial
    % coordinate.
    % The different PDEs vary in the terms they possess and the respective coefficients,
    % because different physical properties are conserved and the materials/
    % material properties are different.
    % '< handle' means, that class inherits from handle, means among others that
    % class methods can write to properties.
    % The class consists of properties and methods. The variables needed
    % to model the reactor are stored in the properties as name value
    % pairs.
    % The methods include among others: the constructor, one method
    % that generates a domain and mesh, one methods that initialize the
    % coefficients, one method that solves the pde system, methods for
    % displaying.
    % The 1 dimensional solution is implemented as solution of a rectangular
    % (2-dimensional) domain with special von neumann BC (gradient 0)
    % at upper and lower domain boundary for every PDE, which results in a quasi
    % infinetely extended domain.
    % 
    %                        upper boundary
    %                            _ _ _ _
    %   left boundary (inlet)   |       |   right boundary (outlet)
    %                            - - - - 
    %                        lower boudnary
    % 
    %   Balance equations for moles, enthalpies are kind of stacked onto
    %   each other (more precisely equations are solved on the same nodes).
    
    properties
        
        % Properties at the inlet
        TInlet                % double: Temperature at inlet in K
        pInlet                % double: pressure at inlet in bar
        speciesNames          % cell column vector: species names
        numSpeciesInlet       % double: number of Species being calculated
        mixInlet              % Instance of MIXTURE: object allowing to calculate thermodynamic properties of inlet reactant mixtures
        yInlet                % double column vector of length numSpeciesInlet: inlet molar fractions / -
        nInlet                % double column vector: initial mole flow in mol/s of all species
        VInlet                % double: Volumetric flow at inlet in m^3/s
        
        % Properties for Weisz modulus
        nRct                  % line vector of double: Reaction order of the species
        rhoBp                 % line vector of double: density at boiling point in mol/m^3
        etaAlkane                % double: viscosity of liquid alkane at 228°C in Pa s
        rhoAlkane                % double: density of liquid alkane at 228°C in kg/m^3
        dCat                    % double: diameter of catalyst particles in m
        epsilonCat            % double: porosity of catalyst
        tauCat                % double: tortuosity of catalyst
        rhoCatP               % double: density of a catalyst particle including porosity kg/m^3
        Psi                   % array of double: Weisz modulus for H2 and CO
        
        % Properties for Mears Criterion
        epsilonBed            % double: porosity of packed bed
        tauBed                % double: tourtuositiy of packed bed
        Bo                    % double: Bodenstein-Zahl
        lhsMears              % array of double: left hand side of mears criterion equation
        rhsMears              % array of double: right hand side of mears criterion equation
        mCrit                 % array of double: contains lhs, rhs and result of Mears Criterion
        
        % Properties for radial Heat transfer criterion
        radHeatLhs
        radHeatLhsMax
        radHeatRhs
        radHeatRhsMax
        radHeatCrit
        
        % Properties relevant for reaction kinetics
        kinFunName            % string: determines which rate law file will be used;
                              %         rate law has to fit the apparent species
        rhoCat                % double: Catalyst density in kg/mol
        massCat               % double: Catalyst mass in kg
        rateFactor            % double: factor to manipulate reaction rate

        % Properties relevant for the reaction zone
        areaCrossRct          % double: Crosssectional area of the reactor in m^2
        dRct                  % double: diameter of the reactor in m
        URct                  % double: circumfence of the reactor volume in m
        wRct                  % double: width of the slit reactor in m
        lengthRct             % double: length of the reactor in m
        areaWallRct           % double: inner wall area in m^2
        volumeRct             % double: reactor volume in m^3
        volRedRct             % double: volumeRct correction in m^3 for internal structures by substraction, only used in kMode=4
        WHSV                  % double: WHSV (Weight Hourly Space velocity) in 1/h
        ratioInlet            % double: molar ration of H2 to CO at inlet
        kRct                  % double: Thermal Transmittance inner Wall in W/m^2/K
        alphaRct              % double: Heat Transfer Coefficient inner Wall in W/m^2/K
        DeltaRctH             % double: reaction Enthalpy in J/mol
                
        % Properties relevant for the wall
        thicknessWall         % double: wall thickness in m
        areaCrossWall         % double: crosssectional area of wall in m^2
        volumeWall            % double: volume of Wall in m^3
        dMid                  % double: diamater in between inner and outer wall at half wall thickness
        dLogMean                % double: logarithmic mean diameter of drct and dcool
        dLogMeanRct             % double: logarithmic mean diameter of drct and dLogMean
        dLogMeanCool            % double: logarithmic mean diameter of dLogMean and dCool
        thicknessWallmeanRct     % double: wall thickness between dRct and dLogMean
        thicknessWallmeanCool    % double: wall thickness between dLogMean and dCool
        solidPhase            % Instance of xMATERIAL containing properties of the solid phase (wall) 
        TWall                 % double: Temperature wall in K

        % Properties relevant for the cooling zone
        dCool                 % double: outer diameter of reactor tube in m
        UCool                 % double: circumfence of outer Wall in m
        areaWallCool          % double: outer Wall area in m^2
        kCool                 % double: Thermal Transmittance outer Wall in W/m^2/K
        rhoWall               % double: density of solidPhase (316L steel) in kg/m^3
        TCool                 % double: temperature the cooling side
        alphaCool             % double: Heat Transfer Coefficient outer Wall in W/m^2/K

        % Properties relevant to thermo-physical data
        numSpeciesGasPhase    % integer: number of species to be considered for properties calculation
                              %          has to be <= numSpeciesInlet
        gasPhase              % Instance of MIXTURE
        D                     % double: Diffusion coefficients in m^2/s
        DFactor
        DFinal                % array double: numSpeciesInlet x location.x Diffusion coefficient in m^2/s
        kFinal                 
        uFinal                % double: Gas velocity in m/s
        TConst                % double: Constant temperature used to calculate MIXTURE properties when constMIXTURESwitch is true in K
        lambdaRctFinal        % double line vector: documentation of thermal conductivity of gasphase at each node / J/(m*K*s)
        aRctFactor
        aWallFactor
        cPRctFinal
        rhoRctFinal
        aRctFinal
        VFinal
        RateFinal
        
        % Properties relevant to PDEtoolbox
        numEquations          % integer: number of PDE to be solved (one for each species+ one for every enthalpy balance)
        meshFactor            % double: factor by which mesh - not reactor dimensions are scaled
        meshLength            % double: length of apparatus
        meshHeight            % double: arbitrary height of mesh
        lengthDomain          % double: length of apparatus
        heightDomain          % double: arbitrary height of mesh
        elementSize           % integer: mesh metric, that influences calculation expense and solution quality
        numNodes              % double: number of nodes in mesh
        diaryOutput           % string: output of diary for one case created by parStudy
        firstResid            % double: first resid reported in the command window / no unit
        resTol                % double: residual tolerance, solver option
        modelPDEsystem        % object: model containing PDEs, mesh, ...
        resultPDE             % object: object containing results of PDE
        reportStats           % char: holds command for PDE option ReportStatistics
        TRctSwitch            % boolean: true, if next to the mole balances the
                              %          enthalpy balance of the reaction zone shall be calculated
        TWallSwitch           % boolean: true, if next to the mole balances and
                              %          enthalpy balance of the reaction zone the enthalpy
                              %          balance of the wall shall be calculated
        kModeSwitch           % boolean: true, if constant kRct shall be assumed
        constMIXTURESwitch    % boolean: true, if MIXTURE properties shall be calculated with constant pressure and temperature
        location              % structure: contains location structure used in coeff_x functions
        locXSort              % array of double: contains sorted coordinates in x direction
        idxSort               % array of double: contains indices of sorted coordinates used to sort xFinal values
        
        % Properties relevant to case
        sDateTime             % string: current date and time in this format: 'yyyy-mm-dd_HH-MM-SS'
        sDate                 % string: current date in this format: 'yyyy-mm-dd'
        computationTime       % double: computation time for the entire program in s
        counter               % double: counts how often coeff functions are called
        caseStatRct           % logical: 0 for error, 1 for success
        caseStatRctMsg        % string: contains information about success or failure of simulation
        
        % Natural constants
        R = 8.3144626181532   % double: universal gas constant in J/(mol*K)
        
        % Results
        xq                    % double array: vector of values for interpolation
        xqMulti               % double array: multiple copies of xq
        yq                    % vector of values for interpolation at y=0.5
        yqMulti               % vector of values for interpolation at different y
        nOut                  % double array: calculated mole flow in mol/s; one column for each step along reactor axis
        nOutMulti             % cell array: calculated mole flow in mol/s; one column for each step along reactor axis, each cell at different y
        nOutMean              % double array: calculated mole flow in mol/s; one column for each step along reactor axis mean
        TRctOut               % double line vector: calculated reactor temperature in K for each step along reactor axis
        TRctOutMulti          % cell array: calculated  reactor temperature in K; one column for each step along reactor axis, each cell at different y
        TRctOutMean           % double line vector: calculated reactor temperature in K for each step along reactor axis mean
        deltaTRctMax          % double: maximum temperature difference of reactor volume in K
        TWallOut              % double line vector: calculated reactor temperature in K for each step along reactor wall
        TWallOutMulti         % cell array: calculated wall temperature in K; one column for each step along reactor axis, each cell at different y
        TWallOutMean          % double line vector: calculated wall temperature in K for each step along reactor axis mean
        deltaTWallMax         % double: maximum temperature difference of wall in K
        w_Hydrocarbons        % double: mass fraction of C1 to C50
        w_Diesel              % double: mass fraction of C10 to C20
        w_DieselOut           % double: mass fraction of C10 to C20 at outlet
        X_H2                  % double: conversion of H2
        X_H2Mean              % double: conversion of H2 calculated from mean solution
        X_CO                  % double: conversion of CO
        X_COMean              % double: conversion of CO calculated from mean solution
        X_TotalMean           % double: conversion of CO and H2
        nOutTot               % column vector: one value of the overall mole flow for each step of the solver
        y_iOut                % double array: calculated mole fraction one column for each solver step
        mOut                  % double array: calculated mass flows of each species one column for each solver step
        mTot                  % double array: total mass flow at each point in the reactor
        w_iOut                % double array: calculated mass fractions of each species one column for each solver step
        slope                 % double array: slope ( ASF plot ln(w_i/nC) ) only at outlet
        S_Carbon_CO_Z         % double: selectivity of CO to each hydrocarbon
        SMolar_CO_C1to4       % double: selectivity of CO to CH4...C4H10
        SCarbon_CO_C1to4      % double: carbon selectivity of CO to CH4...C4H10
        SCarbon_CO_C1to9      % double: carbon selectivity of CO to CH4...C9H20
        SMolar_CO_Diesel      % double: selectivity of CO to Diesel fraction
        SCarbon_CO_Diesel     % double: carbon selectivity of CO to Diesel fraction
        SCarbon_CO_C21to50    % double: carbon selectivity of CO to C21H44...C50H102
        YMolar_CO_C1to4       % double: yield of CO to CH4...C4H10
        YCarbon_CO_C1to4      % double: carbon yield of CO to CH4...C4H10
        YMolar_CO_Diesel      % double: yield of CO to Diesel
        YCarbon_CO_Diesel     % double: carbon yield of CO to Diesel
        STYDiesel             % double: space time yield of Diesel in kg/l/h
        STYHC                 % double: space time yield of hydrocarbons in kg/l/h
        MTYDiesel             % double: mass time yield of Diesel in g/g_cat/h
        MTYHC                 % double: mass time yield of hydrocarbons in g/g_cat/h
        
    end % end of properties
    
    methods
        %% Constructor method of the class
        function obj = stationary_oneDimRct(varargin)
            %ONEDIMRCTWALL Constructs an instance of the class ONEDIM_RCTWALL
            %   The constructor receives information in name value pairs
            %   that is needed to initialize a one-dimensional pde-toolbox
            %   solution. For details regarding the name value pairs
            %   please check out the input parser definition.
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'kineticFunctionName','default');
            addParameter (parser,'rctEnthalpy', -160000);
            addParameter (parser,'massCat', 0.0015);
            addParameter (parser,'numSpeciesInlet', 10);
            addParameter (parser,'numSpeciesGasPhase', 10);
            addParameter (parser,'speciesNames', {});
            addParameter (parser,'meshFactor', 1);
            addParameter (parser,'meshLength', 0.1);
            addParameter (parser,'meshHeight',0.15);
            addParameter (parser,'elementSize',0.01);
            addParameter (parser,'resTol',1e-7);
            addParameter (parser,'Timespan',1:1:50);
            addParameter (parser,'inletMoleFlow', []);
            addParameter (parser,'dRct', 1e-3);
            addParameter (parser,'URct', 16e-3);
%             addParameter (parser,'areaCrossRct', 63.617e-6);
            addParameter (parser,'inletTemperature', 473.15);
            addParameter (parser,'inletPressure', 20);
            addParameter (parser,'wallTemperature', 273 + 150);
            addParameter (parser,'kRct', 100);
            addParameter (parser,'kCool', 100);            
            addParameter (parser,'lengthRct', 0.03144);
            addParameter (parser,'wRct', 8);
            addParameter (parser,'volRedRct', 8);
            addParameter (parser,'rhoCat', 750);
            addParameter (parser,'VelocityCheat', 1);
            addParameter (parser,'reportStatistics', 'off');
            addParameter (parser,'timeStamp',[]);
            addParameter (parser,'dateStamp',[]);
            addParameter (parser,'TRctSwitch',false);
            addParameter (parser,'TWallSwitch',false);
            addParameter (parser,'kModeSwitch',true);
            addParameter (parser,'constMIXTURESwitch',false);
            addParameter (parser,'TConst',273.15);
            addParameter (parser,'TCool',373.15);
            addParameter (parser,'rhoWall',7900);
            addParameter (parser,'thicknessWall',0.001);
            addParameter (parser,'alphaRct',57.5);
            addParameter (parser,'alphaCool',10000);
            addParameter (parser,'WHSV',6.7);
            addParameter (parser,'ratioInlet',2.1);
            addParameter (parser,'rateFactor',1);
            addParameter (parser,'DFactor',1);
            addParameter (parser,'aRctFactor',1);
            addParameter (parser,'aWallFactor',1);
            
            addParameter (parser,'nRct',[0.5 1]);
            addParameter (parser,'etaAlkane',368e-6);
            addParameter (parser,'rhoAlkane',681.2);
            addParameter (parser,'dCat',100e-6);
            addParameter (parser,'epsilonCat',0.5);
            addParameter (parser,'tauCat',2);
            addParameter (parser,'rhoCatP',1830);
            
            addParameter (parser,'epsilonBed',0.5);
            addParameter (parser,'tauBed',2);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Four options for assigning property values:
            % 1a) Name-value-pair in class call.
            % 1b) If no value is provided default values are assigned.
            % 2)  Assign user values.
            % 3)  Assign nothing, because property will be filled by methods.
            
            obj.kinFunName = parser.Results.kineticFunctionName;
            obj.DeltaRctH = parser.Results.rctEnthalpy;
            obj.rhoCat = parser.Results.rhoCat;
            obj.numSpeciesInlet = parser.Results.numSpeciesInlet;
            obj.numSpeciesGasPhase = parser.Results.numSpeciesGasPhase;
            obj.speciesNames = parser.Results.speciesNames;
%             obj.nInlet = parser.Results.inletMoleFlow;
            obj.TInlet = parser.Results.inletTemperature;
            obj.pInlet = parser.Results.inletPressure;
            obj.TWall = parser.Results.wallTemperature;
            obj.kRct = parser.Results.kRct;
            obj.kCool = parser.Results.kCool;
            obj.meshFactor = parser.Results.meshFactor;
            obj.elementSize = parser.Results.elementSize;
            obj.resTol = parser.Results.resTol;
            obj.TRctSwitch = parser.Results.TRctSwitch;
            obj.TWallSwitch = parser.Results.TWallSwitch;
            obj.kModeSwitch = parser.Results.kModeSwitch;
            obj.constMIXTURESwitch = parser.Results.constMIXTURESwitch;
%             obj.TConst = parser.Results.TConst;
%             obj.TCool = parser.Results.TCool;
            obj.TConst = parser.Results.inletTemperature;
            obj.TCool = parser.Results.inletTemperature;
            obj.rhoWall = parser.Results.rhoWall;
            obj.thicknessWall = parser.Results.thicknessWall;
            obj.alphaRct = parser.Results.alphaRct;
            obj.alphaCool = parser.Results.alphaCool;
            obj.WHSV = parser.Results.WHSV;
            obj.ratioInlet = parser.Results.ratioInlet;
            obj.reportStats = parser.Results.reportStatistics;
            obj.rateFactor = parser.Results.rateFactor;
            obj.DFactor = parser.Results.DFactor;
            obj.aRctFactor = parser.Results.aRctFactor;
            obj.aWallFactor = parser.Results.aWallFactor;
            
            obj.nRct = parser.Results.nRct;
            obj.etaAlkane = parser.Results.etaAlkane;
            obj.rhoAlkane = parser.Results.rhoAlkane;
            obj.dCat = parser.Results.dCat;
            obj.epsilonCat = parser.Results.epsilonCat;
            obj.tauCat = parser.Results.tauCat;
            obj.rhoCatP = parser.Results.rhoCatP;
            
            obj.epsilonBed = parser.Results.epsilonBed;
            obj.tauBed = parser.Results.tauBed;
            
            obj.mixInlet = MIXTURE('names', obj.speciesNames(1:obj.numSpeciesInlet),'constMIXTURESwitch',obj.constMIXTURESwitch,'TConst',obj.TConst,'pConst',obj.pInlet);
            obj.gasPhase = MIXTURE('names',obj.speciesNames(1:obj.numSpeciesGasPhase),'constMIXTURESwitch',obj.constMIXTURESwitch,'TConst',obj.TConst,'pConst',obj.pInlet);
            
            obj.rhoBp(1,1) = obj.gasPhase.materials(1,1).rhoBp;
            obj.rhoBp(2,1) = obj.gasPhase.materials(2,1).rhoBp;
            
            obj.wRct = parser.Results.wRct;
            obj.lengthRct = parser.Results.lengthRct;
            obj.volRedRct = parser.Results.volRedRct;
            
%             obj.areaCrossRct = parser.Results.areaCrossRct;
%             obj.dRct = sqrt(4.0*obj.areaCrossRct/pi);
            obj.dRct = parser.Results.dRct;
            
            if obj.kModeSwitch == 4
                % case slit
                obj.areaCrossRct = obj.dRct * obj.wRct;
                obj.areaCrossWall = obj.thicknessWall * obj.wRct * 2;
                obj.URct = 2*obj.wRct;
                obj.UCool = obj.URct;
                obj.areaWallRct = 2 * obj.lengthRct * obj.wRct + 2 * obj.lengthRct * obj.dRct;
                obj.volumeRct = obj.dRct * obj.wRct * obj.lengthRct;
                obj.volumeRct = obj.volumeRct - obj.volRedRct;
                obj.areaWallCool = obj.areaWallRct;
                obj.dCool = obj.dRct + 2*obj.thicknessWall;
                obj.volumeWall = obj.areaWallRct * obj.thicknessWall;
            else
                % case pipe
                obj.areaCrossRct = pi/4*obj.dRct^2;
                obj.URct = pi*obj.dRct;
                obj.areaWallRct = pi * obj.dRct * obj.lengthRct;
                obj.volumeRct = pi/4 * obj.dRct^2 * obj.lengthRct;
                obj.volumeRct = obj.volumeRct - obj.volRedRct;
                obj.areaWallCool = pi * (obj.dRct+2*obj.thicknessWall) * obj.lengthRct;
                obj.dCool = obj.dRct + 2*obj.thicknessWall;
                obj.areaCrossWall = pi/4*(obj.dCool^2 - obj.dRct^2);
                obj.UCool = pi*obj.dCool;
                obj.volumeWall = pi/4 * (obj.dRct+2*obj.thicknessWall)^2 * obj.lengthRct;
                obj.dMid = obj.dRct + obj.thicknessWall;
                obj.dLogMean = (obj.dCool - obj.dRct)/(log(obj.dCool/obj.dRct));
                obj.dLogMeanRct = (obj.dMid - obj.dRct)/(log(obj.dMid/obj.dRct));
                obj.dLogMeanCool = (obj.dCool - obj.dMid)/(log(obj.dCool/obj.dMid));
                obj.thicknessWallmeanRct = obj.dLogMean - obj.dRct;
                obj.thicknessWallmeanCool = obj.dCool - obj.dLogMean;
            end
            
            obj.meshFactor = 1;
            obj.meshHeight = obj.elementSize;           % height of of mesh in m
            obj.meshLength = obj.lengthRct;
            
            obj.lengthDomain = obj.meshLength * obj.meshFactor;
            obj.heightDomain = obj.meshHeight * obj.meshFactor;
            
            % Calculate nInlet from WHSV (Weight Hourly Space velocity) and
            % ratio of H2 and CO
            obj.massCat = obj.volumeRct * obj.rhoCat;
            ratioVector = zeros(obj.mixInlet.N,1);
            obj.nInlet = zeros(obj.mixInlet.N,1);
            ratioVector(1,1) = 1/(1+1/obj.ratioInlet);
            ratioVector(2,1) = 1 - (1/(1+1/obj.ratioInlet));
            nInletTot = obj.WHSV * obj.massCat / (3600 * obj.mixInlet.molarWeight('y',ratioVector));
            obj.nInlet(1,1) = nInletTot * 1/(1+1/obj.ratioInlet);
            obj.nInlet(2,1) = nInletTot * (1 - 1/(1+1/obj.ratioInlet));
            
            % Clarify number of equations to be solved
            obj.numEquations = obj.numSpeciesInlet;
            if obj.TRctSwitch
                obj.numEquations = obj.numEquations + 1;
            end
            if obj.TWallSwitch
                obj.numEquations = obj.numEquations + 1;
            end
            
            % Calculate obj.yInlet
            obj.yInlet = zeros(obj.numSpeciesInlet, 1);
            sumMolarFlow = sum(obj.nInlet);
            for i = 1:1:obj.numSpeciesInlet
                obj.yInlet(i, 1) = obj.nInlet(i, 1)./sumMolarFlow;
            end
            
            % Calculating current volumetric flow
            obj.VInlet = sumMolarFlow * obj.R * obj.TInlet / (obj.pInlet * 1e5);

            obj.counter = 0;

            % Results
            obj.sDateTime = parser.Results.timeStamp;
            obj.sDate = parser.Results.dateStamp;
            
            % wall properties
            obj.solidPhase = xMATERIAL('name','316L');
            
        end

        %% Method for model generation
        function initModel(obj)
            %INITMODEL creates model which contains PDEs, mesh, etc. with
            %   the given number of Equations
            
            % Initialize model
            obj.modelPDEsystem = createpde(obj.numEquations);
            
            % Set solver options
            obj.modelPDEsystem.SolverOptions.ReportStatistics = obj.reportStats;
%             obj.modelPDEsystem.SolverOptions.AbsoluteTolerance = 1e-9;
%             obj.modelPDEsystem.SolverOptions.RelativeTolerance = 1e-9;
            obj.modelPDEsystem.SolverOptions.ResidualTolerance = obj.resTol;
%             obj.modelPDEsystem.SolverOptions.MinStep = 0;
%             obj.modelPDEsystem.SolverOptions.MaxIterations = 500;
            obj.modelPDEsystem.SolverOptions.ResidualNorm = 'energy';
        end
        
        %% Domain and Mesh Generation
        function domainMesh(obj)
            %DOMAINMESH generates domains and mesh for PDE toolbox
            %   calculation representing the 1D reactor, walls and if applicable jacket
            
            % Create geometry matrix, a rectangle in this case
            square = [3 4 0 obj.lengthDomain obj.lengthDomain 0 0 0 obj.heightDomain obj.heightDomain]';
            gd = square;
            sf = 'square';
            ns = char('square')';
            
            % Create regions (faces) from the geometry
            geo = decsg(gd,sf,ns);
            
            % Geometrie dem Modell zuweisen
            obj.modelPDEsystem.Geometry = geometryFromEdges(obj.modelPDEsystem,geo);
            
            % Mesh erzeugen
            warning off MATLAB:subscripting:noSubscriptsSpecified
            obj.modelPDEsystem.Mesh = generateMesh(obj.modelPDEsystem,'Hmax',obj.elementSize);
            
            %Anzahl Knoten als Eigenschaft ablegen
            obj.numNodes = size(obj.modelPDEsystem.Mesh.Nodes,2);
        end
        
        %% Calculate material properties for current state
        % These functions are used to calculate the material properties in
        % dependance of current Temperature, gas mixture and so on.
        % They get called repeatedly by the PDE solver
        
        function testparam = testmethod(~,a)
            testparam = a^2;
        end

        function y = moleFrac(~,nCurr,nSum)
            %MOLEFRAC calculates the molar fraction of every species
            % nCurr represents the current mole flows of the components
            % y for each species is calculated
            % nSum is the total mole flow for each location
            y = nCurr ./ nSum;
        end

        function gasPhaseComposition = gasPhComp(obj,nCurr,nSum)
            % GASPHCOMP ??? normierung der y damit summe wieder 1 ergibt?
            y = moleFrac(obj,nCurr,nSum);
            gasPhaseComposition = y(1:obj.numSpeciesGasPhase,:) ./ sum(y(1:obj.numSpeciesGasPhase,:));
        end
        
        function cPRct = specHeatCap(obj,T,nCurr,nSum)
            % SPECHEATCAP calculates specific heat capacity of the mixture
            % at current Temperature in J/kg/K
            gasPhaseComposition = gasPhComp(obj,nCurr,nSum);
            cPRct = obj.gasPhase.specHeatCap('T', T, 'y', gasPhaseComposition);
        end
        
        function rho_g = densityGas(obj,p,T,nCurr,nSum)
            % DENSITYGAS calculates density of the mixture
            % at current Temperature in kg/m^3
            gasPhaseComposition = gasPhComp(obj,nCurr,nSum);
            rho_g = p*1e5*obj.gasPhase.molarWeight('y', gasPhaseComposition)./(obj.R*T);
        end
       
        function lambdaRct = thermCondGas(obj,pLambda,T,nCurr,nSum)
            % THERMCONDGAS calculates thermal conductivity of the mixture
            % at current Temperature and Pressure in W/m/K
            gasPhaseComposition = gasPhComp(obj,nCurr,nSum);
            lambdaRct = obj.gasPhase.thermCond('p', pLambda, 'T', T, 'y', gasPhaseComposition);
            
%             lambdaRct = real(lambdaRct);
            if ~isreal(lambdaRct)
                lambdaRct = 0;
            end
        end
        
        function V = volumeFlow(obj,p,T,nCurr)
            % VOLUMEFLOW calculates the total volume Flow of the mixture at
            % current pressure and temperature via ideal gas law in m^3/s
            sumMolarFlow = sum(nCurr);
            V = sumMolarFlow * obj.R .* T / (p * 1e5);
        end
        
        function D = diff(obj,pCurr,TCurrRct,nCurr,nSum)
            % DIFF calculates the diffusivity of each component in the
            % mixture at current pressure and temperature in m^2/s
            gasPhaseComposition = gasPhComp(obj,nCurr,nSum);
            D = obj.gasPhase.diffusivityMix('p', pCurr, 'T', TCurrRct, 'y', gasPhaseComposition);
        end
        
        function eta = dynViscGas(obj,pCurr,TCurrRct,nCurr,nSum)
            gasPhaseComposition = gasPhComp(obj,nCurr,nSum);
            eta = obj.gasPhase.dynVisc('p', pCurr, 'T', TCurrRct, 'y', gasPhaseComposition);
        end
        
        %% Methods to evaluate results
        function WeiszModulus(obj)
            %WEISZMODULUS calculates the Weisz modulus for H2 and CO at
            %   each location in the reactor

            % partial pressures in bar
            p_i = obj.y_iOut * obj.pInlet;
            % calculate Henry Constant for CO and H2 at each point in bar
            HiL = ones(2,length(obj.TRctOutMean));
            for i = 1:2
                HiL(i,:) = obj.gasPhase.materials(i).henryConst('T',obj.TRctOutMean);
            end
            % calculate molar fractions in liquid
            xiL = p_i(1:2,:) ./ HiL;
            % calculate concentration in liquid
            ciL = xiL .* obj.rhoAlkane / obj.mixInlet.materials(35).molarWeight;
            % calculate diffusion coefficient in liquid
            DiL = 5.88e-17 * (obj.TRctOutMean * sqrt(obj.mixInlet.materials(35).molarWeight))./(obj.etaAlkane*(obj.rhoBp.^(-0.6)));
            % calculate effective diffusion coefficient in pores
            DiEff = DiL * (obj.epsilonCat/obj.tauCat);
            
%             % % Calculate Rate with xFinal values put out by coeff_x
%             % % functions, doesnt work in parStudy!!
%             % sort location.x and get indices
%             [locXSort,idxSort] = sort(obj.location.x);
%             % write Rate property to variable
%             Rate_i = obj.RateFinal(1:end,:);
%             % sort Rate with the help of indices
%             Rate_i = Rate_i(1:end,idxSort);
            
%             % delete duplicate x locations
%             % round first to 13th digit otherwise unique doesnt work correctly
%             locXSortround = round(locXSort,13);
%             [locXSortUni,idxUni] = unique(locXSortround);
%             % delete duplicates in Rate_i with idx
%             Rate_i = Rate_i(1:end,idxUni);
%             % create interpolated location vector
%             locXSortUniInterp = linspace(locXSortUni(1),locXSortUni(end),length(obj.xq));
%             % interpolate Rate_i to match other arrays
%             Rate_i_interp = ones(2,length(locXSortUniInterp));
%             for i = 1:2
%                 Rate_i_interp(i,:) = interp1(locXSortUni,Rate_i(i,:),locXSortUniInterp);
%             end
%             Rate_i = Rate_i_interp;
            
            % % Calculate Rate with output values
            Rate_i = obj.rateFactor*rateLawKwackViaKirsch(obj.R,obj.TRctOutMean,p_i,obj.numSpeciesInlet);
            
            % calculate Weisz Modulus
            obj.Psi = (obj.nRct' + 1)/2 .* (obj.dCat/6).^2 .* (abs(Rate_i(1:2,:)) * obj.rhoCatP)./(DiEff .* ciL);
        end
        
        function MearsCriterion(obj)
            %MEARSCRITERION calculates the Weisz modulus for H2 and CO at
            %   each location in the reactor
            
            % calculate left hand side of criterion
            obj.lhsMears = obj.lengthRct / obj.dCat;
            
            % calculate Diffusion coefficients with output values
            D_i = obj.diff(obj.pInlet,obj.TRctOutMean,obj.nOut,obj.nOutTot);
            % calculate velocity with output values
            V = obj.volumeFlow(obj.pInlet,obj.TRctOutMean,obj.nOut);
            u = V./(obj.volumeRct/obj.lengthRct);
            % calculate Bodensteinzahl
            obj.Bo = (obj.epsilonBed * D_i(1:2,:) ./ (obj.tauBed * obj.dCat * u) +0.5).^-1;
            % Put conversions into one array
            X_iMean = [obj.X_H2Mean; obj.X_COMean];
            
            % calculate right hand side of criterion
            obj.rhsMears = (8./obj.Bo) .* obj.nRct' .* log(1 ./ (1 - X_iMean));
            
            mCritLogical = obj.lhsMears > obj.rhsMears;
            obj.mCrit = [repmat(obj.lhsMears,1,length(obj.rhsMears));obj.rhsMears;mCritLogical];
        end
        
        function radHeat(obj)
            %RADHEAT calculates the criterion for negligibility of radial
            %heat transfer limitations
            
            lambdaRct = obj.thermCondGas(obj.pInlet,obj.TRctOutMean,obj.nOut,obj.nOutTot);
            p_i = obj.y_iOut * obj.pInlet;
            Rate_i = obj.rateFactor*rateLawKwackViaKirsch(obj.R,obj.TRctOutMean,p_i,obj.numSpeciesInlet);
            
            % calculate left hand side of criterion
            obj.radHeatLhs = ((obj.dRct/2)^2 * abs(obj.DeltaRctH) * obj.rhoCat .* abs(Rate_i(2,:)))./(lambdaRct .* obj.TWallOutMean);
            obj.radHeatLhsMax = max(obj.radHeatLhs);
            
            % Larger one of the two activation energies in rate law
            E_a =  99.5e3;
            % calculate right hand side of criterion,
            obj.radHeatRhs = 0.4*(obj.R * obj.TWallOutMean)/(E_a);
            obj.radHeatRhsMax = max(obj.radHeatRhs);
            
            radHeatCritLogical = obj.radHeatLhs < obj.radHeatRhs;
            obj.radHeatCrit = [obj.radHeatLhs;obj.radHeatRhs;radHeatCritLogical];
        end
        
        %% Build PDE system: set coefficients
        function setCoeff(obj)
            % SETCOEFF sets coefficients of the PDEs

            obj.modelPDEsystem.EquationCoefficients.CoefficientAssignments = specifyCoefficients(obj.modelPDEsystem,...
                'm',0,...
                'd',0,...
                'c',@(location,state)coeff_c(location,state,obj),...
                'a',0,...
                'f',@(location,state)coeff_f(location,state,obj));
        end

        %% Build PDE system: set boundary conditions
        function setBC(obj)
            %SETBC sets the boundary conditions of the PDEs
            
            % BC for all mass balances
            BCInlet = obj.nInlet(:, 1)';
            
            % BC for enthalpy balance reaction zone
            if obj.TRctSwitch
                BCInlet(1, obj.numSpeciesInlet+1) = obj.TInlet; 
            end
            
%             BC for enthalpy balance wall
            if obj.TWallSwitch
                BCInlet(1, obj.numSpeciesInlet+2) = 0;%obj.TInlet; 
            end
            
            X = obj.numEquations;
            H = eye(X);
            H(X,X) = 0;
            Q = zeros(X,X);
            G = zeros(X,1);

            if ~obj.TWallSwitch
                obj.modelPDEsystem.BoundaryConditions.BoundaryConditionAssignments(1) = applyBoundaryCondition(obj.modelPDEsystem,'neumann','Edge',[1 2 3]);
                obj.modelPDEsystem.BoundaryConditions.BoundaryConditionAssignments(2) = applyBoundaryCondition(obj.modelPDEsystem,'dirichlet','Edge',4,'r',BCInlet);
            else
                obj.modelPDEsystem.BoundaryConditions.BoundaryConditionAssignments(1) = applyBoundaryCondition(obj.modelPDEsystem,'mixed', ...
                                                                                                                              'Edge',4,'h',H,'r',BCInlet, ...
                                                                                                                              'q',Q,'g',G);
            end


        end

        %% Build PDE system: set initial conditions - only for initial guess - not for transient calculation
        function setIC(obj)
            % SETIC sets the initial conditions, but not for a transient
            % calculation but as an initial guess

            %% Simple IC
            % IC for all mass balances, derived from inlet conditions
            % ie. the entire reactor is filled with educts
            IC = zeros(obj.numEquations, 1);
            IC(1, 1) = obj.nInlet(1,1);
            IC(2, 1) = obj.nInlet(2,1);

            %% Complex IC
            % Single column vectors from previous simulations are used. For
            % example the gas phase composition from the end of the reactor
            
            % Values from 2021-12-06_16-54-56_rct_55species_528K.mat,
            % rct.resultPDE.NodalSolution(2,:)
%             ICvector = [0.000143524184356123,6.74235905579598e-05,0,0,2.27680581794315e-05,2.43960885897444e-06,1.51086617975392e-07,3.14788703865545e-07,2.88871428985243e-07,2.61839115278109e-07,2.35690217535171e-07,2.11229713466009e-07,1.88752749562725e-07,1.68316218898679e-07,1.49861044970339e-07,1.33272277524189e-07,1.18410353218745e-07,1.21931896108780e-07,1.07669451464738e-07,9.50754905557917e-08,8.39548036571872e-08,7.41350230040832e-08,6.54639496870806e-08,5.78071932033498e-08,5.10460877417599e-08,4.50758519494187e-08,3.98039640840586e-08,3.51487276976045e-08,3.10380046473901e-08,2.74080968111647e-08,2.42027584834143e-08,2.13723242969706e-08,1.88729393706169e-08,1.66658794362966e-08,1.47169504599675e-08,1.29959586693874e-08,1.14762424816710e-08,1.01342592957945e-08,8.94922050323413e-09,7.90276936249735e-09,6.97869642085872e-09,6.16268837298316e-09,5.44210615682885e-09,4.80578914188257e-09,4.24388215708858e-09,3.74768275867436e-09,3.30950640820647e-09,2.92256738348848e-09,2.58087363416356e-09,2.27913394498630e-09,2.01267595679145e-09,1.77737377878581e-09,1.56958406035336e-09,1.38608955504517e-09,1.22404927945422e-09,519.531310403584]';
%             IC = ICvector(1:obj.numSpeciesInlet,1);
%             IC(end) = ICvector(end,1);
            
            % Values from 2021-12-07_09-19-02_rct_SimpleIC.mat,
            % rct.resultPDE.NodalSolution(2,:)
%             ICvector = [0.000166788902496799,7.94670722305565e-05,0,0,8.98544319820834e-06,1.08598281839261e-06,6.93264563631214e-08,1.46441837322504e-07,1.35414726550968e-07,1.23252080492153e-07,1.41456398332431e-07,1.24213467241680e-07,1.09077014539764e-07,9.57891135772403e-08,8.41234942720256e-08,7.38816480781820e-08,6.48894138284239e-08,5.69939840465903e-08,5.00612805831104e-08,4.39736534821989e-08,3.86278632394354e-08,3.39333115698618e-08,2.98104898580906e-08,2.61896185172609e-08,2.30094537455513e-08,2.02162409360318e-08,1.77627967257443e-08,1.56077038502555e-08,1.37146048682347e-08,1.20515826068134e-08,486.912313169468]';
%             IC = ICvector(1:obj.numSpeciesInlet,1);
%             IC(end) = ICvector(end,1);
            


            %% IC for energy balance
            % IC for energy balance reaction zone
            if obj.TRctSwitch
                IC(obj.numSpeciesInlet+1, 1) = obj.TInlet;
            end

            %% IC for wall
            % IC for energy balance wall
            if obj.TWallSwitch
                IC(obj.numSpeciesInlet+2, 1) = obj.TInlet;
            end

            %% Apply simple or complex IC to model
            obj.modelPDEsystem.InitialConditions.InitialConditionAssignments = setInitialConditions(obj.modelPDEsystem,IC);
            
            %% IC from previous results
            % Load IC from a file
%             load('2021-12-07_09-19-02_rct_SimpleIC','rct')
%             obj.modelPDEsystem.InitialConditions.InitialConditionAssignments = setInitialConditions(obj.modelPDEsystem,rct.resultPDE.NodalSolution(end,:)');
%             clear('rct')
            %% NodalSolution IC
            % Use entire NodalSolution from previous simulation
            % Seems to work only for using previous iterations of a time
            % dependent model
            
        end
        
        %% Solve PDE
        function results = solve(obj)
            %SOLVE uses the solvepde solver to solve the PDEs and saves the
            %   solution in the resultPDE-object
%             obj.resultPDE = solvepde(obj.modelPDEsystem);
            results = solvepde(obj.modelPDEsystem);
            obj.resultPDE = results;
        end    
            
        %% Evaluate Solution
        function evaluate(obj)
            %EVALUATE interpolates the solution (derived on 2D domain) to values
            % along the rector axis
            
            % reactor lenght from 0 to obj.lengthRct
            obj.xq = linspace(0,obj.lengthRct)';
            obj.xqMulti = repmat(obj.xq,1,9)';
            % "reactor height" (domain height) from 0 to obj.heightDomain
            obj.yq(1:length(obj.xq),1) = 0.5*obj.heightDomain;
            obj.yqMulti = (0.1:0.1:0.9)'.*ones(9,length(obj.xq)).*obj.heightDomain;

            % Mole balances
            obj.nOut = zeros((size(obj.yqMulti)));
            for j = 1:size(obj.yqMulti,1)
                for i = 1:obj.numSpeciesInlet
                    obj.nOut(i,:) = interpolateSolution(obj.resultPDE,obj.xq,obj.yq,i)';
                    obj.nOutMulti{j}(i,:) = interpolateSolution(obj.resultPDE,obj.xqMulti(j,:),obj.yqMulti(j,:),i)';
                end
            end
            % sum up the interpolated results and devide by number of
            % summands
            obj.nOutMean = sum(cat(3,obj.nOutMulti{:}),3)*1/(size(obj.yqMulti,1));
            
            % Enthalpy balance reaction zone if applicable
            if obj.TRctSwitch
                obj.TRctOut = interpolateSolution(obj.resultPDE,obj.xq,obj.yq,obj.numSpeciesInlet+1)';
                
                for j = 1:size(obj.yqMulti,1)
                    obj.TRctOutMulti{j}(1,:) = interpolateSolution(obj.resultPDE,obj.xqMulti(j,:),obj.yqMulti(j,:),obj.numSpeciesInlet+1)';
                end
                % sum up the interpolated results and devide by number of
                % summands
                obj.TRctOutMean = sum(cat(3,obj.TRctOutMulti{:}),3)*1/(size(obj.yqMulti,1));
            else
                obj.TRctOut = obj.TWall;
            end
            
            % Enthalpy balance wall if applicable
            if obj.TWallSwitch
                obj.TWallOut  = interpolateSolution(obj.resultPDE,obj.xq,obj.yq,obj.numSpeciesInlet+2)';
                
                for j = 1:size(obj.yqMulti,1)
                    obj.TWallOutMulti{j}(1,:) = interpolateSolution(obj.resultPDE,obj.xqMulti(j,:),obj.yqMulti(j,:),obj.numSpeciesInlet+2)';
                end
                % sum up the interpolated results and devide by number of
                % summands
                obj.TWallOutMean = sum(cat(3,obj.TWallOutMulti{:}),3)*1/(size(obj.yqMulti,1));
            end
            
            % Calculate maximum temperature difference
            obj.deltaTRctMax = max(obj.TRctOutMean) - obj.TInlet;
            obj.deltaTWallMax = max(obj.TWallOutMean) - obj.TInlet;
            
            % Calculate total mole flow at each solver step
            obj.nOutTot = sum(obj.nOut);
            
            % Calculate mole fractions at each solver step
            obj.y_iOut = obj.nOut./obj.nOutTot;
            
            % Calculate mass flow of each species (colum vector)for each
            % step (one column vector per step)
            obj.mOut = zeros(obj.numSpeciesInlet, size(obj.nOut,2));
            for i = 1:obj.numSpeciesInlet
                obj.mOut(i,:) = obj.nOut(i,:).*obj.mixInlet.materials(i).molarWeight;
            end
            
            % Calculate mass fractions at each solver step
            obj.mTot = sum(obj.mOut);
            obj.w_iOut = obj.mOut./obj.mTot;
            
            % Mass fraction of hydrocarbons
            obj.w_Hydrocarbons = obj.w_iOut(6:obj.numSpeciesInlet,:);
            obj.w_Hydrocarbons = obj.w_Hydrocarbons./sum(obj.w_Hydrocarbons);
            
            %  Mass fraction of diesel of hydrocarbons
            if obj.numSpeciesInlet < 15
                obj.w_Diesel = 0;
            elseif obj.numSpeciesInlet < 25
                obj.w_Diesel = sum(obj.w_Hydrocarbons(10:end,:));
                obj.w_DieselOut = obj.w_Diesel(1,end);
            else
                obj.w_Diesel = sum(obj.w_Hydrocarbons(10:20,:));
                obj.w_DieselOut = obj.w_Diesel(1,end);
            end
            
            % Calculate slope ( ASF plot ln(w_i/nC) ) only at outlet
            obj.slope = zeros(obj.numSpeciesInlet-5, 1);
            for i=1:obj.numSpeciesInlet-5
                obj.slope(i)=log(obj.w_Hydrocarbons(i,end)/i);
            end
            
            % Calculate conversions
            obj.X_H2 = 1 - obj.nOut(1, end)/obj.nInlet(1);
            obj.X_CO = 1 - obj.nOut(2, end)/obj.nInlet(2);
            
            obj.X_H2Mean = 1 - obj.nOutMean(1, end)/obj.nInlet(1);
            obj.X_COMean = 1 - obj.nOutMean(2, end)/obj.nInlet(2);
            
            obj.X_TotalMean = 1 - sum(obj.nOutMean(1:2,end))/sum(obj.nInlet(1:2));
            
            % Calculate molar selectivity of CO to Methane...Butane
            obj.SMolar_CO_C1to4 = (sum(obj.nOut(6:9,end)) - 0) / (obj.nOut(2,1) - obj.nOut(2,end));
            
            % Calculate carbon selectivity of CO to Methane...Butane
            STemp = zeros(obj.numSpeciesInlet - 5,1);
            % Diesel fraction includes species 15-25 (C10 - C20)
            for i = 6:9
                % Number of carbon atoms is equal to
                % Number of species - 5
                n_C = i - 5;
                STemp(i-5,:) = (obj.nOut(i,end) - 0) / (obj.nOut(2,1) - obj.nOut(2,end)) * n_C;
                obj.SCarbon_CO_C1to4 = sum(STemp);
            end
            
            % Calculate molar selectivity of CO to Diesel fraction (C10 -
            % C(20)
            if obj.numSpeciesInlet < 15
                obj.SMolar_CO_Diesel = 0;
            elseif obj.numSpeciesInlet < 25
                obj.SMolar_CO_Diesel = (sum(obj.nOut(10:end,end)-0)) / (obj.nOut(2,1) - obj.nOut(2,end));
            else
                obj.SMolar_CO_Diesel = (sum(obj.nOut(10:20,end)-0)) / (obj.nOut(2,1) - obj.nOut(2,end));
            end
            
            % Calculate carbon selectivity of CO to each Alkane
            obj.S_Carbon_CO_Z = zeros(obj.numSpeciesInlet - 5,1);
            for i = 6:obj.numSpeciesInlet
                % Number of carbon atoms is equal to
                % Number of species - 5
                n_C = i - 5;
                obj.S_Carbon_CO_Z(i-5,:) = (obj.nOut(i,end) - 0) / (obj.nOut(2,1) - obj.nOut(2,end)) * n_C;
            end
            % Calculate carbon selectivity of CO to C1 - C9
            obj.SCarbon_CO_C1to9 = sum(obj.S_Carbon_CO_Z(1:9,:));
            
            % Calculate carbon selectivity of CO to Diesel fraction (C10 -
            % C(20) and C21 to C50
            if obj.numSpeciesInlet < 15
                obj.SCarbon_CO_Diesel = 0;
                obj.SCarbon_CO_C21to50 = 0;
            elseif obj.numSpeciesInlet < 25
                obj.SCarbon_CO_Diesel = sum(obj.S_Carbon_CO_Z(10:end,:));
                obj.SCarbon_CO_C21to50 = 0;
            else
                obj.SCarbon_CO_Diesel = sum(obj.S_Carbon_CO_Z(10:20,:));
                obj.SCarbon_CO_C21to50 = sum(obj.S_Carbon_CO_Z(21:end,:));
            end
            
            % Calculate yields to corresponding selectivity
            obj.YMolar_CO_C1to4 = obj.X_COMean * obj.SMolar_CO_C1to4;
            obj.YCarbon_CO_C1to4 = obj.X_COMean * obj.SCarbon_CO_C1to4;
            obj.YMolar_CO_Diesel = obj.X_COMean * obj.SMolar_CO_Diesel;
            obj.YCarbon_CO_Diesel = obj.X_COMean * obj.SCarbon_CO_Diesel;
            
            % Calculate space time yield *(3600s/h)/(1000l/m^3) in kg/l/h
            if obj.numSpeciesInlet < 15
                obj.STYDiesel = 0;
            elseif obj.numSpeciesInlet < 25
                obj.STYDiesel = sum(obj.mOut(15:25,end)) *3.6 / obj.volumeRct;
            else
                obj.STYDiesel = sum(obj.mOut(15:25,end)) * 3.6 / obj.volumeRct;
            end
            obj.STYHC = sum(obj.mOut(6:obj.numSpeciesInlet,end)) * 3.6 / obj.volumeRct;
            
            % Calculate mass time yield *(3600s/h) in g/g_cat/h
            if obj.numSpeciesInlet < 15
                obj.MTYDiesel = 0;
            elseif obj.numSpeciesInlet < 25
                obj.MTYDiesel = sum(obj.mOut(15:25,end)) * 3600 / obj.massCat;
            else
                obj.MTYDiesel = sum(obj.mOut(15:25,end)) * 3600 / obj.massCat;
            end
            obj.MTYHC = sum(obj.mOut(6:obj.numSpeciesInlet,end)) * 3600 / obj.massCat;
            
            % calculate Weisz Modulus
            obj.WeiszModulus
            % calculate Mears Criterion
            obj.MearsCriterion
            % criterion for radial heat transfer limitation
            obj.radHeat
        end

        %% Plot solution on 2D domain
        function plotres(obj)
            %PLOTRES plots the solution of the PDE
            figure(2)
            subplot(2,1,1)
            pdeplot(obj.modelPDEsystem,'XYData',obj.resultPDE.NodalSolution(:,1))
            caxis([0 obj.nInlet(1)])
            colormap('jet')
            title('Molar flow n_{H2} in mol/s')
            
            subplot(2,1,2)
            pdeplot(obj.modelPDEsystem,'XYData',obj.resultPDE.NodalSolution(:,obj.numEquations))
            caxis([380 600])
            colormap('jet')
            title('Temperature T_{Rct} in K')
        end
        
        %% Interpolierte 1D Loesung plotten
        function plotInterp(obj, varargin)
            % PLOTINTERP plots the interpolated solution values

            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', true);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            
            % Plots
            fig1 = figure('name', 'Independent variables along reactor axis');
            figure(fig1);
            
            %%% Mole flows
            subplot(3,1,1)
            % plot mole flow and norm x axis by division
            plot(obj.xq./obj.lengthRct,obj.nOut([1 2 5 6 15],:),'Linewidth',1.5)
            % set axis limits
            xlim([0 1])
            ylim([0 obj.nInlet(1)])
            title('(a)    Molar flows in mol/s')
            xlabel('Dimensionless reactor length / -')
            % use latex interpreter for math symbols like dot over n
            % \mathsf changes font back to normal font
            ylabel('$\dot{\mathsf{n}}_\mathsf{i} \mathsf{ ~/~ mol~s}^{\mathsf{-1}}$','Interpreter','latex');
            % putting legend outside of plot changes width of subplot tile
            % therefore get subplot tile size aP before creating legend
            set(gca, 'unit', 'pixel');
            aP = get(gca,'Position');
            lP = legend(obj.speciesNames([1 2 5 6 15],1),'Location','northeastoutside');
            % get width of legend box, set unit to pixel first
            set(lP, 'unit', 'pixel');
            lPwidth = lP.Position(3);
            % set subplot tile size back to what it was before
            set(gca, 'Position', [ aP(1) aP(2)*1.05 aP(3) aP(4)])
            % resetting subplot tile size might move legend outside of
            % figure window so that it gets cut off when saving
            % therefore figure window width gets enlarged by legend width
            % get figure window size
            fP = get(gcf,'Position');
            % add legend width to figure window size
            set(gcf,'Position',[fP(1) fP(2) fP(3)+lPwidth fP(4)]);
            % set line width
            set(gca,'Linewidth',1.5)
            
            %%% Temperature of reaction zone
            subplot(3,1,2)
            if obj.TRctSwitch
                plot(obj.xq./obj.lengthRct,obj.TRctOutMean,'Linewidth',1.5)
                legend('T_{Rct}','Location','northeastoutside')
            else
                plot(obj.xq./obj.lengthRct,obj.TConst*ones(1,length(obj.xq)),'Linewidth',1.5)
                text(0.5*obj.lengthRct, obj.TConst-0.5, 'isothermal calculation','FontSize',16);
                legend('T_{Rct}','Location','northeastoutside')
            end
            xlim([0 1])
            title('(b)    Reactor temperature')
            xlabel('Dimensionless reactor length / -')
            ylabel('T_{Rct} / K')
            % get positions of current subplot tile
            set(gca, 'unit', 'pixel');
            aP2 = get(gca, 'Position');
            % set subplot tile width aP2(3) equal to first tile aP(3)
            % But since aspect ratio is conserved, tile height also
            % increases and title overlaps with tile above. therefore
            % decrease y position to 90%
            set(gca, 'Position', [aP2(1) aP2(2) aP(3) aP2(4)]);
            set(gca,'Linewidth',1.5)
            
            %%% Temperature of reactor wall
            subplot(3,1,3)
            if obj.TWallSwitch
                plot(obj.xq./obj.lengthRct,obj.TWallOutMean,'Linewidth',1.5)
                xlim([0 1])
                title('(c)    Wall temperature')
                xlabel('Dimensionless reactor length / -')
                ylabel('T_{Wall} / K')
                legend('T_{Wall}','Location','northeastoutside')
                % get positions of current subplot tile
                set(gca, 'unit', 'pixel');
                aP3 = get(gca, 'Position');
                % set subplot tile width aP2(3) equal to first tile aP(3)
                set(gca, 'Position', [aP3(1) aP3(2)*0.75 aP(3) aP2(4)]);
                set(gca,'Linewidth',1.5)
            end
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_rct_Results.png');
%                 saveas(gcf, filename);
                pause(1);
                export_fig(filename,'-m4');
            end
            
            %% Physical properties
            fig3 = figure('name', 'Physical properties');
            figure(fig3);
            
            %%% Plot 1
            subplot(2,2,1);
            if size(obj.y_iOut,1) < 25
                plot(obj.xq./obj.lengthRct,obj.gasPhase.molarWeight('y', obj.y_iOut(1:end,:),'switchCheckInput',false),'Linewidth',1.5')
            else
                plot(obj.xq./obj.lengthRct,obj.gasPhase.molarWeight('y', obj.y_iOut(1:25,:),'switchCheckInput',false),'Linewidth',1.5')
            end
            xlim([0 1])
            xlabel('Dimensionless reactor length / -')
            ylabel('M_i / kg mol^{-1}')
            title('(a)    Molar weight')
            legend('gasPhase','Location','northwest')
            set(gca,'Linewidth',1.5)
            
            %%% Plot 2
            if size(obj.y_iOut,1) < 25
                gasPhaseDensity = obj.pInlet*100000*obj.gasPhase.molarWeight('y', obj.y_iOut(1:end,:),'switchCheckInput',false)'/(obj.R*obj.TInlet);
            else
                gasPhaseDensity = obj.pInlet*100000*obj.gasPhase.molarWeight('y', obj.y_iOut(1:25,:),'switchCheckInput',false)'/(obj.R*obj.TInlet);
            end
            subplot(2,2,2);
            plot(obj.xq./obj.lengthRct,gasPhaseDensity,'Linewidth',1.5)
            xlim([0 1])
            xlabel('Dimensionless reactor length / -')
            ylabel([sprintf(strrep('\u03c1','\u','\x')),'_{rct} / kg m^{-3}'])
            title('(b)    Density')
            legend('gasPhase','Location','northwest')
            set(gca,'Linewidth',1.5)
            
            %%% Plot 3
            subplot(2,2,3);
            if obj.numSpeciesGasPhase < 25
                plot(obj.xq./obj.lengthRct,obj.gasPhase.specHeatCap('y', obj.y_iOut(1:obj.numSpeciesGasPhase,:),'T',obj.TRctOut,'switchCheckInput',false)','Linewidth',1.5)
            else
                plot(obj.xq./obj.lengthRct,obj.gasPhase.specHeatCap('y', obj.y_iOut(1:25,:),'T',obj.TRctOut,'switchCheckInput',false)','Linewidth',1.5)
            end
            xlim([0 1])
            xlabel('Dimensionless reactor length / -')
            ylabel('c_{p,rct} / J kg^{-1}K^{-1}')
            title('(c)    Specific heat capacity')
            legend('gasPhase')
            set(gca,'Linewidth',1.5)
            
            % Plot 4
            subplot(2,2,4);
            if obj.numSpeciesGasPhase < 25
                disp('Thermal Conductivity Plot does not work with numSpecies < 25');
            else
                plot(obj.xq./obj.lengthRct,obj.gasPhase.thermCond('y', obj.y_iOut(1:25,:),'T',obj.TRctOut)','Linewidth',1.5)
                xlim([0 1])
                xlabel('Dimensionless reactor length / -')
                ylabel('Thermal conductivity / J/(m*K*s)')
                ylabel([sprintf(strrep('\u03bb','\u','\x')),'_{rct} / J m^{-1}s^{-1}K^{-1}'])
                title('(d)    Thermal conductivity')
                legend('gasPhase')
                set(gca,'Linewidth',1.5)
            end
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_rct_physicalData.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end

        end
        
        %% Result presentation method: Hydrocarbon distribution
        function showResultFT(obj, varargin)
            %SHOWRESULTFT will display hydrocarbon distribution
            %   Mass fraction over carbon number
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', true);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            
            % Plot product mass fractions and diesel fraction
            fig2 = figure('name', 'FT');
            figure(fig2);
            plot(obj.w_Hydrocarbons(:,end),'LineWidth',1.5);
            xlabel('n_C / -','FontSize',20);
            ylabel('w_i / -');
            axis([0 30 0 0.1]);
            ax = gca;
            ax.FontSize = 16;
            title('Mass fraction of diesel in all hydrocarbons');
            text(20, 0.035, cat(2, 'w_{C10-C20} = ', num2str(obj.w_Diesel(end))),'FontSize',16);
            set(gca,'Linewidth',1.5)

            if obj.numSpeciesInlet < 15
                disp('No mass balance for diesel fraction calculated.')
            elseif obj.numSpeciesInlet < 26 
                xverts = [(10:10+length(obj.w_Hydrocarbons(11:end,end))); (11:10+(length(obj.w_Hydrocarbons(11:end,end))+1)); (11:10+(length(obj.w_Hydrocarbons(11:end,end))+1)); (10:10+length(obj.w_Hydrocarbons(11:end,end)))];
                yverts = [zeros(1,1+length(obj.w_Hydrocarbons(11:end,end))); zeros(1,1+length(obj.w_Hydrocarbons(11:end,end))); obj.w_Hydrocarbons(10:end,end)'; obj.w_Hydrocarbons(10:end,end)'];
                patch(xverts,yverts, 'b','LineWidth',0.5);
            else
                xverts = [(10:20); (11:21); (11:21); (10:20)];
                yverts = [zeros(1,11); zeros(1,11); obj.w_Hydrocarbons(11:21,end)'; obj.w_Hydrocarbons(10:20,end)'];
                patch(xverts,yverts, 'b','LineWidth',0.5);
            end
                        
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_rct_ResultsFT.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end
            
            disp(obj.X_CO);
            disp(obj.X_H2);
        end

        %% Result presentation method: Reactor geometry
        function plotReactor(obj, varargin)
            %SHOWREACTOR will display the reactor geometry as PDE geometry
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save', false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            saveSwitch = parser.Results.save;
            
            
            % Plot domain, mesh, reactor 3D model
            fig3 = figure('name', 'Reactor geometry');
            figure(fig3)
            
            % plots the geometry with edge and face labels
            j = subplot(3,1,1);
            pdegplot(obj.modelPDEsystem,'EdgeLabels','on','FaceLabels','on')
            yticks([0 obj.wRct])
            title({'(a)    Geometry plot of edges and faces';''})
            % find all graphic objects in geometry plot
            % write only those with the specific tag into variable
            h = findobj('Tag','PDESubLabel');
            g = findobj('Tag','PDEEdgeLabel');
            % set font size
            % h.Property is a comma sperated list
            % deal overwrites the property of each element
            [h.FontSize] = deal(16);
            [g.FontSize] = deal(16);
            set(gca,'Linewidth',1.5)
            
            % move subplot upwards a little bit
            pos = get(j,'Position');
            set(j,'Position',[pos(1) pos(2)+0.06 pos(3) pos(4)]);
            
            % plots the mesh
            h = subplot(3,1,2);
            pdemesh(obj.modelPDEsystem);
            title({'(b)    Mesh plot';''})
            set(gca,'Linewidth',1.5)
            
            % move subplot upwards a little bit
            pos = get(h,'Position');
            set(h,'Position',[pos(1) pos(2)+0.06 pos(3) pos(4)]);
            
            subplot(3,1,3);
            if obj.kModeSwitch == 4
                % Create geometry of cylinders
                gm = multicuboid([obj.dRct, obj.dRct+2*obj.thicknessWall], [obj.wRct, obj.wRct+2*obj.thicknessWall], obj.lengthRct);
                modelVisualization = createpde;
                modelVisualization.Geometry = gm;
                pdegplot(modelVisualization, 'CellLabels','off','FaceAlpha',0.5)
                a = gca;
                a.ZLim = [0 obj.lengthRct];
                a.XLim = [-(obj.dRct/2+obj.thicknessWall) obj.dRct/2+obj.thicknessWall];
                a.YLim = [-(obj.wRct/2+obj.thicknessWall) obj.wRct/2+obj.thicknessWall];
                a.Box = 'on';
                a.BoxStyle = 'full';
                a.Projection = 'perspective';
                a.View = [70 85];
                camzoom(1.5);
            else
                % Create geometry of cylinders
                gm = multicylinder([obj.dRct/2, obj.dRct/2+obj.thicknessWall], obj.lengthRct);
                modelVisualization = createpde;
                modelVisualization.Geometry = gm;
                pdegplot(modelVisualization, 'CellLabels','off','FaceAlpha',0.5)
                a = gca;
                a.ZLim = [0 obj.lengthRct];
                a.XLim = [-(obj.dRct/2+obj.thicknessWall) obj.dRct/2+obj.thicknessWall];
                a.YLim = [-(obj.dRct/2+obj.thicknessWall) obj.dRct/2+obj.thicknessWall];
                a.Box = 'on';
                a.BoxStyle = 'full';
                a.Projection = 'perspective';
                a.View = [0 0];
            end
            title('(c)    3D-Model reactor')
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_rct_Reactor.png');
%                 saveas(gcf, filename);
                export_fig(filename,'-m4');
            end
        end
        
    end % methods end here
    
end % classdef ends here

