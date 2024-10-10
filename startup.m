function startup(opts)
%STARTUP proxy-mps 
%
%
if nargin < 1
    opts = [];
end

fmmremex = false;
if isfield(opts, 'fmmremex')
    fmmremex = opts.fmmremex;
end

fmmrecompile = false;
if isfield(opts, 'fmmrecompile')
    fmmrecompile = opts.fmmrecompile;
end

fmmremex = or(fmmremex, fmmrecompile);

testfmm = false;
if isfield(opts, 'testfmm')
    testfmm = opts.testfmm;
end


sommerfeldremex = false;
if isfield(opts, 'sommerfeldremex')
    sommerfeldremex = opts.sommerfeldremex;
end

sommerfeldrecompile = false;
if isfield(opts, 'sommerfeldrecompile')
    sommerfeldrecompile = opts.sommerfeldrecompile;
end

sommerfeldremex = or(sommerfeldremex, sommerfeldrecompile);
addpath ./src/


finufftremex = false;
if isfield(opts, 'finufftremex')
    finufftremex = opts.finufftremex;
end

finufftrecompile = false;
if isfield(opts, 'finufftrecompile')
    finufftrecompile = opts.finufftrecompile;
end

finufftremex = or(finufftremex, finufftrecompile);



% run the chunkie part of the startup
if (exist('./chunkie/startup.m', 'file'))
    cd chunkie;
    startup(opts);
    cd ../
else
    msg = "proxymps STARTUP: warning chunkie not found in usual location\n " + ...
        "Check that the submodule was included. ";
    warning('PROXYMPS:flamnotfoundwarn',msg);
end
% 
% run the fmm3dbie part of the startup
if(exist('./fmm3dbie/matlab/startup.m','file'))
    run ./fmm3dbie/matlab/startup.m;
    cd './fmm3dbie'
    icheck = exist(['fmm3dbie_routs.' mexext], 'file');
    if icheck ~=3 || fmmremex || fmmrecompile
        if ismac || isunix
            [status,cmdout] = system('which gfortran');
            if(~status)
                fprintf('------- proxy-mps startup: building fmm3dbie ----- \n');
                fprintf('fortran compiler found at: %s\n',cmdout);
                iffmm = true;
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    [~, result] = system('uname -m');
                    if strcmpi(strtrim(result), 'x86_64')
                        !cp -f make.inc.macos.gnu make.inc;
                    else
                        !cp -f make.inc.macos_arm64.gnu make.inc;
                    end
                end
                if fmmrecompile
                    !make clean
                end
                !make matlab;  
                testfmm = true; % test if new install

            else    
                msg = "PROXY-MPS STARTUP: unable to find a suitable compiler for fmm3dbie. " + ...
                    "See manual install instructions on github";
                warning('PROXYMPSSTARTUP:compilerwarn',msg)
            end
        end
    end
    cd ../
    
    if testfmm
        cd fmm3dbie/matlab/test
        runtests
        close all;
        cd ../../../;
    end
end


% run the fmm3d part of the startup
if(exist('./fmm3dbie/FMM3D/matlab','dir'))
    addpath './fmm3dbie/FMM3D/matlab';
    cd fmm3dbie/FMM3D/;
    icheck = exist(['fmm3d.' mexext], 'file');
    if icheck ~=3 || fmmremex || fmmrecompile
        if ismac || isunix
            [status,cmdout] = system('which gfortran');
            if(~status)
                fprintf('------- proxy-mps startup: building fmm3d ----- \n');
                fprintf('fortran compiler found at: %s\n',cmdout);
                iffmm = true;
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    [~, result] = system('uname -m');
                    if strcmpi(strtrim(result), 'x86_64')
                        !cp -f make.inc.macos.gnu make.inc;
                    else
                        !cp -f make.inc.macos_arm64.gnu make.inc;
                    end
                end
                if fmmrecompile
                    !make clean
                end
                !make matlab;  
                testfmm = true; % test if new install

            else    
                msg = "PROXY-MPS STARTUP: unable to find a suitable compiler for FMM3D. " + ...
                    "See manual install instructions on github";
                warning('PROXYMPSSTARTUP:compilerwarn',msg)
            end
        end
    end
    cd ../../
    
    if testfmm
        cd fmm3dbie/FMM3D/matlab
        runtests
        cd ../../..;
    end
end


% 
% run the sommerfeld part of the startup
if(exist('./sommerfeld/matlab','dir'))
    addpath ./sommerfeld/matlab;
    cd './sommerfeld'
    icheck = exist(['somm2drouts.' mexext], 'file');
    if icheck ~=3 || sommerfeldremex || sommerfeldrecompile
        if ismac || isunix
            [status,cmdout] = system('which gfortran');
            if(~status)
                fprintf('------- proxy-mps startup: building sommerfeld ----- \n');
                fprintf('fortran compiler found at: %s\n',cmdout);
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    [~, result] = system('uname -m');
                    if strcmpi(strtrim(result), 'x86_64')
                        !cp -f make.inc.macos.gnu make.inc;
                    else
                        !cp -f make.inc.macos_arm64.gnu make.inc;
                    end
                end
                if sommerfeldrecompile
                    !make clean
                end
                !make matlab; 
            else    
                msg = "PROXY-MPS STARTUP: unable to find a suitable compiler for sommerfeld. " + ...
                    "See manual install instructions on github";
                warning('PROXYMPSSTARTUP:compilerwarn',msg)
            end
        end
    end
    cd ../
end


% run the finufft part of the startup
if(exist('./finufft/matlab','dir'))
    addpath ./finufft/matlab;
    cd './finufft'
    icheck = exist(['finufft.' mexext], 'file');
    if icheck ~=3 || finufftremex || finufftrecompile
        if ismac || isunix
            [status,cmdout] = system('which gfortran');
            if(~status)
                fprintf('------- proxy-mps startup: building finufft ----- \n');
                fprintf('fortran compiler found at: %s\n',cmdout);
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    !cp -f make.inc.macosx_clang_matlab make.inc;
                end
                if finufftrecompile
                    !make clean
                end
                !make matlab; 
            else    
                msg = "PROXY-MPS STARTUP: unable to find a suitable compiler for finufft. " + ...
                    "See manual install instructions on github";
                warning('PROXYMPSSTARTUP:compilerwarn',msg)
            end
        end
    end
    cd ../
end


if(exist('./surfacefun','dir'))
  cd surfacefun
  run setup.m
  cd ..
end
