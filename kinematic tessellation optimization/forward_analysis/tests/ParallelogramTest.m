classdef ParallelogramTest < matlab.unittest.TestCase
%PARALLELOGRAMTEST  Regression tests for the anisotropic rotating-units code.
%
%   The acceptance criterion is that the SQUARE limit of the generalised code
%   reproduces the existing regular-case code to machine precision, with no
%   special-casing of the regular geometry anywhere in the implementation.
%
%   Run with:   runtests('forward_analysis/tests')

    properties (Constant)
        ThetaDeg = 0:5:90        % hinge angles for the square-limit sweep
        AlphaDeg = 0:5:45        % panel rotations for the general sweep
    end

    methods (TestClassSetup)
        function addFunctionPath(tc)
            here = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            % PathFixture wants a list of folders, not a genpath string
            folders = strsplit(genpath(fullfile(here,'function')), pathsep);
            folders = folders(~cellfun(@isempty, folders));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(folders));
        end
    end

    methods (Test)

        function squareLimitStretches(tc)
        % lambda1 = lambda2 = the regular-case lambda(theta), exactly.
            for th = tc.ThetaDeg*pi/180
                Sq = square_deployment(2,2,1,th);
                M  = parallelogram_metrics(1,1,pi/2,th/2);
                tc.verifyEqual(M.lambda(1), Sq.lambda, 'AbsTol',1e-15, ...
                    sprintf('lambda1 at theta = %.0f deg', th*180/pi));
                tc.verifyEqual(M.lambda(2), Sq.lambda, 'AbsTol',1e-15, ...
                    sprintf('lambda2 at theta = %.0f deg', th*180/pi));
            end
        end

        function squareLimitPoisson(tc)
        % nu = -1 identically on the whole path, not just at the ends.
            for th = tc.ThetaDeg(2:end)*pi/180
                M = parallelogram_metrics(1,1,pi/2,th/2);
                tc.verifyEqual(M.nu, -1, 'AbsTol',1e-14, ...
                    sprintf('nu at theta = %.0f deg', th*180/pi));
            end
        end

        function squareLimitPorosity(tc)
        % p = 1 - 1/(lambda1*lambda2) must collapse onto 1 - lambda^-2.
            for th = tc.ThetaDeg*pi/180
                Sq = square_deployment(2,2,1,th);
                M  = parallelogram_metrics(1,1,pi/2,th/2);
                tc.verifyEqual(M.p, 1 - Sq.lambda^-2, 'AbsTol',1e-15);
            end
            tc.verifyEqual(parallelogram_metrics(1,1,pi/2,pi/4).p, 0.5, 'AbsTol',1e-15);
        end

        function squareLimitIsotropic(tc)
        % mu = 1 exactly, and it is the ONLY isotropic point of the family.
            tc.verifyEqual(parallelogram_metrics(1,1,pi/2).mu, [1 1], 'AbsTol',1e-14);
            for a = [1 1.05 2]
                for g = [50 90 130]*pi/180
                    if a == 1 && abs(g - pi/2) < 1e-12, continue; end
                    tc.verifyGreaterThan(parallelogram_metrics(a,1,g).mu(1), 1 + 1e-6, ...
                        sprintf('a = %g, gamma = %g deg should be anisotropic', a, g*180/pi));
                end
            end
        end

        function anisotropyIdentities(tc)
        % mu1*mu2 = 1 and mu1+mu2 = (a^2+b^2)/(ab sin gamma), from det/trace of M.
            for a = [1 1.3 2.4]
                for b = [0.7 1 1.9]
                    for g = [35 65 90 120]*pi/180
                        M = parallelogram_metrics(a,b,g);
                        tc.verifyEqual(prod(M.mu), 1, 'AbsTol',1e-13);
                        tc.verifyEqual(sum(M.mu), (a^2+b^2)/(a*b*sin(g)), 'RelTol',1e-13);
                        tc.verifyEqual(det(M.M), 1, 'AbsTol',1e-13);
                    end
                end
            end
        end

        function principalFrameIsFixed(tc)
        % The whole point of F being symmetric: the principal directions do not
        % rotate along the deployment path.
            D0 = parallelogram_metrics(1.5,1,65*pi/180,0).dirs;
            for al = tc.AlphaDeg*pi/180
                D = parallelogram_metrics(1.5,1,65*pi/180,al).dirs;
                tc.verifyEqual(abs(D), abs(D0), 'AbsTol',1e-12, ...
                    sprintf('principal frame moved at alpha = %.0f deg', al*180/pi));
            end
        end

        function deformationGradientMapsLattice(tc)
        % F must actually carry the reference lattice onto the deployed one.
            for al = tc.AlphaDeg*pi/180
                M = parallelogram_metrics(1.5,1,65*pi/180,al);
                tc.verifyEqual(M.F*[M.u M.v], [M.s M.t], 'AbsTol',1e-13, ...
                    sprintf('F[u v] ~= [s t] at alpha = %.0f deg', al*180/pi));
            end
        end

        function panelsStayRigidAndHinged(tc)
        % The built sheet must be a mechanism: panel edges keep their lengths
        % and hinged corners stay coincident.
            L0 = [];
            for al = tc.AlphaDeg*pi/180
                S = parallelogram_deployment(4,4,1.5,1,65*pi/180,al);
                L = [];
                for k = 1:numel(S.panels)
                    p = S.panels{k};   q = p([2:end 1],:);
                    L = [L; sqrt(sum((q-p).^2,2))];  %#ok<AGROW>
                end
                if isempty(L0), L0 = L; end
                tc.verifyEqual(L, L0, 'AbsTol',1e-12, 'panels are not rigid');
                tc.verifyLessThan(S.hinge_gap, 1e-12, 'hinges came apart');
            end
        end

        function noOverlapOnTheAdmissiblePath(tc)
        % Exact polygon boolean: while the motion is admissible the union area
        % equals the sum of the panel areas.
            for al = tc.AlphaDeg*pi/180
                S = parallelogram_deployment(3,3,1.5,1,65*pi/180,al);
                tc.verifyLessThan(abs(S.overlap), 1e-10*S.solid_area, ...
                    sprintf('panels overlap at alpha = %.0f deg', al*180/pi));
            end
        end

        function equilateralLimitMatchesTriangleCode(tc)
        % The triangular family's isotropic member must reproduce the existing
        % equilateral rotating-triangles code exactly.
            for th = (0:5:120)*pi/180
                St = triangle_deployment(2,2,1,th);
                M  = rotating_unit_metrics('triangle',1,1,pi/3,th/2);
                tc.verifyEqual(M.lambda(1), St.lambda, 'AbsTol',1e-14);
                tc.verifyEqual(M.lambda(2), St.lambda, 'AbsTol',1e-14);
                tc.verifyEqual(M.p, 1 - St.lambda^-2, 'AbsTol',1e-14);
                if th > 0
                    tc.verifyEqual(M.nu, -1, 'AbsTol',1e-14, ...
                        sprintf('nu at theta = %.0f deg', th*180/pi));
                end
            end
            tc.verifyEqual(rotating_unit_metrics('triangle',1,1,pi/3).mu, ...
                [sqrt(3) sqrt(3)], 'AbsTol',1e-13);
        end

        function unifiedTraceDetLaw(tc)
        % trace(M) = -sum(edge^2)/(2*area) and det(M) = k, for BOTH families.
            for a = [1 1.3 2.4]
                for b = [0.7 1 1.9]
                    for g = [35 60 90 120]*pi/180
                        Mp = rotating_unit_metrics('parallelogram',a,b,g);
                        tc.verifyEqual(trace(Mp.M), -2*(a^2+b^2)/(2*a*b*sin(g)), 'RelTol',1e-12);
                        tc.verifyEqual(det(Mp.M), 1, 'AbsTol',1e-12);

                        Mt = rotating_unit_metrics('triangle',a,b,g);
                        c2 = a^2 + b^2 - 2*a*b*cos(g);
                        tc.verifyEqual(trace(Mt.M), -(a^2+b^2+c2)/(a*b*sin(g)), 'RelTol',1e-12);
                        tc.verifyEqual(det(Mt.M), 3, 'AbsTol',1e-12);
                    end
                end
            end
        end

        function generalTriangleIsAMechanism(tc)
        % A scalene triangle must still close: rigid panels, hinges together,
        % no overlap, right across the sweep.
            L0 = [];
            for al = (0:10:60)*pi/180
                S = triangle_general_deployment(3,3,1.4,0.8,55*pi/180,al);
                L = [];
                for k = 1:numel(S.panels)
                    p = S.panels{k};   q = p([2:end 1],:);
                    L = [L; sqrt(sum((q-p).^2,2))];  %#ok<AGROW>
                end
                if isempty(L0), L0 = L; end
                tc.verifyEqual(L, L0, 'AbsTol',1e-12, 'panels are not rigid');
                tc.verifyLessThan(S.hinge_gap, 1e-12, 'hinges came apart');
                tc.verifyLessThan(abs(S.overlap), 1e-9*S.solid_area, 'panels overlap');
            end
        end

        function equilateralIsTheOnlyIsotropicTriangle(tc)
        % Weitzenboeck: a^2+b^2+c^2 >= 4*sqrt(3)*Area, equality iff equilateral.
            tc.verifyEqual(rotating_unit_metrics('triangle',1,1,pi/3).mu(1), sqrt(3), ...
                'AbsTol',1e-13);
            for a = [1 1.05 2]
                for g = [40 60 100]*pi/180
                    if a == 1 && abs(g - pi/3) < 1e-12, continue; end
                    tc.verifyGreaterThan(rotating_unit_metrics('triangle',a,1,g).mu(1), ...
                        sqrt(3) + 1e-6);
                end
            end
        end

        function squareLimitSheetMatchesRegularCode(tc)
        % End to end: the generalised sheet builder at a = b, gamma = 90 must
        % agree with square_deployment on lambda and stay a valid mechanism.
            for th = tc.ThetaDeg*pi/180
                Sp = parallelogram_deployment(4,4,1,1,pi/2,th/2);
                Sq = square_deployment(4,4,1,th);
                tc.verifyEqual(Sp.lambda(1), Sq.lambda, 'AbsTol',1e-14);
                tc.verifyLessThan(Sp.hinge_gap, 1e-12);
            end
        end
    end
end
