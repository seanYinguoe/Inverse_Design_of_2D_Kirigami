function tests = test_geometry
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
setup_project;
end
function testPanelEdgesStayRigid(testCase)
C = tessellation_deployment(4,4,1,0);
for a = [pi/12, pi/4]
    D = tessellation_deployment(4,4,1,a);
    for k = 1:numel(C)
        for first = 1:4:16
            ix = first:first+3;
            c = C{k}(ix,:); d = D{k}(ix,:);
            testCase.verifyEqual(vecnorm(diff([d;d(1,:)]),2,2), ...
                vecnorm(diff([c;c(1,:)]),2,2),'AbsTol',1e-12);
        end
    end
end
end
function testExpansionMatchesRotatingSquaresLaw(testCase)
C = tessellation_deployment(4,4,1,0); c = vertcat(C{:});
for a = [0, pi/12, pi/4]
    D = tessellation_deployment(4,4,1,a); d = vertcat(D{:});
    testCase.verifyEqual(max(d)-min(d), (cos(a)+sin(a))*(max(c)-min(c)), 'AbsTol',1e-12);
end
end
function testNodesRoundTrip(testCase)
C = tessellation_deployment(2,4,1,pi/9);
testCase.verifyEqual(nodes_to_units(units_to_nodes(C),2,4),C);
end
