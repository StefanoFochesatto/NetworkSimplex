function [graphHist ,bf , xf] = NetworkSimplex(s, e, c, b, x0, maxiter)
% This implementation of network simplex solves the 
% basic minimum cost network flow problem with only positivity contraints. 
% We take in a directed arc set (s, e) where the arc starts at a vertex from s and ends at the corresponding vertex in e.
% We also take in the and the associated cost of each vertex c. Finally we take in 
% x0, an initial basic feasible solution and a vector b, a logical vector
% for identifying the arcs of the feasible spanning tree. 

% First we need to form the correct incidence matrix A for our problem.
G = digraph(s, e, c);
A = -1*(incidence(G));
[m, n] = size(A);

% Set inital BFS
x = x0;

% Pull column indices for B and N matrix
nullIndex = find(~b); 
baseIndex = find(b);
graphHist = {'Basic Spanning Tree', 'Proposed Flow'};



    for i = 1:maxiter
        % Compute the simplex multipliers
        y = NetworkSimpMulti(A, b, c, m);

        % Compute reduced costs
        cHat = reducedCosts(A, b, c, y);


        % Optimiality Check (moved this into here so I can update the graphHist). 
        if cHat >= 0
            xf = x;
            bf = b;
            break 
        end


        % This section of code runs through how the exit variable is chosen
        % First we select the arc which corresponds to the smallest reduced
        % cost and add it to our spanning tree. This must produce a unique
        % cycle. Set the direction (clockwise or counterclockwise about the cycle) 
        % of this new arc (our entering variable) to positive and let arcs
        % in the other direction in the cycle be negative. If we reduce the
        % negative arcs in the cycle by one unit of flow, to continue being a solution
        % we need to increase the flow accross the entering arc. To stay feasible we
        % continue this process until a negative arc has zero flow, and
        % therefore becomes the exiting variable. Finally we update the
        % basic and null indices to reflect our new iteration. 

        [mincost, enterIndex] = min(cHat); % enterIndex is index of edge in N
        enterArc = nullIndex(enterIndex);% Using nullIndex to pull index of edge from A

        % Generating Digraph for proposed flow
        proposedFlow = digraph(s([baseIndex, enterArc]),e([baseIndex, enterArc]), c(s([baseIndex, enterArc])));
   
        graphHist(i,:) = {digraph(s(baseIndex),e(baseIndex), c(s(baseIndex))), proposedFlow};
        
    

        % Pulling enter and exit vertex, for direction purposes
        startVertex = find(~(A(:, enterArc) - 1));
        endVertex = find(~(A(:, enterArc) + 1));

      


        % Find cycle in proposedFlow
        [targetCycle,targetCycleEdges] = allcycles(graph(proposedFlow.adjacency + proposedFlow.adjacency'));
      
        % Checking that cycle is found
        if length(targetCycle)~= 1
            error('cycle not found or solution/ iteration was not BF')
        end
        
        % Pulling -1* cycle incidence matrix. 
        Cycle = -1*proposedFlow.incidence;
        Cycle = Cycle(:, targetCycleEdges{1});

        % Pull arc index from full incidence matrix A
        % I have to do this whole thing because the Cycle()
        % reorders the columns. 
        arcIndex = [];
        for j = 1:length(targetCycle{1})
            arcIndex = [arcIndex, find(all(Cycle(:,j) - A == 0))];
        end



        % Determine direction of arcs with respect to proposed additional arc. 
        % This is a mess... maybe skip to the next bit. 
        arcDirection = 1;
        ArcIndex = [find(all( A(:, enterArc) - Cycle == 0))];
        currentVertex = endVertex;
        temp = Cycle(currentVertex, :);
        temp(ArcIndex(1)) = 0;
        ArcIndex(2) = find(temp);
        if  Cycle(currentVertex,ArcIndex(1)) == Cycle(currentVertex,ArcIndex(2))
            arcDirection = [arcDirection, -1*arcDirection(1)];
            else
                arcDirection = [arcDirection, arcDirection(1)];
        end
        CurrentArc = Cycle(:,ArcIndex(2));
        CurrentArc(currentVertex) = 0;
        currentVertex = find(CurrentArc);

        
        for j = 3:length(targetCycle{1})
            temp = Cycle(currentVertex, :);
            temp(ArcIndex(j - 1)) = 0;
            ArcIndex(j) = find(temp);
            if  Cycle(currentVertex,ArcIndex(j - 1)) == Cycle(currentVertex,ArcIndex(j))
                arcDirection = [arcDirection, -1*arcDirection(j - 1)];
            else
                arcDirection = [arcDirection, arcDirection(j - 1)];
        end
        CurrentArc = Cycle(:,ArcIndex(j));
        CurrentArc(currentVertex) = 0;
        currentVertex = find(CurrentArc);

        end


        % Pull arc that has smallest flow in the opposite
        % direction of the entry arc. 
        [minExitFlow, minExitIndex] = min(x(arcIndex(find(arcDirection(ArcIndex) - 1))));
        exitArc = arcIndex(minExitIndex);


        % Apply minimum flow to entry arc, subtract minimum flow from all
        % arcs of opposite direction.
        x(enterArc) = minExitFlow;
        x(arcIndex(find(arcDirection(ArcIndex) - 1))) = x(arcIndex(find(arcDirection(ArcIndex) - 1))) - minExitFlow;

        % Add proposed arc to the basic set and remove the minimum flow arc
        % that points in the opposite direction. 
        b(enterArc) = 1;
        b(exitArc) = 0;

        % Update null matrix, base matrix, and current iterate. 
        nullIndex = find(~b);
        baseIndex = find(b);
        xf = x;
        bf = b;

    end

end





% Awful Spanning Tree Traversal Algorithm
function y = NetworkSimpMulti(A, b, c, m)
    y = zeros(m, 1);
    TraversalHistory = ones(m, 1);
    treeIndex = find(b);
    B = A(:, treeIndex); % Form incidence matrix for x

    i = 1;
    TraversalHistory(1) = 0;
   while any(TraversalHistory, 'all')
          incidentV = B(i, :);
          costIndex = treeIndex(~~incidentV);
          ci = c(costIndex);
          Q = B(:,~~incidentV);
          direction = Q(i,:); 
          Q(i,:) = zeros(size(Q, 2),1)';

          [multiplierIndex, incidentEdgeIndex] = find(Q);
          
          y(multiplierIndex) = y(i) - direction.*ci;
          possibleVertices = intersect([1, find(y)'], find(TraversalHistory));
          i = possibleVertices(1);
          TraversalHistory(i) = 0;
   end
   
end





% Computing reduced cost. 
function cHat = reducedCosts(A, b, c, y)
    nullIndex = find(~b);
    N = A(:,nullIndex);
    n = size(N, 2);
    cHat = zeros(n, 1);
    for i = 1:n
        cHat(i) = c(nullIndex(i)) - y'*N(:,i); 
    end
end

% Example input
% s = [1,1,2,2,3,3,4,5,5,5,6,7]
% e = [3,4,4,5,4,6,6,4,6,7,8,8]
% c = [5,8,12,2,4,6,7,8,5,7,8,3]
% b = [1,0,0,1,1,0,1,0,1,0,1,1]
% x = [10,0,0,15,10,0,10,0,15,0,25,0]


