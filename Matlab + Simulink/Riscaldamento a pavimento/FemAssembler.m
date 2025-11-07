classdef FemAssembler<FemAssemblerAbstract
    %FEM analysis in Cartesian 2D coordinates for the
    %Transport-Diffusion-Reaction problem.
    %Matrices and constant terms vectors are calculated excluding known
    %(Dirichlet) nodes

    %2022 Paolo Bardella, paolo.bardella@polito.it
    methods(Static)
        function [D,T,R]=BuildStiffness(Me, MergeMatrices)
            %BUILDSTIFFNESS creates the stiffness matrix
            %Input:
            %   Me:             a mesh2D object
            %   MergeMatrices:  optional flag, default:true. If true, the
            %                   matrix D contains diffusion, transport and
            %                   reaction terms. If false, the 3 contributions
            %                   are separately returned
            %Output:
            %   D, T, R:        Stiffnes matrices with diffusion, transport and
            %                   reaction terms
            if nargin<1
                error([mfilename ':BuildStiffness:NotEnoughInputs'],'At least 1 input arguments is required');
            end
            if nargin==1
                MergeMatrices=true;
            end
            if ~MergeMatrices && nargout~=3
                error([mfilename ':BuildStiffness:InvalidInput1'],'3 output arguments (D,T,R) are require when Split=true');
            elseif MergeMatrices && nargout~=1
                error([mfilename ':BuildStiffness:InvalidInput2'],'1 output argument (D) is require when Split=false');
            end
            %for clarity, call some properties of Me with shorter names
            %Indices of vertices of each triangle; 3 cols. matrix
            V = Me.Triangles.Vertices;
            Twin=Me.Nodes.TwinNode;
            %Area of each triangle; column vector
            Areas = Me.Triangles.Areas;
            %Numbering of Dirichlet(<0) and unknown (>0) nodes
            Dof = Me.Nodes.Dof;
            %number of internal nodes: we know that the numDof unknown nodes are
            %numbered from 1 to numDof in Me.Nodes.Dof; its maximum is therefore the
            %number of unknown (degrees of freedom)
            numDof = max(Dof);
            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            row = zeros(Me.MatrixContributions, 1,FemAssembler.SparseIndexFormat);
            col = zeros(Me.MatrixContributions, 1,FemAssembler.SparseIndexFormat);
            d = zeros(Me.MatrixContributions, 1);
            t = zeros(Me.MatrixContributions, 1);
            r = zeros(Me.MatrixContributions, 1);
            pos = 1;  %we start from the element in position 1, we'll increase this index
            %evaluate the coefficient in front of the Laplace operator

            mu=Me.mu;
            IncludeDiffusion= nnz(mu)>0;
            %evaluate the coefficient in front of the divergence operator
            %NB. beta is a two-columns matrix, indicating the speed in the x
            %and y directions
            beta=Me.beta;
            rho=Me.rho;
            IncludeTransport=nnz(beta)>0;
            %evaluate the advection/reaction coefficient
            sigma=Me.sigma;
            IncludeReaction= nnz(sigma)>0;
            %dx and dy for all the triangles
            Dx=Me.Triangles.Dx;
            Dy=Me.Triangles.Dy;
            %main loop on each triangle
            for e = 1:size(V, 1)
                Dxe=Dx(e,:);
                Dye=Dy(e,:);
                %for each vertex of this triangle
                for ni = 1:3
                    %look at the "unknown" numbering: if the node is positive, it
                    %corresponds to a degree of freedom of the problem
                    Veni=V(e,ni);
                    ii = Dof(Veni);
                    %is it unknown?
                    if ii > 0
                        %yes it is! second loop on the vertices
                        for nj = 1:3
                            jj = Dof(V(e, nj));
                            %%is it unknown as well?
                            if jj > 0
                                %add the contribution to the stiffness matrix
                                row(pos) =ii;
                                col(pos) = jj;
                                if IncludeTransport
                                    t(pos)=(beta(e,1)*Dye(nj)-beta(e,2)*Dxe(nj))*1/6*rho(e);
                                end
                                if IncludeDiffusion
                                    d(pos)=mu(e)*(Dye(ni)*Dye(nj)+Dxe(ni)*Dxe(nj))/(4.0*Areas(e));
                                end
                                if IncludeReaction
                                    r(pos)=sigma(e)*Areas(e)*(1+(ni==nj))/12;
                                end
                                row(pos)=ii;
                                col(pos)=jj;
                                pos = pos + 1;

                                %Non sparse solution:
                                %D(ii,jj)=D(ii,jj) + mu(e)*(Dy(e,ni)*Dy(e,nj)+Dx(e,ni)*Dx(e,nj))/(4.0*Area) ;
                                %Periodic B.C.
                                gem=Twin(Veni);
                                if gem>0
                                    row(pos)=gem;
                                    if nj~=ni
                                        col(pos)=jj;
                                    else
                                        col(pos)=gem;
                                    end
                                    if IncludeTransport
                                        t(pos)=t(pos-1);
                                    end
                                    if IncludeDiffusion
                                        d(pos)=d(pos-1);
                                    end
                                    if IncludeReaction
                                        r(pos)=r(pos-1);
                                    end
                                    pos=pos+1;
                                end
                                %end periodic B.C.
                            end
                        end
                    end
                end
            end
            %assemble the stiffness matrix D
            pos=pos-1;
            if MergeMatrices
                if length(d)<pos || length(t) <pos || length(r)<pos
                    D = sparse(row(1:pos), col(1:pos), d(1:pos), numDof, numDof)+...
                        sparse(row(1:pos), col(1:pos), r(1:pos), numDof, numDof)+...
                        sparse(row(1:pos), col(1:pos), t(1:pos), numDof, numDof);
                else
                    D = sparse(row(1:pos), col(1:pos), d(1:pos)+t(1:pos)+r(1:pos), numDof, numDof);
                end
            else
                D = sparse(row(1:pos), col(1:pos), d(1:pos), numDof, numDof);
                R = sparse(row(1:pos), col(1:pos), r(1:pos), numDof, numDof);
                T = sparse(row(1:pos), col(1:pos), t(1:pos), numDof, numDof);
            end
        end

        function b=BuildForce(Me,f,bIn)
            %Assemble the vector b of the Transport-Diffusion-Reaction problem
            %Input:
            %   Me     :a Mesh2D object
            %   f      :MATLAB function of (x,y) which returns the values of the
            %           external source.
            %   bIn    :vector to add to the calculatd b vector
            %
            %Output:
            %   b      :constant terms vector
            if nargin<2
                error([mfilename ':BuildForce:NotEnoughInputs'],'Not enough input arguments');
            end
            %for clarity, call some properties of Me with shorter names
            V=Me.Triangles.Vertices;
            Areas=Me.Triangles.Areas;
            CenterOfMass=Me.Triangles.CenterOfMass;
            Dof=Me.Nodes.Dof;

            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);

            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            b = zeros(numDof,1);
            %we evaluate the external force in the center of mass of each triangle
            if isa(f,'function_handle')
                force = f(CenterOfMass.X,CenterOfMass.Y);
            else
                force = f*ones(size(CenterOfMass.X));
            end
            AreaForceDiv3=Areas.*force/3.0;
            %main loop on each triangle
            for e=1:size(V,1)
                %for each vertex of this triangle
                for ni=1:3
                    %look at the "unknown" numbering: if the node is positive, it
                    %corresponds to a degree of freedom of the problem
                    ii = Dof(V(e,ni));
                    %is it unknown?
                    if ii > 0
                        %build the constant terms vector adding the external
                        %contribution
                        b(ii) = b(ii)+ AreaForceDiv3(e);
                    end
                end
            end
            if nargin>2 %add bIn to b
                b=b+bIn;
            end
        end

        function [D,b]=AddStabilization(Me, DIn, bConstIn, StabCoeff)
            %Assemble the matrix D and the vector b of the Transport-Diffusion-Reaction
            %problem with non homogeneous Dirichlet B.C.s with artificial diffusion
            %Input:
            %   Me            :a Mesh2D object
            %   DIn, bIn      :Optional terms to add to D and bConst respectively
            %   StabCoeff     :Non negative number which scales the correction terms
            %Output:
            %   D, bConst     :stabilization terms
            if nargin<4
                StabCoeff=1;
            end
            if StabCoeff<0 
                error([mfilename ':AddStabilization:InvalidStabCoeff'],'StabCoeff must be a non negative number');
            end
            %for clarity, call some properties of Me with shorter names
            V=Me.Triangles.Vertices;
            Areas=Me.Triangles.Areas;
            Dof=Me.Nodes.Dof;
            Twin=Me.Nodes.TwinNode;
            Dx=Me.Triangles.Dx;
            Dy=Me.Triangles.Dy;
            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);
            numDirichletNodes=-min(Dof);
            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            bval = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            brow = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            bpos=1;
            Drow = zeros(Me.MatrixContributions,1);
            Dcol = zeros(Me.MatrixContributions,1);
            Dval=zeros(Me.MatrixContributions,1);
            dpos=1;

            %evaluate the coefficient in front of the Laplace operator
            mu=Me.mu;
            %evaluate the coefficient in front of the divergence operator
            %NB. beta is a two-coefficient vector, indicating the speed in the x
            %and y directions
            Peclet=Me.Triangles.Peclet;

            %main loop on each triangle
            for e=1:size(V,1)
                if Peclet(e)>1
                    muStab=mu(e)*Peclet(e)*StabCoeff;
                    %for each vertex of this triangle
                    for ni=1:3
                        %look at the "unknown" numbering: if the node is positive, it
                        %corresponds to a degree of freedom of the problem
                        Veni=V(e,ni);
                        ii = Dof(Veni);
                        %is it unknown?
                        if ii > 0
                            %yes it is! second loop on the vertices
                            for nj=1:3
                                dtmp=muStab*(Dy(e,ni)*Dy(e,nj)+Dx(e,ni)*Dx(e,nj))/(4.0*Areas(e));
                                jj = Dof(V(e,nj));
                                %%is it unknown as well?
                                if jj > 0
                                    Drow(dpos)=ii;
                                    Dcol(dpos)=jj;
                                    Dval(dpos)=dtmp;
                                    dpos=dpos+1;
                                    gem=Twin(Veni);
                                    if gem>0
                                        Drow(dpos)=gem;
                                        Dcol(dpos)=jj;
                                        Dval(dpos)=dtmp;
                                        dpos=dpos+1;
                                    end
                                else
                                    val=Me.BC.DirichletNodes(-jj,2);
                                    bval(bpos)=-dtmp*val;
                                    brow(bpos)=ii;
                                    bpos=bpos+1;
                                    gem=Twin(Veni);
                                    if gem>0
                                        bval(bpos)=-dtmp*val;
                                        brow(bpos)=gem;
                                        bpos=bpos+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            %assemble the correction matrix DStab
            %We have entries just for the triangles with Peclet >0,
            %so row and col contain zeros at the end if there are triangles with
            %Peclet <=1
            D=sparse(Drow(1:dpos-1), Dcol(1:dpos-1), Dval(1:dpos-1), numDof, numDof);
            b=sparse(brow(1:bpos-1), ones(bpos-1,1), bval(1:bpos-1), numDof, 1);
            if nargin>1 && ~isempty(DIn)
                D=D+DIn;
            end
            if nargin>2 && ~isempty(bConstIn)
                b=b+bConstIn;
            end
        end

        function [M, MVar]=BuildMass(Me)
            %BUILDMASS calculates the mass matrix
            %Input:
            %   Me:     a mesh2D object
            %Output:
            %   M:      mass matrix
            %   MVar:   cell array of vectors containing the correction for
            %           time dependent Dirichlet B.C.s associated to edges with
            %           group ID >=1

            %for clarity, call some properties of Me with shorter names
            V = Me.Triangles.Vertices;
            Dof = Me.Nodes.Dof;
            Areas=Me.Triangles.Areas;
            Twin=Me.Nodes.TwinNode;
            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);
            numGroups=numel(Me.BC.Groups);
            numDirichletNodes=-min(Dof);
            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            row = zeros(Me.MatrixContributions, 1);
            col = zeros(Me.MatrixContributions, 1);
            m = zeros(Me.MatrixContributions, 1);
            pos = 1; %we start from the element in position 1, we'll increase this index
            %everytime we add an entry

            rowvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            colvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            mvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            posvar = 1; %we start from the element in position 1, we'll increase this index
            %everytime we add an entry

            %evaluate the coefficient in front of the time derivative
            rho = Me.rho;
            %main loop on each triangle
            for e = 1:size(V, 1)
                %for each vertex of this triangle
                for ni = 1:3
                    %look at the "unknown" numbering: if the node is positive, it
                    %corresponds to a degree of freedom of the problem
                    Veni=V(e, ni);
                    ii = Dof(Veni);
                    %is it unknown?
                    if ii > 0
                        %yes it is! second loop on the vertices
                        for nj = 1:3
                            jj = Dof(V(e, nj));
                            %%is it unknown as well?
                            entry = 1/12 * Areas(e) * rho(e) * ((ii == jj) + 1);
                            if jj > 0
                                %add the contribution to the mass matrix
                                row(pos) = ii;
                                col(pos) = jj;
                                m(pos) = entry;
                                pos = pos + 1;
                                gem=Twin(Veni);
                                if gem>0
                                    row(pos)=gem;
                                    if nj~=ni
                                        col(pos)=jj;
                                    else
                                        col(pos)=gem;
                                    end
                                    m(pos)=entry;                                    
                                    pos=pos+1;
                                end                                
                            else
                                group=Me.BC.DirichletNodes(-jj,3);
                                if group~=0
                                    mvar(posvar)=entry;
                                    rowvar(posvar)=ii;
                                    colvar(posvar)=group;
                                    posvar = posvar + 1;
                                    gem=Twin(Veni);
                                    if gem>0
                                        mvar(posvar)=entry;
                                        rowvar(posvar)=gem;
                                        colvar(posvar)=group;
                                        posvar = posvar + 1;
                                    end      
                                end
                            end
                        end
                    end
                end
            end
            %assemble the mass matrix M
            pos=pos-1;
            M = sparse(row(1:pos), col(1:pos), m(1:pos), numDof, numDof);
            posvar=posvar-1;
            MVar = sparse(rowvar(1:posvar), colvar(1:posvar), mvar(1:posvar), numDof, numGroups-1);
        end

        function [M, MVar]=BuildMassLumping(Me)
            %BUILDMASSLUMPING calculates the mass matrix using the mass lumping
            %approach
            %Input:
            %   Me:     a mesh2D object
            %Output:
            %   M:      mass matrix
            %   MVar:   cell array of vectors containing the correction for
            %           time dependent Dirichlet B.C.s associated to edges with
            %           group ID >=1

            V = Me.Triangles.Vertices;
            Dof = Me.Nodes.Dof;
            Areas=Me.Triangles.Areas;
            Twin=Me.Nodes.TwinNode;
            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);
            numGroups=numel(Me.BC.Groups);

            %vectors preallocation
            M = zeros(numDof,1);
            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            rowvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            colvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            mvar = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            posvar = 1; %we start from the element in position 1, we'll increase this index

            %evaluate the coefficient in front of the time derivative
            rho = Me.rho;
            %main loop on each triangle
            for e = 1:size(V, 1)
                %for each vertex of this triangle
                for ni = 1:3
                    %look at the "unknown" numbering: if the node is positive, it
                    %corresponds to a degree of freedom of the problem
                    Veni=V(e, ni);
                    ii = Dof(Veni);
                    %is it unknown?
                    if ii > 0
                        %yes it is! second loop on the vertices
                        for nj = 1:3
                            jj = Dof(V(e, nj));
                            %%is it unknown as well?
                            entry = Areas(e) * rho(e) / 3;
                            if jj > 0
                                %add the contribution to the mass matrix
                                M(ii) = M(ii) + entry;
                                gem=Twin(Veni);
                                if gem>0
                                    M(gem) = M(gem) + entry;
                                end   
                            else
                                group=Me.BC.DirichletNodes(-jj,3);
                                if group~=0
                                    mvar(posvar)=entry;
                                    rowvar(posvar)=ii;
                                    colvar(posvar)=group;
                                    posvar = posvar + 1;
                                    gem=Twin(Veni);
                                    if gem>0
                                        mvar(posvar)=entry;
                                        rowvar(posvar)=gem;
                                        colvar(posvar)=group;
                                        posvar = posvar + 1;
                                    end  
                                end
                            end

                        end
                    end
                end
            end
            %assemble the mass matrix MVar
            posvar=posvar-1;
            MVar = sparse(rowvar(1:posvar), colvar(1:posvar), mvar(1:posvar), numDof, numGroups-1);
        end

        function [bConst,bVar]=BuildDirichlet(Me,bConstIn,bVarIn)
            %BUILDDIRICHLET Assembles the vector b with non-homogeneous Dirichlet B.C.s
            %Input:
            %   Me:                 a Mesh2D object
            %   bConstIn, bVarIn:   optional terms to add to bConst and bVar
            %
            %Output:
            %   bConst:             constant terms vector for edges in group 0
            %   bVar:               array of cells with constant terms vector for
            %                       edges in the other groups
            if nargin<1
                error([mfilename ':BuildDirichlet:NotEnoughInputs'],'At least 1 input arguments is required');
            end
            %for clarity, call some properties of Me with shorter names
            V=Me.Triangles.Vertices;
            %Area of each triangle; column vector
            Areas = Me.Triangles.Areas;
            Dof=Me.Nodes.Dof;
            Dx=Me.Triangles.Dx;
            Dy=Me.Triangles.Dy;
            Twin=Me.Nodes.TwinNode;
            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);
            numDirichletNodes=-min(Dof);
            numGroups=numel(Me.BC.Groups);
            %evaluate the coefficient in front of the Laplace operator
            mu=Me.mu;
            %evaluate the coefficient in front of the divergence operator
            %NB. beta is a two-columns matrix, indicating the speed in the x
            %and y directions
            beta=Me.beta;
            rho=Me.rho;
            %evaluate the advection/reaction coefficient
            sigma=Me.sigma;
            row = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1,FemAssembler.SparseIndexFormat);
            col = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1,FemAssembler.SparseIndexFormat);
            b = zeros(numDirichletNodes*FemAssembler.DirichletBCScaling, 1);
            pos = 1;  %we start from the element in position 1, we'll increase this index
            %main loop on each triangle
            for e=1:size(V,1)
                %for each vertex of this triangle
                for ni=1:3
                    Veni=V(e,ni);
                    %look at the "unknown" numbering: if the node is positive, it
                    %corresponds to a degree of freedom of the problem
                    ii = Dof(Veni);
                    %is it unknown?
                    if ii > 0
                        %yes it is! second loop on the vertices
                        for nj=1:3

                            jj = Dof(V(e,nj));
                            %%is it unknown as well?
                            if jj < 0
                                transport=(beta(e,1)*Dy(e,nj)-beta(e,2)*Dx(e,nj))*1/6*rho(e);
                                diffusion=mu(e)*(Dy(e,ni)*Dy(e,nj)+Dx(e,ni)*Dx(e,nj))/(4.0*Areas(e));
                                reaction=sigma(e)*Areas(e)*(1+(ni==nj))/12;
                                val=Me.BC.DirichletNodes(-jj,2);
                                group=Me.BC.DirichletNodes(-jj,3);
                                entry=-(transport+diffusion+reaction)*val;

                                b(pos)=entry;
                                row(pos)=ii;
                                col(pos)=group+1;
                                pos = pos + 1;
                                gem=Twin(Veni);
                                if gem>0
                                    b(pos)=entry;
                                    row(pos)=gem;
                                    col(pos)=group+1;
                                    pos = pos + 1;
                                end


                            end
                        end
                    end
                end
            end

            pos=pos-1;
            b = sparse(row(1:pos), col(1:pos), b(1:pos), numDof, numGroups);
            bConst=b(:,1);
            bVar=b(:,2:end);
            if nargin>1
                bConst=bConst+bConstIn;
            end
            if nargin>2
                bVar=bVar+bVarIn;
            end
        end



        function [bConst,bVar]=BuildNeumann(Me,bConstIn, bVarIn)
            %BUILDNEUMANN Assembles the vector b with non-homogeneous Neumann B.C.s
            %Input:
            %   Me:                 a Mesh2D object
            %   bConstIn, bVarIn:   optional terms to add to bConst and bVar
            %
            %Output:
            %   bConst:             constant terms vector for edges in group 0
            %   bVar:               array of cells with constant terms vector for
            %                       edges in the other groups
            if nargin<1
                error([mfilename ':BuildNeumann:NotEnoughInputs'],'At least 1 input arguments is required');
            end
            %for clarity, call some properties of Me with shorter names
            Dof=Me.Nodes.Dof;
            Edges=Me.Edges;
            Neumann=Me.BC.NeumannEdges;
            Twin=Me.Nodes.TwinNode;

            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numNeumannEdges=size(Neumann,1);
            numDof = max(Dof);
            numGroups=numel(Me.BC.Groups);
            brow = zeros(numNeumannEdges*2, 1,FemAssembler.SparseIndexFormat);
            bcol = zeros(numNeumannEdges*2, 1,FemAssembler.SparseIndexFormat);
            bval = zeros(numNeumannEdges*2, 1);
            bpos = 1;  %we start from the element in position 1, we'll increase this index

            for k=1:numNeumannEdges
                g=Me.BC.NeumannEdges(k,2);
                if g==0
                    continue;
                end
                L=Me.BC.NeumannEdges(k,3);
                group=Me.BC.NeumannEdges(k,4);
                for n=1:2
                    Node=Edges(Neumann(k),n);
                    DofNode=Dof(Node);
                    if DofNode>0
                        entry=g/2*L;
                        bval(bpos)=entry;
                        brow(bpos)=DofNode;
                        bcol(bpos)=group+1;
                        bpos = bpos + 1;
                        gem=Twin(Node);
                        if gem>0
                            bval(bpos)=entry;
                            brow(bpos)=gem;
                            bcol(bpos)=group+1;
                            bpos = bpos + 1;
                        end
                    end
                end
            end
            %assemble the stiffness matrix D from the
            bpos=bpos-1;
            b = sparse(brow(1:bpos), bcol(1:bpos), bval(1:bpos), numDof, numGroups);
            bConst=b(:,1);
            bVar=b(:,2:end);
            if nargin>1
                bConst=bConst+bConstIn;
            end
            if nargin>2
                bVar=bVar+bVarIn;
            end
        end

        function [D,bConst,bVar]=BuildRobin(Me, DIn, bConstIn, bVarIn)
            %BUILDROBIN Assembles the vector b with Robin B.C.s
            %Input:
            %   Me:                 a Mesh2D object
            %   DIn,bConstIn, bVarIn:   optional terms to add to D,bConst,bVar
            %
            %Output:
            %   D:                  effect of the Robin B.C. con the stiffness matrix
            %   bConst:             constant terms vector for edges in group 0
            %   bVar:               array of cells with constant terms vector for
            %                       edges in the other groups
            if nargin<1
                error([mfilename ':BuildRobin:NotEnoughInputs'],'At least 1 input arguments is required');
            end
            if nargout<2
                error([mfilename ':BuildRobin:InvalidOutputParametersNumber'],'When calling BuildRobin you must specify at least 2 OUTPUT parameters');
            end
            Dof=Me.Nodes.Dof;
            Edges=Me.Edges;
            Robin=Me.BC.RobinEdges;
            Twin=Me.Nodes.TwinNode;
            %number of internal nodes: we know that the N unknown nodes are numbered from
            %1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
            %(degrees of freedom)
            numDof = max(Dof);
            numGroups=numel(Me.BC.Groups);
            numRobinEdges=size(Robin,1);
            %vectors preallocation: instead of allocating the (sparse) diffusion matrix,
            %we save the rows, columns and values corresponding to each contribution;
            %at the end, we'll call sparse(...) to obtain the diffusion matrix
            %NOTE: the +10 and +20 are reasonable numbers for periodic BC.
            brow = zeros(numRobinEdges*2+10,1);
            bcol = zeros(numRobinEdges*2+10,1);
            bval = zeros(numRobinEdges*2+10,1);
            bpos=1;  %we start from the element in position 1, we'll increase this index
            %everytime we add an entry
            drow = zeros(numRobinEdges*4+20,1);
            dcol = zeros(numRobinEdges*4+20,1);
            dval = zeros(numRobinEdges*4+20,1);
            dpos=1;  %we start from the element in position 1, we'll increase this index
            %everytime we add an entry


            for k=1:numRobinEdges
                Node1=Edges(Robin(k,1),1);
                Node2=Edges(Robin(k,1),2);
                ii1=Dof(Node1);
                ii2=Dof(Node2);
                h=Robin(k,2);
                g=Robin(k,3);
                L=Robin(k,4);
                groupEdge=Robin(k,5);
                if ii1>0 && ii2<0 %ii1 is unknown, ii2 is known
                    DirichletValue=Me.BC.DirichletNodes(-ii2,2);                    
                    group=Me.BC.DirichletNodes(-ii2,3);
                    entry=-h*L/6*DirichletValue;
                    bval(bpos)=entry;
                    brow(bpos)=ii1;
                    bcol(bpos)=group+1;
                    bpos = bpos + 1;
                    entry=g/2*L;
                    bval(bpos)=entry;
                    brow(bpos)=ii1;
                    bcol(bpos)=groupEdge+1;
                    bpos = bpos + 1;
                    drow(dpos)=ii1;
                    dcol(dpos)=ii1;
                    dval(dpos)=h*L/3;
                    dpos=dpos+1;
                    %D(ii1,ii1)=D(ii1,ii1)+h*L/3;
                elseif ii1<0 && ii2>0 %ii1 is known, ii2 is unknown
                    DirichletValue=Me.BC.DirichletNodes(-ii1,2);
                    entry=-h*L/6*DirichletValue;
                    group=Me.BC.DirichletNodes(-ii1,3);

                    bval(bpos)=entry;
                    brow(bpos)=ii2;
                    bcol(bpos)=group+1;
                    bpos = bpos + 1;
                    entry=g/2*L;
                    bval(bpos)=entry;
                    brow(bpos)=ii2;
                    bcol(bpos)=groupEdge+1;
                    bpos = bpos + 1;


                    drow(dpos)=ii2;
                    dcol(dpos)=ii2;
                    dval(dpos)=h*L/3;
                    dpos=dpos+1;
                    %D(ii2,ii2)=D(ii2,ii2)+h*L/3;
                else  %both are unknwon
                    entry=g/2*L;
                    bval(bpos)=entry;
                    brow(bpos)=ii1;
                    bcol(bpos)=groupEdge+1;
                    bpos=bpos+1;
                    bval(bpos)=entry;
                    brow(bpos)=ii2;
                    bcol(bpos)=groupEdge+1;
                    bpos=bpos+1;
                    drow(dpos:dpos+3)=[ii1;ii2;ii1;ii2];
                    dcol(dpos:dpos+3)=[ii1;ii2;ii2;ii1];
                    dval(dpos:dpos+3)=[2;2;1;1]*h*L/6;
                    dpos=dpos+4;
                    gem1=Twin(Node1);
                    if gem1>0
                        bval(bpos)=entry;
                        brow(bpos)=gem1;
                        bcol(bpos)=groupEdge+1;
                        bpos=bpos+1;
                        drow(dpos:dpos+1)=[gem1,gem1];
                        dcol(dpos:dpos+1)=[gem1,ii2];
                        dval(dpos:dpos+1)=[2,1]*h*L/6;
                        dpos=dpos+2;
                    end
                    gem2=Twin(Node2);
                    if gem2>0
                        bval(bpos)=entry;
                        brow(bpos)=gem2;
                        bcol(bpos)=groupEdge+1;
                        bpos=bpos+1;
                        drow(dpos:dpos+1)=[gem2,gem2];
                        dcol(dpos:dpos+1)=[gem2,ii1];
                        dval(dpos:dpos+1)=[2,1]*h*L/6;
                        dpos=dpos+2;
                    end
  
                    %D(ii1,ii1)=D(ii1,ii1)+h*L/3;
                    %D(ii2,ii2)=D(ii2,ii2)+h*L/3;
                    %D(ii1,ii2)=D(ii1,ii2)+h*L/6;
                    %D(ii2,ii1)=D(ii2,ii1)+h*L/6;
                end
            end
            bpos=bpos-1;
            dpos=dpos-1;
            %assemble the stiffness matrix D from the
            b = sparse(brow(1:bpos), bcol(1:bpos), bval(1:bpos), numDof, numGroups);
            bConst=b(:,1);
            bVar=b(:,2:end);
            D = sparse(drow(1:dpos), dcol(1:dpos), dval(1:dpos), numDof, numDof);
            if nargin>1
                D=D+DIn;
            end
            if nargin>2
                bConst=bConst+bConstIn;
            end
            if nargin>3
                bVar=bVar+bVarIn;
            end

        end
    end
    %%%%%%%%%%% static parameters used to fine tune the memory allocation
    methods (Access = private, Static = true)
        function val = DirichletBCScaling(newval)
            persistent DirichletBCScalingValue;
            if nargin >= 1
                DirichletBCScalingValue = newval;
            end
            if isempty(DirichletBCScalingValue)
                val=10;
            else
                val = DirichletBCScalingValue;
            end
        end

    end
    methods(Static)
        function setDirichletBCScaling(newval)
            FemAssembler.DirichletBCScaling(newval);
        end
        function v = getDirichletBCScaling()
            v = FemAssembler.DirichletBCScaling;
        end
    end

end




