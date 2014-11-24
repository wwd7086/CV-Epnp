%Estimating rotations and translations
%Alyssa Novelia (20130501)

%Implementation of the four methods used to estimate rotations and
%translations using tracking data from four material points on the object.

%Details on the method can be found at http://rotations.berkeley.edu/?page_id=1379
%Worked example for noiseless case can be found at
%http://rotations.berkeley.edu/?page_id=1524.

clear; clc; close all;

%%  Preliminaries
%In this section, the equation and figure numbers pertain to the worked
%example webpage (http://rotations.berkeley.edu/?page_id=1524)

%Define book's dimension (equation 2)
l = 8;  %length
w = 6;  %width
b = 1;  %thickness

%Choose 4 material points (landmarks) in reference configuration to track
%(figure 1)
r1 = 0.5*[-l; w; -b];
r2 = 0.5*[l; w; -b];
r3 = 0.5*[-l; -w; -b];
r4 = 0.5*[-l; w; b];

%Parameterize rotation using a set of 3-2-1 Euler angles: L1, L2, L3 are
%matrix components of rotation tensors
L3 = [cos(pi/4) sin(pi/4) 0; -sin(pi/4) cos(pi/4) 0; 0 0 1];
L2 = [cos(pi/6) 0 -sin(pi/6); 0 1 0; sin(pi/6) 0 cos(pi/6)];
L1 = [1 0 0; 0 cos(pi/4) sin(pi/4); 0 -sin(pi/4) cos(pi/4)];

%Take the product of the matrices and transpose it to obtain the rotation matrix 
Q = (L1*L2*L3)';

%Define translation: see equation 4 on example webpage
d = [1; 1; -10];

%Calculate the positions of landmarks in the present configuration: 
p1 = Q*r1 + d;
p2 = Q*r2 + d;
p3 = Q*r3 + d;
p4 = Q*r4 + d;

%Noise will be modeled using a random number generator within the interval
%defined below.
noise_mag = 0;  %set to 0 for noiseless case
noise_int = noise_mag*[-1 1];

%In this example we will only add noise to the present configuration's
%measurements. The measured positions in the present configuration are as
%follows
m1 = p1 + noise_int(1) + (noise_int(2)-noise_int(1)).*rand(3,1);
m2 = p2 + noise_int(1) + (noise_int(2)-noise_int(1)).*rand(3,1);
m3 = p3 + noise_int(1) + (noise_int(2)-noise_int(1)).*rand(3,1);
m4 = p4 + noise_int(1) + (noise_int(2)-noise_int(1)).*rand(3,1);

f = figure; 
set(f, 'Name', 'Visualizing estimated rotation error. Blue: estimated rotation, yellow: actual rotation');

%%  Method 1: Naive method
%From this section onwards, the equation and figure numbers pertain to the
%details webpage (http://rotations.berkeley.edu/?page_id=1379).

%Create augmented M and R matrices (denoted as m' and r' in equation 8).
augM = [m1 m2 m3 m4]; augM = [ones(1,4); augM];
augR = [r1 r2 r3 r4]; augR = [ones(1,4); augR];

%Solve for the augmented solution (equation 9)
augX = augM*augR^(-1);

%Extract rotation matrix and translation vector
R_method1 = augX(2:end, 2:end);
d_method1 = augX(2:end,1);

%Plot comparison of results from the 4 methods here
s1 = subplot(2,2,1); title('Naive method');
set(s1, 'CameraPosition', [14 8 75]);
drawbook(Q,d,l,w,b,'y');
drawbook(R_method1,d_method1,l,w,b,'b');
xlabel('E_1'); ylabel('E_2'); zlabel('E_3');
axis equal;

%%  Method 2: Triad method
%There are 24 possible ways to construct vectors V_i and v_i in equation 12.
%One way to see this is to imagine a cuboid where 4 of its vertices are
%r/m_1,2,3,4. We can calculate the remainder of the vertices of the cuboid
%from r/m_1,2,3,4.

%The vertex on each face of the cuboid is a possible origin for the triad
%formed by {V_i} or {v_i}. The triad will have two vectors {V_1,2 or v_1,2}
%parallel to that particular face of the cuboid and one vector {V_3 or v_3}
%that is normal to the face of the cuboid. The vectors parallel to the face
%of the cuboid will be based on the two edges that meet at that particular
%vertex.

%Since there are 6 faces of the cube, and for each face there are 4
%vertices, there are 6x4 = 24 possible triads to be constructed.

%Construct a cell array to store the vertices on each face. There are 6
%faces and 4 vertices on each face, hence we have a 6x4 cell array.
Rvertices = {r1 r2 r2+r4-r1 r4;
            r2 r2+r3-r1 r2+r3+r4-2*r1 r2+r4-r1;
            r2+r3-r1 r3 r3+r4-r1 r2+r3+r4-2*r1;
            r1 r3 r3+r4-r1 r4;
            r4 r2+r4-r1 r2+r3+r4-2*r1 r3+r4-r1;
            r1 r2 r2+r3-r1 r3};

Mvertices = {m1 m2 m2+m4-m1 m4;
            m2 m2+m3-m1 m2+m3+m4-2*m1 m2+m4-m1;
            m2+m3-m1 m3 m3+m4-m1 m2+m3+m4-2*m1;
            m1 m3 m3+m4-m1 m4;
            m4 m2+m4-m1 m2+m3+m4-2*m1 m3+m4-m1;
            m1 m2 m2+m3-m1 m3};

%Append all the position vectors together
R = [r1 r2 r3 r4];
M = [m1 m2 m3 m4];

%Prepare variables to store error, calculated R matrices and d vectors
err_list = [];
R_list = zeros(3,3,size(Rvertices,1)*size(Rvertices,2));
d_list = [];

for iii = 1:size(Rvertices,1)   %loop over faces
    for jjj = 1:size(Rvertices,2)   %loop over vertices

        %Choose an origin
        pivot = jjj;
        %Choose two other vertices that will define the two edges on the
        %plane
        next = jjj+1; if next > 4, next = next - 4; end
        prev = jjj-1; if prev < 1, prev = prev + 4; end

        %Define r/m_1,2,3 used to calculate V/v_1,2,3 (equation 12 and 13)
        Rpivot = Rvertices{iii,pivot};
        Rnext = Rvertices{iii,next};
        Rprev = Rvertices{iii,prev};

        Mpivot = Mvertices{iii,pivot};
        Mnext = Mvertices{iii,next};
        Mprev = Mvertices{iii,prev};

        %Calculate V_1,2,3 and v_1,2,3 (equation 12 and 13)
        V1 = (Rnext-Rpivot); V1 = V1/norm(V1);
        V2 = cross(V1,(Rprev-Rpivot)); V2 = V2/norm(V2);
        V3 = cross(V1,V2); V3 = V3/norm(V3);

        v1 = (Mnext-Mpivot); v1 = v1/norm(v1);
        v2 = cross(v1,(Mprev-Mpivot)); v2 = v2/norm(v2);
        v3 = cross(v1,v2); v3 = v3/norm(v3);

        %Append V_1,2,3 and v_1,2,3 (equation 14)
        V = [V1 V2 V3];
        v = [v1 v2 v3];

        %Calculate rotation matrix (equation 15)
        R2 = v*V';
        
        %Calculate translation (equation 16 and 17)
        d2 = (m1 + m2 + m3)./3 - R2*((r1 + r2 + r3)./3);

        %Calculate error (equation 18)
        err_square = norm(M - R2*R + [d2 d2 d2 d2]);

        %Append errors, calculated R and d for each iteration
        err_list = [err_list; err_square];
        d_list = [d_list d2];
        R_list(:,:,size(d_list,2)) = R2;
    end
end

%Find the rotation and translation with minimum error and assign the value
[dummy idx] = min(err_list);
R_method2 = R_list(:,:,idx);
d_method2 = d_list(:,idx);

%Plot comparison of results from 4 methods here
s2 = subplot(2,2,2); title('TRIAD method');
set(s2, 'CameraPosition', [14 8 75]);
drawbook(Q,d,l,w,b,'y');
drawbook(R_method2,d_method2,l,w,b,'b');
xlabel('E_1'); ylabel('E_2'); zlabel('E_3');
axis equal;


%%  Method 3: Singular value decomposition method
%Sum the appended matrices R and M across the columns, and divide by the
%column size of the matrices to calculate the mean positions (equation 19)
meanR = sum(R,2)./size(R,2);   
meanM = sum(M,2)./size(M,2);

%Calculate and append position vectors relative to mean positions (equation
%21)
C = M - repmat(meanM,1,4);
D = R - repmat(meanR,1,4);

%Calculate S (equation 20)
S = C*D';

%Note that the third matrix output by the SVD function is the transpose
%of R2 (equation 22)
[R1 LAM R2T] = svd(S);
R2 = R2T';  %calculate R2

%Calculate the rotation matrix and translation vector (equation 23 and 24)
R_method3 = R1*[1 0 0; 0 1 0; 0 0 det(R1*R2)]*R2;
d_method3 = meanM - R_method3*meanR;

%Plot comparison of results from the 4 methods here
s3 = subplot(2,2,3); title('SVD method');
set(s3, 'CameraPosition', [14 8 75]);
drawbook(Q,d,l,w,b,'y');
drawbook(R_method3,d_method3,l,w,b,'b');
xlabel('E_1'); ylabel('E_2'); zlabel('E_3');
axis equal;

%%  Method 4: q-method
%Assign weights to each position vector and calculate weighted position
%vectors- if all position vectors have equal weight, we expect the q-method
%to yield the same result as the SVD method (equation 25)
alpha = [1 1 1 1];
weightR = alpha*R'; weightR = weightR'./sum(alpha);
weightM = alpha*M'; weightM = weightM'./sum(alpha);

%Define B matrix (equation 32)
B = zeros(3,3);
for iii = 1:size(B,1)   %loop over rows, iii = i
    for mmm = 1:size(B,1)   %loop over columns, mmm = m
        B1 = zeros(size(alpha));
        %Calculate the term in bracket, call that B1
        for kkk = 1:size(alpha,2)   %to over alpha to calculate summation, kkk = k
            B1(kkk) = alpha(kkk)*M(iii,kkk)*R(mmm,kkk);
        end
        B1 = sum(B1,2)./sum(alpha);
        B(iii,mmm) = B1 - weightM(iii)*weightR(mmm);
    end
end

%Define z matrix (equation 31)
z2 = zeros(3,size(alpha,2));
for kkk = 1:size(alpha,2)   %loop over alpha to calculate summation, kkk = k
    z2(:,kkk) = cross(alpha(kkk)*M(:,kkk), R(:,kkk));
end
z2 = sum(z2,2)./sum(alpha);
z = cross(weightM,weightR) - z2;

%Define augmented matrix K (equation 30)
K = [B + B' - trace(B)*eye(size(B)) z];
K = [K; z' trace(B)];

%Find maximum eigenvalue of K, and corresponding eigenvector (equation 34)
[dum idx] = max(eig(K));
[eigV lam] = eig(K);
q = eigV(:,idx);    %where q = [e1 e2 e3 e0]

%Reconstruct rotation matrix using Euler-Rodriguez parameters (equation 28)
R_method4 = (q(4)^2 - q(1)^2 - q(2)^2 - q(3)^2)*eye(3,3) + ...
            2*q(1:3)*q(1:3)' + ...
            2*[0 -q(4)*q(3) q(4)*q(2); q(4)*q(3) 0 -q(4)*q(1); -q(4)*q(2) q(4)*q(1) 0];

%Calculate translation (equation 34)
d_method4 = weightM - R_method4*weightR;    

%Plot comparison of results from the 4 methods here
s4 = subplot(2,2,4); title('q-method');
set(s4, 'CameraPosition', [14 8 75]);
drawbook(Q,d,l,w,b,'y');
drawbook(R_method4,d_method4,l,w,b,'b');
xlabel('E_1'); ylabel('E_2'); zlabel('E_3');
axis equal;