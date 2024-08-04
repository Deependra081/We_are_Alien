// Parameters
hcutc = 1973;
mcsq = 940 * 10^6;
r0 = 0.131349;
alpha = 1.44;
D = 0.755501;
rmin = 0.05;
rmax = 1;

// Input: Number of intervals
n = 750
s = (rmax - rmin) / n;

// Compute r and r_dash
for i = 1:n
    r(1, i) = rmin + s * (i - 1);
    r_dash(1, i) = (r(1, i) - r0) / r0;
end

// Compute potential matrix
factor1 = -(hcutc^2) / (2 * mcsq * s^2);
V_matrix = zeros(n, n);
for i = 1:n
    V_matrix(i, i) = D * (exp(-2 * alpha * r_dash(1, i)) - exp(-alpha * r_dash(1, i)));
end

// Plot potential curve
plot(r, diag(V_matrix)/6.5,'b','linewidth',2);
title('Plot of first four Energy eigen values and Probability density for eigen functions');
xlabel('Internuclear separation (r, Angstrom)');
ylabel('Energy (MeV) and scaled wavefunctions');

// Construct kinetic energy matrix
K = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i == j
            K(i, j) = -2;
        elseif i == (j - 1) || i == (j + 1)
            K(i, j) = 1;
        else
            K(i, j) = 0;
        end
    end
end

// Compute total Hamiltonian matrix
T_matrix = factor1 * K;
H_matrix = V_matrix + T_matrix;

// Compute eigenvalues
eval = spec(H_matrix);
[a , b ]= spec ( H_matrix ) ;
// Print lowest vibrational energy
printf("The Lowest vibrational energy: %f MeV", eval(i, 1) / 10^6);

// Plot energy eigenvalues
for i = 1:4
    eval_num = eval(i);
    eval_num_vec = eval_num * ones(1, n);
    plot(r, eval_num_vec,'g');
end

// Compute probability densities
normalisation = sum(a(:, 1) .* conj(a(:, 1)));
psisq1 = (1 / normalisation) * (a(:, 1) .* conj(a(:, 1)));
psisq2 = (1 / normalisation) * (a(:, 2) .* conj(a(:, 2)));
psisq3 = (1 / normalisation) * (a(:, 3) .* conj(a(:, 3)));
psisq4 = (1 / normalisation) * (a(:, 4) .* conj(a(:, 4)));

// Plot scaled wavefunctions
plot(r(1, :), 15 * psisq1 + eval(1), 'r');
plot(r(1, :), 15 * psisq2 + eval(2), 'r');
plot(r(1, :), 15 * psisq3 + eval(3), 'r');
plot(r(1, :), 15 * psisq4 + eval(4), 'r');
legend(['Morse Potential Curve',(['Energy Eigen values','','','']),'Energy density']);
