function atomdensities

enrichment = .032;
M_U235 = 235.0439299; % g/mol
M_U238 = 238.0507882; % g/mol
M_O = 15.9994; % g/mol
M_H1 = 1.00782503207; % g/mol

N = .60221415; % cm^2/b/mol

rho_H2O = .72; % g/cm^3
rho_UO2 = 11; % g/cm^3

M_U = (enrichment/M_U235 + (1-enrichment)/M_U238)^(-1);

rho_U235 = enrichment * M_U/(M_U + 2*M_O) * rho_UO2;
rho_U238 = (1-enrichment)* M_U/(M_U + 2*M_O) * rho_UO2;
rho_O = 2*M_O / (M_U + 2*M_O) * rho_UO2;

N_H2O = rho_H2O/(2*M_H1 + M_O) * N
N_U235 = rho_U235 / M_U235 * N
N_U238 = rho_U238 / M_U238 * N
N_O = rho_O / M_O * N

save atomdensities.mat N_H2O N_U235 N_U238 N_O rho_UO2 rho_H2O
