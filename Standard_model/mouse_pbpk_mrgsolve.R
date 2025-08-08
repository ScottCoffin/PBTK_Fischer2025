library(mrgsolve)
library(dplyr)

# Allow PFAS to be selected at runtime (default PFBS)
args <- commandArgs(trailingOnly = TRUE)
pfas <- if (length(args) > 0) args[1] else "PFBS"

# Locate parameter table relative to script location
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_dir <- if (length(script_path) > 0) dirname(script_path) else "."
params <- read.csv(file.path(script_dir, "pfas_parameters.csv"))
pfas_params <- subset(params, PFAS == pfas)

# mrgsolve model code translated from generic PBPK structure
MousePBPK.code <- '
$PARAM @annotated
BW   : 0.03 : kg Body weight
QCC  : 13.5 : L/h/kg^0.75 Cardiac output
QLC  : 0.25 : Fraction blood flow to liver
QKC  : 0.21 : Fraction blood flow to kidney
P_app        : 1e-6 : Permeability (cm/s)
K_FABP_exp   : 0 : log10 FABP binding
K_SA_exp     : 0 : log10 serum albumin binding
K_Glob_exp   : 0 : log10 globulin binding
K_SP_exp     : 0 : log10 structural protein binding
K_ML_exp     : 0 : log10 membrane lipid binding

$CMT AST ASI APlasma ALiver AKidney ARest Afeces Aurine

$MAIN
// Derived flows (L/h)
double QC = QCC*pow(BW,0.75);
double QL = QLC*QC;
double QK = QKC*QC;
double QR = QC-QL-QK;
// Simple absorption rate constant (h^-1)
double Kabs = P_app * 3600; // convert cm/s -> cm/h

$ODE
dxdt_AST    = -Kabs*AST;
dxdt_ASI    = Kabs*AST - Kabs*ASI;
// simplistic blood flow-limited distribution
// volumes (L) are scaled to BW assuming density 1 kg/L
double VPlas = 0.05*BW;
double VL    = 0.04*BW;
double VK    = 0.01*BW;
double VRest = BW - VPlas - VL - VK;

dxdt_APlasma = Kabs*ASI + QL*(ALiver/VL - APlasma/VPlas) +
                QK*(AKidney/VK - APlasma/VPlas) +
                QR*(ARest/VRest - APlasma/VPlas);
dxdt_ALiver  = QL*(APlasma/VPlas - ALiver/VL);
dxdt_AKidney = QK*(APlasma/VPlas - AKidney/VK);
dxdt_ARest   = QR*(APlasma/VPlas - ARest/VRest);
dxdt_Afeces  = 0.0;
dxdt_Aurine  = 0.0;

$TABLE
capture Plasma = APlasma/VPlas;
'

mod <- mcode("mousepbpk", MousePBPK.code)

# Set chemical-specific parameters from table
pars <- c(P_app = pfas_params$P_app,
          K_FABP_exp = pfas_params$K_FABP_exp,
          K_SA_exp = pfas_params$K_SA_exp,
          K_Glob_exp = pfas_params$K_Glob_exp,
          K_SP_exp = pfas_params$K_SP_exp,
          K_ML_exp = pfas_params$K_ML_exp)

mod <- param(mod, pars)

# Simple single oral dose of 1 mg/kg
ex <- ev(amt = 1, cmt = "AST")

tgrid <- tgrid(0, 24*5, 1)

out <- mod %>% mrgsim_d(data = ex, tgrid = tgrid)
print(head(out))
