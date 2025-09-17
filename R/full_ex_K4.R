#tree1: {1} vs {2,3,4}|D -> {2} vs{3,4}|D -> {3} vs {4}|D
#tree2: {1,4} vs {2,3}|E -> {1} vs {4}| D OR {2} vs {3}| D
#tree3: {1,2} vs {3,4}|D -> {1} vs {2}| E OR {3} vs {4}| E
#tree4: {1,2} vs {3,4}|D -> {1} vs {2}| D OR {3} vs {4}| D


## ===========================================================
## 0) Setup & switches
## ===========================================================
library(nimble)
set.seed(4242)

## ---- Choose gate type ----
useSoftmaxGate <- 0   # 0 = Dirichlet (constant mixture), 1 = softmax covariate gate

## Problem size (feel free to bump once it works)
N <- 250    # persons
I <- 8      # items
R <- 4      # categories
Tn <- 4     # number of trees
J  <- 3     # nodes per tree

## Person & item covariates (only used if useSoftmaxGate == 1)
if (useSoftmaxGate == 1) {
        P <- 2
        Q <- 1
        X_person <- scale(cbind(rnorm(N), rbinom(N, 1, .5))) # N x P
        W_item   <- scale(matrix(rnorm(I*Q), nrow = I, ncol = Q)) # I x Q
} else {
        P <- 0; Q <- 0
        X_person <- matrix(0, N, 0)
        W_item   <- matrix(0, I, 0)
}

ilogit <- function(x) 1/(1+exp(-x))

## ===========================================================
## 1) Define trees (indicator matrices) & node types
## ===========================================================
## Conventions:
## I[t, r, j] in {1,0,-1}: 1 = success/right, 0 = fail/left, -1 = node not on path
## traitIdx[t, j]: 1 = Direction (theta_dir), 2 = Extremity (theta_ext)
I_arr <- array(NA_integer_, dim = c(Tn, R, J))
traitIdx <- matrix(NA_integer_, nrow = Tn, ncol = J)

## Tree 1: D -> D -> D
I_arr[1,,] <- rbind(
        c(0, -1, -1),  # cat1
        c(1,  0, -1),  # cat2
        c(1,  1,  0),  # cat3
        c(1,  1,  1)   # cat4
)
traitIdx[1,] <- c(1,1,1)   # D,D,D

## Tree 2: E at root; then D on extremes & D on middle
I_arr[2,,] <- rbind(
        c(1, 0, -1),   # cat1: extremes then 1
        c(0, -1, 0),   # cat2: middle then 2
        c(0, -1, 1),   # cat3: middle then 3
        c(1, 1, -1)    # cat4: extremes then 4
)
traitIdx[2,] <- c(2,1,1)   # E,D,D

## Tree 3: D at root; then E within sides
I_arr[3,,] <- rbind(
        c(0, 1, -1),   # cat1: low side then extreme => 1
        c(0, 0, -1),   # cat2: low side then less extreme => 2
        c(1, -1, 0),   # cat3: high side then less extreme => 3
        c(1, -1, 1)    # cat4: high side then extreme => 4
)
traitIdx[3,] <- c(1,2,2)   # D,E,E

## Tree 4: D at root; then D within sides
I_arr[4,,] <- rbind(
        c(0, 0, -1),   # cat1: low side then 1
        c(0, 1, -1),   # cat2: low side then 2
        c(1, -1, 0),   # cat3: high side then 3
        c(1, -1, 1)    # cat4: high side then 4
)
traitIdx[4,] <- c(1,1,1)   # D,D,D

## ===========================================================
## 2) True parameters & simulation
## ===========================================================
## Person traits
theta_dir_true <- rnorm(N, 0, 1)
theta_ext_true <- rnorm(N, 0, 1)

## Item node difficulties: b[t, j, i]
## Keep priors modest to retain separation between trees
b_true <- array(NA_real_, dim = c(Tn, J, I))
for (t in 1:Tn) for (j in 1:J) for (i in 1:I) {
        b_true[t, j, i] <- rnorm(1, 0, 0.7)
}
## (Optional: tie shared D-nodes across trees to help identifiability)
## e.g., Tree 1 node3 and Tree 4 node3 are both D on {3 vs 4}
## b_true[4,3,] <- b_true[1,3,]

## Gate parameters
if (useSoftmaxGate == 1) {
        ## Multinomial logit (softmax); set tree 1 as baseline (gamma = 0) for identifiability
        gamma0_true <- c(0.0, 0.4, -0.2, 0.1)           # length Tn; gamma0_true[1] ignored (baseline 0)
        gamma0_true[1] <- 0
        gamma_p_true   <- matrix(0, nrow = Tn, ncol = P) # Tn x P
        gamma_i_true   <- matrix(0, nrow = Tn, ncol = Q) # Tn x Q
        if (P > 0) {
                gamma_p_true[2,] <- c( 0.7, -0.4)
                gamma_p_true[3,] <- c(-0.5,  0.3)
                gamma_p_true[4,] <- c( 0.2,  0.2)
        }
        if (Q > 0) {
                gamma_i_true[2,] <- c( 0.5)
                gamma_i_true[3,] <- c(-0.2)
                gamma_i_true[4,] <- c( 0.0)
        }
} else {
        ## Constant mixture over trees
        omega_true <- c(0.35, 0.25, 0.20, 0.20)  # length Tn; sums to 1
}

## Simulate responses
y <- matrix(NA_integer_, N, I)

for (n in 1:N) for (i in 1:I) {

        ## Node probabilities for each tree/node
        p_tj <- matrix(NA_real_, nrow = Tn, ncol = J)
        for (t in 1:Tn) for (j in 1:J) {
                eta <- - b_true[t, j, i]
                if (traitIdx[t, j] == 1) eta <- eta + theta_dir_true[n]
                if (traitIdx[t, j] == 2) eta <- eta + theta_ext_true[n]
                p_tj[t, j] <- ilogit(eta)
        }

        ## Category probs for each tree via path product
        pi_t_r <- matrix(NA_real_, nrow = Tn, ncol = R)
        for (t in 1:Tn) for (r in 1:R) {
                prod <- 1.0
                for (j in 1:J) {
                        if (I_arr[t, r, j] == 1)  prod <- prod * p_tj[t, j]
                        if (I_arr[t, r, j] == 0)  prod <- prod * (1 - p_tj[t, j])
                        ## if -1: node not used; multiply by 1
                }
                pi_t_r[t, r] <- prod
        }

        ## Gate: mixture over trees
        if (useSoftmaxGate == 1) {
                eta_t <- rep(0, Tn)
                eta_t[1] <- 0
                for (t in 2:Tn) {
                        add_p <- if (P > 0) sum(gamma_p_true[t,] * X_person[n,]) else 0
                        add_i <- if (Q > 0) sum(gamma_i_true[t,] * W_item[i,])   else 0
                        eta_t[t] <- gamma0_true[t] + add_p + add_i
                }
                eeta <- exp(eta_t - max(eta_t))
                w_t  <- eeta / sum(eeta)
        } else {
                w_t <- omega_true
        }

        ## Mixture over trees and draw Y
        mix <- as.numeric(w_t %*% pi_t_r)
        y[n, i] <- sample.int(R, size = 1, prob = mix)
}

## ===========================================================
## 3) NIMBLE model (collapsed mixture over trees)
## ===========================================================
modelCode <- nimbleCode({

        ## Person traits
        for (n in 1:N) {
                theta_dir[n] ~ dnorm(0, 1)
                theta_ext[n] ~ dnorm(0, 1)
        }

        ## Item node difficulties
        for (t in 1:Tn) for (j in 1:J) for (i in 1:I) {
                b[t, j, i] ~ dnorm(0, 1/(1.5^2))   # sd ~ 1.5
        }

        ## Gate parameters
        if (useSoftmaxGate == 1) {
                ## Multinomial logit; tree 1 is baseline: gamma0[1]=0, gamma_p[1,]=0, gamma_i[1,]=0 (not sampled)
                for (t in 2:Tn) {
                        gamma0[t] ~ dnorm(0, 1/(2^2))
                        for (p in 1:P) gamma_p[t, p] ~ dnorm(0, 1/(2^2))
                        for (q in 1:Q) gamma_i[t, q] ~ dnorm(0, 1/(2^2))
                }
        } else {
                ## Constant mixture weights
                omega[1:Tn] ~ ddirch(alpha_dirch[1:Tn])
        }

        ## Likelihood (collapsed mixture)
        for (n in 1:N) for (i in 1:I) {

                ## Node probabilities for all trees/nodes
                for (t in 1:Tn) for (j in 1:J) {
                        eta_node[n, i, t, j] <- - b[t, j, i] +
                                (traitIdx[t, j] == 1) * theta_dir[n] +
                                (traitIdx[t, j] == 2) * theta_ext[n]
                        p_tj[n, i, t, j] <- ilogit(eta_node[n, i, t, j])
                }

                ## Per-tree, per-category path products
                for (t in 1:Tn) for (r in 1:R) {
                        comp[n, i, t, r] <- 1
                        for (j in 1:J) {
                                comp[n, i, t, r] <- comp[n, i, t, r] *
                                        pow(p_tj[n, i, t, j],      I_arr[t, r, j] == 1) *
                                        pow(1 - p_tj[n, i, t, j],  I_arr[t, r, j] == 0)
                        }
                }

                ## Gate probabilities over trees
                if (useSoftmaxGate == 1) {
                        ## Tree 1 baseline => eta_gate[n,i,1] = 0
                        eta_gate[n, i, 1] <- 0
                        for (t in 2:Tn) {
                                ## person contribution
                                gate_p[n, i, t] <- 0
                                for (p in 1:P) gate_p[n, i, t] <- gate_p[n, i, t] + gamma_p[t, p] * X_person[n, p]
                                ## item contribution
                                gate_i[n, i, t] <- 0
                                for (q in 1:Q) gate_i[n, i, t] <- gate_i[n, i, t] + gamma_i[t, q] * W_item[i, q]
                                eta_gate[n, i, t] <- gamma0[t] + gate_p[n, i, t] + gate_i[n, i, t]
                        }
                        ## softmax
                        sum_e[n, i] <- 0
                        for (t in 1:Tn) {
                                e_eta[n, i, t] <- exp(eta_gate[n, i, t])
                                sum_e[n, i] <- sum_e[n, i] + e_eta[n, i, t]
                        }
                        for (t in 1:Tn) gateprob[n, i, t] <- e_eta[n, i, t] / sum_e[n, i]

                } else {
                        for (t in 1:Tn) gateprob[n, i, t] <- omega[t]
                }

                ## Mixture over trees → category probabilities
                for (r in 1:R) {
                        mixtureProb[n, i, r] <- 0
                        for (t in 1:Tn) mixtureProb[n, i, r] <- mixtureProb[n, i, r] + gateprob[n, i, t] * comp[n, i, t, r]
                }

                ## Observation
                y[n, i] ~ dcat(mixtureProb[n, i, 1:R])
        }
})

## ===========================================================
## 4) Build, compile, run MCMC
## ===========================================================
constants <- list(N = N, I = I, R = R, Tn = Tn, J = J,
                  I_arr = I_arr, traitIdx = traitIdx,
                  useSoftmaxGate = useSoftmaxGate,
                  P = P, Q = Q)
dataList <- list(y = y)
if (useSoftmaxGate == 1) {
        dataList$X_person <- X_person
        dataList$W_item   <- W_item
} else {
        ## Dirichlet prior concentration (can be all 1's)
        dataList$alpha_dirch <- rep(1, Tn)
}
## Initial values
inits <- list(
        theta_dir = rnorm(N, 0, 0.5),
        theta_ext = rnorm(N, 0, 0.5),
        b = array(rnorm(Tn*J*I, 0, 0.5), dim = c(Tn, J, I))
)
if (useSoftmaxGate == 1) {
        inits$gamma0 <- c(0, rnorm(Tn-1, 0, 0.2))
        if (P > 0) inits$gamma_p <- matrix(0, nrow = Tn, ncol = P)
        if (Q > 0) inits$gamma_i <- matrix(0, nrow = Tn, ncol = Q)
} else {
        inits$omega <- rep(1/Tn, Tn)
}

Rmodel <- nimbleModel(modelCode, data = dataList, constants = constants, inits = inits, check = FALSE)
Cmodel <- compileNimble(Rmodel)

## Monitors
mon <- c("b", "theta_dir", "theta_ext")
if (useSoftmaxGate == 1) {
        mon <- c(mon, "gamma0", "gamma_p", "gamma_i")
} else {
        mon <- c(mon, "omega")
}

conf <- configureMCMC(Rmodel, monitors = mon)
## For better mixing on gate params (softmax), use slice samplers
if (useSoftmaxGate == 1) {
        conf$removeSamplers(c("gamma0"))
        for (t in 2:Tn) conf$addSampler(target = paste0("gamma0[", t, "]"), type = "slice")
        if (P > 0) {
                for (t in 2:Tn) for (p in 1:P)
                        conf$addSampler(target = paste0("gamma_p[", t, ", ", p, "]"), type = "slice")
        }
        if (Q > 0) {
                for (t in 2:Tn) for (q in 1:Q)
                        conf$addSampler(target = paste0("gamma_i[", t, ", ", q, "]"), type = "slice")
        }
}
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, niter = 8000, nburnin = 2000, thin = 10,
                   nchains = 2, samplesAsCodaMCMC = TRUE)

## ===========================================================
## 5) Sanity checks
## ===========================================================
library(coda)
sm <- summary(samples)
print(sm$statistics[1:min(30, nrow(sm$statistics)), , drop=FALSE])

## --- Gate recovery ---
if (useSoftmaxGate == 1) {
        ## posterior means
        ms <- lapply(samples, as.matrix)
        M  <- as.matrix(samples)
        ## Pull gamma0 (t=2..Tn), gamma_p, gamma_i posterior means
        g0_hat <- colMeans(M[, grep("^gamma0\\[", colnames(M)), drop=FALSE])
        gp_hat <- if (P>0) colMeans(M[, grep("^gamma_p\\[", colnames(M)), drop=FALSE]) else numeric(0)
        gi_hat <- if (Q>0) colMeans(M[, grep("^gamma_i\\[", colnames(M)), drop=FALSE]) else numeric(0)

        cat("\nTRUE gamma0 (t=2..Tn):   ", paste(round(gamma0_true[2:Tn], 3), collapse=", "), "\n")
        cat("POST mean gamma0 (t=2..Tn):", paste(round(g0_hat, 3), collapse=", "), "\n")

        if (P>0) {
                cat("\nTRUE gamma_p (rows t=2..Tn):\n"); print(round(gamma_p_true[2:Tn,,drop=FALSE], 3))
                cat("POST mean gamma_p (vectorized t-major):\n"); print(round(matrix(gp_hat, nrow=Tn, byrow=FALSE), 3))
        }
        if (Q>0) {
                cat("\nTRUE gamma_i (rows t=2..Tn):\n"); print(round(gamma_i_true[2:Tn,,drop=FALSE], 3))
                cat("POST mean gamma_i (vectorized t-major):\n"); print(round(matrix(gi_hat, nrow=Tn, byrow=FALSE), 3))
        }

} else {
        M  <- as.matrix(samples)
        omega_hat <- colMeans(M[, grep("^omega\\[", colnames(M)), drop=FALSE])
        cat("\nTRUE omega:    ", paste(round(omega_true, 3), collapse=", "), "\n")
        cat("POST mean omega:", paste(round(omega_hat, 3), collapse=", "), "\n")
}

## --- Tree separability diagnostic (average KL between trees per item) ---
## Using the true parameters from simulation to see if trees were separable
kl <- function(p, q) sum(p * log((p + 1e-12)/(q + 1e-12)))
avgKL <- numeric(I)
for (i in 1:I) {
        val <- 0
        for (n in 1:N) {
                ## recompute p_tj under true params
                p_tj <- matrix(NA, nrow=Tn, ncol=J)
                for (t in 1:Tn) for (j in 1:J) {
                        eta <- - b_true[t, j, i]
                        if (traitIdx[t, j] == 1) eta <- eta + theta_dir_true[n]
                        if (traitIdx[t, j] == 2) eta <- eta + theta_ext_true[n]
                        p_tj[t, j] <- ilogit(eta)
                }
                ## per-tree category probs
                pi_t_r <- matrix(NA, nrow=Tn, ncol=R)
                for (t in 1:Tn) for (r in 1:R) {
                        prod <- 1
                        for (j in 1:J) {
                                if (I_arr[t, r, j] == 1) prod <- prod * p_tj[t, j]
                                if (I_arr[t, r, j] == 0) prod <- prod * (1 - p_tj[t, j])
                        }
                        pi_t_r[t, r] <- prod
                }
                ## average pairwise KL for this (n,i)
                for (t1 in 1:(Tn-1)) for (t2 in (t1+1):Tn)
                        val <- val + 0.5 * (kl(pi_t_r[t1,], pi_t_r[t2,]) + kl(pi_t_r[t2,], pi_t_r[t1,]))
        }
        avgKL[i] <- val / N
}
cat("\nAvg pairwise KL between trees per item:\n")
print(round(avgKL, 3))
