setwd("E:/RStudio/thesis/Code/FCIR")
source("generate_without_sparsity.r")
source("estimate_FCIR.r")
source("estimate_penalized_FCIR.r")

# 1. Define the parameter grid you want to sweep across
Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# Global fixed parameters
L = 3          
K = 2          
B_reps = 1000  
seed = 42      

# Record the total start time
total_start = Sys.time()

# 2. Start the nested loops
for (n in Ns) {
  for (p in Ps) {
    
    # NEW FILENAMES for dense experiments to avoid overwriting existing sparse results
    data_filename = paste0("../../Simulation_Results/FCIR_data_Dense_N", n, "_P", p, ".Rdata")
    est_filename  = paste0("../../Simulation_Results/FCIR_estimates_Dense_N", n, "_P", p, ".Rdata")
    est_filename_pen = paste0("../../Simulation_Results/FCIR_estimates_penalized_Dense_N", n, "_P", p, ".Rdata")
    
    if (!dir.exists(dirname(data_filename))) {
      dir.create(dirname(data_filename), recursive = TRUE)
    }
    
    if (!file.exists(data_filename)) {
      print(paste("Generating DENSE data ( N =", n, ", P =", p, ")..."))
      generate_fully_dense_fcir_data(N = n, P = p, L = L, K = K, B_reps = B_reps, 
                                     seed = seed, filename = data_filename)
    } else {
      print(paste("DENSE Data already exists for N =", n, ", P =", p, "- Skipping generation."))
    }
    
    print("Loading DENSE data...")
    load(data_filename) 
    
    # ---------------- Unpenalized Estimation ----------------
    if (!file.exists(est_filename)) {
      print(paste("Starting 1000 Monte Carlo Unpenalized estimations for Dense Data N =", n, "P =", p, "..."))
      est_beta_0  = matrix(NA, nrow = B_reps, ncol = L)
      est_B_mat   = array(NA, dim = c(L, K, B_reps))
      est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
      est_A_mat   = array(NA, dim = c(L, K, B_reps))
      
      for (b in 1:B_reps) {
        Y_b = Y[, , b] 
        result = estimate_unpenalized_FCIR(Y = Y_b, X = X, Tr = Tr)
        
        est_beta_0[b, ]  = result$beta_0
        est_B_mat[,, b]  = result$B_mat
        est_alpha_0[b, ] = result$alpha_0
        est_A_mat[,, b]  = result$A_mat
      }
      
      # Save unpenalized results
      save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat,
           beta_0, B_mat, alpha_0, A_mat,
           N = n, P = p, L, K, B_reps, seed,
           file = est_filename)
    } else {
      print(paste("Unpenalized Dense estimates already exist for N =", n, ", P =", p, "- Skipping."))
    }
    
    # ---------------- Penalized Estimation (Disabled temporarily) ----------------
    # if (!file.exists(est_filename_pen)) {
    #   print(paste("Starting 1000 Monte Carlo Penalized estimations for Dense Data N =", n, "P =", p, "..."))
    #   est_beta_0  = matrix(NA, nrow = B_reps, ncol = L)
    #   est_B_mat   = array(NA, dim = c(L, K, B_reps))
    #   est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
    #   est_A_mat   = array(NA, dim = c(L, K, B_reps))
    #   
    #   for (b in 1:B_reps) {
    #     Y_b = Y[, , b] 
    #     # Penalized (alpha=1 for Lasso by default)
    #     result_pen = estimate_penalized_FCIR(Y = Y_b, X = X, Tr = Tr)
    #     
    #     est_beta_0[b, ]  = result_pen$beta_0
    #     est_B_mat[,, b]  = result_pen$B_mat
    #     est_alpha_0[b, ] = result_pen$alpha_0
    #     est_A_mat[,, b]  = result_pen$A_mat
    #   }
    #   
    #   save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat,
    #        beta_0, B_mat, alpha_0, A_mat,
    #        N = n, P = p, L, K, B_reps, seed,
    #        file = est_filename_pen)
    # } else {
    #   print(paste("Penalized Dense estimates already exist for N =", n, ", P =", p, "- Skipping."))
    # }
    
  }
}

print(paste("Total execution time (Dense Experiment):", Sys.time() - total_start))