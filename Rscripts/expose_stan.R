library(cmdstanr)
library(Rcpp)
library(RcppEigen)



expose_cmdstanr_functions <- function(model_path, include_paths = NULL,
                                     expose_to_global_env = FALSE) {
  required_pkgs <- c("Rcpp", "RcppEigen", "cmdstanr")
  found_pkgs <- required_pkgs %in% rownames(installed.packages())
  if (!all(found_pkgs)) {
    stop(
      "The following required packages are missing: ",
      paste0(required_packages[!found_pkgs], collapse = ", "),
      "."
    )
  }
  if (cmdstanr::cmdstan_version() < "2.26.0") {
    stop("Please install CmdStan version 2.26 or newer.", call. = FALSE)
  }
  get_cmdstan_flags <- function(flag_name) {
    cmdstan_path <- cmdstanr::cmdstan_path()
    flags <- processx::run(
      "make", 
      args = c(paste0("print-", flag_name)),
      wd = cmdstan_path
    )$stdout
    flags <- gsub(
      pattern = paste0(flag_name, " ="),
      replacement = "", x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = " stan/", replacement = paste0(" ", cmdstan_path, "/stan/"),
      x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = "-I lib/", replacement = paste0("-I ", cmdstan_path, "/lib/"),
      x = flags, fixed = TRUE
    )
    flags <- gsub(
      pattern = "-I src", replacement = paste0("-I ", cmdstan_path, "/src"),
      x = flags, fixed = TRUE
    )
    gsub("\n", "", flags)
  }
  temp_stan_file <- tempfile(pattern = "model-", fileext = ".stan")
  temp_cpp_file <- paste0(tools::file_path_sans_ext(temp_stan_file), ".cpp")
  file.copy(model_path, temp_stan_file, overwrite = TRUE)
  if (isTRUE(.Platform$OS.type == "windows")) {
    stanc3 <- "./bin/stanc.exe"
  } else {
    stanc3 <- "./bin/stanc"
  }
  processx::run(
    stanc3,
    args = c(
      temp_stan_file,
      "--standalone-functions",
      paste0("--include-paths=", include_paths),
      paste0("--o=",temp_cpp_file)
    ),
    wd = cmdstanr::cmdstan_path()
  )
  code <- paste(readLines(temp_cpp_file), collapse = "\n")
  code <- paste(
    "// [[Rcpp::depends(RcppEigen)]]",
    "#include <stan/math/prim/fun/Eigen.hpp>",
    "#include <RcppCommon.h>
    #include <boost/random/additive_combine.hpp>
    #include <iostream>

    namespace Rcpp {
      SEXP wrap(boost::ecuyer1988 RNG);
      SEXP wrap(boost::ecuyer1988& RNG);
      SEXP wrap(std::ostream stream);
      template <> boost::ecuyer1988 as(SEXP ptr_RNG);
      template <> boost::ecuyer1988& as(SEXP ptr_RNG);
      template <> std::ostream* as(SEXP ptr_stream);
      namespace traits {
        template <> class Exporter<boost::ecuyer1988&>;
        template <> struct input_parameter<boost::ecuyer1988&>;
      }
    }

    #include <Rcpp.h>

    namespace Rcpp {
      SEXP wrap(boost::ecuyer1988 RNG){
        boost::ecuyer1988* ptr_RNG = &RNG;
        Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
        return Xptr_RNG;
      }

      SEXP wrap(boost::ecuyer1988& RNG){
        boost::ecuyer1988* ptr_RNG = &RNG;
        Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
        return Xptr_RNG;
      }

      SEXP wrap(std::ostream stream) {
        std::ostream* ptr_stream = &stream;
        Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
        return Xptr_stream;
      }

      template <> boost::ecuyer1988 as(SEXP ptr_RNG) {
        Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
        boost::ecuyer1988& RNG = *ptr;
        return RNG;
      }

      template <> boost::ecuyer1988& as(SEXP ptr_RNG) {
        Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
        boost::ecuyer1988& RNG = *ptr;
        return RNG;
      }

      template <> std::ostream* as(SEXP ptr_stream) {
        Rcpp::XPtr<std::ostream> ptr(ptr_stream);
        return ptr;
      }

      namespace traits {
        template <> class Exporter<boost::ecuyer1988&> {
        public:
          Exporter( SEXP x ) : t(Rcpp::as<boost::ecuyer1988&>(x)) {}
          inline boost::ecuyer1988& get() { return t ; }
        private:
          boost::ecuyer1988& t ;
        } ;

        template <>
        struct input_parameter<boost::ecuyer1988&> {
          typedef
          typename Rcpp::ConstReferenceInputParameter<boost::ecuyer1988&> type ;
          //typedef typename boost::ecuyer1988& type ;
        };
      }
    }

    RcppExport SEXP get_stream_() {
      std::ostream* pstream(&Rcpp::Rcout);
      Rcpp::XPtr<std::ostream> ptr(pstream, false);
      return ptr;
    }

    RcppExport SEXP get_rng_(SEXP seed) {
      int seed_ = Rcpp::as<int>(seed);
      boost::ecuyer1988* rng = new boost::ecuyer1988(seed_);
      Rcpp::XPtr<boost::ecuyer1988> ptr(rng, true);
      return ptr;
    }
    ",
    "#include <RcppEigen.h>",
    code,
    sep = "\n"
  )
  code <- gsub("// [[stan::function]]",
               "// [[Rcpp::export]]", code, fixed = TRUE)
  code <- gsub(
    "stan::math::accumulator<double>& lp_accum__, std::ostream* pstream__ = nullptr){",
    "std::ostream* pstream__ = nullptr){\nstan::math::accumulator<double> lp_accum__;",
    code,
    fixed = TRUE
  )
  code <- gsub("__ = nullptr", "__ = 0", code, fixed = TRUE)

  get_stream <- function() {
    return(.Call('get_stream_'))
  }
  get_rng <- function(seed=0L) {
    if (!identical(seed, 0L)) {
      if (length(seed) != 1)
        stop("Seed must be a length-1 integer vector.")
    }
    return(.Call('get_rng_', seed))
  }
  if (expose_to_global_env) {
    env = globalenv()
  } else {
    env = new.env()
  }
  compiled <- withr::with_makevars(
    c(
      USE_CXX14 = 1,
      PKG_CPPFLAGS = "",
      PKG_CXXFLAGS = get_cmdstan_flags("CXXFLAGS"),
      PKG_LIBS = paste0(
        get_cmdstan_flags("LDLIBS"),
        get_cmdstan_flags("LIBSUNDIALS"),
        get_cmdstan_flags("TBB_TARGETS"),
        get_cmdstan_flags("LDFLAGS_TBB")
      )
    ),
    Rcpp::sourceCpp(code = code, env = env)
  )
  for (x in compiled$functions) {
    FUN <- get(x, envir = env)
    args <- formals(FUN)
    args$pstream__ <- get_stream()
    if ("lp__" %in% names(args)) args$lp__ <- 0
    if ("base_rng__" %in% names(args)) args$base_rng__ <- get_rng()
    formals(FUN) <- args
    assign(x, FUN, envir = env)
  }
  assign("stan_rng__", get_rng, envir = env)
  if (expose_to_global_env) {
    invisible(NULL)
  } else {
    return(env)
  }
}

udfs <- expose_cmdstanr_functions("stan_models/Linear.neutral.FO.stan")
udes <- expose_cmdstanr_functions("stan_models/Linear.neutral.TAA4.1.stan")
print(udfs$chi_source(12))

initial_conds <- c(0, 0, 1187133, 6083180)
initial_cond1 <- c( 0, 0, 106.3357, 92984.5643)
parms <- c(0.027, 0.009, 0.026, 5.3)
parms2 <- c(0.01, 0.01, 0.01, 0.01, 67)
parms3 <- c(0.01, 0.01, 0.01, 0.01, 88)
#ageatbmt <- c(45, 52, 55)
solve_time1 <- c(59, 69, 76, 88, 95, 102, 108, 109, 113, 119, 122, 123, 124, 141, 156, 158, 183, 212, 217, 219, 235, 261, 270, 289, 291, 306, 442, 524, 563, 566, 731)
solve_ageatBMT <- c(41, 41, 41, 41, 74, 74, 80, 74, 99, 41, 80, 41, 42, 51, 58, 65, 65, 58, 70, 65, 74, 51, 89, 51, 101, 89, 64, 90, 87, 90, 87)

#combine solve_ ageatbmt and solve_time1 in 2d array

solve_time <- t(rbind(solve_time1, solve_ageatBMT))



ageatbmt1 <- rep(41, length(x = solve_time1))
solve_time2 <- c(88, 95, 102, 108, 109, 113, 119, 122, 123, 124, 141, 156, 158, 183, 212, 217, 219, 235, 261, 270, 289, 291, 306, 442, 524, 563, 566, 731)
ageatbmt2 <- rep(67, length(x = solve_time2))
solve_time3 <- c(102, 108, 109, 113, 119, 122, 123, 124, 141, 156, 158, 183, 212, 217, 219, 235, 261, 270, 289, 291, 306, 442, 524, 563, 566, 731)
ageatbmt3 <- rep(88, length(x = solve_time3))
# initial_gen <- udfs$initial_cond_generator(initial_conds, ageatbmt1[1], parms)
# print(initial_gen)
# initial_gen1 <- udfs$initial_cond_generator(initial_cond1, ageatbmt2[1], parms2)
# print(initial_gen1)
# initial_gen2 <- udfs$initial_cond_generator(initial_cond1, ageatbmt3[1], parms3)
# print(initial_gen2)
# intial_gen_T2 <- udes$initial_cond_generator(initial_conds, ageatbmt1[1], parms)

mysol <- udes$ode_solver(initial_conds, solve_time1, solve_ageatBMT, parms)
# mysol_T2 <- udes$ode_pred(initial_conds, solve_time1, ageatbmt1[1], parms)
# mysol1 <- udfs$ode_pred(initial_conds, solve_time2, ageatbmt2[1], parms2)
# mysol2 <- udfs$ode_pred(initial_conds, solve_time3, ageatbmt3[1], parms3)
# print(mysol)
print(mysol)
# print(mysol1)
# print(mysol2)

# #initialize vectors to store solutions
# X1 <- numeric(length(mysol)/4)
# X2 <- numeric(length(mysol)/4)
# X3 <- numeric(length(mysol)/4)
# X4 <- numeric(length(mysol)/4)

# #separate four solutions from mysol

# for(i in 1:length(mysol)/4     )
# { X1[i] <- mysol[4*i-3]
#   X2[i] <- mysol[4*i-2]
#   X3[i] <- mysol[4*i-1]
#   X4[i] <- mysol[4*i]
# }

# print(X1)
# print(X2)
# print(X3)
# print(X4)
# #save solution to csv file

write.csv(mysol, "mysol.csv")
# write.csv(mysol_T2, "mysol_T2.csv")

#Calculate Mz_total and Nfd



#for (i in 1:3) {
 # myinit <- udfs$initial_cond_generator(initial_conds, ageatbmt[i], parms)
  #print(ageatbmt[i])
  #print(myinit)
#}

# for (i in 1:3) {
#  mysol <- udfs$prediction_generator(initial_conds, solve_time[i], ageatbmt[i], parms)
#   print(solve_time[i])
#   print(mysol)
# }

#myfin <- udfs$ode_solver(initial_conds, solve_time, ageatbmt, parms)
#print(myfin)

print("Done")




