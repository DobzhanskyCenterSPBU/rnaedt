f <- function(a, b, m) {
    t <- c(0:(m - 1));
    if (sum(is.na(log(a + t) - log(b - t - 1))) > 0){
        print(c(a, b, m))
        stop("Nan")
    }
    return(cumsum(c(lbeta(a, b), log((a + t) / (b - t - 1)))))
}

L_A <- function(i, a_edt, b_edt, a_0, b_0, edt_coverage, edt_ref_count, edt_edt_count, f, p) {
    n <- edt_coverage[i];
    n_a <- edt_ref_count[i];
    n_g <- edt_edt_count[i];
    n_ct <- n - n_a - n_g;
    m1 <- matrix(c(1:(n_a + 1)), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c(0:n_g), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_edt, b_edt + n_a + n_g, n_a + n_g) - lbeta(a_edt, b_edt)
    m_edt <- y[m]
    m1 <- matrix(c(0:n_a), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c((n_g + 1):1), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_0 + n - n_a - n_g, b_0 + n_a + n_g, n_a + n_g) - lbeta(a_0, b_0) - c(-n_g:n_a) * log(3)
    m_n <- y[m]
    m_ca <- matrix(p[[n_a + 1]], nrow = n_a + 1, ncol = n_g + 1);
    m_cg <- t(matrix(p[[n_g + 1]], nrow = n_g + 1, ncol = n_a + 1));
    m_const <- matrix(n_ct * log(2) - (n - n_a) * log(3), nrow = n_a + 1, ncol = n_g + 1)
    l <- m_edt + m_n + m_ca + m_cg + m_const;
    l_max <- max(l);
    l <- l - l_max;
    if (log(sum(exp(l))) + l_max == Inf | is.na(log(sum(exp(l))) + l_max)) {
        print(c(a_edt, b_edt, a_0, b_0, n, n_a, n_g))
        stop("Inf occured")
    }
    return(log(sum(exp(l))) + l_max)
}

#' @keywords internal
#'
#' @importFrom parallel mclapply
L <- function(a_edt, b_edt, a_0, b_0, k, noise, edt_number, edt_coverage, edt_ref_count, edt_edt_count, p) {
    cat(a_edt, b_edt, a_0, b_0, "\n")
    noise_counts <- noise$noise_counts
    reference_counts <- noise$reference_counts
    l1 <- lbeta(a_0 + noise_counts, b_0 + reference_counts) - lbeta(a_0, b_0)
    x <- sum(unlist(mclapply(
        1:edt_number,
        L_A,
        a_edt = a_edt,
        b_edt = b_edt,
        a_0 = a_0,
        b_0 = b_0,
        edt_coverage = edt_coverage,
        edt_ref_count = edt_ref_count,
        edt_edt_count = edt_edt_count,
        f = f,
        p = p,
        mc.cores = k)))
    if ((-sum(l1)-x == Inf) | is.na(-sum(l1)-x)) {
        print(c(a_edt, b_edt, a_0, b_0, -sum(l1), x))
        stop("Inf occured")
    }
    return(-sum(l1)-x)
}

I_1 <- function(i, a_edt, b_edt, a_0, b_0, edt_coverage, edt_ref_count, edt_edt_count, f, p){
    n <- edt_coverage[i]
    n_a <- edt_ref_count[i]
    n_g <- edt_edt_count[i]
    m1 <- matrix(c(1:(n_a + 1)), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c(0:n_g), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_edt + 1, b_edt + n_a + n_g, n_a + n_g) - lbeta(a_edt, b_edt)
    m_edt <- y[m]
    m1 <- matrix(c(0:n_a), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c((n_g + 1):1), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_0 + n - n_a - n_g, b_0 + n_a + n_g, n_a + n_g) - lbeta(a_0, b_0) - c(-n_g:n_a) * log(3)
    m_n <- y[m]
    m_ca <- matrix(p[[n_a + 1]], nrow = n_a + 1, ncol = n_g + 1);
    m_cg <- t(matrix(p[[n_g + 1]], nrow = n_g + 1, ncol = n_a + 1));
    l <- m_edt + m_n + m_ca + m_cg;
    l_max <- max(l);
    l <- l - l_max;
    if (log(sum(exp(l))) + l_max == Inf | is.na(log(sum(exp(l))) + l_max)) {
        print(c(a_edt, b_edt, a_0, b_0, n, n_a, n_g))
        stop("Inf occured")
    }
    return(log(sum(exp(l))) + l_max)
}

#' @keywords internal
#'
#' @importFrom stats pbeta
J_in <- function(i, t, a_edt, b_edt, a_0, b_0, edt_coverage, edt_ref_count, edt_edt_count, f, p) {
    n <- edt_coverage[i]
    n_a <- edt_ref_count[i]
    n_g <- edt_edt_count[i]
    m1 <- matrix(c(1:(n_a + 1)), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c(0:n_g), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_edt, b_edt + n_a + n_g, n_a + n_g) - lbeta(a_edt, b_edt)
    tmp <- c(0:(n_a + n_g))
    y <- y + pbeta(t, a_edt + tmp, b_edt + n_a + n_g - tmp, log.p = TRUE)
    if (sum(is.na(y)) | sum(y == Inf) | sum(y == -Inf)){
        print(c(t, n, n_a, n_g))
        print("Inf occured")
    }
    m_edt <- y[m]
    m1 <- matrix(c(0:n_a), nrow = n_a + 1, ncol = n_g + 1)
    m2 <- t(matrix(c((n_g + 1):1), nrow = n_g + 1, ncol = n_a + 1))
    m <- m1 + m2
    y <- f(a_0 + n - n_a - n_g, b_0 + n_a + n_g, n_a + n_g) - lbeta(a_0, b_0) - c(-n_g:n_a) * log(3)
    m_n <- y[m]
    m_ca <- matrix(p[[n_a + 1]], nrow = n_a + 1, ncol = n_g + 1);
    m_cg <- t(matrix(p[[n_g + 1]], nrow = n_g + 1, ncol = n_a + 1));
    l <- m_edt + m_n + m_ca + m_cg;
    l_max <- max(l);
    l <- l - l_max;
    if (log(sum(exp(l))) + l_max == Inf | is.na(log(sum(exp(l))) + l_max)){
        print(c(a_edt, b_edt, a_0, b_0, n, n_a, n_g))
        print("Inf occured")
    }
    return(log(sum(exp(l))) + l_max)
}


#' @title Performs RNA editing analysis
#'
#' @author Olga Dudkina
#'
#' @description Processes data supplied from \code{bamrd} and estimates RNA editing probability
#'
#' @param input A list of the data frames returned by \code{bamrd} function
#' @param strand_specific A \code{logical(1)} indicating the presence either of the strand-specific (\code{TRUE}) RNA-seq experiment, or non-strand-specific (\code{FALSE}) one
#' @param threads An \code{integer(1)} number of threads
#'
#' @return Data frame with results of RNA editing analysis
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate filter select arrange
#' @importFrom rlang .data
#' @importFrom parallel makeCluster parLapply stopCluster mclapply
#' @importFrom MASS fitdistr
#' @importFrom stats coef
#' @importFrom stats4 mle
#'
#' @export
#'
#' @example
#' rnaedt
#'
rnaedt <- function(input, strand_specific = TRUE, threads = 20) {
    all_sites <- input$pileup_table
    if (strand_specific) {
        all_sites <- all_sites %>%
            filter(.data$coverage > 5) %>%
            mutate(
                reference_counts = ifelse(.data$reference == 'A', .data$A, ifelse(.data$reference == 'C', .data$C, ifelse(.data$reference == 'T', .data$T, .data$G))),
                noise_counts = .data$coverage - .data$reference_counts,
                fraction = .data$noise_counts / .data$coverage,
                is_editing = ifelse((.data$strand == '+' & .data$reference == 'A') | (.data$strand == '-' & .data$reference == 'T'), TRUE, FALSE),
                is_noise = ifelse(.data$is_editing == TRUE, FALSE, TRUE)) %>%
            filter(.data$fraction < 0.95)
    } else {
        all_sites <- all_sites %>%
            filter(.data$coverage > 5) %>%
            mutate(
                reference_counts = ifelse(.data$reference == 'A', .data$A, ifelse(.data$reference == 'C', .data$C, ifelse(.data$reference == 'T', .data$T, .data$G))),
                noise_counts = .data$coverage - .data$reference_counts,
                fraction = .data$noise_counts / .data$coverage,
                is_editing = ifelse(.data$reference == 'A' | .data$reference == 'T', TRUE, FALSE),
                is_noise = ifelse(.data$is_editing == TRUE, FALSE, TRUE)) %>%
            filter(.data$fraction < 0.95)
    }

    cat("length of RNA sequence", nrow(all_sites), "\n")

    # combinations table calculation
    M_0 <- max(
        all_sites$A[all_sites$is_editing == T & all_sites$reference == 'A'],
        all_sites$G[all_sites$is_editing == T & all_sites$reference == 'A'],
        all_sites$T[all_sites$is_editing == T & all_sites$reference == 'T'],
        all_sites$C[all_sites$is_editing == T & all_sites$reference == 'T'])

    k <- threads
    cl <- makeCluster(k)
    system.time(p <- parLapply(cl, c(0:M_0), function(x) lchoose(x, c(0:x))))
    stopCluster(cl)

    # choose start points
    noise <- filter(all_sites, .data$is_noise == TRUE)
    fractions_for_models <- select(.data$noise, .data$fraction)

    # moments
    mu <- mean(noise$fraction)
    var <- var(noise$fraction)
    a_mom <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    b_mom <- a_mom * (1 / mu - 1)

    # fit distribution
    coeffs <- tryCatch(
        {
            m <- fitdistr(fractions_for_models$fraction, 'beta', start = list(shape1 = a_mom, shape2 = b_mom), lower = c(0, 0))
            c(m$estimate[[1]], m$estimate[[2]])
        },
        error = function(err) {
            cat("Can't fit noise distribution. Using moments instead.")
            return(c(a_mom, b_mom))
        }) # END tryCatch

    a_0_mass <- coeffs[1]
    b_0_mass <- coeffs[2]
    cat("MASS noise","a_0", a_0_mass, "b_0", b_0_mass, "\n")

    editing <-
        filter(all_sites, .data$is_editing == TRUE) %>%
        mutate(edt_count = ifelse(.data$reference == 'A', .data$G, .data$C), edt_fraction = .data$edt_count / (.data$edt_count + .data$reference_counts))
    fractions_for_models <- select(editing, .data$edt_fraction)

    mu <- mean(editing$edt_fraction)
    var <- var(editing$edt_fraction)
    a_mom <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    b_mom <- a_mom * (1 / mu - 1)

    # fit distribution
    coeffs <- tryCatch(
        {
            m <- fitdistr(fractions_for_models$edt_fraction, 'beta', start = list(shape1 = a_mom, shape2 = b_mom), lower = c(0, 0))
            c(m$estimate[[1]], m$estimate[[2]])
        },
        error = function(err) {
            cat("Can't fit edt distribution. Using moments instead. \n")
            return(c(a_mom, b_mom))
        }) # END tryCatch

    a_edt_mass <- coeffs[1]
    b_edt_mass <- coeffs[2]
    cat("MASS edt", "a_edt =", a_edt_mass, "b_edt =", b_edt_mass, "\n")

    edt_number <- nrow(editing)
    edt_coverage <- editing$coverage
    edt_ref_count <- editing$reference_counts
    edt_edt_count <- editing$edt_count

    ############################# PARALLEL CALCULATIONS #################

    k <- threads
    cat("number of edt sites = ", edt_number, "\n")
    cat("One iteration of L function with parallel\n")
    cat("system.time", system.time(y <- L(a_edt_mass, b_edt_mass, a_0_mass, b_0_mass)), "\n")
    cat("L = ", y, "\n")

    cat("MLE with start MASS\n")
    cat(
        "system.time",
        system.time(
            p2 <- mle(
                L,
                start = list(
                    a_edt = a_edt_mass,
                    b_edt = b_edt_mass,
                    a_0 = a_0_mass,
                    b_0 = b_0_mass,
                    k = k,
                    noise = noise,
                    edt_number = edt_number,
                    edt_coverage = edt_coverage,
                    edt_ref_count = edt_ref_count,
                    edt_edt_count = edt_edt_count,
                    p = p),
                lower = c(0.0000001, 0.0000001, 0.0000001, 0.0000001))),
        "\n")

    cat("MLE+MASS coef\n")
    cat(coef(p2), "\n")
    coef(p2)

    ### Does this work?
    a_edt <- coef(p2)[1]
    b_edt <- coef(p2)[2]
    a_0 <- coef(p2)[3]
    b_0 <- coef(p2)[4]

    J_t <- double()
    J_1 <- double()
    I <- double()
    p_edt <- double()
    PEP <- double()
    t_0 <- 0.15

    I <- unlist(mclapply(
        1:edt_number,
        I_1,
        a_edt = a_edt,
        b_edt = b_edt,
        a_0 = a_0,
        b_0 = b_0,
        edt_coverage = edt_coverage,
        edt_ref_count = edt_ref_count,
        edt_edt_count = edt_edt_count,
        f = f,
        p = p,
        mc.cores = k))
    J_1 <- unlist(mclapply(
        1:edt_number,
        L_A,
        a_edt = a_edt,
        b_edt = b_edt,
        a_0 = a_0,
        b_0 = b_0,
        edt_coverage = edt_coverage,
        edt_ref_count = edt_ref_count,
        edt_edt_count = edt_edt_count,
        f = f,
        p = p,
        mc.cores = k))
    J_coef <- (edt_coverage - edt_ref_count - edt_edt_count) * log(2) + lfactorial(edt_coverage) - (edt_coverage - edt_ref_count) * log(3) -
        lfactorial(edt_ref_count) - lfactorial(edt_edt_count) - lfactorial(edt_coverage - edt_ref_count - edt_edt_count)
    p_edt <- exp(I - J_1)

    J_t <- unlist(mclapply(
        1:edt_number,
        J_in,
        t = t_0,
        a_edt = a_edt,
        b_edt = b_edt,
        a_0 = a_0,
        b_0 = b_0,
        edt_coverage = edt_coverage,
        edt_ref_count = edt_ref_count,
        edt_edt_count = edt_edt_count,
        f = f,
        p = p,
        mc.cores = k))

    PEP <- exp(J_t - J_1)
    J_t <- J_coef + J_t
    J_t_max <- max(J_t)
    J_t <- J_t - J_t_max
    J_1 <- J_coef + J_1
    J_1_max <- max(J_1)
    J_1 <- J_1 - J_1_max

    editing <- editing %>%
        mutate(p_edt = p_edt, PEP = PEP, J_1 = J_1, J_t = J_t) %>%
        arrange(PEP) %>%
        mutate(bv = exp(log(cumsum(exp(.data$J_t))) + .data$J_t_max - log(cumsum(exp(.data$J_1))) - .data$J_1_max))

    editing
}
