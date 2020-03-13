#' @export
calcSigma <- function(X1,X2,l=12) {
    Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
    for (i in 1:nrow(Sigma)) {
        for (j in 1:ncol(Sigma)) {
            Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
        }
    }
    return(Sigma)
}

#' @importFrom MASS mvrnorm
#' @export
fitGP <- function(in_df, x_range, n_samples, l = 12){
    x_range <- sort(unique(c(x_range, in_df$x)))
    
    mean_df <- in_df %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(y = mean(y, na.rm = TRUE)) %>%
        dplyr::ungroup()
    sigma <- calcSigma(x_range,x_range, l = l)
    # Calculate the covariance matrices
    # using the same temperature_range values as above
    k_xx <- calcSigma(mean_df$x, mean_df$x, l = l)
    k_xxs <- calcSigma(mean_df$x, x_range, l = l)
    k_xsx <- calcSigma(x_range,mean_df$x, l = l)
    k_xsxs <- calcSigma(x_range, x_range, l = l)
    # get the standard deviation of the noise
    sigma_n <- sapply(sort(unique(in_df$x)), function(xi){
        sd(filter(in_df, x == xi)$y, na.rm = TRUE)
    })
    # Recalculate the mean and covariance functions
    f_star_bar <- k_xsx %*% solve(k_xx + sigma_n^2*diag(1, ncol(k_xx))) %*% mean_df$y
    cov_f_star <- k_xsxs - k_xsx %*% solve(k_xx + sigma_n^2*diag(1, ncol(k_xx))) %*% k_xxs
    
    # Recalulate the sample functions
    values <- matrix(rep(0,length(x_range) * n_samples), ncol = n_samples)
    for (i in seq_len(n_samples)) {
        values[,i] <- mvrnorm(1, f_star_bar, cov_f_star)
    }
    values <- cbind(x = x_range, as.data.frame(values))
    values <- gather(values, key, value, -x)
    
    mean_vals <- values %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(value = mean(value)) %>%
        dplyr::ungroup()
    
    gg <- ggplot(values,aes(x = x,y = value)) +
        geom_line(aes(group = key), colour="grey80", alpha = 0.25) +
        geom_line(colour="grey20",data = mean_vals) +
        geom_errorbar(data = mean_df, aes(x = x, y = NULL,
                                          ymin = y - 2 * sigma_n,
                                          ymax = y + 2 * sigma_n), width=0.2) +
        geom_point(data = mean_df,aes(x = x ,y = y)) +
        theme_bw() +
        xlab("input, x")
    
    residuals <- (mean_vals %>%
                      left_join(in_df, by = "x") %>%
                      mutate(res = value - y))$res
    
    return(list("mean values" = mean_vals,
                "full values" = values,
                "plot" = gg,
                "residuals" = residuals))
}

#' @export
evalGP <- function(gp_model, in_df, model_type = "H0"){
    out_tab <- tibble(
        nObs = nrow(in_df),
        rss = sum(gp_model$residuals^2, na.rm = TRUE)
    )
    colnames(out_tab) <- paste0(colnames(out_tab), model_type)
    return(out_tab)
}

#' @export
fitH0ModelTppPhGP <- function(in_df,
                              x_range = seq(37, 67, len = 50),
                              n_samples = 500, l = 5){
    out_df <- in_df %>%
        mutate(x = temperature, y = rel_value) %>% 
        group_by(gene_name, mod_sequence) %>%
        do({
            h0_model = fitGP(., x_range = x_range, 
                             n_samples = n_samples, l = l)
            evalGP(h0_model, ., model_type = "H0")
        }) %>%
        ungroup
    
    return(out_df)
}

#' @export
fitH1ModelTppPhGP <- function(in_df,
                              x_range = seq(37, 67, len = 50),
                              n_samples = 500, l = 5){
    out_df <- in_df %>%
        mutate(x = temperature, y = rel_value) %>% 
        group_by(gene_name, mod_sequence, var) %>%
        do({
            #print(.$gene_name[1])
            h1_model = fitGP(., x_range = x_range, 
                             n_samples = n_samples, l = l)
            evalGP(h1_model, ., model_type = "H1")
        }) %>%
        group_by(gene_name, mod_sequence) %>% 
        summarize(nObsH1 = sum(nObsH1),
                  rssH1 = sum(rssH1)) %>% 
        ungroup
    
    return(out_df)
}

#' @export
fitH0ModelTppGP <- function(in_df,
                              x_range = seq(37, 67, len = 50),
                              n_samples = 500, l = 5){
    out_df <- in_df %>%
        mutate(x = temperature, y = rel_value) %>% 
        group_by(gene_name) %>%
        do({
            h0_model = fitGP(., x_range = x_range, 
                             n_samples = n_samples, l = l)
            evalGP(h0_model, ., model_type = "H0")
        }) %>%
        ungroup
    
    return(out_df)
}

#' @export
fitH1ModelTppGP <- function(in_df,
                              x_range = seq(37, 67, len = 50),
                              n_samples = 500, l = 5){
    out_df <- in_df %>%
        mutate(x = temperature, y = rel_value) %>% 
        group_by(gene_name, var) %>%
        do({
            #print(.$gene_name[1])
            h1_model = fitGP(., x_range = x_range, 
                             n_samples = n_samples, l = l)
            evalGP(h1_model, ., model_type = "H1")
        }) %>%
        group_by(gene_name) %>% 
        summarize(nObsH1 = sum(nObsH1),
                  rssH1 = sum(rssH1)) %>% 
        ungroup
    
    return(out_df)
}