#' @title gPlot
#' @description
#' Plot \eqn{g(\bold{\beta^{(t+1)}},\theta^{(t+1)},\sigma^{(t+1)}|\bold{\beta^{(t)}},\theta^{(t)},
#' \sigma^{(t)})} over iterations at a specified spike parameter.
#' \eqn{g(\bold{\beta},\theta,\sigma|\bold{\beta^{(t)}},\theta^{(t)},\sigma^{(t)})} is
#' the expected value of the log likelihood function with respect to the latent
#' variables \eqn{X} and \eqn{\bold{\gamma}}, conditioning on current
#' estimated parameters.
#'
#' @param fit_obj Fitted object from function 'MUSS'.
#' @param spike_param The value of the spike parameter that needed to be specified,
#' which must be in \code{spike_params} designated in `MUSS` function.
#' Each \code{spike_param} corresponds to a unique gPlot.
#' If default \code{spike_params} is used in `MUSS`, one can check the values of
#' \code{spike_params} by calling `$spike_params`
#' from the returned object.
#' @export

gPlot = function(fit_obj, spike_param){

  spike_index = which(fit_obj$spike_params == spike_param)
  g_vals = fit_obj$g_List[[spike_index]]
  iter_num = fit_obj$iter_nums[spike_index]

  theta_path = fit_obj$theta_path
  ggplot()+
    geom_line(aes(x = 1:iter_num, y = g_vals, colour = "g"),size = 1)+
    geom_point(aes(x = 1:iter_num, y = g_vals, colour = "g"),size = 3)+
    xlab("iteration")+ylab("g")+
    labs(title ="g over iterations")+
    theme(plot.title = element_text(size=13,hjust=0.5))+
    scale_color_manual("",values = c("g" = "lightblue3"))
}

