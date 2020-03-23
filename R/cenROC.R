#' Estimation of time-dependent ROC curve for right censored survival data
#'
#' @description  This function computes the time-dependent ROC curve for the right censored survival data using the cumulative sensitivity and dynamic specificity definitions.
#'  The ROC curves can be either empirical (non-smoothed) or smoothed with/wtihout boundary correction. It also calculates the time-dependent area under the ROC curve (AUC).
#' @usage cenROC(Y, M, censor, t, U = NULL, h = NULL, bw = "NR", len = 151, method = "tra",
#'     ktype = "normal", ktype1 = "normal", plot = TRUE)
#' @param Y The numeric vector of event-times or observed times.
#' @param M The numeric vector of marker values for which the time-dependent ROC curves is computed.
#' @param censor The censoring indicator, \code{1} if event, \code{0} otherwise.
#' @param t A scaler time point at which we want to compute the time-dependent ROC curve.
#' @param U The vector of grid points where the ROC curve is estimated. The default is a sequence of \code{151} numbers between \code{0} and \code{1}.
#' @param len The length of the grid points. The default is \code{151}.
#' @param bw A character string specifying the bandwidth estimation method. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and the cross-validation "\code{CV}". The default is the "\code{NR}" normal reference method. It is also possible to use a numeric value.
#' @param h A scaler for the bandwidth of Beran's weight calculaions. The defualt is using the method of Sheather and Jones (1991).
#' @param method The method of ROC curve estimation. The possible options are "\code{emp}" emperical metod; "\code{untra}" smooth without boundary correction and "\code{tra}" is smooth ROC curve estimation with boundary correction. The default is the "\code{tra}" smooth ROC curve estimate with boundary correction.
#' @param ktype A character string giving the type kernel distribution to be used for smoothing the ROC curve: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @param ktype1 A character string specifying the desired kernel needed for Beran weight calculation. The possible options are "\code{normal}", "\code{epanechnikov}", "\code{tricube}", "\code{boxcar}", "\code{triangular}", or "\code{quartic}". The defaults is "\code{normal}" kernel density.
#' @param plot The logical parameter to see the ROC curve plot. The default is \code{TRUE}.
#' @details The empirical (non-smoothed) ROC estimate and the smoothed ROC estimate with/without boundary correction can be obtained using this function.
#' The smoothed ROC curve estimators require selecting two bandwidth parametrs: one for Beran’s weight calculation and one for smoothing the ROC curve.
#' For the latter, three data-driven methods: the normal reference "\code{NR}", the plug-in "\code{PI}" and the cross-validation "\code{CV}" were implemented.
#' To select the bandwidth parameter needed for Beran’s weight calculation, by default, the plug-in method of Sheather and Jones (1991) is used but it is also possible to use numeric value.
#' The time-dependent AUC can be computed either using the direct method (i.e. \eqn{\widehat{AUC}_{t}=1-{n^{-1}\sum_{i=1}^n\hat{\mathcal{W}}_i\hat{Z}_{i}}}) or using the numerical integration (i.e. \eqn{\widehat{AUC}_{t}=\int_{0}^1\widehat{ROC}_{t}(u)du}).
#' The details about these methods can be found in the paper of Beyene and El Ghouch (2020).
#' @return Returns the following items:
#' @return    \code{ROC      } The vector of estimated ROC values. These will be numeric numbers between \code{0} and \code{1}.
#' @return    \code{U        } The vector of grid points used.
#' @return    \code{AUC      } The estimated area under the ROC curve at a given \code{t} using direct method.
#' @return    \code{AUC1     } The estimated area under the ROC curve at a given \code{t} using numerical integration method.
#' @return    \code{bw       } The computed value of bandwidth. For the empirical method this is always \code{1}.
#' @importFrom stats pnorm qnorm quantile approx bw.SJ
#' @importFrom condSURV Beran
#' @importFrom graphics abline legend segments text 
#' @examples library(cenROC)
#'
#' data(mayo)
#' data <- mayo[ ,c( "time","censor","mayoscore5" )]
#' t <- 365*6
#'
#' resu <- cenROC(Y=data$time, M=data$mayoscore5, censor=data$censor, t=t, U=NULL,
#'          len=151, bw="PI", method="tra", plot=TRUE)
#'
#'  resu$AUC
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data.
#' @references Sheather, S. J. and Jones, M. C. (1991). A Reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society}. Series B (Methodological) 53(3): 683–690.
#' @export

cenROC <- function(Y, M, censor, t, U = NULL, h = NULL, bw = "NR", len = 151, method = "tra", ktype = "normal", ktype1 = "normal", plot = TRUE)
  {
    if (is.null(U)) {U <- seq(0, 1, length.out = len)}
    if (!is.vector(Y, mode = "numeric") |
        !is.vector(M, mode = "numeric") |
        !is.vector(censor, mode = "numeric"))
      print("Error! all numeric vectors Y, M and censor should be specified")
    else{
      Dt <- Csurv(Y = Y, M = M, censor = censor, t = t, h = h, kernel = ktype1)$positive
      estim <- RocFun(U = U, D = Dt, M = M, method = method, bw = bw, ktype = ktype)
      ROC <- estim$roc
      AUC <- 1 - estim$auc
      AUC1 <- integ(U, ROC, method = "trap")
    }
    if (plot == "TRUE") {
      plot(c(0, U, 1), c(0, (ROC), 1), type = "l", lwd = 2, col.lab = "blue", col = "blue",
           xlab = "False positive rate", ylab = "True positive rate")
      abline(coef = c(0,1))
    }

    return(list(AUC = AUC, AUC1 = AUC1, ROC = ROC, U = U, bw = estim$bw))
  }
